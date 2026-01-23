#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# cibuildwheel macOS repair helper
# Usage: cibuildwheel passes {wheel} and {dest_dir} — some setups pass only the wheel path as $1.
WHEEL="$1"
DEST_DIR="${2:-$(dirname "$WHEEL")}"   # fallback if {dest_dir} isn't provided

echo "Repairing wheel: $WHEEL -> $DEST_DIR"

# helper function: print to stderr
err() { echo "$@" >&2; }

TMP="$(mktemp -d)"
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT

# unzip into tmp and work there
pushd "$TMP" >/dev/null
unzip -q "${GITHUB_WORKSPACE:-$PWD}/$WHEEL" || unzip -q "$WHEEL" || { err "Failed to unzip $WHEEL"; exit 2; }

# Find package root(s) that contain a .dylibs folder or .so files.
# We'll handle each package we find; common layout is <pkg>/.dylibs
PKG_CANDIDATES=()
while IFS= read -r -d '' p; do PKG_CANDIDATES+=("$p"); done < <(find . -type d -name ".dylibs" -print0 2>/dev/null)

# If none found, try heuristics: directories that contain .so files at top-level
if [ "${#PKG_CANDIDATES[@]}" -eq 0 ]; then
  while IFS= read -r -d '' dir; do
    PKG_CANDIDATES+=("$dir")
  done < <(find . -maxdepth 2 -type f -name "*.so" -exec dirname {} \; -print0 2>/dev/null)
fi

# Deduplicate candidates (dirname may appear multiple times)
PKG_UNIQ=()
for p in "${PKG_CANDIDATES[@]}"; do
  dir="$(dirname "$p")"
  case " ${PKG_UNIQ[*]} " in
    *" $dir "*) ;;
    *) PKG_UNIQ+=("$dir");;
  esac
done

# If still empty, fallback to "uedge" (your package) if present
if [ "${#PKG_UNIQ[@]}" -eq 0 ] && [ -d "./uedge" ]; then
  PKG_UNIQ+=("./uedge")
fi

if [ "${#PKG_UNIQ[@]}" -eq 0 ]; then
  echo "No package directories detected; nothing to repair." >&2
else
  echo "Package roots to process: ${PKG_UNIQ[*]}"
fi

# For each package root:
for pkg in "${PKG_UNIQ[@]}"; do
  # ensure .dylibs exists
  DDIR="${pkg}/.dylibs"
  mkdir -p "$DDIR"

  # Move any top-level dylibs into .dylibs (common in some wheels)
  # also move any direct matches (filename) that live at top-level
  find . -maxdepth 2 -type f -name '*.dylib' -print0 2>/dev/null | while IFS= read -r -d '' lib; do
    # skip if already under this package's .dylibs
    case "$lib" in
      ./"${DDIR}"/*) continue ;;
    esac
    mv -f "$lib" "${DDIR}/" || true
  done

  # Collect the set of bundled dylib basenames
  bundled=()
  while IFS= read -r -d '' f; do bundled+=("$(basename "$f")"); done < <(find "$DDIR" -maxdepth 1 -type f -name '*.dylib' -print0 2>/dev/null || true)

  # Step A: set each dylib's ID to @rpath/<name>
  for d in "${bundled[@]}"; do
    full="${DDIR}/${d}"
    if [ -f "$full" ]; then
      newid="@rpath/${d}"
      echo "Setting LC_ID_DYLIB: $full -> $newid"
      install_name_tool -id "$newid" "$full" || { err "install_name_tool -id failed for $full"; }
    fi
  done

  # Step B: add LC_RPATH '@loader_path/.dylibs' to every consumer binary and bundled dylib
  # and simultaneously rewrite any absolute build-path references to @rpath/<name>
  # For more robust resolution we will make binaries use @rpath/<name> and have rpath point to loader_path/.dylibs
  find "$pkg" -type f \( -name '*.so' -o -name '*.dylib' \) -print0 | while IFS= read -r -d '' bin; do
    # add rpath (ignore error if already present)
    echo "Ensuring LC_RPATH @loader_path/.dylibs on $bin"
    install_name_tool -add_rpath "@loader_path/.dylibs" "$bin" 2>/dev/null || true

    # iterate deps and rewrite those that match bundled basenames or absolute build paths
    otool -L "$bin" 2>/dev/null | sed -n '2,$p' | awk '{print $1}' | while read -r dep; do
      [ -z "$dep" ] && continue
      base="$(basename "$dep")"
      # If dep references a bundled file by name or is an absolute path to a bundled file, rewrite to @rpath/<base>
      if printf '%s\n' "${bundled[@]}" | grep -xq "$base"; then
        # get current dep string and change to @rpath/<base> if needed
        if [ "$dep" != "@rpath/${base}" ]; then
          echo "  Changing $dep -> @rpath/${base} in $bin"
          install_name_tool -change "$dep" "@rpath/${base}" "$bin" 2>/dev/null || err "warning: change failed for $bin ($dep)"
        fi
      else
        # If dep is an absolute build path (e.g. /DLC/...), and a bundled copy by basename exists somewhere in wheel, rewrite
        if [[ "$dep" == /* ]]; then
          # find any matching basename in wheel (brute force)
          if find . -type f -name "$base" | grep -q .; then
            echo "  Absolute dep $dep appears bundled somewhere; changing -> @rpath/${base}"
            install_name_tool -change "$dep" "@rpath/${base}" "$bin" 2>/dev/null || err "warning: change failed for $bin ($dep)"
          fi
        fi
      fi
    done
  done
done

# Step C: ad-hoc codesign any modified mach-o files to avoid AMFI kills
echo "Codesigning modified Mach-O files (ad-hoc) ..."
# find all Mach-O objects (.so and any .dylib under .)
find . -type f \( -name '*.so' -o -name '*.dylib' \) -print0 | while IFS= read -r -d '' m; do
  echo "  codesign $m"
  codesign --force --sign - --timestamp=none "$m" 2>/dev/null || err "codesign failed (continuing): $m"
done

# Step D: strict post-repair validation
echo "Post-repair validation..."
fail=0

# helper: get LC_ID_DYLIB (last line of otool -D)
for d in $(find . -path '*/.dylibs/*.dylib' -type f 2>/dev/null || true); do
  id="$(otool -D "$d" 2>/dev/null | tail -n 1 | tr -d '\r')"
  if [ -z "$id" ]; then
    err "ERROR: couldn't read LC_ID_DYLIB for $d"; fail=1; continue
  fi
  # id must be @rpath/<name> for packaged dylibs
  case "$id" in
    @rpath/*) ;;
    *)
      err "ERROR: Dylib id is not @rpath/ for $d -> $id"; fail=1;;
  esac
done

# validate dependencies for all shipped binaries
find . -type f \( -name '*.so' -o -path '*/.dylibs/*.dylib' \) -print0 | while IFS= read -r -d '' bin; do
  otool -L "$bin" 2>/dev/null | sed -n '2,$p' | awk '{print $1}' | while read -r dep; do
    [ -z "$dep" ] && continue
    # Reject raw build-machine absolute prefixes (/DLC/...)
    if [[ "$dep" == /DLC/* ]]; then
      err "ERROR: build-path dependency remains: $bin -> $dep"; fail=1
    fi
    # Reject unexpected absolute non-system libs (except /usr /System /Library)
    if [[ "$dep" == /* ]] && [[ "$dep" != /usr/* ]] && [[ "$dep" != /System/* ]] && [[ "$dep" != /Library/* ]]; then
      err "ERROR: unexpected absolute dependency: $bin -> $dep"; fail=1
    fi
    # @rpath deps are OK (we set dylib ids to @rpath and add rpaths)
  done
done

if [ "$fail" -ne 0 ]; then
  err "Post-repair validation FAILED. Inspect output above."
  exit 10
fi

# Repack wheel and move back to destination
OUTNAME="$(basename "$WHEEL")"
# create cleaned wheel in TMP; zip contents (preserve all files)
zip -q -r "/tmp/${OUTNAME}" . || { err "zip failed"; exit 3; }
# move back to DEST_DIR (cibuildwheel expects the repaired wheel at the same path)
mkdir -p "$DEST_DIR"
mv -f "/tmp/${OUTNAME}" "${GITHUB_WORKSPACE:-$PWD}/${WHEEL}" || mv -f "/tmp/${OUTNAME}" "${DEST_DIR}/${OUTNAME}"

popd >/dev/null
echo "Repaired $WHEEL -> placed at ${DEST_DIR}/${OUTNAME}"
exit 0

