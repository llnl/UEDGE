#!/usr/bin/env bash
set -euo pipefail

# tools/tls_audit.sh <target-dir> <package-dir-or-name>
# - target-dir: usually build/lib.macosx-*/ or build/bdist... (pass '.' to scan defaults)
# - package-dir-or-name: package directory inside build (e.g. "uedge") or package name
#
# The script:
# - finds all .so files under build paths and (if present) inside dist/*.whl
# - for each .so, runs `otool -L` and `nm -m` (if available) and writes per-extension reports
# - does NOT fail the build on missing files; it prints warnings instead
#
# Dependency: macOS command-line tools: otool, nm, unzip
# NOTE: uses /usr/bin/nm (classic nm) when available; falls back to whatever `nm` in PATH.

BUILD_ROOT="${1:-build}"
PKGNAME="${2:-uedge}"
OUTDIR="${BUILD_ROOT}/linkmaps"
TMPDIR_BASE="$(mktemp -d "${PWD}/.tls_audit_tmp.XXXX")"

cleanup() {
  rm -rf "$TMPDIR_BASE"
}
trap cleanup EXIT

mkdir -p "$OUTDIR"

echoblue() { printf "\033[1;34m%s\033[0m\n" "$1"; }
echored()  { printf "\033[1;31m%s\033[0m\n" "$1"; }
echoyellow(){ printf "\033[1;33m%s\033[0m\n" "$1"; }

echoblue "== TLS audit: scanning for .so files under ${BUILD_ROOT} and wheels in dist/ ..."

# collect candidate .so locations (generalized)
declare -a SO_CANDIDATES=()

# 1) build/lib.*/*.so (flexible)
while IFS= read -r -d '' sfile; do SO_CANDIDATES+=("$sfile"); done < <(find "${BUILD_ROOT}" -type f -name '*.so' -print0 2>/dev/null)

# 2) src/package .so files if built in-tree
if [ -d "build/lib" ]; then
  while IFS= read -r -d '' sfile; do SO_CANDIDATES+=("$sfile"); done < <(find build/lib -type f -name '*.so' -print0 2>/dev/null)
fi

# 3) dist/*.whl (inspect wheels)
shopt -s nullglob
for whl in dist/*.whl; do
  echoyellow "== Inspecting wheel: $whl"
  wtmp="${TMPDIR_BASE}/wheel-$(basename "$whl" .whl)"
  mkdir -p "$wtmp"
  unzip -q -o "$whl" -d "$wtmp"
  while IFS= read -r -d '' sfile; do
    # prefix with 'wheel:' so user knows origin
    SO_CANDIDATES+=("wheel:${wtmp}/${sfile}")
  done < <(find "$wtmp" -type f -name '*.so' -print0)
done
shopt -u nullglob

if [ ${#SO_CANDIDATES[@]} -eq 0 ]; then
  echoyellow "== TLS audit: no .so files found under ${BUILD_ROOT} or inside dist/*.whl (warning, not fatal)."
  exit 0
fi

echoblue "== TLS audit: found ${#SO_CANDIDATES[@]} candidate .so files."

# helper to run commands but write to per-ext map file
run_for_so() {
  local src="$1"
  local label="$2"
  local mapfile="${OUTDIR}/${label}.map"
  echo "----" > "$mapfile"
  echo "# Source: $src" >> "$mapfile"
  echo "# Report generated: $(date -u)" >> "$mapfile"
  echo "" >> "$mapfile"

  # determine actual filesystem path
  if [[ "$src" == wheel:* ]]; then
    local realpath="${src#wheel:}"
  else
    local realpath="$src"
  fi

  if [ ! -e "$realpath" ]; then
    echo "WARNING: file disappeared: $realpath" | tee -a "$mapfile"
    return 0
  fi

  {
    echo "=== otool -L output ==="
    if command -v otool >/dev/null 2>&1; then
      otool -L "$realpath" 2>/dev/null || echo "(otool -L failed)"
    else
      echo "(otool not found)"
    fi
    echo ""
    echo "=== nm -m (checking for emutls symbols) ==="
    if command -v nm >/dev/null 2>&1; then
      # prefer /usr/bin/nm if available
      if [ -x /usr/bin/nm ]; then
        /usr/bin/nm -m "$realpath" 2>/dev/null || echo "(nm failed)"
      else
        nm -m "$realpath" 2>/dev/null || echo "(nm failed)"
      fi
    else
      echo "(nm not found)"
    fi
    echo ""
    echo "=== greps (emutls presence) ==="
    # try to detect emutls descriptors/definitions
    if command -v nm >/dev/null 2>&1; then
      if ( ( [ -x /usr/bin/nm ] && /usr/bin/nm -m "$realpath" 2>/dev/null ) || nm -m "$realpath" 2>/dev/null ); then
        if ( ( [ -x /usr/bin/nm ] && /usr/bin/nm -m "$realpath" 2>/dev/null ) | egrep -q "___emutls_v\\." ); then
          echo "[TLS] emutls_v symbols present"
        else
          echo "[TLS] no emutls_v symbols found"
        fi
      else
        echo "(nm listing failed; can't check emutls)"
      fi
    else
      echo "(nm missing; can't check emutls)"
    fi
  } >> "$mapfile"

  printf "wrote %s\n" "$mapfile"
}

# iterate candidates
count=0
for candidate in "${SO_CANDIDATES[@]}"; do
  # construct a label safe for filenames
  # source path -> label by replacing slashes with '__' and stripping leading dirs
  if [[ "$candidate" == wheel:* ]]; then
    path="${candidate#wheel:}"
    label="wheel__$(basename "$path")"
  else
    path="$candidate"
    # prefer relative path from project root if possible
    rel="$(realpath --relative-to="${PWD}" "$path" 2>/dev/null || echo "$path")"
    label="${rel//\//__}"
  fi
  # limit filename length a bit
  label="$(echo "$label" | sed 's/[^A-Za-z0-9_.-]/_/g' | cut -c1-200)"
  run_for_so "$candidate" "$label"
  count=$((count+1))
done

echoblue "== TLS audit done. Reports written to ${OUTDIR} (count=${count})."
exit 0

