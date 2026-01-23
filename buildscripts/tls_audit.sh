#!/usr/bin/env bash
set -euo pipefail

# tools/tls_audit.sh <target-dir> <package-dir-or-name> [--fail-on-dup]
#
# - target-dir: directory to scan for built .so files (default: build)
# - package-dir-or-name: package directory inside build (for labeling) (default: uedge)
# - --fail-on-dup: optional flag; when present, script will exit non-zero
#                  if more than one extension defines any ___emutls_v. symbols
#
# Also supports environment variable FAIL_ON_DUP=1 to enable the same strict behaviour.
#
# The script:
# - discovers .so files under build and inside dist/*.whl
# - writes per-extension reports to build/linkmaps/<label>.map
# - by default prints warnings but exits 0
# - if --fail-on-dup or FAIL_ON_DUP=1 is set, it fails with code 2 when
#   multiple extensions define emutls descriptors (___emutls_v.)
#
# Dependencies: macOS CLI tools: otool, nm, unzip (script prefers /usr/bin/nm)
# NOTE: Use TAB characters for Makefile recipe lines that call this script.

BUILD_ROOT="${1:-build}"
PKGNAME="${2:-uedge}"
shift 2 || true

# parse optional flags
FAIL_ON_DUP_FLAG=0
for arg in "$@"; do
  case "$arg" in
    --fail-on-dup) FAIL_ON_DUP_FLAG=1 ;;
    *) echo "WARNING: unknown arg '$arg' (ignored)";;
  esac
done

# environment var also enables fail mode
if [ "${FAIL_ON_DUP:-0}" = "1" ] || [ "${FAIL_ON_DUP:-0}" = "true" ]; then
  FAIL_ON_DUP_FLAG=1
fi

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
declare -a SO_CANDIDATES=()

# flexible discovery: look under BUILD_ROOT (if exists)
if [ -d "${BUILD_ROOT}" ]; then
  while IFS= read -r -d '' sfile; do SO_CANDIDATES+=("$sfile"); done < <(find "${BUILD_ROOT}" -type f -name '*.so' -print0 2>/dev/null)
fi

# Also check common build/lib.* locations (compatibility)
while IFS= read -r -d '' sfile; do SO_CANDIDATES+=("$sfile"); done < <(find build -type f -name '*.so' -print0 2>/dev/null || true)

# inspect wheels in dist/*.whl (wheel-aware)
shopt -s nullglob
for whl in dist/*.whl; do
  echoyellow "== Inspecting wheel: $whl"
  wtmp="${TMPDIR_BASE}/wheel-$(basename "$whl" .whl)"
  mkdir -p "$wtmp"
  unzip -q -o "$whl" -d "$wtmp"
  # find .so inside wheel tree and prefix with wheel:
  while IFS= read -r -d '' sfile; do
    SO_CANDIDATES+=("wheel:${wtmp}/${sfile}")
  done < <(find "$wtmp" -type f -name '*.so' -print0)
done
shopt -u nullglob

if [ ${#SO_CANDIDATES[@]} -eq 0 ]; then
  echoyellow "== TLS audit: no .so files found under ${BUILD_ROOT} or inside dist/*.whl (warning, not fatal)."
  exit 0
fi

echoblue "== TLS audit: found ${#SO_CANDIDATES[@]} candidate .so files."

# helper: run per-so report
run_for_so() {
  local src="$1"
  local label="$2"
  local mapfile="${OUTDIR}/${label}.map"
  echo "----" > "$mapfile"
  echo "# Source: $src" >> "$mapfile"
  echo "# Report generated: $(date -u)" >> "$mapfile"
  echo "" >> "$mapfile"

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
    echo "=== nm -m (symbol listing) ==="
    if command -v nm >/dev/null 2>&1; then
      if [ -x /usr/bin/nm ]; then
        /usr/bin/nm -m "$realpath" 2>/dev/null || echo "(nm failed)"
      else
        nm -m "$realpath" 2>/dev/null || echo "(nm failed)"
      fi
    else
      echo "(nm not found)"
    fi
    echo ""
    echo "=== emutls detection ==="
    # check for any ___emutls_v.<name> lines that are NOT marked (undefined)
    if command -v nm >/dev/null 2>&1; then
      if [ -x /usr/bin/nm ]; then
        nm -m "$realpath" 2>/dev/null | egrep "___emutls_v\\." || true
      else
        nm -m "$realpath" 2>/dev/null | egrep "___emutls_v\\." || true
      fi
    else
      echo "(nm missing; cannot detect emutls symbols)"
    fi
  } >> "$mapfile"

  printf "wrote %s\n" "$mapfile"
}

# iterate candidates, collect list of modules that define emutls symbols
declare -a MODULES_WITH_EMUTLS=()
declare -A EMUTLS_SYMBOL_MAP  # symbol -> comma-separated module list (for extra diagnostics)

count=0
for candidate in "${SO_CANDIDATES[@]}"; do
  if [[ "$candidate" == wheel:* ]]; then
    path="${candidate#wheel:}"
    label="wheel__$(basename "$path")"
    realpath="$path"
  else
    path="$candidate"
    realpath="$path"
    rel="$(realpath --relative-to="${PWD}" "$path" 2>/dev/null || echo "$path")"
    label="${rel//\//__}"
  fi

  label="$(echo "$label" | sed 's/[^A-Za-z0-9_.-]/_/g' | cut -c1-200)"
  run_for_so "$candidate" "$label"

  # check if emutls defs exist (not "(undefined)")
  if command -v nm >/dev/null 2>&1; then
    nm_cmd_output=""
    if [ -x /usr/bin/nm ]; then
      nm_cmd_output="$(/usr/bin/nm -m "$realpath" 2>/dev/null || true)"
    else
      nm_cmd_output="$(nm -m "$realpath" 2>/dev/null || true)"
    fi
    # find defined emutls lines (exclude lines with "(undefined)")
    # Example matching lines:
    #   ___emutls_v.__foo (undefined) external
    #   00000000000abcd ( __TEXT,__text ) non-external ... ___emutls_get_address
    # We'll consider a symbol "defined" when the nm output line does NOT include "(undefined)"
    while IFS= read -r line; do
      # filter only emutls_v symbols
      if echo "$line" | egrep -q "___emutls_v\\."; then
        if ! echo "$line" | egrep -q "\\(undefined\\)"; then
          MODULES_WITH_EMUTLS+=("$label")
          # extract symbol names (the symbol token containing ___emutls_v.<name>)
          # simplistic: find occurrences in the line
          syms="$(echo "$line" | egrep -o "___emutls_v\\.[A-Za-z0-9_]+" || true)"
          for s in $syms; do
            prev="${EMUTLS_SYMBOL_MAP[$s]:-}"
            if [ -z "$prev" ]; then
              EMUTLS_SYMBOL_MAP[$s]="$label"
            else
              EMUTLS_SYMBOL_MAP[$s]="${prev},${label}"
            fi
          done
        fi
      fi
    done <<< "$nm_cmd_output"
  fi

  count=$((count+1))
done

echoblue "== TLS audit done. Reports written to ${OUTDIR} (count=${count})."

# Report summary about emutls definitions
unique_modules=()
declare -A seen_mod
for m in "${MODULES_WITH_EMUTLS[@]}"; do
  if [ -z "${seen_mod[$m]:-}" ]; then
    seen_mod[$m]=1
    unique_modules+=("$m")
  fi
done

if [ ${#unique_modules[@]} -gt 0 ]; then
  echoyellow "== Found emutls definitions in the following extension modules:"
  for m in "${unique_modules[@]}"; do
    printf " - %s\n" "$m"
  done

  echoyellow "== emutls symbol mapping (symbol -> modules):"
  for k in "${!EMUTLS_SYMBOL_MAP[@]}"; do
    printf " %s -> %s\n" "$k" "${EMUTLS_SYMBOL_MAP[$k]}"
  done
else
  echoblue "== No emutls definitions detected in any extension."
fi

if [ "$FAIL_ON_DUP_FLAG" -eq 1 ]; then
  if [ ${#unique_modules[@]} -gt 1 ]; then
    echored "ERROR: more than one extension defines ___emutls_v symbols (count=${#unique_modules[@]})."
    echored "This usually causes TLS runtime collisions (libgomp / libgcc TLS)."
    echored "Modules: ${unique_modules[*]}"
    echored "See ${OUTDIR} for per-extension reports."
    exit 2
  fi
  # optional: also fail if same emutls symbol shows up in more than 1 module
  # (i.e. EMUTLS_SYMBOL_MAP contains a symbol mapped to multiple modules)
  for k in "${!EMUTLS_SYMBOL_MAP[@]}"; do
    IFS=',' read -r -a mods <<< "${EMUTLS_SYMBOL_MAP[$k]}"
    if [ ${#mods[@]} -gt 1 ]; then
      echored "ERROR: emutls symbol $k defined in multiple modules: ${mods[*]}"
      echored "This is a concrete collision of the same descriptor symbol across modules."
      echored "Failing due to --fail-on-dup."
      exit 2
    fi
  done
fi

# success
exit 0

