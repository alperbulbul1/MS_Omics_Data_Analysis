#!/usr/bin/env bash
# configure.sh — point the release at your own copy of the data and your own interpreters.
#
# Every script in this release was written against absolute paths on the authors' machine. Rather
# than rewrite 300 scripts by hand (and risk changing behaviour), the release ships them with the
# two machine-specific prefixes replaced by placeholders:
#
#     __MS_GEO_ROOT__   the project root that holds the data directories
#     __PYTHON_BIN__    the Python interpreter with the environment from env/requirements.txt
#
# This script substitutes both, in place, across scripts/. It is idempotent: running it again with
# different values re-substitutes from the current values, and running it with the same values is a
# no-op. Verify the result with --check before running anything expensive.
#
# Usage:
#   ./configure.sh /path/to/MS_GEO_data /path/to/python
#   ./configure.sh --check
#
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"
STATE="$HERE/.configured"

# grep exits 1 when it matches nothing, which under `set -e -o pipefail` would abort the script
# on the success case; swallow that here so the count is what decides.
placeholders_left() {
  { grep -rl '__MS_GEO_ROOT__\|__PYTHON_BIN__' "$SCRIPTS" 2>/dev/null || true; } | wc -l | tr -d ' '
}

if [[ "${1:-}" == "--check" ]]; then
  n=$(placeholders_left)
  if [[ -f "$STATE" ]]; then
    echo "configured with:"; cat "$STATE"
  fi
  if [[ "$n" -gt 0 ]]; then
    echo "NOT configured: $n file(s) still contain placeholders. Run ./configure.sh ROOT PYTHON"
    exit 1
  fi
  echo "OK: no placeholders remain in scripts/"
  # a configured tree must not contain the authors' own paths either
  # the authors' home prefix, split so this guard does not match itself
  AUTHOR_PREFIX='/Users/'"alperbulbul"
  if grep -rq "$AUTHOR_PREFIX" "$SCRIPTS" 2>/dev/null; then
    echo "WARNING: author-specific paths are still present:"
    grep -rl "$AUTHOR_PREFIX" "$SCRIPTS" | head
    exit 1
  fi
  echo "OK: no author-specific absolute paths in scripts/"
  exit 0
fi

if [[ $# -lt 2 ]]; then
  echo "usage: $0 /path/to/MS_GEO_data /path/to/python" >&2
  echo "       $0 --check" >&2
  exit 2
fi

ROOT="$(cd "$1" && pwd)"           # must exist; the data lives here
PY="$2"
[[ -x "$PY" ]] || { echo "not executable: $PY" >&2; exit 2; }

# current values, so re-running with new values works
OLD_ROOT="__MS_GEO_ROOT__"; OLD_PY="__PYTHON_BIN__"
if [[ -f "$STATE" ]]; then
  OLD_ROOT="$(sed -n 's/^root=//p' "$STATE")"
  OLD_PY="$(sed -n 's/^python=//p' "$STATE")"
fi

n=0
while IFS= read -r -d '' f; do
  # BSD and GNU sed differ on -i; write to a temp file instead
  tmp="$f.tmp.$$"
  sed -e "s|${OLD_ROOT}|${ROOT}|g" -e "s|${OLD_PY}|${PY}|g" "$f" > "$tmp" && mv "$tmp" "$f"
  n=$((n+1))
done < <(find "$SCRIPTS" -type f \( -name '*.R' -o -name '*.py' -o -name '*.sh' -o -name '*.ipynb' \) -print0)

printf 'root=%s\npython=%s\n' "$ROOT" "$PY" > "$STATE"
echo "rewrote $n file(s)"
echo "  root   -> $ROOT"
echo "  python -> $PY"
echo "run ./configure.sh --check to verify"
