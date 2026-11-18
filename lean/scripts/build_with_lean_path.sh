#!/usr/bin/env bash
# Build the Lean development with plain `lean` against a Mathlib v4.35.0 build, via LEAN_PATH.
# Default: the repository's own Lake environment, lean/.lake/packages (create it with
# `cd lean && lake exe cache get`, or use `lake build` directly).  Override: set MATHLIB_PACKAGES to
# the .lake/packages directory of an existing Mathlib v4.35.0 build.  Outputs go to lean/build/ only.
set -euo pipefail
cd "$(dirname "$0")/.."
PKGS="${MATHLIB_PACKAGES:-$PWD/.lake/packages}"
if [ ! -d "$PKGS/mathlib" ]; then
  echo "no Mathlib build found at $PKGS" >&2
  echo "run 'lake exe cache get' in lean/, or set MATHLIB_PACKAGES to an existing .lake/packages directory" >&2
  exit 1
fi
OUT=build/lib
mkdir -p "$OUT/SteklovTriangles"
LP="$PWD/$OUT"
for d in "$PKGS"/*; do LP="$LP:$d/.lake/build/lib/lean"; done
export LEAN_PATH="$LP"
for m in AreaCorollary Reformulation CosBound; do
  echo "== SteklovTriangles.$m"
  lean -o "$OUT/SteklovTriangles/$m.olean" "SteklovTriangles/$m.lean"
done
lean -o "$OUT/SteklovTriangles.olean" SteklovTriangles.lean
echo "== axioms"
lean scripts/print_axioms.lean
