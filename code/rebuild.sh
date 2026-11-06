#!/bin/bash
# Rebuild the Fortran bindings of the GMRES copy (msspec_gmres/src).
# Run from any directory. Needs the venv at code/_venv (see BUILD.md).
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC="$ROOT/msspec_gmres/src"
VENV="$ROOT/_venv"

source "$VENV/bin/activate"
# BUILDDIR must be absolute (meson backend, numpy >= 2); see BUILD.md section 4.
echo "[1/2] Building phagen..."
cd "$SRC/msspec/phagen/fortran"
make all BUILDDIR="$SRC/msspec/phagen/fortran/build" VERBOSE=1

echo "[2/2] Building spec (phd_se/mi/ce, eig_mi/pw, comp_curves)..."
cd "$SRC/msspec/spec/fortran"
make all BUILDDIR="$SRC/msspec/spec/fortran/build" VERBOSE=1

echo "Done. The tests load the copy with PYTHONPATH, so no install is needed:"
echo "  cd $ROOT && PYTHONPATH=msspec_gmres/src python tests/test_tree_gmres.py"
echo "Then run the release gate: $ROOT/check_release.sh"
