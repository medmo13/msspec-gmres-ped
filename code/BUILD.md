# MsSpec GMRES copy: build, test and release gates

Folder: `code/`. Copy of the code: `msspec_gmres/` (the only code
this study uses). Venv: `code/_venv`. Verified toolchain: Python 3.14.4,
gfortran 15.2.0, numpy 2.5.4, scipy 1.18.1, ase 3.29.0. Exact hashes and
the BLAS link are in `BUILD_MANIFEST.json`.

## 1. Venv and Python dependencies

```bash
cd code
python3 -m venv _venv
source _venv/bin/activate
pip install --upgrade pip setuptools wheel
pip install numpy scipy ase
```

System needs: `gfortran` and the LAPACK/BLAS runtime (`liblapack.so.3`,
`libopenblas.so.0` as linked in the manifest).

## 2. Build the Fortran bindings

```bash
cd code
./rebuild.sh
```

`rebuild.sh` runs `make all` in `msspec_gmres/src/msspec/phagen/fortran` and in
`msspec_gmres/src/msspec/spec/fortran`, with the venv active. The bindings
are written next to their sources as `.so` files (`libphagen.so`,
`_phd_se_noso_nosp_nosym.so`, `_phd_mi_noso_nosp_nosym.so`, `_phd_gm_noso_nosp_nosym.so` (GMRES),
`_phd_ce_noso_nosp_nosym.so`, `_eig_mi.so`, `_eig_pw.so`, `_comp_curves.so`).
After a source change, rebuild and then regenerate the manifest (section 6).

## 3. Patches in the copy

All patches are in `msspec_gmres/src`. Each one is needed to build with the
current toolchain.

1. Absolute `BUILDDIR`. The numpy >= 2 f2py build uses the meson backend, so
   a relative `-I./build` fails with `Include dir ./build does not exist`.
   `options.mk` sets `BUILDDIR ?= $(abspath build)`, so no command-line
   override is needed. Passing `BUILDDIR=$PWD/build` still works.
2. `msspec/spec/fortran/eig/common/plotfd.f`: `USE DIM_MOD` is present. Without
   it: `Error: Variable 'natm' cannot appear...`.
3. `msspec/VERSION` holds `4.2.0`. The copy has no `.git`, so `git describe`
   gives nothing and the package version cannot be derived.
4. `msspec/calcio.py` (line 673): the removed `np.chararray(..., unicode=True)`
   is replaced by `np.full(..., dtype='U64')`. Without it `run_spec()` fails
   with `module 'numpy' has no attribute 'chararray'`.
5. Makefiles use `$(sort $(wildcard ...))` for every source list in
   `msspec/spec/fortran/*.mk`, so the link order does not depend on the file
   system.
6. GMRES changes (this study): `MS_TREE_GMRES` and `MS_TREE_SOLVE` in
   `msspec/spec/fortran/phd_gm_noso_nosp_nosym/axb_mat_ms_la.f`. See
   `PLAN.md` for the design and the open confirmation on the in-cycle stop.

## 4. Load the copy

The tests load the copy with `PYTHONPATH`; nothing is installed:

```bash
cd code
source _venv/bin/activate
PYTHONPATH=msspec_gmres/src python -c "from msspec.calculator import MSSPEC; print('ok')"
```

## 5. Tests and release gates

Run from `code` with the venv active:

```bash
PYTHONPATH=msspec_gmres/src python tests/test_tree_gmres.py   # C4 regression, expects PASS
PYTHONPATH=msspec_gmres/src python tests/material_check.py    # Ag, Ni, Cu, writes tests/material_check_report.json
```

- `test_tree_gmres.py` checks tree-GMRES against dense LU, that the iterate dump
  does not change the solution, that the last iterate of each RHS equals the
  returned solution, that forced failure recovers through dense LU, and that
  disabling dense LU gives exit 99 with `none_available`.
- `material_check.py` gates: the direct (`inversion`) result is byte-identical to
  the original build; each tree result is within 1e-4 relative L2 of LU. The
  tree-to-tree difference is reported, not gated (see PLAN.md section 11).

## 6. Manifest

```bash
python tests/write_build_manifest.py
```

Writes `BUILD_MANIFEST.json`: source and binary SHA-256, toolchain, flags,
Python and package versions, BLAS link, CPU, kernel and governor. Regenerate it
after any source or binary change.

## 7. Files

- (the plan and status files are not part of the deposit)
- `code/rebuild.sh`: rebuild the bindings.
- `code/tests/`: regression, material check and manifest scripts.
- `code/BUILD_MANIFEST.json`: frozen build record.
- `code/msspec_gmres/`: the GMRES copy.
- `code/_venv/`: working venv (git-ignored).
