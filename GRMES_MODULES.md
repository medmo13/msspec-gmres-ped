# GMRES code: modules and how they work

Written 2026-10-24 from the source at commit `128c866`. Updated 2026-10-31 for the v3 build (no dense fallback; iteration-0 dump; stop switch). This file describes the code, not the results. For the results and the to-do list, read [GRMES_STATUS.md](GRMES_STATUS.md).

All solver code lives in one Fortran file:
`code/msspec_gmres/src/msspec/spec/fortran/phd_gm_noso_nosp_nosym/axb_mat_ms_la.f`.
The Python side is `code/msspec_gmres/src/msspec/calculator.py` and `parameters.py`.

## 1. The problem in one paragraph

PED solves the multiple-scattering system `(I - G T) X = B` for each energy and each emitter.
`G T` is the scattering operator. The code calls the full matrix `A = I - G T`.
The exact method stores `A` and factors it with dense LU. This costs memory and time that grow fast with the cluster size.
The GMRES path never stores `A`. It solves the same system with an iterative method. A block tree factorisation is the preconditioner.

## 2. Module map

| Routine | Role | Called by |
|---|---|---|
| `INV_MAT_MS` | Entry point for one energy. Chooses the solver. A failed tree solve stops the run. | the spectroscopy driver |
| `MS_TREE_SERIES_BATCH` | Builds the tree, the block factors and the pair blocks. Starts the solve. | `INV_MAT_MS` |
| `MS_TREE_GMRES` | Restarted GMRES on `M^-1 A x = M^-1 b`. Writes the solution. | `MS_TREE_SERIES_BATCH` |
| `MS_TREE_SOLVE` | Applies `M^-1` to one or more vectors. | `MS_TREE_SERIES_BATCH`, `MS_TREE_GMRES` |
| `MS_TREE_PAIR_BLOCK` | Computes one exact pair block `A(parent, child)` or `A(child, parent)`. | `MS_TREE_SERIES_BATCH` |
| `MS_MATRIX_ACTION_BLOCK` | Computes `Y = A X` for several vectors, matrix-free. | `MS_TREE_GMRES`, `MS_TREE_SERIES_BATCH` |
| `MS_AXIAL_PAIR_ACTION` | Pair action by rotation, axial translation and inverse rotation. | `MS_MATRIX_ACTION_BLOCK` |
| `MS_TREE_ITER_DUMP` | Writes each iterate to a text file when the dump switch is on. | `MS_TREE_GMRES` |
| `MS_MATRIX_ACTION` | Single-vector version of the matrix action. Used in test code paths. | `INV_MAT_MS` (test only) |

## 3. Call chain

```
calculator.py  (algorithm='gmres', spectroscopy='PED')
   sets environment variables from linear_solver_* parameters
   calls the spectroscopy run inside those variables
        |
        v
INV_MAT_MS(JE)                                 one energy
   reads MSSPEC_TEST_TREE_SERIES, MSSPEC_TREE_SERIES_ORDER
   calls MS_TREE_SERIES_BATCH(JE, N, NRHS, B, ORDER, INFO)
        |
        +-- build tree (Prim), diagonal LU, pair blocks EP and EC
        +-- x0 = M^-1 b                         MS_TREE_SOLVE
        +-- if MSSPEC_TREE_GMRES=1:
        |      MS_TREE_GMRES(...)              restarted GMRES
        |         loop per cycle:
        |            MS_MATRIX_ACTION_BLOCK     w = A v
        |            MS_TREE_SOLVE              w = M^-1 w
        |            Arnoldi and Givens
        |         update x, then true residual  MS_MATRIX_ACTION_BLOCK + MS_TREE_SOLVE
        |      else: finite correction series   (the order study, not the GMRES path)
        |
        +-- if INFO ~= 0 (solver failed):
        |      status line "none_available", STOP 99 (method failure; no dense fallback)
        |
        v
   right-hand sides X are returned to the spectroscopy code
```

## 4. Module details

### 4.1 Python wrapper (`calculator.py`, `parameters.py`)

- `calculator.py` line 204: `algorithm='gmres'` is accepted only for `spectroscopy='PED'`. Any other spectroscopy raises a `ValueError`.
- `calculator.py` lines 472 to 483: when `algorithm == 'gmres'`, the wrapper builds an environment dictionary and runs the spectroscopy with it:
  - `MSSPEC_TEST_TREE_SERIES = 1`: use the tree path in `INV_MAT_MS`.
  - `MSSPEC_TREE_SERIES_ORDER` = `linear_solver_max_iterations`.
  - `MSSPEC_TREE_TOL` = `linear_solver_tolerance`.
  - `MSSPEC_TREE_GMRES = 1`: use GMRES, not the correction series.
  - `MSSPEC_TREE_GMRES_RESTART` = `linear_solver_restart`.
- `parameters.py` lines 1472 to 1480 define the user knobs and their defaults: tolerance, maximum iterations (default 12) and restart (default 6).
- Any other algorithm runs the spectroscopy without these variables. The dense LU path is then used as before.

The wrapper passes the knobs through the environment. Fortran reads them with `GET_ENVIRONMENT_VARIABLE`. This keeps the Fortran interface unchanged.

### 4.2 `INV_MAT_MS`: entry point and failure stop

- Reads the tree switches. If `MSSPEC_TEST_TREE_SERIES = 1`, it calls `MS_TREE_SERIES_BATCH` with the right-hand sides saved in `IN`.
- The tree order is read from `MSSPEC_TREE_SERIES_ORDER`. The default is 2. Values outside 1 to 30 stop the run with `STOP 98`.
- A copy of the right-hand sides is kept in `BSAVE`, because the tree solver overwrites `IN`.
- If the tree solver fails (`TREE_INFO ~= 0`), it prints `MSSPEC_SOLVER_STATUS je,tree_info,status= ... none_available` and stops with `STOP 99`. There is no dense fallback (removed in v2; the study counts this as a method failure).
- If the tree path is off (`algorithm='inversion'`), the dense LU path (`ZGETRF`) runs as before. That is the reference.

### 4.3 `MS_TREE_SERIES_BATCH`: tree, blocks and start point

This routine has four jobs.

1. **Read switches.** It reads `MSSPEC_TREE_GMRES` (only the value `1` turns GMRES on), `MSSPEC_TREE_GMRES_RESTART` (clamped to 2 to 12), `MSSPEC_TREE_TOL` (clamped to 1e-12 to 1e-1, default 1e-2), and the diagnostic and check switches.
2. **Build the tree.** Prim's algorithm chooses a minimum-distance spanning tree over the atoms. It starts from atom 1. Each step adds the closest unvisited atom to the tree. The routine stores each atom's `PARENT`, its visit order in `VISIT`, its block size `Q` and its offset `START`.
3. **Factorise the tree.** For each atom, the diagonal block of `A` is factorised with `ZGETRF`. For each tree edge, the pair blocks are computed by `MS_TREE_PAIR_BLOCK`:
   - `EP(:,:,child) = A(parent, child)`,
   - `EC(:,:,child) = A(child, parent)`.
   Leaves are eliminated towards the root. Each diagonal Schur block is factorised. The root block is factorised last.
4. **Start and dispatch.** It computes `x0 = M^-1 b` with `MS_TREE_SOLVE`. If `MSSPEC_TREE_GMRES = 1`, it calls `MS_TREE_GMRES`. Otherwise it runs the correction series: `x <- M^-1 (b - (A - M) x)`, with `ORDER` terms. The series is used only for the order study.

`M` is the block tree matrix. It keeps the atom diagonal blocks and the exact pair blocks for tree edges only. It drops all other pairs. `M` is the preconditioner. The full `A` is never stored on the GMRES path.

### 4.4 `MS_TREE_GMRES`: the restarted solver

This routine solves `M^-1 A x = M^-1 b` for all right-hand sides together.

Setup:

- `MMAX = MAX(2, MIN(12, RESTART))` is the Arnoldi steps per cycle.
- `MAXIT` is the total number of Arnoldi steps allowed per right-hand side. The wrapper sets it to `linear_solver_max_iterations` (default 12). With the default restart of 6, that is two cycles at most.
- `BNORMS` stores the norm of each `b`. `ACTIVE`, `INNER`, `ITERS` and `BROKE` track the state of each right-hand side.

Each outer cycle:

1. Compute the preconditioned residual `BPRE = M^-1 (b - A x)` and its norm `BETA`. If `BETA <= TOL * ||b||`, that right-hand side is done.
2. Inner Arnoldi loop, for `J = 1 .. MMAX`:
   - `W = A V(:,J)` with `MS_MATRIX_ACTION_BLOCK`, then `W = M^-1 W` with `MS_TREE_SOLVE`.
   - Modified Gram–Schmidt against the previous basis vectors gives the new column of `H`.
   - Complex Givens rotations reduce the new column. They keep the rotated right-hand side `G`.
   - The estimate of the preconditioned residual norm is `|G(J+1)|`.
3. **Early stop.** If `|G(J+1)| <= TOL * ||b||`, that right-hand side leaves the inner loop at once. This is the in-cycle stop. It does not form the solution yet.
4. **Cycle end.** For each right-hand side, the code solves the small triangular system `H y = G`, then updates `X <- X + V y`.
5. **True residual.** It computes `AX = M^-1 (A X)`, then `R = BPRE - AX`, and its norm `RNORMS`. Only this norm decides convergence. If `RNORMS <= TOL * ||b||`, the right-hand side is done. Otherwise the next cycle starts.

Failure handling:

- A non-finite norm, a zero diagonal in the triangular solve, or a failed `MS_TREE_SOLVE` gives an error code. The codes are negative, one offset per right-hand side: -150 to -200 for the GMRES steps, and -40 to -60 for the tree solves.
- If `ITERS` reaches `MAXIT` before convergence, the routine returns an error. `INV_MAT_MS` then stops with `none_available` (section 4.2).
- If the Arnoldi vector is numerically zero (breakdown), the cycle ends for that right-hand side. The true residual step then decides.

The solver prints `TREE_GMRES iterations,rhs=` for each right-hand side. These lines show the work per solve.

### 4.5 `MS_TREE_SOLVE`: applying `M^-1`

- Input: a block vector `B` with `NRHS` columns.
- Sweep 1 (leaves to root): for `K = NA .. 2`, it takes the child in `VISIT(K)`, solves with the child's factored Schur block, and subtracts the child's contribution from the parent's right-hand side through `EC`.
- Sweep 2 (root to leaves): it solves the root block. Then, for each child, it updates the child's right-hand side from its parent's solution through `EP` and solves the child's block.
- Output: `B` holds `M^-1 B`.

Each call costs one pass over the tree. It does not touch the full `A`.

### 4.6 `MS_TREE_PAIR_BLOCK`: exact pair blocks

- Computes one directed pair block `A(parent, child)` or `A(child, parent)` exactly.
- It uses the same expressions as the production operator: Gaunt coefficients, spherical Hankel functions, spherical harmonics, the attenuation factor `exp(-Im(k) r)`, and the t-matrix.
- The pairs are the tree edges only, so the cost is proportional to the number of tree edges, not to `N^2`.

### 4.7 `MS_MATRIX_ACTION_BLOCK`: matrix-free `A`

- Computes `Y = (I - G T) X` for `NRHS` vectors at once. The matrix is never stored.
- The loop runs over ordered atom pairs. The geometry of a pair (distance, angles) is evaluated once and reused for all right-hand sides.
- Threads: `MSSPEC_TREE_ACTION_THREADS` sets the number of OpenMP threads. The value is capped at the maximum the OpenMP runtime gives.
- Geometry cache: `MSSPEC_TREE_GEOMETRY_CACHE=1` stores the pair geometry for later calls. It is used only when the cache size stays below 512 MiB (`CACHELIMIT`). The cache is valid for one energy and one cluster size.
- Gaunt table: the real Gaunt table is kept between calls. It is rebuilt when the largest angular momentum changes.
- Axial path (default): for each pair, the action is evaluated in the rotated frame by `MS_AXIAL_PAIR_ACTION`. The direct path is used when `MSSPEC_TREE_AXIAL_TRANSLATION=0`.
- Profile: `MSSPEC_TREE_ACTION_PROFILE=1` prints timings of the action.

### 4.8 `MS_AXIAL_PAIR_ACTION`: rotate, translate, rotate back

- Computes the action of one directed pair in three steps:
  1. Rotate the source vector so that the pair axis is along `z`. The rotation uses Wigner-d matrices. The code calls `DJMN2` for these.
  2. Apply the axial translation. For a `z`-directed translation, only equal `m` values couple. This makes the step cheap.
  3. Rotate the result back to the lab frame.
- The default variant is set by `MSSPEC_TREE_AXIAL_VARIANT` (default 1, clamped to 0 to 3). The validate switch `MSSPEC_TREE_AXIAL_VALIDATE=1` compares the axial result with the direct result. Use it only for checks.
- The file header calls this routine a prototype.

### 4.9 `MS_TREE_ITER_DUMP`: iterate dump

- Runs only when `MSSPEC_TREE_ITER_DUMP=1`.
- Writes `msspec_tree_iter.dat`. Iteration 0 is the preconditioner alone `M^-1 b`. Then one block per Arnoldi step, with the iterate, the absolute preconditioned residual and the step index (D20, D21).
- At the end of each right-hand side, a `FINAL` block is written. Its content equals the returned solution.
- The test `code/tests/test_tree_gmres.py` checks that the dump is bit-identical between runs and that the last iterate equals the result.

## 5. Switches and defaults

Values are read from the source (`axb_mat_ms_la.f`). The wrapper sets the first four for the GMRES path.

| Variable | Default in Fortran | Set by wrapper | Meaning |
|---|---|---|---|
| `MSSPEC_TEST_TREE_SERIES` | off | `1` | Use the tree path in `INV_MAT_MS`. |
| `MSSPEC_TREE_SERIES_ORDER` | `2` | `linear_solver_max_iterations` (12) | Series order, or GMRES maximum steps per right-hand side. Range 1 to 30. |
| `MSSPEC_TREE_GMRES` | off | `1` | `1` selects GMRES. Otherwise the correction series runs. |
| `MSSPEC_TREE_GMRES_RESTART` | `6` | `linear_solver_restart` (6) | Arnoldi steps per cycle. Clamped to 2 to 12. |
| `MSSPEC_TREE_TOL` | `1e-2` | `linear_solver_tolerance` | Relative tolerance. Clamped to 1e-12 to 1e-1. |
| `MSSPEC_TREE_DIAGNOSTICS` | off | no | Extra norms of retained and omitted blocks. Costs one more pass. |
| `MSSPEC_TREE_VERIFY` | off | no | One extra full matrix action as a check. |
| `MSSPEC_TREE_ITER_DUMP` | off | no | Writes `msspec_tree_iter.dat`. |
| `MSSPEC_TREE_STOP_AFTER` | unset (-1) | no | Study only (E16, D22). `k=0` returns `M^-1 b`; `k=1..restart` returns the iterate after exactly k Arnoldi steps, with no residual check. Other values stop with `none_available`. |
| `MSSPEC_TREE_ACTION_THREADS` | `1` | no | OpenMP threads for the matrix action. |
| `MSSPEC_TREE_ACTION_PROFILE` | off | no | Prints action timings. |
| `MSSPEC_TREE_GEOMETRY_CACHE` | off | no | Caches pair geometry (used below 512 MiB). |
| `MSSPEC_TREE_AXIAL_TRANSLATION` | on | no | `0` forces the direct Gaunt action. |
| `MSSPEC_TREE_AXIAL_VALIDATE` | off | no | Compares the axial and direct results. |
| `MSSPEC_TREE_AXIAL_VARIANT` | `1` | no | Internal variant of the axial action, 0 to 3. |
| `MSSPEC_TREE_AXIAL_FULL_LOOP_REFERENCE` | off | no | Uses the full loop bounds (a control for checks). |

Two other switches, `MSSPEC_CHECK_MV` and `MSSPEC_PACK_RHS`, are read in the same file. They are not part of the GMRES path and are not described here.

## 6. Status lines and exit behaviour

- `MSSPEC_SOLVER_STATUS je,tree_info,status= ... none_available`: the tree solve failed (always the case on failure now). The run stops with exit code 99.
- `MSSPEC_TIMING`: prints timings of the tree build and the solve.
- `TREE_GMRES iterations,rhs=`: the number of Arnoldi steps for each right-hand side.
- `TREE_GMRES stop_after=k`: printed when the stop switch is set.

## 7. How the pieces fit

- Accuracy depends on the tree. The tree keeps the strongest links, so `M` is close to `A`. The GMRES method then needs only a few steps.
- Memory depends on the tree and the pair blocks. The dense matrix is never built on the GMRES path. The only stored blocks are the atom diagonal factors and the tree-edge pair blocks.
- Time depends on the number of GMRES steps. Each step costs one matrix action and one tree solve.
- Safety comes from three checks. The in-cycle early stop only skips work. The true residual decides convergence. A failure is never hidden: the run stops and is counted as a method failure.

## 8. What was not verified here

- These descriptions come from reading the source. They were not checked by a new run in this step.
- The tests last passed on the current build (see [GRMES_STATUS.md](GRMES_STATUS.md), section 2.2). Run `code/tests/test_tree_gmres.py` again after any change to this file.
