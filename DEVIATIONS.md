# Deviations and decisions

Each entry gives the date, what changed from the plan or the protocol, and why.
Nothing here changes a result. Entries are added, never removed.

## D1. In-cycle early stop changes the iteration path (decided 2026-10-19)

- Plan: the new tree-GMRES should match the original tree output to 1e-12.
- What happens: the in-cycle stop (C1) ends a cycle once the estimate meets tau.
  The new build uses 10 to 11 iterations where the original uses 12, so the
  two tree results differ by up to 2.4e-5 (relative).
- Accepted gate: each tree result within 1e-4 relative L2 of dense LU. Measured:
  original 1.4e-6 to 2.6e-6, new 6.1e-6 to 8.5e-6 (Ag, Ni, Cu).
- Decision: confirmed by the author on 2026-10-19 in chat ("fix all"). The
  strict 1e-12 comparison is reported, not gated.

## D2. Tag `v1-gmres` created (2026-10-19)

- Annotated local tag on commit `e56f1b8`, the frozen state (code and manifest).
  Created on the author's explicit instruction after an earlier attempt was denied
  by the permission classifier. Nothing is pushed.

## D3. Benchmark inputs: located and partly copied (2026-10-19)

- The two CSVs needed by E2b and E9 exist, not under `paper_v2/data/` at the
  top of `msspec_algo/` but under `msspec_algo/order_sweep/paper_v2/data/`
  (identical copy in `msspec_ped_package/`). Copied unchanged to
  `data/inputs/paper_v2/data/`:
  - `random_emitters.csv` sha256 b313cad5ad480ec4c0e2de84368efd2584ad97c3ac505c3dedadf7b2c81d506c
  - `ml21_selection.csv` sha256 e26eb44467d91b1b28a0636486d54ccfc5b67e9df154370532c8056fb5a988db
- `order_sweep/solver_order_bench/run_order_bench.py` is NOT copied. It names
  another iterative method and compares it with GMRES, which this project must not do. The
  runners for this study are written new in `runners/`, following its
  structure only where useful.
- E2a input `Materials/stereographic_comparison/protocol.json` is used as in D4.
- Status: E2b and E9 inputs are in place. The RND part of E3 uses `random_emitters.csv`
  and is unblocked. Runners not written yet.

## D4. E2a input file names another iterative method (2026-10-19)

- `Materials/stereographic_comparison/protocol.json` sets `"solver": "tree-preconditioned another iterative method"`.
- This project uses only GMRES and must not name another iterative method. The file
  is not copied into `article_gmres/`. The E2a settings in `PROTOCOL.md` section 2
  and 5 are taken from this document, not from that file's solver field.
- The file's other settings (100 eV, 5,269 directions, 8 threads, dense LU for the
  reference) are used as stated in `PROTOCOL.md`.

## D5. Pilot not run (2026-10-19)

- The pilot is part of Phase 1 execution. It has not run, because D3 blocks part
  of the benchmark set and the pilot's cases come from MAT-B and E8.
- The pilot may set only timeouts, the run budget and the memory preflight
  constant. Solver settings in `PROTOCOL.md` section 2 are already fixed.

## D6. Project moved to `msspec_GRMES` (2026-10-22)

- On the author's instruction the project was copied from
  `msspec_algo/article_gmres/` to `/home/mohammed/Desktop/REAL-WORK/msspec_GRMES`
  (git history and tag `v1-gmres` kept). All work continues only there.
- The old folder `msspec_algo/article_gmres/` is left untouched; it is frozen at
  commit `9989acd` and is not used again.
- Only the venv scripts held the old path; they were fixed in the new copy.
  Checked in the new location: imports, `check_release.sh` PASS, C4 test PASS
  (same numbers as before the move).
- Inputs read from the sibling project (`Materials/stereographic_comparison/protocol.json`,
  `Materials/` cluster files) stay read-only references under `msspec_algo/`.

## D7. Amendment A1: sets E13 and E14 added (2026-10-23)

- Added to `PROTOCOL.md` before any reported run, after the pilot size series.
- Reason: the pilot shows that LU memory limits the size series near 150 Cu atoms
  at 300 eV. E14 measures GMRES beyond that limit; E13 answers the expected
  question on thread dependence of the speed-up.
- The buried emitter for E2b is the per-material entry of
  `data/inputs/paper_v2/data/random_emitters.csv` (emitter index, depth, lateral
  offset), as PLAN.md section 5 already states. No choice is left open.

## D8. Project moved to `MsSpec_GRMES` (2026-10-23)

- On the author's instruction the project was copied from
  `/home/mohammed/Desktop/REAL-WORK/msspec_GRMES` to
  `/home/mohammed/Desktop/REAL-WORK/MsSpec_GRMES` (git history and tag
  `v1-gmres` kept). All work continues only there.
- The duplicate `msspec_GRMES/` was then deleted on the author's request
  (2026-10-23), after checking that it had no uncommitted work and that its
  last commit `56f1207` is in this repository. `msspec_algo/article_gmres/`
  is left untouched and is not used again.
- Venv scripts fixed for the new path. Checked in the new location:
  `check_release.sh` PASS, C4 test PASS with the same numbers as before.

## D9. Series orders 1 to 8 on the material platform (2026-10-24)

- `PROTOCOL.md` section 2 and `GRMES_STATUS.md` section 4 allow series orders 1 to 3 only.
  The material platform (`platform/`) runs orders 1 to 8, on the author's request.
- Orders 4 to 8 are exploratory. They are reported as such and are not part of
  the pre-registered E11 result.
- Reason: the author wants to see how the series behaves at higher order on
  materials with large spectral radius.

## D10. Definition and threshold of rho for the material platform (2026-10-24)

- rho is the spectral radius rho(G0T) of the multiple-scattering operator, computed
  per material with MsSpec `EIG` and `algorithm='power'` (`platform/rho_probe.py`).
  The series converges only if rho < 1 (`msspec_ML/presentation/main.tex`).
- The high group was first set at rho >= 3 (the author's "rho >> 1"). The author
  then approved lowering it to rho >= 1.5, and then to rho >= 1.2 (author's
  choice, 2026-10-24). The threshold is recorded with each row.
- Reason: the scan (`platform/rho_scan.csv`, 35 rows, 33 finished) found no
  material with rho >= 3 and only 3 with rho >= 1.5 among clusters of at most 150
  atoms. The author stopped the scan, kept the 150-atom limit, and then resumed
  the scan (same order, same rules) to reach 50 materials at rho >= 1.2.
- Absorber for each run is read from `config.json` (`absorber_idx`), not from the
  flag column of `cluster.xyz`, which is not in one layout across runs. The
  absorber position is checked against `absorber_xyz`.
- Failed scan rows stay in `rho_scan.csv` with rc not 0. A run is retried only
  when it has no rc 0 row.

## D11. High-rho group taken from the bo_rho set (2026-10-24)

- Author's choice (2026-10-24): the high-rho group is the 50 highest-rho cases of
  `msspec_SE/msspec_standard/material_platform/high_rho/` (97 cases, indexed in
  its `index.csv`; read only), not cases from the prod3000 scan. The list is
  `platform/high_group_top50.csv`. Its rho range is 2.11 to 4.57. The low-rho
  group (rho < 1) still comes from prod3000.
- The first list (before this change) was the 50 divergent cases in
  `bo_divergent_manifest.csv`. It is replaced by the top 50 of `high_rho/`. The
  two lists differ: the manifest reaches down to rho 2.01, while the top 50 starts
  at rho 2.11, so the manifest's lower-rho cases are not in the new list.
- Differences from the prod3000 legs, recorded here:
  - Clusters are hemispherical (`hemispherical_cluster`, `high_rho/<tag>/run.py`),
    not prod3000 clusters. They have no surface orientation, so orientation
    coverage in the high group comes from crystal structure and cluster size only.
  - Settings: lmaxt 19 (prod3000 legs use lmaxt 8) and Debye model off
    (prod3000 legs use it on).
  - Materials: 21 distinct materials in the 50 cases.
  - rho is taken from `bo_rho/bo_results.csv`. The per-run `results.dat` is not
    present in `bo_rho/runs/`, so these values are not checked run by run. A
    spot check agrees: SiC 10.6 eV (2.04316) and W 6.0 eV 5 planes (2.73917,
    row 91 of `bo_results.csv`).
  - The run folders hold `run.py` and `log.txt` only (no `cluster.xyz`,
    `config.json`), so the platform rebuilds each cluster from `run.py`.
- Atom cap, decided by the author (2026-10-24): the 150-atom cap does not apply to
  the high group. The high group has 50 cases; 18 have at most 150 atoms and 32 have
  more (up to 298 atoms). All 50 are kept. A case whose dense LU matrix does not fit
  the memory preflight is recorded as a preflight failure, not dropped. The 150-atom
  cap still applies to the low group (prod3000).

## D12. Material selection for the platform (2026-10-24)

- `platform/select.py` writes `platform/materials_100.csv`: the 50 high-rho cases
  (`platform/high_group_top50.csv`, D11) and up to 50 low-rho cases from
  `platform/rho_scan.csv` (rc 0, rho < 1, at most 150 atoms, one run per distinct
  (name, prototype), lowest ke first).
- Shortfalls at the time of writing:
  - Low group: 46 of 50 (the scan stopped at 91 of 589 and was resumed; it is
    running again in the background, see `platform/work/rho_scan_ke3.log`).
  - High group: 50 cases but only 21 distinct materials. The original goal of 100
    different materials is therefore not met. This is recorded, not changed.
- `select.py` is rerun after the scan finishes; the table is rewritten each time.

## D13. Build v2: dense fallback code deleted (2026-10-30)

- Author's decision (2026-10-30): the dense LU fallback is removed from the study.
  Source change in `axb_mat_ms_la.f`: the fallback branch, the dense-fallback
  switches and `MS_DENSE_ASSEMBLE` are deleted. A failed tree solve prints
  `MSSPEC_SOLVER_STATUS ... status=,none_available` and stops (`STOP 99`). The
  status label was renamed from `fallback=` to `status=`.
- `git diff --stat v1-gmres` for that file: 6 insertions, 93 deletions.
- Tests (`code/tests/test_tree_gmres.py`): `forced_fail` now expects a method
  failure (non-zero exit, `none_available`, no `dense_lu`). `forced_nodense` is
  removed, because its settings are the same as `forced_fail` and no fallback
  exists to disable.
- Checks run on 2026-10-30:
  - `rebuild.sh` log (`platform/work/rebuild_v2.log`): no error lines. The `.so`
    is newer than the edited source.
  - Built `_phd_mi_noso_nosp_nosym.so`: strings `ms_dense_assemble`,
    `MSSPEC_FALLBACK_DENSE` and `dense_lu` count 0; `MSSPEC_SOLVER_STATUS` count 1.
  - `test_tree_gmres.py --case` lu, tree, tree_dump, forced_fail: all PASS.
    `tree` vs LU rel L2 = 4.18e-3 (tolerance 1e-2, the pre-registered settings).
    This is above the 1e-4 accuracy gate of D1, which was set for the other
    settings (tolerance 1e-5, restart 6, cap 12 in `material_check.py`).
  - `check_release.sh`: PASS.
  - `write_build_manifest.py` rerun; `source_tree_sha256` is now
    `1e7eab4cb23304be472cc0120f86651611a55b76fa37476ce3abdd20e5a41916`.
- Tag: `v2-gmres` is not created. It needs the author's yes. `v1-gmres` (commit
  `e56f1b8`) stays as the frozen record.
- Update (2026-10-30, same day, on the author's yes): the v2 build is committed as
  `3fe9a7c` ("GMRES v2: dense fallback code removed") and the local annotated tag
  `v2-gmres` is on that commit. Nothing is pushed. Checked with `git tag -n1` and
  `git rev-parse v2-gmres^{commit}`. The `v1-gmres` tag object is `2a28855`, and it
  peels to `e56f1b8`, which matches D2.
- Not checked: `material_check.py` not rerun (its `orig_tree` branch is unresolved,
  Step 0). The full `test_tree_gmres.py` without `--case` not run. The per-case
  logs do not hold the solver status line (the worker output is not printed), so
  `none_available` was confirmed only by the test's own assertion.
- Open: `code/tests/test_tree_gmres.py` and `code/tests/material_check.py` call
  `tempfile.TemporaryDirectory()` without `dir=`. The runs here used `TMPDIR`,
  which points into `platform/work/`. Both need `dir=` in code (project rule).

## D14 (2026-10-31): amendment A2 written before any reported run

- PROTOCOL.md amendment A2 (A2.1 to A2.12) replaces the dense fallback, the set
  list for the main evidence, the materials, the metrics and the memory rule.
- Reason: the v2 build removed the fallback (D13); the rho scan and the pilot gave
  new facts (D15 to D17). Changes are made before any reported run, so no result
  is affected. The original PROTOCOL text stays and is marked as superseded.

## D15 (2026-10-31): the pilot did run; correction of D5

- D5 said the pilot was not run. That was wrong. The pilot ran on 2026-10-23:
  `pilot/small_compare.csv` (4 rows: radius 1.0, 1.5, 2.0, 2.5 a0 with 9, 34, 70 and
  143 atoms) and `pilot/small_compare.log`. Not reported.
- The radius 3.0 runs FAILED (a spin-polarisation message on stderr). The cause is
  not diagnosed. E14 is blocked until it is (A2.12).
- Memory from `pilot/small_compare.csv` (column `lu_rss_mib`): 13,643 MiB at 143
  atoms. Fitted 0.667 MiB per atom^2 from that one point; a predicted 8 GiB limit
  near 111 atoms (computed 2026-10-31). The 70-atom point gives 0.699, so the model
  is approximate. Extrapolations to 150 atoms are estimates, not measurements.

## D16 (2026-10-31): D11 superseded (high-group materials)

- D11 used `platform/high_group_top50.csv` (from `high_rho/`, lmaxt 19, hemispherical
  clusters, Debye off). Recounted 2026-10-31: 32 of the 50 cases are above 150 atoms
  (maximum 298). For those no dense LU reference can be made on this machine (A2.9).
- The high group and the low group also differ in four settings, so a group
  difference could not be interpreted. D11 is superseded by A2.4 and A2.5. The text
  of D11 stays as written.

## D17 (2026-10-31): D12 counts are stale

- Recounted 2026-10-31 from `platform/rho_scan.csv` (602 rows; columns include rc,
  rho, natoms, name, prototype):
  - rc 0: 600 rows; rc 1: 2 rows (Co 80 atoms, Bi 91 atoms; cause not established).
  - rho at least 1.0: 200 rows, 138 distinct materials (name, prototype), 14 to 150 atoms.
  - rho at least 1.2: 99 rows, 79 distinct, 23 to 150 atoms.
  - rho at least 1.5: 25 rows, 21 distinct, 23 to 148 atoms.
  - rho below 1: 400 rows, 258 distinct.
  - Maximum rho: 3.21.
- These numbers replace the scan counts in D10 and D12. The scan is complete for
  the purpose of selection; no scan process is running.

## D18 (2026-10-31): test files and results after the temp-folder fix

- `code/tests/test_tree_gmres.py`: `tempfile.TemporaryDirectory` now uses
  `dir=platform/work/` in `main()` and `run_one()` (project rule, never `/tmp`).
- `code/tests/material_check.py` already had `dir=WORK` (line 108). The earlier
  note saying it did not was wrong; `material_check.py` is v1-only and is not rerun
  (author instruction, 2026-10-30).
- Results, from `platform/work/test_v2b_*.log` (read 2026-10-31):
  - `--case lu`, `tree`, `tree_dump`, `forced_fail`: all PASS.
  - `tree` and `tree_dump`: rel L2 vs LU = 0.00418127.
  - Full run: tree vs LU 0.004181274, dump bit-identical True, dump rhs 6,
    last iterate vs returned 0.0, forced_fail exit 99, PASS.
- Not checked: `material_check.py` (v1, not rerun); the status line of each
  per-case run (worker output is not printed); LU bit-reproducibility (step 4.2).

## D19 (2026-10-31): corrections of D15 and of amendment A2.9

- D15 said `pilot/small_compare.csv` has 4 rows (radius 1.0 to 2.5). It has 3 data
  rows: radius 1.5, 2.0 and 2.5 a0 (34, 70, 143 atoms). Checked with Python `csv`
  this turn. The 9-atom point (radius 1.0) is in `pilot/small_compare_r1.0.txt`, a
  separate run of the same script, with LU RSS 225 and 226 MiB and GMRES RSS 199
  MiB. It is not in the CSV.
- The pilot CSV values are the mean of two repeats (for example LU RSS 976.6 MiB at
  34 atoms, from 976 and 977 MiB in the log). Per-atom^2 values: 0.845 (34 atoms),
  0.699 (70), 0.667 (143) MiB.
- The radius 3.0 explanation has two sources that disagree. `pilot/small_compare.log`
  shows FAILED for LU and for GMRES, both with the spin-polarisation message ("Enable
  spin polarization in the constructor of your Calculator"). `small_compare_r1.0.txt`
  says the radius 3.0 run "was stopped on purpose" because LU would need about 35
  GiB. The stop reason fits LU only; GMRES also failed with the spin message, so the
  cause is not established. Both statements are recorded. The 35 GiB figure is an
  extrapolation and is not supported by any recorded atom count for radius 3.0.
- A2.9 said cases of "about 110 to 150 atoms" use the +3 GiB exception. That is
  wrong. PLAN.md section 8 makes the +3 GiB rule a free-RAM condition (prediction +
  3 GiB below free RAM, swap off), not an atom range. A2.9 is corrected in
  `PROTOCOL.md` with the rule as written in PLAN.md. The fitted 8 GiB limit is at
  about 111 atoms and 30 GiB at about 215 atoms (both from the 0.667 fit, computed
  2026-10-31).
- The design file's statement that LU fits in 30 GiB "only up to about 150 atoms"
  is withdrawn. It was not consistent with the fit (about 215 atoms) and it was not
  computed.

## D20 (2026-10-31): iteration-0 dump added to the GMRES build (author approved)

- Change: `MS_TREE_GMRES` in `axb_mat_ms_la.f` writes an `ITERATE` block with
  iteration 0 when `MSSPEC_TREE_ITER_DUMP=1`. The block holds the preconditioner
  alone, `M**(-1) b`, before the first Arnoldi step (needed for E16). Its residual is
  the preconditioned residual `|BPRE - M**(-1) A BPRE|` (first version divided by
  `|BPRE|`, a unit error corrected in D21), computed with one
  extra action and tree solve that run only under the dump. `AX` is scratch there and
  is overwritten in the first cycle, so the solver path is unchanged.
- No change to the solver numerics or to the default (dump off) path.
- Source hash (`source_tree_sha256` in `code/BUILD_MANIFEST.json`): `3eef5dec906c98a5…`
  (was `1e7eab4cb23304be…` at HEAD, the v2 hash). Written by
  `tests/write_build_manifest.py` this turn.
- Checks this turn (outputs in `platform/work/iter0_*.log`): `rebuild.sh` rc 0, the
  `.so` is newer than the source; cases `lu`, `tree`, `tree_dump`, `forced_fail` each
  PASS; full run PASS (`tree_vs_lu` 0.004181, `dump_bit_identical` True,
  `last_iterate_vs_returned` 6.97e-18, `forced_fail_exit` 99); the full run also
  checks that the first block of each RHS is iteration 0; `check_release.sh` PASS.
- Naming: the build is not tagged. Whether it becomes `v3-gmres` is the author's call.
- Follow-up done the same day: `PROTOCOL.md` (A2.1, A2.11) and `GRMES_STATUS.md`
  (§2.1) now cite `3eef5dec...`. The tag `v2-gmres` still points at the old commit
  `3fe9a7c`; the new source is not committed or tagged. `manuscript/ARTICLE_DESIGN.md`
  has no hash reference. (The hash in this entry and in the three documents was
  replaced by D21.)

## D21 (2026-10-31): unit error in the iteration-0 residual, new hash

- Found in a review run: in D20 the iteration-0 residual was divided by `|BPRE|`
  (1.062 in the test case), while the later ITERATE blocks hold the absolute
  `|G(J+1)|`. The sequence then looked like 0.240, 0.250 (a growth), which GMRES with
  x0 = 0 cannot show. Fix: iteration 0 is now absolute. Test case, RHS 1:
  0.255, 0.250, 0.0608, 0.0152, 0.00401 (read from the dump, not typed).
- `test_tree_gmres.py` now asserts that the residual sequence of each RHS does not
  grow (tolerance 1e-6 relative), so this error cannot return unseen.
- New source hash `6926ee97...` (replaces `3eef5dec...` in D20, PROTOCOL A2.1, A2.11,
  A3 header, and GRMES_STATUS §2.1). Checks: rebuild rc 0, `.so` newer than the source,
  full test PASS (`tree_vs_lu` 0.004181, `dump_bit_identical` True), `check_release.sh`
  PASS.
- Found in the same review: the dump writes vectors, not observables. Theta scans of
  an iterate need a switch that stops GMRES after k steps. Done in D22.

## D22 (2026-10-31): stop switch for E16 (author approved), new hash

- Change: `MSSPEC_TREE_STOP_AFTER=k` in `MS_TREE_GMRES`. k = 0 returns the
  preconditioner alone, `M**(-1) b`. k = 1 to the restart length returns the iterate
  after exactly k Arnoldi steps: no early stop, no residual check, no convergence
  error. The normal pipeline then computes the pattern from that iterate. Unset or -1:
  the normal solve (unchanged). A value that is not an integer, or above the restart
  length, gives status -210, `STOP 99`, `none_available` (never a silent normal solve).
  The solver prints `TREE_GMRES stop_after=k`; the iterations printed must equal k.
  A breakdown before step k would return fewer steps, so the runner must compare.
- Study use only (E16, A3.2). Reported accuracy of the method never uses it.
- New source hash `b64d4c64fd404d4c...` (replaces `6926ee97...`).
- Checks (output in `platform/work/review/full3.log`): rebuild rc 0, `.so` newer than
  the source; the full test PASS with the default path unchanged (`tree_vs_lu`
  0.004181, `dump_bit_identical` True); `stop0` and `stop2` returned vectors equal to
  the dumped iterates 0 and 2 (max relative difference 0.0); observable error against
  LU on the test case: k = 0: 0.4444, k = 2: 0.07235, converged: 0.004181 (relative
  L2); values 99 and `abc` are refused; `check_release.sh` PASS.
- Not checked: k = 1, 3, 4 (the code path is the same loop; E16 small test covers
  them); behaviour with several energies and emitters in one run.


## D23 (2026-10-31): radius 3.0 diagnosis (E14 blocker)

- Test: one GMRES job, same pilot worker (`pilot/small_compare.py`), Cu hemisphere at
  radius 3.0 a0, 300 eV theta scan, current build (hash `b64d4c64...`), full stdout and
  stderr kept in `platform/work/review/r30/` (script `platform/work/review/r30_run.py`).
- Result: return code 0, 174.9 s wall, 259 atoms, peak RSS 919 MiB (+150 MiB children),
  `gmres.npy` written. `TREE_GMRES iterations,rhs=` shows 12 for rhs 1 to 6 (checked
  lines 1 to 6 only). stderr holds only the known warnings; the "spin polarization"
  text is the harmless warning for the ignored `integrals` parameter
  (`parameters.py` lines 1663 and 1678), not an error.
- So the GMRES solver is not blocked at 259 atoms. The 2026-10-23 pilot failure of
  GMRES at 3.0 is NOT reproduced; its true cause is not established. The pilot kept only
  the last 300 characters of stderr, so the real message was lost.
- The LU failure at 3.0 is explained by memory: 0.667 MiB per atom^2 times 259^2 is
  about 44.7 GiB, above the 30 GiB of this machine. This corrects the "about 35 GiB"
  figure in `small_compare_r1.0.txt` (D19).
- Open hypothesis (not tested): the 2026-10-23 pilot used the v1 build, which could
  fall back to dense LU after a failed tree solve, and that fallback would need 44.7 GiB.
  Testing it needs the v1 build; not done.
- Effect: E14 is unblocked as far as the solver is concerned. The accuracy reference at
  this size is converged GMRES (A1), as already pre-registered. The pilot must not be
  cited for the radius 3.0 outcome.
- Not checked: whether the run converged at every rhs (the status line was not read for
  all 6 beyond iterations=12, equal to the cap); the answer is the return code 0 with no
  fallback in the build; load average was 4 (another job ran), so the time is not a
  benchmark.

## D24 (2026-10-31): P100 selection, runner, evaluator; series orders 2 and 3 fail with the Debye model

- New files, no solver source changed (source hash `b64d4c64` unchanged): `platform/select2.py`
  and `platform/materials_v2.csv` (50 high, 50 low, seed 20261030), `platform/check_materials.py`
  (prints PASS), `code/runners/run_set.py` and `worker.py`, `analysis/evaluate.py` with
  `analysis/test_evaluate.py` (PASS, 7 groups). Only smoke runs were made
  (`platform/work/runs/`); none is a reported result.
- R_chi definition (author has not yet approved it): chi = (I - I0) / I0, R_chi = sum (chi - chi_ref)^2 / sum (chi^2 + chi_ref^2), points with I0 <= 0 left out. It is a supporting metric, not a pass criterion.
- Not implemented: the untimed warm-up of A2.8. Times from the runner are not benchmark times until it is added.
- Consequence of rule A2.9: 16 of the 100 cases (15 high, 1 low) predict LU >= 8 GiB. They are
  `lu_not_run` unless swap is off and free RAM allows one attempt. For them GMRES has no reference.
- Finding: with the Debye model on (as in prod3000), series orders 2 and 3 return all-NaN
  intensities (rc 0, L43, 12 atoms, rho 0.49) in this build. Order 1 is finite. With the Debye
  model off, order 2 is finite. Raising the stack limit did not help. `MSSPEC_DEBYE_QD` is not
  read anywhere in `code/msspec_gmres/src`; the prod3000 note says orders >= 2 need a fixed
  correlated-displacement routine that this copy does not have. Cause in the Fortran not
  established. The runner now records such a run as `method_failure`, reason `nonfinite_result`.
  Open decision for the author: (a) report orders 2 and 3 as failing under Debye, (b) run the series
  with the Debye model off for all methods (new settings set, needs a protocol amendment), or
  (c) port the fixed routine (changes the source hash).
- Not checked: series orders 4 to 8; whether the NaN also occurs on high-rho cases; GMRES with
  Debye on at several cases (only the smoke cases were run).

## D25 (2026-10-31): series expansion runs on the FINAL dev tree

- Reason (D24): our `code/msspec_gmres/src` cannot run series orders 2 and 3 with the Debye
  model on (all-NaN). The dev tree reads `MSSPEC_DEBYE_QD` in `phd_se_noso_nosp_nosym/do_main.f`
  (the cluster's dummy cell gives an infinite Debye wave vector otherwise).
- Decision (author, 2026-10-31): the series methods (`series1` to `series8`) use
  `/home/mohammed/Desktop/REAL-WORK/msspec_ML/FINAL/07_msspec_dev/msspec_python3/src`, read-only.
  LU, GMRES and k-iteration runs stay on `code/msspec_gmres` (hash `b64d4c64`, unchanged).
- Binary: `_phd_se_noso_nosp_nosym.so`, sha256 starts `e009c24cbdbaec58`. `run_set.py` refuses
  to start series runs if it differs (`SERIES_SO_SHA`). Each row has a `series_src` column.
- Check done: the series Fortran source and the `.so` are identical (`diff -rq`, `cmp`) to
  `msspec_ML/msspec_dev/msspec_python3`. On smoke case L43 (12 atoms, Debye on, 15 theta points)
  orders 1, 2, 3 are finite: eps_2 0.067, 0.094, 0.094; eps_inf 0.241, 0.353, 0.366 (smoke only, not a result).
- Consequence: the series and the LU reference come from different MsSpec code bases (same
  physics settings). The paper must say so. Not checked: the dev tree's LU against ours;
  orders 4 to 8; high-rho cases; that the dev series equals our series where ours is finite.
- The series source is outside the `v3-gmres` hash. If the dev tree changes, the `.so` hash gate stops the run.

## D26 (2026-10-31): E16 emitters changed to a diverse set (author approved)

- Change: PROTOCOL.md A3.2b replaces "one surface emitter per stratum" with four
  cases: H11 (BaTiO3, 1s, coordination 2, radial 0.80), H10 (NaS, 2p, coordination 5,
  radial 0.73), H08 (ReMg, 4f7/2, coordination 12, radial 0.48), L45 (InTe, 3d5/2,
  coordination 6, radial 0.10). Source: `platform/materials_v2.csv` (read-only check
  in this turn, not written to any file).
- Reason: the small test used one emitter on a 12-atom cluster with 4 nearest
  neighbours (a boundary atom, not a deep emitter). The author asked for diverse
  emitters.
- Kinetic energies: each case keeps its own energy from its run folder (269.64, 44.15,
  60.65 and 434.53 eV). The author has not confirmed one energy for all four; this is
  an assumption to confirm before the run.
- The coordination and radial measures are proxies, not a measured depth below a surface.
- Not changed: the map grid (5,269 directions versus a product grid) is still open.
- Not checked: whether the four emitters give visibly different convergence in the small
  test; the depth proxy against a true depth.

## D27 (2026-10-31): E16 full run settings (author: per-emitter energy, product grid, go)

- Energy: each emitter at its own kinetic energy from `platform/materials_v2.csv`
  (H08 60.65 eV, H10 44.15 eV, H11 269.64 eV, L45 434.53 eV).
- Grid: product grid, theta -70 to 70 in 71 steps (2 degrees), phi 0 to 350 in 36 steps
  (10 degrees), 2,556 directions. The author chose "product grid" but not the ranges;
  these ranges are my choice, fixed here before the run. The small test used theta -70 to 70 in 101 steps at phi 0.
- Methods: LU, k0 to k4 (`MSSPEC_TREE_STOP_AFTER`), converged GMRES, one repetition.
- Dense LU for H11 (118 atoms, predicted 9.07 GiB) is refused by rule A2.9 (swap on).
  H11 runs its GMRES and k-step panels now, without an LU reference, until the author decides.
- Not checked: run time per case (not a benchmark, A2.8 warm-up not in place).

## D28 (2026-10-31): swap-off condition removed from the dense LU rule (author's decision)

- Removed: the A2.9 condition that dense LU above 8 GiB runs only with swap off. The
  runner's `lu_decision` no longer refuses on swap; PROTOCOL.md A2.9 and the runner
  docstring are amended to match.
- Kept: the prediction (0.667 MiB x atoms^2) must be below MemAvailable minus 3 GiB; the run
  is killed when MemAvailable falls below 2 GiB; a refusal is written as `lu_not_run`.
- Reason (author): the swap rule does not help the research question and blocks cases that
  fit in RAM (H11, 118 atoms, predicted 9.07 GiB, 30 GiB machine).
- Residual risk, stated: with swap on, the kernel may page instead of failing. The
  MemAvailable and 2 GiB floor checks are the only protection.
- Not checked: the effect on memory and time of an LU run that partly pages.

### Correction to D27 (2026-10-31, added after the failed E16 run)

- The E16 grid in D27 was wrong for a map. Theta from -70 to 70 with phi from 0 to 350
  counts each direction twice: (theta = -t, phi = p) is the same direction as
  (theta = t, phi = p + 180). A map needs theta from 0 (surface normal) to 70.
- The failure message "DIMENSION OF NPH_M TOO SMALL, SHOULD BE AT LEAST 12496" comes from
  `treat_phd.f` lines 66 to 70 (the stereographic-scan branch). That branch recomputes the
  number of phi points as (PHI1 - PHI0) (NTHETA - 1) / (THETA1 - THETA0) + 1 = 350 x 70 / 140 + 1
  = 176, and 71 x 176 = 12496 (recomputed with Python). So the limit was hit by the grid
  I chose, not by a fixed size of the build. NPH_M is passed at run time
  (`memalloc/allocation.f` line 112), not compiled in. My earlier explanation to the
  author (a compiled array, options "raise NPH_M" or "fewer phi") was wrong.
- No E16 result exists. The grid must be chosen again before a new run.

## D29 (2026-10-31): map scans fixed (run-time NPH_M, theta from 0); E16 map small test

- Fix, runner only (solver source and hash `b64d4c64` unchanged): `code/runners/worker.py`
  passes `malloc={'NPH_M': 200000}` to `get_theta_phi_scan`, the run-time override used by
  the stereo examples in `msspec_gpu/` (e.g. `stereo_Cu.py`). `run_set.py` refuses a map
  scan with negative theta.
- E16 grid (replaces the D27 grid): theta 0 to 70 in 36 points (2 degrees), phi 0 to 360 in
  181 points (2 degrees), `--scan thetaphi:0:70:36:0:360:181`. MsSpec chooses its own phi
  points on each theta ring: the result has 3,493 directions (36 theta values), not 36 x 181.
- Small test (A3.6), L45 (InTe, 51 atoms, 434.53 eV), `platform/work/runs/e16_maptest/`:
  LU, k0 and converged GMRES all rc 0, 3,493 of 3,493 values finite. Against LU:
  converged GMRES eps_2 0.00313, eps_inf 0.02433 (pass); k0 eps_2 0.48509, eps_inf 8.71432.
  Smoke test only, not a result.
- Not checked: the high-rho emitters (H08, H10, H11) on this grid; LU memory at 118 atoms.
- Update (same day, after the author pointed to prod3000): the E16 grids are aligned with
  earlier work instead of the 0 to 70 grid above.
  - Theta scan: theta -90 to 89 in 1 degree steps (180 points) at the case's phi, equal to
    prod3000 `THETA = np.arange(-90., 90., 1.)` (checked: `np.linspace(-90, 89, 180)` is equal).
    `--scan theta:-90:89:180`. Negative theta is correct here: a theta scan at fixed phi
    crosses the normal.
  - Map: theta 0 to 90 in 46 points, phi 0 to 360 in 181 points (2 degrees), as
    `msspec_gpu/.../stereo_Cu.py`. `--scan thetaphi:0:90:46:0:360:181`. The result has
    5,269 directions, the count named earlier for the E16 map.
  - Small tests, L45: theta scan LU, k0, GMRES rc 0, 180 of 180 finite; GMRES eps_2 0.00324,
    eps_inf 0.02908 (pass); k0 eps_2 0.34640. Map: rc 0, 5,269 of 5,269 finite; GMRES eps_2
    0.00342, eps_inf 0.06263 (pass); k0 eps_2 0.49288.
  - Our LU theta scan of L45 against the saved prod3000 inversion scan of the same run
    folder (`scan_inv.npz`): same theta values, relative L2 difference 5.98e-07, max relative
    difference 6.71e-06. So the runner reproduces the prod3000 settings for this case.

## D30 (2026-10-31): Amendment A4 written (P100 scan, series source, divergence, R_chi, order)

Written after the author's review of the whole plan ("the plan is not being followed").
The decisions are in PROTOCOL.md Amendment A4 (A4.1 to A4.8) and are not repeated here.
No P100 run existed at the time (0 rows). PROTOCOL.md line 144 now says that
`select2.py` is written. GRMES_STATUS.md section 3 was rewritten as the step list.

## D31 (2026-10-31): tolerance study on the E16 cases (author approved 1e-6 and cap 30)

- Scope: the four E16 cases (A3.2b), GMRES at convergence tolerance 1e-2, 1e-3 and 1e-6,
  and GMRES stopped after 0, 1, 2 and 3 steps (k0 to k3), against dense LU.
- The author approved tolerance 1e-6, which is outside PROTOCOL section 2 (1e-2) and the
  E5 list (1e-2, 1e-3, 1e-4), and a step cap of 30 for these runs, as E5 uses.
- The main P100 runs stay at tolerance 1e-2 and cap 12. Only this study uses the other settings.
- Reason: to show how the error falls with more steps and with a tighter tolerance.

## D32 (2026-11-01): supplementary GMRES sets on all P100 cases (A5)

- The author asked for GMRES at tolerance 1e-3 and 1e-6 and at 0 to 3 steps on all 100
  P100 cases; the primary P100 run had only tolerance 1e-2. PROTOCOL.md Amendment A5 adds
  the sets `p100_tol3`, `p100_tol6`, `p100_steps`. The primary result is not changed.
- D31 said the main P100 runs stay at tolerance 1e-2 and cap 12. They do; the new sets are
  supplementary and are reported separately.
- The LU reference of `p100` is reused (copied), so LU is not repeated three times.

## D33 (2026-11-01): build manifest regenerated after the release-script change

- `code/check_release.sh` gained a block that scans `manuscript/` for forbidden names (paper v2).
  The manifest hashes that file, so the source hash changed from `b64d4c64fd40` to `8a6fb984ffd5`.
- A per-file comparison against the old manifest this turn showed that only `check_release.sh`
  differs (256 other source files and all 7 binaries are identical). No solver code changed.
  The P100 rows keep `b64d4c64fd40`; E3 rows carry `8a6fb984ffd5`.

## D34 (2026-11-01): E3 run on the P100 cases, not on MAT-S, MAT-B and RND

- PROTOCOL.md E3 names MAT-S, MAT-B and RND with map scans. This run used the 100 P100 cases
  with the theta scan (theta -90 to 89, 180 points), the same scan as the rest of P100, and
  the iterate dump switched on (`code/runners/e3_run.py`). The sets: `runs/p100_e3` (LU, fixed
  steps 0 to 6, converged GMRES at tol 1e-2, cap 12, restart 6) and `runs/p100_e3_dump`.
- Reason: E3 needs the error of the cross section against the residual, and the P100 cases already have
  the dense-LU reference used by the paper. The pre-registered E3 set is still to be run.
- Result: 97 cases with a reference (3 have no LU), 92 converged. Analysis in
  `analysis/e3_convergence.py` and `analysis/plot_e3.py`; macros in `analysis/make_results_numbers.py`.

## D34 (2026-11-01): iteration cap raised to 67 for a rerun of the 5 clusters without a result

- The author asked to rerun the 5 clusters that gave no GMRES result at tolerance 1e-2, cap 12, with cap 67.
- The cap was limited to 30 in `parameters.py` (linear_solver_max_iterations) and in `axb_mat_ms_la.f` line 345;
  both limits are now 67. No array is sized by the cap. Rebuilt; regression case `tree` PASS; manifest regenerated
  (source hash `73f86005`).
- Set `runs/p100_cap67` (H13, H17, H23, H24, H34; tolerance 1e-2, restart 6, cap 67); LU copied from `p100`.
  All 5 converged (11 to 22 steps). Primary p100 result (cap 12) is unchanged; the paper reports both.
- New measures eps_rel, R (theta-adapted Pendry), R2, Pearson r: `analysis/extra_metrics.py`.

## D35 (2026-11-01): rerun of the 29 clusters without a result at tolerance 1e-6, cap 67

- The author asked to rerun the 29 clusters that gave no GMRES result at tolerance 1e-6, cap 30, with cap 67.
- Set `runs/p100_tol6_cap67` (tolerance 1e-6, restart 6, cap 67, build `73f86005`); LU copied from `p100`.
- Result: 28 of 29 converged (28 to 62 steps); one strong scatterer (rho 2.29, 144 atoms, H24) did not reach 1e-6 in 67 steps.
  Together with the cap-30 runs, 96 of 97 clusters have a result at 1e-6. Primary `p100_tol6` is unchanged.

## D36 (2026-11-02): map runs of BaTiO3 (118 atoms) need NDIM_M raised at run time

- The map run `runs/maps_e16` (4 emitters, theta 0 to 90, 5,269 directions) failed for BaTiO3 (H11, 118 atoms) in every
  method with "DIMENSION OF THE ARRAYS TOO SMALL IN THE INCLUDE FILE FOR THE TREAT_PHD SUBROUTINE, NDIM_M BE AT LEAST 116564"
  (`treat_phd.f` line 749; the default NDIM_M is 100000, `calculator.py` line 406). The other three emitters ran.
- Fix, runner only: `code/runners/worker.py` passes `NDIM_M: 400000` next to `NPH_M: 200000` in the map branch (same
  run-time `malloc` override as D29; the solver source and hash `73f86005` are unchanged). The array size has no effect on
  the numbers. The H11 rows of the failed run were moved to `platform/work/failed_H11_ndim/` and H11 was run again.
- Not checked: that the other three maps are identical with the larger NDIM_M (they were not rerun).

## D37 (2026-11-02): energy series E17 run as energy scans (`get_energy_scan`), 4 clusters, 10 to 100 eV in 5 eV steps

- PROTOCOL A3.3 fixed a grid of 6 energies (10, 20, 50, 100, 150, 200 eV) as separate theta scans. The author asked (2026-11-02)
  for an energy scan with `MSSPEC.get_energy_scan()` and then for 4 high-rho samples from 10 to 100 eV in 5 eV steps
  (19 energies). Changed: one `get_energy_scan` per (cluster, method), `--scan energy:10:100:19:0` (polar angle 0, the
  case's phi), methods LU, GMRES (tol 1e-2, cap 12, restart 6) and series order 3. The 150 and 200 eV points are dropped.
- Clusters (all rho >= 1.2 at their own production energy; sizes, orientation classes and emitter depth proxies from
  `platform/select_e17.py`): H42 Ru hcp 0001 polar (26 atoms), H17 PbB rocksalt 011 low-index (48), H28 TiO2 rutile m122
  high-index (85), H48 Hf hcp 10m1m1 semipolar (106). rho at each energy is not recomputed in this run.
- Runner only (`run_set.py`: `energy:` scan type and `iters_list` column; `worker.py`: `get_energy_scan`); solver source and
  hash `73f86005` unchanged. `get_energy_scan` fails in its plot-title code when `theta` is a list, so one polar angle per scan.
  Each energy has 10 right-hand sides; `iters_list` holds 10 counts per energy in order.
- Small test (Ru, 5 energies, `runs/e17_scan_test`, load gate relaxed to 4 for that test only): GMRES 0.32 %, series 76.4 % against LU.
- Cap 67 for the energy scans: a Fortran STOP at the cap ends the whole process, so the cap-12 scans of PbB and TiO2 (10 eV
  needs 18 and 19 steps) were lost; those two, and the timing set `e17_timing` and `e17_big`, use cap 67 (`runs/e17_scan_cap67`).
- 2026-11-04, fifth cluster for the energy scan: B01 LiF rocksalt m201 (277 atoms, from `platform/materials_big.csv`), the largest
  cluster that dense LU can run here at 10 to 100 eV. The runner's LU memory model (0.667 MiB per atom^2, for lmax 8 at production
  energy) predicts 51 GB; the measured LU peak in the energy scans is 0.12 MiB per atom^2 (Hf, 106 atoms: 1336 MiB; LiF: 7974 MiB
  of about 20 GB free). New option `--lu-mib-per-atom2` (default 0.667, used 0.15 for this set); the RAM watchdog (stop below 2 GB
  free) stays. Sets: `e17_big_test` (3 energies, not cited), `e17_big` (19 energies, warm-up + 1 timed run per method, fixed order),
  `e17_timing` (four smaller clusters, warm-up + 3 timed runs, rotated order). Solver hash `73f86005` unchanged.

## D38 (2026-11-04): default iteration cap raised from 12 to 67 (code and runner)
Author's decision. `linear_solver_max_iterations` default in `parameters.py` and `--cap` default in `code/runners/run_set.py` are now 67 (the parameter maximum). Reason: the paper reports 67 as the maximum and 12 or 30 only as first-run caps; a default of 12 stopped strong scatterers before they converged. Results already saved were run with the cap in their `row.json` (12, 30 or 67) and are unchanged. Runs started without `--cap` from now on use 67. Tolerance 1e-2 and restart 6 are unchanged.

## D39 (2026-11-06): dense LU run on the 7 clusters that the LU memory estimate had excluded (JCTC review)
Reason: the runner's LU memory model (0.667 MiB per atom^2) was a prediction, not a measurement (D37 already measured 0.12). It
marked 3 test clusters `lu_not_run` (H19, H37, H49: 148 to 149 atoms, rho 1.22 to 1.30) and was used to call the 4 large clusters
(B01 to B04, 262 to 277 atoms) too large for dense LU. Measured: dense LU ran on all 7 on the 31 GB machine (peak 2.4 GiB for the
three 148-atom clusters, 7.0 to 7.8 GiB for the four large ones; wall 46 s and 236 to 481 s).
- Sets: `runs/big_lu` (B01 to B04, `theta:-90:89:180`) and `runs/p100_lu_excluded` (H19, H37, H49, same scan), `--lu-mib-per-atom2 0.12`,
  build `88affa55a477`. The runner's default `--expect-hash` (`b64d4c64`) was stale and was passed on the command line.
- The GMRES and series runs compared with these references are the saved ones (`runs/big_ref`, `runs/big`, `runs/p100`); no solver
  setting changed. The 97-case statistics of the test set are unchanged; the 3 extra references are reported separately
  (`analysis/lu_new_refs.py`, `analysis/out/lu_new_refs.csv`).
- Paper: the claims that these clusters are "too large for dense LU" are removed; the large-cluster figure now uses dense LU as reference.

## D40 (2026-11-06): option to switch the tree preconditioner off (ablation M = I) (JCTC review)
Reason: the review asked what the tree preconditioner contributes beyond restarted GMRES itself. The diagonal blocks of A are the
identity, so "no tree" is GMRES with M = I.
- Code: `MS_TREE_SOLVE` returns at once when `MSSPEC_TREE_NOPRECOND=1` (read once, default off, so the default path is unchanged:
  regression case `tree` PASS). Runner method `nopre` (same tol, cap, restart as `gmres`, same threads). New build hash `4aa9cef1b73e`
  (`code/BUILD_MANIFEST.json` regenerated). The tree set-up is still done in a `nopre` run (0.1 to 0.2 s), so its wall time includes it.
- Settings unchanged: restart 6, tol 1e-2, cap 67; `PROTOCOL.md` unchanged.
- Small tests (not cited as results): `runs/nopre_test` (BaTiO3 118 atoms, InTe 51 atoms), `runs/nopre_test_tol2` and `nopre_test_tol6`
  (three strong scatterers). Steps with the tree against without, at tol 1e-2: 8-9 vs 9-11, 6-7 vs 7-9, 14-15 vs 22-24, 11-13 vs 14-16,
  13-14 vs 15-16; at tol 1e-6: 59-62 vs no convergence in 67, 38-41 vs 48-51, 44-45 vs 50-52. Errors against dense LU are alike.
  So on these clusters the tree lowers the step count by about 15 to 35 %, and most of the gain over the series comes from GMRES itself.
- Pending the author's "go": the full set (97 clusters, tol 1e-2 and 1e-6).
- Full ablation run (2026-11-06, author's go): `runs/nopre_tol2` and `runs/nopre_tol6` (100 clusters each, build `4aa9cef1b73e`, restart 6, cap 67, scan
  `theta:-90:89:180`, start load at most 1.48). Analysis `analysis/ablation_nopre.py` -> `analysis/out/ablation_nopre.csv`, `figures/ablation_nopre.png`.

## D41 (2026-11-06): cap 100 for the 1e-6 cases that did not converge at 67 (author's request)
`MSSPEC_TREE_SERIES_ORDER` guard in `axb_mat_ms_la.f` and the parameter limit `linear_solver_max_iterations` raised from 67 to 100 (arrays depend on the
restart, not on the cap); the default stays 67. Build hash `9cf29d2d`; regression case `tree` PASS. Sets `p100_tol6_cap100` (tree GMRES) and `nopre_tol6_cap100`
(no tree), tol 1e-6, restart 6, scan `theta:-90:89:180`. Results: tree converges on H23 (62 steps, was done at 67), H24 (144 atoms, rho 2.29: 79-85 steps,
eps2 2.1e-6 against dense LU) and H49 (148 atoms, rho 1.26: 32-33 steps; its earlier 'failure' was a first run at cap 30 that was never rerun at 67).
Without the tree: H23 87-95 steps, H49 38 steps (cap-67 run: 'ok' too), H24 does not converge in 100. So at tol 1e-6 and cap 100 the tree converges on
all 100 clusters, and GMRES without the tree on 99.

## D42 (2026-11-06): cap 200 for the one no-tree case that did not converge in 100 (author's request)
`MSSPEC_TREE_SERIES_ORDER` guard and the parameter limit `linear_solver_max_iterations` raised from 100 to 200 (default stays 67). Build hash `2198ac08`;
regression case `tree` PASS. Set `nopre_tol6_cap200` (H24: 144 atoms, rho 2.29, tol 1e-6, restart 6, no tree). It converges in 108-116 steps (wall 129 s,
rss 268 MiB); eps2 2.0e-6 and eps_inf 1.1e-5 against dense LU. With the tree: 79-85 steps, 99 s, eps2 2.1e-6. So at tol 1e-6 with cap 200, GMRES
converges on 100/100 clusters both with and without the tree.

## D43 (2026-11-06): series of order 4 on a subset (author's request, JCTC review)
`PROTOCOL.md` limits the series to orders 1 to 3 (A.11). The review asked for a higher order. Order 4 was run with the series binary of D25 (`--series-src` default),
`theta:-90:89:180`, set `runs/series4`, build gate `2198ac08`. It is much slower than orders 1 to 3 (4 to 111 times order 3; 3.9 s to 482 s), so the run was stopped by
the author after the first 13 clusters of the test set (H01 to H13, all strong scatterers, rho 1.21 to 2.81, 34 to 141 atoms); `runs/series4/H14` is an unfinished case
with no row. Order 5 on a 107-atom cluster (H13) had not finished after 134 min (one core) and was stopped by PID (`analysis/out/series5_stopped.csv`); no result.
Order 4 errors against dense LU (eps2) are 13 to 29,288 %, not systematically below order 3 (lower on 9 of 13). Analysis: `analysis/series4_subset.py`
-> `analysis/out/series4_subset.csv`, `figures/series4_subset.png`; paper Section 4.2 states it as a subset. Orders 1 to 3 stay the comparison on all 100 clusters.

Note on D42 (added 2026-11-06): the errors there ("eps2 2.0e-6, eps_inf 1.1e-5" without tree; "eps2 2.1e-6" with tree) are fractions, not percent. In percent
(the unit of the paper) they are 2.0e-4 %, 1.1e-3 % and 2.1e-4 %.

## D44 (2026-11-06): keep the finite-value guards of the solver in the build (JCTC review)
Finding: `-ffast-math` (`GFORTRAN_FFLAGS` in `options.mk`) implies `-ffinite-math-only`, which lets gfortran fold `IEEE_IS_FINITE` to true. A test program with the same flags
does not fire its guard on a NaN (it does at `-O2`). In `axb_mat_ms_la.f` (four guards, lines 834, 881, 919, 1005 and the NaN test in the Arnoldi step) the object
file had 16 floating-point comparisons; with `-fno-finite-math-only` it has 24 (26 at `-O2` without `-ffast-math`).
Change: `phd_gm_noso_nosp_nosym.mk` adds `-fno-finite-math-only` for this one object only (target-specific `FFLAGS`); all other files keep `-O2 -ffast-math`. No source change, so
`source_tree_sha256` (`2198ac08`) is unchanged; `_phd_gm_noso_nosp_nosym.so` changed (`55f863fa...` -> `91934e73...`, `code/BUILD_MANIFEST.json`).
Checks: `tests/test_tree_gmres.py` PASS; GMRES and no-tree runs of 4 clusters (`runs/fm44_test`, 8 runs) give cross sections bit-identical to the saved runs
(`runs/p100`, `runs/nopre_tol2`). Not tested: a run that really produces a NaN (none occurred in the 729 saved GMRES and no-tree rows, `n_nonfinite` 0 in all).
The solver settings (restart 6, tolerances, cap) are unchanged. The saved runs were made with the earlier binary; the results are identical on the 8 re-run cases only.

## D45 (2026-11-06): H24 at tol 1e-3 rerun with cap 100 (author's request)
Set `p100_tol3_cap100` (H24: 144 atoms, rho 2.29, tol 1e-3, restart 6, cap 100, build 2198ac08, scan `theta:-90:89:180`). The first run (cap 30) had no result; the rerun
converges in 33-36 steps (wall 55.8 s). Paper Section 7.3 and the cap sentence in Section 5.3 now state the caps of each set and the two cap-100 reruns of H24 (macros `CapH*`);
the main-table statistics are unchanged (still without these two runs).
