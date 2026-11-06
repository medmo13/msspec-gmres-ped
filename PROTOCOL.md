# Pre-registered protocol: tree-preconditioned GMRES for PED

Status: written 2026-10-19, before any benchmark run. Each change after this
date goes in `DEVIATIONS.md` with the reason and the date. No result from the
pilot changes a setting below; the pilot may set only the timeouts, the run
budget and the memory preflight constant.

## 1. Build and platform

- Build: the frozen copy `code/msspec_gmres`, recorded in `code/BUILD_MANIFEST.json`
  (source and binary hashes). Tag `v1-gmres` marks the frozen v1 state. From A2 the
  build is v2 (tag `v2-gmres`, fallback code removed); see A2.1.
- Machine: one workstation, Intel i7-13700H, 20 logical CPUs, 30 GiB RAM.
- Threads: 8 for both methods. Dense LU uses 8 BLAS threads. GMRES uses 8
  OpenMP action threads and 1 BLAS thread.
- One process at a time. A fresh process for every run.
- Timing runs use a clean slot: no other load, governor recorded, load average
  recorded before and after the run.

## 2. Solver settings (fixed here, not tuned on results)

- Restart length m = 6.
- Tolerance tau = 1e-2 (default), with 1e-3 and 1e-4 only in E5.
- Iteration cap = 12 (default), 30 only in E5.
- Series orders for E11: 1, 2 and 3 only, plain series, no renormalised variant.
- Dense LU fallback: superseded by A2. The v2 build has no dense fallback. A failed
  solve stops with `STOP 99` and `MSSPEC_SOLVER_STATUS ... none_available`, and the
  study counts it as a method failure (A2.2).
- The C3 structured status is out of scope for v1 and v2.

## 3. Metrics

Computed by one script (`analysis/evaluate.py`) from the saved arrays.

- eps_2: relative L2 error of the intensity against the same-session dense LU.
- eps_inf: worst-direction relative error. The denominator is floored at 1 % of
  the pattern maximum (primary). The unfloored value goes in the supporting
  information.
- R_chi: reliability factor between modulation functions (maps only). Always reported.
- Pass gate: eps_2 < 10 % and eps_inf < 10 %.
- Cost: complete-process wall time, MS-stage solve time, peak RSS, iterations,
  number of operator actions, restart cycles.
- Every failure is listed. No case is replaced.

## 4. Statistics

- Per case: median over repeats.
- Across cases: median, with a percentile bootstrap 95 % CI (10,000 resamples,
  fixed seed).
- Speed: two-sided Wilcoxon signed-rank test on log time ratios.

## 5. Benchmark sets

| ID | Systems | Scan | Methods, repeats | Question |
|---|---|---|---|---|
| E2a MAT-S | 29 materials, 100 eV, surface emitter | map, 5,269 directions | LU, GMRES; 3 repeats | Q1, Q4 |
| E2b MAT-B | same clusters, sub-surface emitter (seeded) | map | LU, GMRES; 3 repeats | Q1, Q4 |
| E3 Iterations | MAT-S, MAT-B, RND | map, iterate dump (C2) | 1 run each | Q2 |
| E4 Restart | MAT-B, m in {2, 4, 6, 8, 12} | map | 1 repeat per m, clean slot | Q3 |
| E5 Tolerance | MAT-S, MAT-B, tau in {1e-2, 1e-3, 1e-4}, cap 30 | map | 1 repeat, clean slot | Q3 |
| E6 Preconditioner | 5 smallest MAT-B clusters: none vs tree | map, timeout 20x LU time | 1 repeat | Q4 |
| E7 Action | Cu size series: axial vs direct Gaunt action | theta scan | 3 repeats (1 at the largest) | Q4 |
| E8 Size | Cu hemispheres, 3 to 7 planes, 300 eV | theta scan and one map per size | LU, GMRES; 3 repeats (1 at the largest) | Q4 |
| E9 RND-21 | 21 stratified cases | map and theta scan | LU, GMRES; 3 repeats | Q1, Q4 |
| E10 Energy | 3 materials (light, transition metal, oxide), 6 energies 50 to 800 eV, 4 planes | map | LU, GMRES; 2 repeats | Q1, Q2 |
| E11 Series | MAT-B (and MAT-S if time): series orders 1, 2, 3 | map | 3 repeats, clean slot | Q5 |
| E12 Robustness | regression test and count of method failures over all runs (A2.2) | - | - | Q6 |

Questions: Q1 accuracy against LU; Q2 convergence; Q3 restart and tolerance;
Q4 cost and memory; Q5 comparison with series orders 1 to 3; Q6 robustness.

From A2, the sets E2a, E2b, E9 and E11 are not main sets; their evidence moves to
the 100-material set P100 (A2.3). The sections above stay as written, with the
superseded parts marked.

## 6. Execution order and gates

1. Pilot (at most one day, not reported): 3 MAT-B cases and 2 E8 cases, to set
   timeouts, run budget and the memory preflight constant only.
2. E2a, E2b, E3. Gate: no unexplained failure. If eps_inf fails often at tau = 1e-2,
   the default is not changed; the result is reported and E5 answers it.
3. E4 and E5, clean-slot timing.
4. E7 and E8, largest sizes last.
5. E9 and E10.
6. E6 and E11, then the E12 summary.

Runner rules: one resumable runner per set. Raw CSV rows are written at once.
LU and GMRES arrays are saved as `.npy` for the metrics. Each row records load
average, governor, start and end time, return code and solver status lines.
A memory preflight runs before LU. Each run has a timeout. The log shows a
progress counter.

## 7. Reporting

Every number in the paper comes from `analysis/numbers.tex`. Captions state the
scan type, the timing scope (complete process or MS stage) and the repeats.
Drawn figures say "schematic". Each limit is stated next to its claim.

## Amendment A1 (2026-10-23, before any reported run)

Added after the pilot size series (`pilot/small_compare.*`, not reported), which
showed that dense LU needs 13.6 GiB at 143 Cu atoms (300 eV) and would need
about 35 GiB at the next size. No reported run exists yet, so these sets are
pre-registered, not deviations.

| ID | Systems | Scan | Methods, repeats | Question |
|---|---|---|---|---|
| E13 Threads | 3 sizes of E8 (smallest, middle, largest with LU); threads 1, 2, 4, 8 | theta scan | LU, GMRES; 3 repeats, clean slot | Q4 |
| E14 Beyond LU | Cu hemispheres above the LU memory limit (sizes fixed by the author, see `manuscript/ARTICLE_DESIGN.md` section 11) | theta scan and one map per size | GMRES only, tau = 1e-2; reference: converged GMRES tau = 1e-4; check against tau = 1e-5 on the smallest size; 1 repeat | Q7 |

Q7: how large a cluster a workstation can solve, and how accurate it is there.
In E14 every caption says that the reference is converged GMRES, not LU.

## Amendment A2 (2026-10-31, before any reported run)

Added after the v2 build (no dense fallback), the rho scan (`platform/rho_scan.csv`)
and the pilot memory check. The earlier text stays; the items below replace the
parts they name. No reported run exists yet.

- **A2.1 Build.** The study uses the v2 code with the iteration-0 dump added (D20,
  2026-10-31, corrected by D21, stop switch added by D22): source hash `b64d4c64...` (`code/BUILD_MANIFEST.json`). The solver
  numerics are unchanged. The tag `v2-gmres` (local only) still points at commit
  `3fe9a7c`, which has the old source (hash `1e7eab4c...`); it is not the build used
  for reported runs until the new commit is made and tagged. `v1-gmres` is the frozen
  record of the fallback build and is not used for reported runs.
- **A2.2 Method failure, not recovery.** The v2 build has no dense fallback. A failed
  GMRES solve stops with `STOP 99` and `MSSPEC_SOLVER_STATUS ... none_available`. The
  run is a method failure: it is written as an error row (`error=1`, `error_reason`,
  return code, status line), counted as a fail in P1, and never replaced by LU.
  E12 is the count of method failures over all runs, with the reasons listed.
- **A2.3 Main sets.** E2a, E2b, E9 and E11 are not main sets. Their evidence moves
  to P100, the 100-material set in A2.4 to A2.6.
- **A2.4 One settings set.** All 100 cases use the prod3000 settings of the rho scan:
  lmax 8, surface emitter, and the `use_debye_model` value recorded in each case's
  `config.json`. Every run row records lmax and the Debye flag. A case whose value
  differs from the others is excluded and listed. Stratum (high or low rho) is a label
  and a covariate, not a second settings set.
- **A2.5 Strata and selection.** From the rc 0 rows of `platform/rho_scan.csv`
  (600 rows; the two rc 1 rows, Co 80 atoms and Bi 91 atoms, are listed as excluded).
  High stratum: rho at least 1.2 and at most 150 atoms (99 rows, 79 distinct
  materials); choose 50, stratified by prototype and orientation class, with a seed
  fixed in the code and recorded. Low stratum: rho below 1 and at most 150 atoms
  (50 chosen, stratified the same way). The selection script is `platform/select2.py`
  (written 2026-10-31, D24; `select.py` is kept unchanged). This replaces the high group of D11.
- **A2.6 Metrics and statistics.** eps_inf floor: a single global value,
  1 % of the pattern maximum, as in section 3. Bootstrap: percentile, 10,000
  resamples, resampling by material, seed 20261030. Wilcoxon on log time ratios uses
  converged pairs. Wilson 95 % intervals for every pass and failure rate.
- **A2.7 Endpoints.** P1: GMRES pass rate (fails include method failures); P2:
  method-failure rate with reasons; P3: series orders 1 to 3, pass and divergence
  rates per stratum. S1: median eps_2 and eps_inf over converged runs with the
  bootstrap CI. S2: time and memory ratios on the timing subset. S3: descriptive
  regression on log rho and log atoms. The three P endpoints are the family for
  multiple comparisons; the rest are descriptive.
- **A2.8 Timing.** Complete-process wall time, with the MS-stage time recorded
  beside it. Untimed warm-up for each method, alternating method order, clean slot,
  load average and governor recorded. Timing subset and timeout constants are set
  from the pilot step (A2.10) and written here with their date.
- **A2.9 Memory.** Dense LU needs about 16 n^2 bytes. The pilot CSV gives
  `lu_rss_mib` = 977, 3,425 and 13,643 MiB at 34, 70 and 143 atoms, which is
  0.845, 0.699 and 0.667 MiB per atom^2. The fit is not exact, so it is checked
  against E8 before any prediction is quoted. Under the 0.667 fit, the 8 GiB
  preflight limit is at about 111 atoms and 30 GiB at about 215 atoms. The rule of
  `PLAN.md` section 8 applies unchanged: run LU when the prediction is below 8 GiB.
  Above it, one LU attempt, alone, is allowed when the prediction + 3 GiB is below
  the free RAM at that moment (MemAvailable), and it is stopped when free RAM falls
  below 2 GiB. The swap-off condition was removed on 2026-10-31 (D28, author's decision). Every attempt writes a row (status, peak RSS, time to stop,
  free RAM before the run). If LU gives no result, the case is "LU not run" (not a
  GMRES failure), is listed, and takes the converged-GMRES reference only where
  A1 (E14) applies.
- **A2.10 Pilot constants.** The pilot is a time-boxed step (about one day, not
  reported). It sets timeouts (20 times the LU time), the timing subset size and
  the preflight constant only. Solver settings do not change.
- **A2.11 Load and runner rules.** Moved here from `PLAN.md`: a run is refused when
  the 1-minute load average is above 1.5 at start. The runner refuses a GMRES run
  whose build hash is not the hash `b64d4c64...` of the current build (D20 to D22; the build
  has no fallback switch), and it records that hash in the row.
- **A2.12 Series and extensions.** Series orders 1, 2 and 3 are pre-registered.
  Orders 4 to 8 are exploratory and labelled as such in the SI. E14 (A1) waits until
  the radius 3.0 failure is diagnosed. E15 (comparison with a measured pattern) is
  not in the plan; the author decides (`manuscript/ARTICLE_DESIGN.md` section 11).
  Note 2026-10-31 (D23): the radius 3.0 diagnosis is done; GMRES runs at 3.0 a0 (259 atoms),
  so E14 is no longer waiting for it. Its sizes still need the author's choice.

## Amendment A3 (approved by the author 2026-10-31; in force)

Proposed by the author on 2026-10-31 and approved by the author in chat the same day.
Approval decisions: A3 approved; A3.1 MemAvailable definition confirmed; A3.2 code
change approved (implemented as D20 to D22, build hash `b64d4c64...`); A3.3 grid confirmed;
A3.4 orders 1 to 3 only. No A3 run may start before the small test of A3.6 is done and
the author says "go". Where A3 and A2 differ, A3 governs the new sets only.

- **A3.1 LU reference up to the memory limit.** Clusters above 143 Cu atoms are run
  with dense LU as reference while memory allows. Using the fit of A2.9 (0.667 MiB
  per atom²): 14 GiB at about 147 atoms, 15 GiB at about 152, 16 GiB at about 157.
  The rule of A2.9 applies: one LU attempt, alone, when prediction + 3 GiB is below
  the free RAM at that moment. "Free RAM" is defined here as **MemAvailable** (the
  kernel's estimate of memory usable without swapping), not the "free" column, which
  is much smaller because of the page cache. At the time of drafting, MemAvailable
  was 18 GiB of 30 GiB. Each LU run records MemAvailable before and after. Clusters
  above about 165 atoms (predicted 18 GiB or more) get no LU result here; they belong
  to E14 (A1).
- **A3.2 Convergence figures (set E16).** Theta scans and stereographic maps of
  the observable for: dense LU (reference); the preconditioner alone (no GMRES
  steps); GMRES after 1, 2, 3 and 4 steps; plus the converged GMRES result. Cases:
  one surface emitter per stratum. Each panel shows the error against LU.
  Method: one full run per panel with `MSSPEC_TREE_STOP_AFTER=k`, k = 0 (the
  preconditioner alone, `M**(-1) b`) and k = 1 to 4, plus one normal converged run
  and one LU run; the pattern of each run comes from the normal pipeline. This is
  implemented (D20 to D22). The runner checks that the iterations printed equal k.
  E16 is unblocked for its small test; the full run waits for the author's "go"
  (A3.6).
- **A3.2b Diverse emitters for E16 (amendment, 2026-10-31, approved by the author
  the same day).** Replaces "one surface emitter per stratum" in A3.2. Four cases
  from `platform/materials_v2.csv`, chosen to span coordination number (nearest
  neighbours, 2 to 12) and radial position (0 = centre, 1 = edge of the cluster):
  H11 (BaTiO3, 1s, coordination 2, radial 0.80), H10 (NaS, 2p, coordination 5,
  radial 0.73), H08 (ReMg, 4f7/2, coordination 12, radial 0.48), L45 (InTe, 3d5/2,
  coordination 6, radial 0.10). Each case keeps its own kinetic energy from its run
  folder (269.64, 44.15, 60.65 and 434.53 eV). The two measures are geometric
  proxies, not a measured depth below a surface. All four cases use the single
  prod3000 settings set (A2.4). The panels and the k-steps are as in A3.2.
- **A3.3 Energy series, high rho (set E17).** High-rho materials (rho ≥ 1.2 at
  the energy used) at kinetic energies from 10 to 200 eV, on a grid fixed here
  (10, 20, 50, 100, 150, 200 eV; confirmed by the author 2026-10-31). For each case and
  energy: rho is recomputed at that energy with the EIG/power method; dense LU is the
  reference; GMRES is run as in A2; the classical MsSpec series expansion is run at
  the orders in A3.4 on the same settings. Series divergence is measured, not assumed:
  the claim that the series does not converge for rho > 1 is tested in this set.
- **A3.4 Series expansion to compare.** The comparison uses the classical MsSpec
  series expansion, not the tree-series path (`MSSPEC_TREE_GMRES=0`, which is a
  study-only mode). The option is `algorithm='expansion'` (Rehr-Albers series), with the
  order set by `scattering_order` (default 3, range 1 to 10) and the basis cut-off by
  `RA_cutoff` (default 1, range 0 to 8). `renormalization_mode` stays at its default
  (None, the plain series). Checked 2026-10-31 in `msspec/parameters.py` of the project
  build and of `msspec_acceleration/msspec_standard`. The reference scripts
  (`msspec_acceleration/test_materials/*_expansion.py`) do not set the order, so they use
  the default 3. Only orders 1 to 3 are run in A3 (author decision 2026-10-31).
  Orders 4 to 8 are not run under A3.
- **A3.5 Materials and settings.** Materials come from the high-rho rows of
  `platform/rho_scan.csv`, but rho must be recomputed at 10 to 200 eV, so the scan
  is not reused as is. The same cap, tolerance and restart as A2 apply.
- **A3.6 Test before full run.** Each new set starts with a small test of about 100
  samples, and the full run waits for the author's "go" (REAL-WORK CLAUDE.md,
  Experiments and long runs).

## Amendment A4 (2026-10-31, approved by the author; before any reported P100 run)

Written after a review of the whole plan (D30). Where A4 and A2/A3 differ, A4 governs.

- **A4.1 P100 scan and methods.** Theta scan from -90 to 89 degrees in 1 degree steps
  (180 points) at each case's phi, equal to prod3000 `THETA = np.arange(-90., 90., 1.)`
  (`--scan theta:-90:89:180`). Methods: `lu, gmres, series1, series2, series3`, one
  repetition (accuracy is deterministic).
- **A4.2 Code sources.** LU and GMRES: `code/msspec_gmres`, build hash `b64d4c64...`.
  Series: the FINAL dev tree
  `/home/mohammed/Desktop/REAL-WORK/msspec_ML/FINAL/07_msspec_dev/msspec_python3/src`,
  `.so` sha256 starting `e009c24cbdbaec58` (D25). The runner refuses any other `.so`.
- **A4.3 Divergence.** A series run has diverged when its result is non-finite or its
  eps_2 is 1 or more. P3 reports, per order and per stratum, the pass rate and the
  divergence rate, each with a Wilson 95 % interval.
- **A4.4 R_chi.** Supporting metric only (definition in D24). It is never used for
  pass or fail.
- **A4.5 LU gate.** As D28: LU runs when the prediction plus 3 GiB is below MemAvailable
  and is stopped if MemAvailable falls below 2 GiB. Otherwise the row is `lu_not_run`;
  it is listed and never counted as a GMRES failure.
- **A4.6 E16 grids.** The theta scan of A4.1, and the map `thetaphi:0:90:46:0:360:181`
  (5,269 directions) with run-time `malloc NPH_M 200000` (D29). Maps start at theta 0.
- **A4.7 Order of work.** P100 accuracy, then the timing pilot and timing subset, then
  E16, then E17.
- **A4.8 LU consistency check.** Each P100 LU theta scan is compared with the prod3000
  `scan_inv.npz` of the same case. A relative L2 difference above 1e-4 is flagged and
  listed; it is not a failure.

## Amendment A5 (2026-11-01, approved by the author; the primary P100 run of A4 was already finished)

- **A5.1 Supplementary GMRES settings on all P100 cases.** The author asked that GMRES at
  convergence tolerance 1e-3 and 1e-6 (step cap 30), and GMRES stopped after 0, 1, 2 and 3
  steps (`k0` to `k3`, tolerance 1e-2, cap 12), are run on all 100 P100 cases. These are
  supplementary sets `p100_tol3`, `p100_tol6` and `p100_steps`. The primary endpoints P1 to
  P3 stay as A4 defines them (tolerance 1e-2, cap 12, set `p100`).
- **A5.2 Reference.** Each supplementary set uses the dense-LU run of `p100` as its
  reference (the `lu_0` folder is copied unchanged into the case folder). Same build
  (`b64d4c64`), same scan. Cases with `lu_not_run` stay without a reference.
- **A5.3 Record.** Wall time, peak memory and load are recorded on every row, as in A4.
  A GMRES run that does not reach the tolerance within the cap is a method failure
  (`none_available`) and counts as a fail, never replaced.
