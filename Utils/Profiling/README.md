# TEOBResumS profiling & optimization (eccentric BBH/BNS)

Tooling and findings for speeding up eccentric/BNS waveform generation, where the
full EOB dynamics ODE must be integrated (no post-adiabatic shortcut).

## Build targets (`C/Makefile`)

| target | flags | use |
|--------|-------|-----|
| `make` | `-g`, no `-O` | default; **unoptimized** — do not use for timing |
| `make opt` | `+ -O2` | realistic timing, baselines, regression checks |
| `make prof` | `+ -O2 -pg -DPROFILING=1` | gprof + PROF section timers / RHS counters |

## Instrumentation

`C/src/TEOBResumSProf.c` (+ macros in `TEOBResumS.h`) adds `clock_gettime` section
timers and RHS-call counters behind the `PROFILING` macro. When `PROFILING=0` the
`PROF_*` macros expand to nothing, so the default build is **bit-identical** to the
un-instrumented code (verified: all data files match a pristine `HEAD` build).

Timed sections: `ODE_loop`, `RHS`, `flux`, `Fphi_ecc`, `metric+rc`, `wav_store`,
`NQC`, `interp_spline`, `compute_hpc`. Counters: `rhs_calls_total`,
`rhs_calls_storage` (per-accepted-step), `rhs_calls_stepper` (= total − storage).
Output goes to `prof_timers.txt` in the run's `output_dir` and is echoed to stdout
as `# PROF ...` lines.

> The fine-grained RHS/flux timers fire ~1e5–1e6×; their own `clock_gettime`
> overhead perturbs those wall-times (and `prof_timer_start` shows up in gprof).
> Trust the **RHS-call counts** and the **gprof flat profile** for the true
> per-function breakdown; read the section wall-times as relative indicators.

## Scripts

- `profile.py` — builds `prof`, runs every config in `parfiles/`, prints a section-
  timer + RHS-count table and per-config gprof top functions. `--csv`, `--config`,
  `--no-build`, `--no-gprof`.
- `bench_py.py` — in-process (`EOBRunPy`) stepper timing + accuracy for sub-second
  BBH, where the C-exe fork/exec + I/O overhead swamps the computation. Repeats
  the identical call and takes the min. `--reps`, `--bns`.
- `stepper_grid.py` — random parameter-space sweep of stepper mismatch vs rkf45
  (BBH + BNS), to validate a stepper before promoting it. `--n`, `--nbns`, `--seed`.
- `compare_revs.py` — before/after regression report across two git revisions
  (walltime speedup + mismatch of new vs base). Worktree-based; safe cleanup.
- `freeze_baselines.sh [DEST]` — builds `opt`, regenerates all configs, stores data
  files + an `md5` manifest for bit-identical regression checks.

## Benchmark configs (`parfiles/`)

`ecc_bbh_q1_e01` (q1 nonspin e=0.1), `ecc_bbh_q3_e03` (q3 e=0.3, stresses
non-circular flux), `ecc_bbh_q5_s09s-04_e01` (spinning q5 e=0.1),
`ecc_bns_e005` (BNS + tides, e=0.05), `ecc_bns_long` (low-freq BNS, output OFF,
multi-second — the reliable timing benchmark), `qc_bbh_control` (quasi-circular
Giotto regression control).
