# Profiling TEOBResumS

Tools to time the code and see where the time goes, mostly for Dali, where the whole dynamics is
evolved with the ODE solver.

## Builds (`C/Makefile`)

| target | flags | for |
|--------|-------|-----|
| `make` | `-g` | the default, not for timings |
| `make opt` | `-O2` | timings and output checks |
| `make opt-lto` | `-O2 -flto` | as `opt`, a few % faster, same output |
| `make prof` | `-O2 -pg -DPROFILING=1` | gprof, and the timers below |

`-pg` needs gcc, Apple's clang does not support it.

## Timers

With `PROFILING=1`, `C/src/TEOBResumSProf.c` times the sections `ODE_loop`, `RHS`, `flux`,
`Fphi_ecc`, `metric+rc`, `wav_store`, `NQC`, `interp_spline` and `compute_hpc`, and counts the
r.h.s. calls (all of them, those that store a step, and those the stepper makes). The results go to
`prof_timers.txt` in the output directory and to stdout. Without it the macros are empty and the
output is the same.

The r.h.s. and flux timers run so often that their own cost shows in what they measure, so compare
the call counts and the gprof profile, and read those times as rough.

## Scripts

All take parfiles on the command line, e.g. those in `C/par`.

- `profile.py PARFILES`: builds `prof`, runs the parfiles, prints timers, call counts and the top
  of gprof. `--no-build`, `--no-gprof`, `--csv`.
- `compare_revs.py --base REV [--new REV] PARFILES`: time of two git revisions, and the mismatch
  of the new (2,2) against the old one.
- `bench_py.py`: times the ODE steppers through `EOBRunPy`, for systems that take well under a
  second, and gives their mismatch against rkf45. `--reps`, `--bns`.
- `stepper_grid.py`: mismatch of the steppers against rkf45 at random points in parameter space.
  `--n`, `--nbns`, `--seed`.
