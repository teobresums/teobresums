#!/usr/bin/env python3
"""
In-process ODE-stepper benchmark for TEOBResumS via EOBRun_module.

For sub-second (BBH) systems the C-executable path is dominated by fork/exec +
parfile I/O + output-file writing, which swamps the ~0.1 s computation with noisy
fixed overhead. Calling EOBRunPy in-process instead measures pure computation,
returns the modes in memory (no disk), and lets us repeat the identical call many
times and take the min - which is what makes sub-second timing reliable.

Usage:
  python3 bench_py.py                 # default eccentric-BBH config
  python3 bench_py.py --reps 200      # timing repetitions per stepper
  python3 bench_py.py --bns           # eccentric-BNS-like config (with tides)

Reports, per stepper, the min/median wall time and speedup vs rkf45, plus the
accuracy of each stepper's (2,2) mode against the rkf45 baseline (time+phase
maximized mismatch).
"""

import argparse
import os
import sys
import time
import numpy as np

# Prefer the freshly-built in-place module in Python/ over any pip-installed copy
# in site-packages (which may lag the current C sources / lack new parameters).
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "Python"))
import EOBRun_module

STEPPERS = ["rkf45", "rk8pd", "rkck", "msadams"]


def modes_to_k(modes):
    return [int(x[0] * (x[0] - 1) / 2 + x[1] - 2) for x in modes]


def base_pars(bns=False):
    p = dict(
        M=40.0, q=3.0, chi1=0.3, chi2=0.0,
        distance=400.0, inclination=0.5, coalescence_angle=0.0,
        use_geometric_units="no", initial_frequency=20.0,
        domain=0, interp_uniform_grid="yes", srate_interp=4096.0,
        ecc=0.1, ecc_freq=2,
        use_mode_lm=modes_to_k([[2, 2]]), arg_out="no",
    )
    if bns:
        p.update(M=2.7, q=1.0, chi1=-0.05, chi2=-0.05,
                 initial_frequency=30.0, ecc=0.05,
                 LambdaAl2=400.0, LambdaBl2=400.0,
                 use_lambda234_fits="YAGI13", use_tidal="TEOBRESUM")
    return p


def time_stepper(pars, stepper, reps):
    p = dict(pars); p["ode_stepper"] = stepper
    # warm up (lazy init / first-touch allocation)
    for _ in range(3):
        EOBRun_module.EOBRunPy(p)
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        EOBRun_module.EOBRunPy(p)
        ts.append(time.perf_counter() - t0)
    ts = np.array(ts)
    return ts.min(), np.median(ts)


def get_h22(pars, stepper):
    p = dict(pars); p["ode_stepper"] = stepper
    p["use_mode_lm"] = modes_to_k([[2, 2]])
    t, hp, hc = EOBRun_module.EOBRunPy(p)
    return np.asarray(t), np.asarray(hp) - 1j * np.asarray(hc)


def flat_mismatch(h0, h):
    """Time-and-phase-maximized mismatch (flat PSD). The time maximization is the
    IFFT of the cross-spectrum; the |.| gives the phase maximization. A naive
    fixed-index phase-at-merger comparison is NOT reliable here: near the steep
    merger a sub-sample time offset between two steppers' output grids shows up as
    a large spurious phase difference, which made rk8pd look inaccurate when it is
    only mis-aligned by a fraction of a sample. Use this maximized match instead."""
    n = max(len(h0), len(h))
    N = 1 << int(np.ceil(np.log2(2 * n)))
    H0 = np.fft.fft(h0, N); H = np.fft.fft(h, N)
    z = np.fft.ifft(H0 * np.conj(H))
    ov = np.max(np.abs(z)) / np.sqrt(np.vdot(h0, h0).real * np.vdot(h, h).real)
    return 1.0 - ov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--bns", action="store_true")
    args = ap.parse_args()

    pars = base_pars(bns=args.bns)
    label = "ecc BNS-like" if args.bns else "ecc BBH"
    print(f"config: {label}  (reps={args.reps}, in-process EOBRunPy)\n")

    _, h0 = get_h22(pars, "rkf45")            # accuracy baseline (current default)
    tmin_ref = None
    print(f"{'stepper':<9}{'min_s':>9}{'med_s':>9}{'speedup':>9}{'mismatch':>12}")
    for st in STEPPERS:
        try:
            mn, md = time_stepper(pars, st, args.reps)
        except Exception as e:
            print(f"{st:<9} FAILED: {e}")
            continue
        if tmin_ref is None and st == "rkf45":
            tmin_ref = mn
        _, h = get_h22(pars, st)
        mm = flat_mismatch(h0, h)
        sp = tmin_ref / mn if tmin_ref else 1.0
        print(f"{st:<9}{mn:>9.4f}{md:>9.4f}{sp:>8.2f}x{mm:>12.2e}")


if __name__ == "__main__":
    main()
