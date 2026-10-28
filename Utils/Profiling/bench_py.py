#!/usr/bin/env python3
"""
Times the ODE steppers through EOBRunPy, for systems that take well under a second, where the
executable's start-up and file output would dominate. Repeat calls and keep the minimum.
Prints, per stepper, time, speed-up and mismatch of the (2,2) against rkf45.

Usage:
  python3 bench_py.py [--reps 200] [--bns]
"""

import argparse
import os
import sys
import time
import numpy as np

# the module built in Python/, not an installed one
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
    # warm-up call
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
    """Flat-PSD mismatch, maximized over time shift and phase"""
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
