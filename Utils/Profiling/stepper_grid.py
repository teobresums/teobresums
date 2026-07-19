#!/usr/bin/env python3
"""
Broad stepper-accuracy sweep for TEOBResumS (eccentric BBH + BNS).

Randomly samples the (q, chi1, chi2, ecc, anomaly) parameter space and, for each
point, generates the (2,2) mode with the candidate steppers and the rkf45
baseline, then reports the time-and-phase-maximized mismatch (flat PSD) of each
candidate vs rkf45. Purpose: confirm a faster stepper (rk8pd) stays within an
acceptable mismatch everywhere before promoting it to the default, and surface
any parameter regions where it degrades or the solver fails.

Usage:
  python3 stepper_grid.py --n 60            # 60 random BBH configs
  python3 stepper_grid.py --n 60 --nbns 15  # add 15 random BNS configs
  python3 stepper_grid.py --candidates rk8pd,msadams --seed 42
"""

import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "Python"))
import EOBRun_module


def modes_to_k(modes):
    return [int(x[0] * (x[0] - 1) / 2 + x[1] - 2) for x in modes]


K22 = modes_to_k([[2, 2]])


def h22(pars, stepper):
    p = dict(pars); p["ode_stepper"] = stepper; p["use_mode_lm"] = K22
    t, hp, hc = EOBRun_module.EOBRunPy(p)
    return np.asarray(hp) - 1j * np.asarray(hc)


def mismatch(h0, h):
    n = max(len(h0), len(h))
    N = 1 << int(np.ceil(np.log2(2 * n)))
    H0 = np.fft.fft(h0, N); H = np.fft.fft(h, N)
    z = np.fft.ifft(H0 * np.conj(H))
    ov = np.max(np.abs(z)) / np.sqrt(np.vdot(h0, h0).real * np.vdot(h, h).real)
    return 1.0 - ov


def sample_bbh(rng):
    q = rng.uniform(1.0, 8.0)
    return dict(
        M=50.0, q=q,
        chi1=rng.uniform(-0.8, 0.9), chi2=rng.uniform(-0.8, 0.9),
        distance=400.0, inclination=rng.uniform(0, np.pi), coalescence_angle=0.0,
        use_geometric_units="no", initial_frequency=20.0,
        domain=0, interp_uniform_grid="yes", srate_interp=4096.0,
        ecc=rng.uniform(0.0, 0.5), anomaly=rng.uniform(0.0, 2 * np.pi), ecc_freq=2,
        arg_out="no",
    )


def sample_bns(rng):
    lam = rng.uniform(100.0, 1500.0)
    return dict(
        M=2.7, q=rng.uniform(1.0, 1.8),
        chi1=rng.uniform(-0.1, 0.1), chi2=rng.uniform(-0.1, 0.1),
        distance=100.0, inclination=rng.uniform(0, np.pi), coalescence_angle=0.0,
        use_geometric_units="no", initial_frequency=30.0,
        domain=0, interp_uniform_grid="yes", srate_interp=8192.0,
        ecc=rng.uniform(0.0, 0.2), anomaly=rng.uniform(0.0, 2 * np.pi), ecc_freq=2,
        use_tidal="TEOBRESUM", use_lambda234_fits="YAGI13",
        LambdaAl2=lam, LambdaBl2=lam, arg_out="no",
    )


def run(sampler, n, candidates, rng, label):
    results = {c: [] for c in candidates}
    fails = {c: 0 for c in candidates}
    base_fail = 0
    worst = {c: (0.0, None) for c in candidates}
    for i in range(n):
        pars = sampler(rng)
        try:
            h0 = h22(pars, "rkf45")
        except Exception:
            base_fail += 1
            continue
        if not np.all(np.isfinite(h0)) or len(h0) < 8:
            base_fail += 1
            continue
        for c in candidates:
            try:
                h = h22(pars, c)
                mm = mismatch(h0, h)
                if not np.isfinite(mm):
                    raise ValueError("non-finite mismatch")
            except Exception:
                fails[c] += 1
                continue
            results[c].append(mm)
            if mm > worst[c][0]:
                key = {k: round(pars[k], 3) for k in ("q", "chi1", "chi2", "ecc", "anomaly")}
                worst[c] = (mm, key)
    print(f"\n=== {label}: {n} random configs ({n - base_fail} baseline-valid) ===")
    if base_fail:
        print(f"  (rkf45 baseline failed/invalid on {base_fail} configs, skipped)")
    print(f"{'stepper':<9}{'n_ok':>6}{'n_fail':>7}{'median':>11}{'p90':>11}{'max':>11}")
    for c in candidates:
        a = np.array(results[c])
        if len(a) == 0:
            print(f"{c:<9}{0:>6}{fails[c]:>7}  (all failed)")
            continue
        print(f"{c:<9}{len(a):>6}{fails[c]:>7}{np.median(a):>11.2e}"
              f"{np.percentile(a, 90):>11.2e}{a.max():>11.2e}")
    for c in candidates:
        if worst[c][1]:
            print(f"  worst {c}: mismatch={worst[c][0]:.2e} at {worst[c][1]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--nbns", type=int, default=0)
    ap.add_argument("--candidates", default="rk8pd,msadams")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    cands = args.candidates.split(",")
    rng = np.random.default_rng(args.seed)
    print(f"stepper accuracy sweep (mismatch vs rkf45, flat PSD, time+phase maximized)")
    run(sample_bbh, args.n, cands, rng, "ecc BBH")
    if args.nbns:
        run(sample_bns, args.nbns, cands, rng, "ecc BNS")


if __name__ == "__main__":
    main()
