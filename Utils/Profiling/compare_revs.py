#!/usr/bin/env python3
"""
Before/after regression report across two git revisions.

For each benchmark parfile it builds both revisions (opt), measures wall time
(min of N, file output OFF), and generates the (2,2) mode (file output ON) to
compute the time-and-phase-maximized mismatch of the NEW output vs the BASE.
The result is one table:

  config | t_base | t_new | speedup | mismatch(new vs base) | stepper_new

This both documents the speedup and re-proves that the output change stays
within the accuracy gate. Bit-identical changes show mismatch ~0; the stepper
default change shows the (small) stepper mismatch vs the old rkf45 output.

Usage:
  python3 compare_revs.py --base REV --new REV
  python3 compare_revs.py --reps 4

Note: base REV must have the `opt` Makefile target (i.e. be at or after the
profiling-infrastructure commit); otherwise pass an -O2-capable rev.

SAFETY: temporary git worktrees are created under the system temp dir and are
removed ONLY via `git worktree remove` (which refuses to touch the main
worktree). This script never calls rm/rmtree on a path itself.
"""

import argparse
import os
import re
import subprocess
import tempfile
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PARDIR = os.path.join(HERE, "parfiles")
TMPROOT = tempfile.gettempdir()


def sh(cmd, cwd=None, check=True):
    return subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)


def output_dir_of(parfile):
    m = re.search(r'^output_dir\s*=\s*"\./([^"]+?)/?"', open(parfile).read(), re.M)
    return m.group(1).rstrip("/") if m else None


def build_worktree(rev):
    wt = tempfile.mkdtemp(prefix=f"teob_cmp_{rev[:8]}_")
    sh(["git", "-C", REPO, "worktree", "add", "-q", "--detach", wt, rev])
    sh(["make", "opt"], cwd=os.path.join(wt, "C"))
    return wt


def cleanup_worktree(wt):
    # Git-managed removal only; refuses to remove the main worktree, and we
    # additionally require the path to live under the temp dir with our prefix.
    if wt and wt.startswith(TMPROOT) and "teob_cmp_" in os.path.basename(wt):
        sh(["git", "-C", REPO, "worktree", "remove", "--force", wt], check=False)


def run_parfile(wt, parfile, overrides):
    cdir = os.path.join(wt, "C")
    tmp = os.path.join(cdir, "_cmp.par")
    txt = open(parfile).read()
    for k, v in overrides.items():
        txt = re.sub(rf'^{k}\s*=.*$', f'{k} = {v}', txt, flags=re.M)
        if not re.search(rf'^{k}\s*=', txt, re.M):
            txt += f'\n{k} = {v}\n'
    open(tmp, "w").write(txt)
    odir = output_dir_of(tmp)
    subprocess.run(["./TEOBResumS.x", "-p", "_cmp.par"], cwd=cdir, capture_output=True)
    return os.path.join(cdir, odir) if odir else None


def timeit(wt, parfile, reps):
    off = {"output_hpc": '"no"', "output_multipoles": '"no"', "output_dynamics": '"no"'}
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        run_parfile(wt, parfile, off)
        best = min(best, time.perf_counter() - t0)
    return best


def load_h22(odir):
    f = os.path.join(odir, "hlm_interp_l2_m2.txt")
    if not odir or not os.path.exists(f):
        return None
    d = np.loadtxt(f)
    return d[:, 1] * np.exp(-1j * d[:, 2])


def mismatch(h0, h):
    if h0 is None or h is None:
        return float("nan")
    n = max(len(h0), len(h)); N = 1 << int(np.ceil(np.log2(2 * n)))
    z = np.fft.ifft(np.fft.fft(h0, N) * np.conj(np.fft.fft(h, N)))
    return 1 - np.max(np.abs(z)) / np.sqrt(np.vdot(h0, h0).real * np.vdot(h, h).real)


def stepper_of(odir):
    p = os.path.join(odir, "params.txt") if odir else None
    if p and os.path.exists(p):
        m = re.search(r'^ode_stepper\s*=\s*"([^"]+)"', open(p).read(), re.M)
        if m:
            return m.group(1)
    return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--new", default="HEAD")
    ap.add_argument("--reps", type=int, default=3)
    args = ap.parse_args()
    base_sha = sh(["git", "-C", REPO, "rev-parse", "--short", args.base]).stdout.strip()
    new_sha = sh(["git", "-C", REPO, "rev-parse", "--short", args.new]).stdout.strip()

    print(f"Regression report: base={base_sha}  new={new_sha}  (reps={args.reps})")
    wt_b = wt_n = None
    on = {"output_multipoles": '"yes"'}
    try:
        wt_b = build_worktree(args.base)
        wt_n = build_worktree(args.new)
        print(f"\n{'config':<26}{'t_base':>9}{'t_new':>9}{'speedup':>9}{'mismatch':>12}  stepper_new")
        print("-" * 82)
        for pf in sorted(os.listdir(PARDIR)):
            if not pf.endswith(".par"):
                continue
            parfile = os.path.join(PARDIR, pf)
            name = pf[:-4]
            tb = timeit(wt_b, parfile, args.reps)
            tn = timeit(wt_n, parfile, args.reps)
            od_b = run_parfile(wt_b, parfile, on)
            od_n = run_parfile(wt_n, parfile, on)
            mm = mismatch(load_h22(od_b), load_h22(od_n))
            print(f"{name:<26}{tb:>9.3f}{tn:>9.3f}{tb/tn:>8.2f}x{mm:>12.2e}  {stepper_of(od_n)}")
    finally:
        cleanup_worktree(wt_b)
        cleanup_worktree(wt_n)


if __name__ == "__main__":
    main()
