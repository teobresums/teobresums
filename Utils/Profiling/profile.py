#!/usr/bin/env python3
"""
Profiling driver for TEOBResumS eccentric/BNS performance work.

Builds the C executable with the `prof` target (-O2 -pg -DPROFILING=1), runs a
suite of representative eccentric-BBH / eccentric-BNS / quasi-circular configs,
and collects for each:

  * wall time
  * the PROF section timers + RHS-call counters (from prof_timers.txt / stdout)
  * the gprof flat profile (top functions), if gprof is available

Outputs a summary table to stdout and, with --csv, a machine-readable CSV.

Usage:
  python3 profile.py                 # build (prof) + run all configs + table
  python3 profile.py --no-build      # skip rebuild, reuse existing TEOBResumS.x
  python3 profile.py --csv out.csv   # also write a CSV summary
  python3 profile.py --config NAME   # run a single config by basename

Note on numbers: the fine-grained RHS/flux section timers are perturbed by their
own clock_gettime overhead (they fire ~1e5-1e6 times). Trust the RHS-call
*counts* and the gprof flat profile for the true per-function breakdown; treat
the section wall-times as indicative of relative section cost only.
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CDIR = REPO / "C"
EXE = CDIR / "TEOBResumS.x"
PARDIR = HERE / "parfiles"


def build(target="prof"):
    print(f"[profile] make {target} ...", flush=True)
    subprocess.run(["make", target], cwd=CDIR, check=True,
                   stdout=subprocess.DEVNULL)


def output_dir_of(parfile):
    txt = parfile.read_text()
    m = re.search(r'^output_dir\s*=\s*"\./([^"]+)/?"', txt, re.M)
    return CDIR / m.group(1).rstrip("/") if m else None


def parse_prof_stdout(stdout):
    """Return (timers, counts) from the '# PROF ...' lines echoed to stdout."""
    timers, counts = {}, {}
    for line in stdout.splitlines():
        m = re.match(r"# PROF (\d+) (\S+)\s+(\d+)\s+(\S+)\s+(\S+)", line)
        if m:
            timers[m.group(2)] = dict(calls=int(m.group(3)),
                                      avg=float(m.group(4)),
                                      total=float(m.group(5)))
            continue
        m = re.match(r"# PROF (rhs_calls_\w+)\s+(\d+)", line)
        if m:
            counts[m.group(1)] = int(m.group(2))
    return timers, counts


def gprof_top(n=12):
    gmon = CDIR / "gmon.out"
    if not gmon.exists():
        return []
    try:
        out = subprocess.run(["gprof", "-p", "-b", str(EXE), str(gmon)],
                             cwd=CDIR, capture_output=True, text=True).stdout
    except FileNotFoundError:
        return []
    rows = []
    for line in out.splitlines():
        m = re.match(r"\s*([\d.]+)\s+[\d.]+\s+[\d.]+\s+(\d+)?\s*\S*\s*\S*\s*(\w+)\s*$", line)
        if m and m.group(1) != "0.00":
            rows.append((float(m.group(1)), m.group(3)))
        if len(rows) >= n:
            break
    return rows


def run_config(parfile, do_gprof=True):
    odir = output_dir_of(parfile)
    if odir and odir.exists():
        subprocess.run(["rm", "-rf", str(odir)])
    t0 = time.perf_counter()
    proc = subprocess.run([str(EXE), "-p", str(parfile)],
                          cwd=CDIR, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        print(f"  ! {parfile.name} FAILED (rc={proc.returncode})")
        print(proc.stdout[-500:], proc.stderr[-500:])
        return None
    timers, counts = parse_prof_stdout(proc.stdout)
    gp = gprof_top() if do_gprof else []
    return dict(name=parfile.stem, wall=wall, timers=timers,
                counts=counts, gprof=gp, odir=odir)


def print_summary(results):
    print("\n" + "=" * 96)
    print("SUMMARY  (wall = full run; ODE/RHS/flux/Fphi = PROF section timers, seconds)")
    print("=" * 96)
    hdr = f"{'config':<26}{'wall':>8}{'ODE':>8}{'RHS':>8}{'flux':>8}{'Fphi':>8}" \
          f"{'wavst':>8}{'rhs_step':>10}{'rhs_tot':>10}"
    print(hdr)
    print("-" * 96)
    for r in results:
        if r is None:
            continue
        t = r["timers"]
        c = r["counts"]
        def tot(k): return t.get(k, {}).get("total", 0.0)
        print(f"{r['name']:<26}{r['wall']:>8.3f}{tot('ODE_loop'):>8.3f}"
              f"{tot('RHS'):>8.3f}{tot('flux'):>8.3f}{tot('Fphi_ecc'):>8.3f}"
              f"{tot('wav_store'):>8.3f}"
              f"{c.get('rhs_calls_stepper', 0):>10}"
              f"{c.get('rhs_calls_total', 0):>10}")
    print("-" * 96)
    for r in results:
        if r is None or not r["gprof"]:
            continue
        print(f"\n[{r['name']}] gprof top self-time functions:")
        for pct, fn in r["gprof"]:
            print(f"    {pct:5.1f}%  {fn}")


def write_csv(results, path):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["config", "wall_s", "ode_s", "rhs_s", "flux_s",
                    "fphi_ecc_s", "wav_store_s", "rhs_stepper", "rhs_total"])
        for r in results:
            if r is None:
                continue
            t, c = r["timers"], r["counts"]
            def tot(k): return t.get(k, {}).get("total", 0.0)
            w.writerow([r["name"], f"{r['wall']:.4f}", f"{tot('ODE_loop'):.4f}",
                        f"{tot('RHS'):.4f}", f"{tot('flux'):.4f}",
                        f"{tot('Fphi_ecc'):.4f}", f"{tot('wav_store'):.4f}",
                        c.get("rhs_calls_stepper", 0), c.get("rhs_calls_total", 0)])
    print(f"[profile] wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--config", help="run only this config (basename, no .par)")
    ap.add_argument("--csv", help="write CSV summary to this path")
    ap.add_argument("--no-gprof", action="store_true")
    args = ap.parse_args()

    if not args.no_build:
        build("prof")
    elif not EXE.exists():
        sys.exit("TEOBResumS.x not found; run without --no-build first")

    pars = sorted(PARDIR.glob("*.par"))
    if args.config:
        pars = [p for p in pars if p.stem == args.config]
        if not pars:
            sys.exit(f"no config named {args.config}")

    results = []
    for p in pars:
        print(f"[profile] running {p.name} ...", flush=True)
        results.append(run_config(p, do_gprof=not args.no_gprof))

    print_summary(results)
    if args.csv:
        write_csv(results, args.csv)


if __name__ == "__main__":
    main()
