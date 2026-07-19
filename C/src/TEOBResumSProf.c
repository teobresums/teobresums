/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2026 See AUTHORS file
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 */

/** \file TEOBResumSProf.c
 *  \brief Lightweight portable profiling: section timers + RHS-call counters.
 *
 * Enabled with -DPROFILING=1 (see the `prof` Makefile target). When disabled
 * the PROF_* macros in TEOBResumS.h expand to nothing, so the timed code is
 * bit-identical to the un-instrumented build. Timers are indexed by the PROF_*
 * enum for O(1), strcmp-free start/stop on the hot RHS/flux paths.
 *
 * Note: fine-grained timers around the r.h.s. fire ~1e5-1e6 times per waveform,
 * so their clock_gettime overhead perturbs the very quantity they measure. The
 * RHS-call counters (prof_rhs_calls, prof_rhs_store_calls) are the robust metric
 * for comparing ODE steppers; treat the RHS/flux wall times as indicative only.
 */

#include "TEOBResumS.h"

/* clock_gettime requires _POSIX_C_SOURCE >= 199309L; TEOBResumS.h defines
 * _GNU_SOURCE (which implies it), and the `prof` target passes the flag too. */

long prof_rhs_calls       = 0;
long prof_rhs_store_calls = 0;

#if (PROFILING)

typedef struct {
  double t_start;    /* wall time of the last start() */
  double t_total;    /* accumulated wall time */
  long   n;          /* number of start/stop pairs */
  const char *name;
} prof_timer_t;

static prof_timer_t prof_timers[PROF_NTIMERS];

static double prof_now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

void prof_init(void)
{
  static const char *names[PROF_NTIMERS] = {
    "ODE_loop", "RHS", "flux", "Fphi_ecc",
    "metric+rc", "wav_store", "NQC", "interp_spline", "compute_hpc"
  };
  for (int i = 0; i < PROF_NTIMERS; i++) {
    prof_timers[i].t_start = 0.0;
    prof_timers[i].t_total = 0.0;
    prof_timers[i].n       = 0;
    prof_timers[i].name    = names[i];
  }
  prof_rhs_calls       = 0;
  prof_rhs_store_calls = 0;
}

void prof_timer_start(int i)
{
  prof_timers[i].t_start = prof_now();
  prof_timers[i].n++;
}

void prof_timer_stop(int i)
{
  prof_timers[i].t_total += prof_now() - prof_timers[i].t_start;
}

void prof_timer_output(void)
{
  char fname[STRLEN];
  strcpy(fname, EOBPars->output_dir);
  strcat(fname, "/prof_timers.txt");
  FILE *fp = fopen(fname, "w");
  if (!fp) {
    /* output_dir may not exist when no file output was requested; fall back
     * to stdout only rather than aborting a profiling run. */
    fprintf(stderr, "[prof] could not open %s, printing to stdout only\n", fname);
  } else {
    fprintf(fp, "# Index Name Calls Avg-time Tot-time\n");
    for (int i = 0; i < PROF_NTIMERS; i++)
      fprintf(fp, "%02d %-14s %10ld %.6e %.6e\n", i, prof_timers[i].name,
              prof_timers[i].n,
              prof_timers[i].n ? prof_timers[i].t_total / prof_timers[i].n : 0.0,
              prof_timers[i].t_total);
    fprintf(fp, "# rhs_calls_total       %ld\n", prof_rhs_calls);
    fprintf(fp, "# rhs_calls_storage     %ld\n", prof_rhs_store_calls);
    fprintf(fp, "# rhs_calls_stepper     %ld\n", prof_rhs_calls - prof_rhs_store_calls);
    fclose(fp);
  }

  /* Also echo to stdout so Python-driven runs can scrape without a file. */
  printf("# PROF Index Name Calls Avg-time Tot-time\n");
  for (int i = 0; i < PROF_NTIMERS; i++)
    printf("# PROF %02d %-14s %10ld %.6e %.6e\n", i, prof_timers[i].name,
           prof_timers[i].n,
           prof_timers[i].n ? prof_timers[i].t_total / prof_timers[i].n : 0.0,
           prof_timers[i].t_total);
  printf("# PROF rhs_calls_total   %ld\n", prof_rhs_calls);
  printf("# PROF rhs_calls_storage %ld\n", prof_rhs_store_calls);
  printf("# PROF rhs_calls_stepper %ld\n", prof_rhs_calls - prof_rhs_store_calls);
}

#else

/* Un-instrumented build: symbols exist but are never called (PROF_* macros are
 * empty), so they carry no runtime cost. */
void prof_init(void)         {}
void prof_timer_start(int i) { (void)i; }
void prof_timer_stop(int i)  { (void)i; }
void prof_timer_output(void) {}

#endif /* PROFILING */
