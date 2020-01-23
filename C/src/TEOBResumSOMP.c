/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 See AUTHORS file
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

#include "TEOBResumS.h"

#ifdef _OPENMP

void openmp_pr_info()
{
  int id = omp_get_thread_num();
  int nt = omp_get_num_threads();
  int nm = omp_get_max_threads();
  printf("OMP: Thread %d of %d (max threads %d)\n", id,nt,nm);
}

/* Rough timers system */
int timeron = 0; /* screen timers? */
struct timer
{
  double ti, te, tt;
  char name[STRLEN];
  int n;
};
#define NTIMERS (42)
struct timer timers[NTIMERS];
int ntimer = 0;

void openmp_timer_start(char *name)
{
  if (!USETIMERS) return;
  int n;
  for (n=0; n<ntimer; n++)
    if (STREQUAL(name, timers[n].name)) break;
  if (n==ntimer) {
    if (n>=NTIMERS) errorexit("Not enough memory for timers (increase NTIMERS)");
    ntimer++;
    strcpy (timers[n].name,name);
  }
  timers[n].n++;
  if (timeron) printf("> [%s] (%06d) time-mark\n",timers[n].name, timers[n].n);
  timers[n].ti = omp_get_wtime();
}

void openmp_timer_stop(char *name)
{
  if (!USETIMERS) return;
  double te = omp_get_wtime();
  int n;
  for (n=0; n<ntimer; n++)
    if (STREQUAL(name, timers[n].name)) break;
  timers[n].te = te - timers[n].ti;
  timers[n].tt += timers[n].te;
  if (timeron) printf("< [%s] (%06d) time-elapsed = %e\n",timers[n].name,timers[n].n,timers[n].te);
}

void openmp_timer_output()
{
  if (!USETIMERS) return;
  FILE *fp;
  char fname[STRLEN];
  strcpy(fname,EOBPars->output_dir);
  strcat(fname,"/openmp_timers.txt");
  fp = fopen(fname,"w");
  if (!fp) errorexit("failed to open file");
  fprintf(fp,"# Index Name Avg-time Tot-time\n");
  for (int n=0; n<ntimer; n++) 
    fprintf(fp,"%02d %-20s %.4e %.4e\n", n, timers[n].name,
	    timers[n].tt/timers[n].n, timers[n].tt); 
  fclose(fp); 
}

void openmp_init(const int verbose) 
{
  const int nt  = EOBPars->openmp_threads;
  if (nt==0) {
    if (verbose) printf("OMP: Use automatic settings\n");
  } else {
    omp_set_dynamic(0); // disable dynamic teams  
    if (nt>0) {
      if (verbose) printf("OMP: Set num_threads = %d\n", nt);
      omp_set_num_threads(nt);
    } else { 
      if (verbose) printf("OMP: Set num_threads = 1\n");
      omp_set_num_threads(1);   
    }
  }
  if (USETIMERS) {
    timeron = EOBPars->openmp_timeron;
    openmp_timer_start("main");
  }
}

void openmp_free() 
{
  if (USETIMERS) {
    openmp_timer_stop("main");
    openmp_timer_output();
  }
}

#endif
