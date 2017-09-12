/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include "TEOBResumS.h"

int TEOBResumS(
	       Waveform **hpp,       /** (h+,hx) return array  **/
	       Waveform_lm **hlm      /** hlm return array    **/
	       )
{

  int size = par_get_i("size");

  double q = par_get_d("q");
  double chi1 = par_get_d("chi1");
  double chi2 = par_get_d("chi2");

  /** Alloc memory for dynamical vars */
  Dynamics *dyn;
  Dynamics_alloc(&dyn, size);
      
  /** Set parameters */
  dyn->dt            = par_get_d("dt");
  dyn->t1            = par_get_d("ode_tmax");
  dyn->MOmg_prev     = 0.;
  dyn->t_stop        = 0.;
  dyn->Omg           = 0.;
  dyn->Omg_orb       = 0.;
  dyn->A             = 0.;
  dyn->ddotr         = 0.;
  dyn->stop_flag     = false;
  dyn->MOmgpeak_flag = false;
  // ...

  /* Compute light-ring (if needed) */
  if (par_get_i("use_tidal")) {
    dyn->rLR = AdiabLR(&params);
    dyn->rLSO = 6.0; 
    par_set_d("rLR", dyn->rLR);
  }
    
  /** Computing the initial conditions */
  gsl_odeiv2_system sys = {rhs, NULL , EOB_EVOLVE_VARS, &params};
  if (par_get_i("use_spins")) {
    sys = {s_RHS, NULL, EOB_EVOLVE_VARS, &params};
    s_initial(dyn->y0, &params);
  } else {
    sys     = {rhs, NULL, EOB_EVOLVE_VARS, &params};
    initial(dyn->y0, &params);
  }
    
  /** Initial conditions: t, r, phi, prstar, pphi */
  dyn->t = 0.0;
  dyn->y[EOB_EVOLVE_RAD]    = dyn->y0[EOB_ID_RAD];
  dyn->y[EOB_EVOLVE_PHI]    = 0.;
  dyn->y[EOB_EVOLVE_PRSTAR] = dyn->y0[EOB_ID_PRSTAR];
  dyn->y[EOB_EVOLVE_PPH]    = dyn->y0[EOB_ID_PPH];

  /** Final BH */
  if (!(par_get_i("use_tidal"))) {
    double BH_final_mass = HealyBBHFitRemnant(chi1, chi2, q);
    par_set_d("BH_final_mass", BH_final_mass);
  }
    
  /** Initialize ODE system solver */
  const double ode_abstol = par_get_d("ode_abstol");
  const double ode_reltol = par_get_d("ode_relstol");

  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, dyn->dt, ode_abstol, ode_reltol);    
  gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, EOB_EVOLVE_NVARS);
  gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (ode_abstol, ode_reltol);
  gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (EOB_EVOLVE_NVARS);
    
  /** Sove ODE */
  const int solver_scheme = par_get_i("solver_scheme");
  int STATUS = OK;
  int iter = 0;
  int k;
  while (dyn->stop_flag == false) {
    
    if (solver_scheme == 0) {
      /* Uniform timestepping */
      dyn->ti = dyn->t + dyn->dy;
      STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;
      }
    } else if (solver_scheme == 1) {
      /* Adaptive timestepping until LSO ... */
      if (y[EOB_EVOLVE_RAD]>dyn->rLSO) {
	STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t1, &dyn->dt, dyn->y);
	if (STATUS != GSL_SUCCESS) {
	  printf ("ODE solver failed. Error = %d\n", STATUS);
	  return STATUS;
	}
      } else {
	/* ... uniform afterwards */
	dyn->ti = dyn->t + dyn->dt;
	STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
	if (STATUS != GSL_SUCCESS) {
	  printf ("ODE solver failed. Error = %d\n", STATUS);
	  return STATUS;
	}
      }
    } else {
      /* Adaptive timestepping */
      STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t1, &dyn->dt, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;	
      }
    }
    
    /** Unpack data */
    dyn->r      = y[EOB_EVOLVE_RAD];
    dyn->phi    = y[EOB_EVOLVE_PHI];
    dyn->prstar = y[EOB_EVOLVE_PRSTAR];
    dyn->pphi   = y[EOB_EVOLVE_PPHI];
        
    /** Checking whether the dynamics produces NaN values
	this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
	printf ("ODE solver return NaN radius.\n");
	return ERROR;      
    }

    /** Waveform computation*/
    // fixme: call and routine:
    s_waveform(dyn->t, dyn->y, &params, dyn->Omg, dyn->Omg_orb, dyn->A, dyn->ddotr,
	       h_form);

    /** Update size and push arrays (if needed) */
    iter++;
    if (iter>size) {

      size = iter;
      par_set_i("size", size);
      
      Waveform_push (&hpp, size);
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);

    }
    
    /** Append dynamics and waveform to vectors */
    dyn->time[iter] = t;
    dyn->data[EOB_RAD][iter]  = dyn->r;
    dyn->data[EOB_MOMG][iter] = dyn->Omg;
    dyn->data[EOB_PPHI][iter] = dyn->pphi;
    dyn->data[EOB_PRSTAR] = dyn->prstar;
    dyn->data[EOB_OMGORB] = dyn->Omg_orb;
    dyn->data[EOB_DDOTR] = dyn->ddotr;

    //hlm_rad_vec.push_back(h_form[lm].dat[0]);
    //hlm_phase_vec.push_back(h_form[lm].dat[1]);
    for (k = 0; k < KMAX; k++) {
      hlm->real[k][iter] =0,; // amplitude
      hlm->imag[k][iter] = 0.; //phase
    }
    
    /** Check when to break the computation
	find peak of omega curve and continue for delta_t=10. afterwards */
    //MOmg = Omg; //NOTE: was MOmg = Omg_orb; before!!! (only for the spinning case)
    if (params.flags.spin==1) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }
    
    if (dyn->MOmgpeak_flag==false) {
      if (dyn->MOmg < dyn->MOmg_prev) {	  
	dyn->MOmgpeak_flag = true;
	dyn->t_stop     = dyn->t + 4.*dt;
      } else {
	dyn->MOmg_prev = dyn->MOmg;
      }
    } else {
      if (t >= t_stop) {
	stop_flag = true;
      }
    }
  
  } /* end time iteration */
  
  /** Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);
  


  // SB stops here waiting for improved NQC & Ringdown ...............
  //
  // ------------------------------












    /** To compute the NQC corrections, a precise determination of the time at which Omega_orb peaks is required,
	in order to construct a grid which passes from that t_peak.
	The first step is to find the t_peak on the grid. */
    
    double t_max_grid = t_vec[0];
    double omg_max_grid = 0.;
    int index_max = 0.;
    
    /* Find max and index */
    for (int index=0; index<Omg_orb_vec.size(); index++)
      {
	if(Omg_orb_vec[index] > omg_max_grid)
	  {
	    index_max    = index;
	    omg_max_grid = Omg_orb_vec[index_max];
	    t_max_grid   = t_vec[index_max];
	  }
      }
    cout << "t_max_grid:\t" <<  t_max_grid << "\nomg_max_grid:\t" <<  omg_max_grid << endl;
    
    /* Then take a few points around the peak and analytically interpolate these points
       to find a better approximation for t_peak. */
    double x1 = t_vec[index_max - 1];
    double x2 = t_vec[index_max    ];
    double x3 = t_vec[index_max + 1];
    double y1 = Omg_orb_vec[index_max - 1];
    double y2 = Omg_orb_vec[index_max    ];
    double y3 = Omg_orb_vec[index_max + 1];
    
    double c1 = (pow(x3,2)*(y1 - y2) + pow(x1,2)*(y2 - y3) +
		 pow(x2,2)*(-y1 + y3))/((x1 - x2)*(x1 - x3)*(x2 - x3));
    double c2 = (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))/
      ((x1 - x2)*(x1 - x3)*(x2 - x3));
    
    double t_max = (-c1)/(2.*c2);
    double omega_max = interp1d (3, t_max, 3, &Omg_orb_vec[index_max - 1], &t_vec[index_max - 1]);
    cout << "t_max_interp:\t" << t_max << "\nomg_max_interp:\t" <<  omega_max << endl;
    
    /** Interpolate quantities on a new grid of spacing dt and passing by t_max */

    // HERE THERE ARE PROBLEMS:

    int grid_length = (int)((t_vec.back()-t_vec[0])/dt + 1);        
    int grid_length2 = t_vec.size(); // why using this gives SEG FAULT !?!? CHECKME
    printf(" %d %d\n",grid_length2,grid_length ); // they are different ?
    
    int N_before = int((t_max - t_vec[0])/dt +1);
    vector<double> t_vecg(grid_length);
    for (int k = N_before; k >= 0 ; k--)
     {
         t_vecg[k] = t_max - (N_before-k)*dt;
     }
     for (int k = N_before+1; k < grid_length; k++)
    {
         t_vecg[k] = t_max + (k-N_before)*dt;
    }


     // print out
     cout << "t_max_interp:\t" << t_max << "\tt_max_grid:\t" << t_vecg[N_before] << endl;
     for(int k = 0; k < grid_length; k++)
     {
     cout << "i:\t" << k << "\tt_vecg[i]:\t" << t_vecg[k] << endl;
     }

     // seg fault right after here
     //exit(1);
    printf("***********************************");

    // alloc memory
    vector<double> r_vecg         = r_vec;
    vector<double> MOmg_vecg      = MOmg_vec;
    vector<double> pph_vecg       = pph_vec;
    vector<double> prstar_vecg    = prstar_vec;
    vector<double> ddotr_vecg     = ddotr_vec;
    vector<double> OmgOrb_vecg    = Omg_orb_vec;
    
    vector<double> hlm_phase_vecg = hlm_phase_vec;
    vector<double> hlm_rad_vecg   = hlm_rad_vec;
    std::vector<vector<double> > hlm_ampl_g(35);
    std::vector<vector<double> > hlm_phase_g(35);
    for (int m=35; m--; )
      {
	hlm_ampl_g[m] = hlm_phase_vec;
	hlm_phase_g[m] = hlm_phase_vec;
      }

    // 4th order interp all fields on tg
    const int order = 4;
    double* tt  = &t_vec[0];
    double* tg = &t_vecg[0];
    double wg[4];
    double* datap; // this is because some idiots should live in one of the other 11D
    
    printf("***********************************");

    for(int k = 0; k < grid_length; k++)
      {
	// NN & weights
	int ix = find_point_bisection(tg[k], grid_length, tt, order/2);	
	baryc_weights(4, &tt[ix], wg);

	// dynamics
	if (params.flags.tidal == 0)
	  {
	    datap = &r_vec[ix];       r_vecg[k]      = baryc_f_weights(tg[k], order, datap, &tt[ix], wg); 
	    datap = &MOmg_vec[ix];    MOmg_vecg[k]   = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);  
	    datap = &pph_vec[ix];     pph_vecg[k]    = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);  
	    datap = &prstar_vec[ix];  prstar_vecg[k] = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);  
	    datap = &ddotr_vec[ix];   ddotr_vecg[k]  = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);  
	    datap = &Omg_orb_vec[ix]; OmgOrb_vecg[k] = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);  
	  }
	
	// wave 
	datap = &hlm_phase_vec[ix]; hlm_phase_vecg[k] = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);
	datap = &hlm_rad_vec[ix];   hlm_rad_vecg[k]   = baryc_f_weights(tg[k], order, datap, &tt[ix], wg);
	for (int m=35; m--; )
	  {
	    double *amplitude = &hlm_ampl[m][0];
	    double *phase     = &hlm_phase[m][0];
	    hlm_ampl_g[m][k]   = baryc_f_weights(tg[k], order, &amplitude[ix], &tt[ix], wg);
	    hlm_phase_g[m][k]  = baryc_f_weights(tg[k], order, &phase[ix]    , &tt[ix], wg);
	  }
      }
    
    /** NQCs corrections */
    /** NOTE THAT IF YOU REMOVE PARAMS.SPIN==TRUE EVERYTHING IS FUCKED UP FOR SOME REASON */
    
#if (DEBUG)
    std::FILE* waveform_preNQC = std::fopen("waveform_preNQC.dat", "w");
    for (int j=0;j<hlm_ampl_g[1].size();j++)
      std::fprintf(waveform_preNQC, "%f\t%e\t%e\n", (double)j*dt, hlm_ampl_g[1][j], hlm_phase_g[1][j]);
    std::fclose(waveform_preNQC);
    
    std::FILE* waveform_NQC = std::fopen("waveform_nqc.dat", "w");
#endif
    
    if (params.flags.tidal==0 && params.flags.spin==1)
      {
	
	vector<vector<gsl_complex> > nqc = find_a1a2a3(t_vecg,r_vecg,MOmg_vecg,pph_vecg,prstar_vecg,hlm_phase_g,OmgOrb_vecg,hlm_ampl_g,ddotr_vecg,&params);
	
        for (int k=35; k--; )
	  {
	    if (k==1)
	    {
	      for (int i=0; i<grid_length; i++ )
		{
		  hlm_ampl_g[k][i]  = hlm_ampl_g[k][i]  * nqc[k][i].dat[0];
		  hlm_phase_g[k][i] = hlm_phase_g[k][i] + nqc[k][i].dat[1];
#if (DEBUG)
		  std::fprintf(waveform_NQC, "%f\t%e\t%e\n", t_vecg[i], nqc[k][i].dat[0], nqc[k][i].dat[1]);
#endif	
		}
	    }
	  }
	
      }

#if (DEBUG)
      std::fclose(waveform_NQC);
      std::FILE* waveform_postNQC   = std::fopen("waveform_postNQC.dat", "w");
      for (int j=0;j<hlm_ampl_g[1].size();j++)
	  std::fprintf(waveform_postNQC, "%f\t%e\t%e\n", (double)j*dt, hlm_ampl_g[1][j], hlm_phase_g[1][j]);
      std::fclose(waveform_postNQC);
#endif

      /** Define a time vector for each multipole
	  These will be cut by the ringdown, where
	  each multipole has its own starting time */
      
    /** Ringdown attachment */
    if (params.flags.tidal==0)
    {

      vector<vector<double> > t_g(35);
      for (int k=35; k--; )
	{
	  t_g[k] = t_vecg;
	}
      
      ringdown(params,t_g,OmgOrb_vecg,hlm_ampl_g,hlm_phase_g);
    }
















    // h+ hx computation


    /** All multipoles will now have size N+Nringdown */
    /** Multipole for which no ringdown model is available will be filled with 0s */
    /** We pick the index 1 since it is the 22 mode and it is always computed */
    
    int N = hlm_ampl_g[1].size();

    /** Allocate hplus and hcross */
    Waveform *hplus_out = (Waveform *)malloc(sizeof(Waveform));
    if (hplus_out == NULL)
    {
        printf("ERROR allocating hplus.\n");
        exit(-1);
    }
    hplus_out->data = (double *)malloc(N*sizeof(double));
    if (hplus_out->data == NULL)
    {
        printf("ERROR allocating hplus->data.\n");
        exit(-1);
    }
    hplus_out->length = N;
    Waveform *hcross_out = (Waveform *)malloc(sizeof(Waveform));
    if (hcross_out == NULL)
    {
        printf("ERROR allocating hcross.\n");
        exit(-1);
    }
    hcross_out->data = (double *)malloc(N*sizeof(double));
    if (hcross_out->data == NULL)
    {
        printf("ERROR allocating hcross->data.\n");
        exit(-1);
    }
    hcross_out->length = N;
    
    /** Set them to zero initially */
    
    memset(hplus_out->data, 0, N*sizeof(double));
    memset(hcross_out->data, 0, N*sizeof(double));
    
    /** Spherical harmonics projection **/
    /** construct hplus and hcross **/
    /** h22 = 1/R * (nu*M)*G/c^2 h_code_output */

    double mtot_m = 1.;
    double amplitude_prefactor = 1.;    
    if (!(params.flags.geometric_units))
    {
      mtot_m = (m1+m2)*MSUN_M;
      amplitude_prefactor = params.nu*mtot_m/(distance*MPC_M);
    } 

    if (params.flags.multipoles == 1)
    {
        for (i=0; i<N; i++)
        {
            hplus_out->data[i]   = hlm_ampl_g[lm][i];
            hcross_out->data[i]  = hlm_phase_g[lm][i];
        }

    }
    else
    {
        for (int k=35; k--; )
        {
            if (k==1)
            {
                double Y_real, Y_imag;
                spinsphericalharm(&Y_real, &Y_imag, -2, L[k], M[k], polarisation, inclination);
                
                /** there is a MINUS SIGN in the phase h = A exp(-i phase) **/
                for (i=0; i<N; i++)
                {
                    double Aki = hlm_ampl_g[k][i]*amplitude_prefactor;
                    double cosPhi = cos(hlm_phase_g[k][i]);
                    double sinPhi = -sin(hlm_phase_g[k][i]);
                    hplus_out->data[i] += Aki*(cosPhi*Y_real - sinPhi*Y_imag);
                    hcross_out->data[i] -= Aki*(cosPhi*Y_imag + sinPhi*Y_real);
                }
            }
        }
    }
    *hplus = hplus_out;
    *hcross= hcross_out;











    /** Free memory for dynamical vars */
    Dynamics_free(dyn);

    return OK;
}
