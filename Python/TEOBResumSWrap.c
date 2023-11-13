/* Python wrapper for TEOBResumS
 * See: 
 * https://scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
 * and tests in CWrap/ 
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include "lib/TEOBResumS.h"

/*
 * Utilities
 */

/* Util from python cookbook example:
   Set pointer to 1D Carray from PyArray 
   Assumes PyArray is contiguous in memory  */
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
  /* int n = arrayin->dimensions[0]; */
  return (double *) arrayin->data;
}

int SetOptionalVariables(PyObject* dict){
  /* Options */
  if ( PyDict_GetItemString(dict, "r0") != NULL ) {
    EOBPars->r0 = PyFloat_AsDouble(PyDict_GetItemString(dict, "r0"));
  }
  if ( PyDict_GetItemString(dict, "initial_frequency") != NULL ) { 
    EOBPars->initial_frequency = PyFloat_AsDouble(PyDict_GetItemString(dict, "initial_frequency"));
  }
  if ( PyDict_GetItemString(dict, "use_geometric_units") != NULL ) { 
    EOBPars->use_geometric_units = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_geometric_units")));
  }
  if ( PyDict_GetItemString(dict, "use_spins") != NULL ) { 
    EOBPars->use_spins = (int) PyLong_AsLong(PyDict_GetItemString(dict, "use_spins"));
  }
  if ( PyDict_GetItemString(dict, "use_tidal") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_tidal"));
    for(EOBPars->use_tidal=0; EOBPars->use_tidal<=TIDES_NOPT; EOBPars->use_tidal++){
      if (EOBPars->use_tidal == TIDES_NOPT) EOBPars->use_tidal = TIDES_OFF;
      if (STREQUAL(val,tides_opt[EOBPars->use_tidal])) break;
    }
  }
  if ( PyDict_GetItemString(dict, "use_tidal_gravitomagnetic") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_tidal_gravitomagnetic"));
    for(EOBPars->use_tidal_gravitomagnetic=0; EOBPars->use_tidal_gravitomagnetic<=TIDES_NOPT; EOBPars->use_tidal_gravitomagnetic++){
      if (EOBPars->use_tidal_gravitomagnetic == TIDES_GM_NOPT) EOBPars->use_tidal_gravitomagnetic = TIDES_GM_OFF;
      if (STREQUAL(val,tides_gravitomagnetic_opt[EOBPars->use_tidal_gravitomagnetic])) break;
    }  
  }
  if ( PyDict_GetItemString(dict, "use_lambda234_fits") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_lambda234_fits"));
    for(EOBPars->use_lambda234_fits=0; EOBPars->use_lambda234_fits<=Lambda234_fits_NOPT; EOBPars->use_lambda234_fits++){
      if (EOBPars->use_lambda234_fits == Lambda234_fits_NOPT) EOBPars->use_lambda234_fits = Lambda234_fits_YAGI13;
      if (STREQUAL(val,use_lambda234_fits_opt[EOBPars->use_lambda234_fits])) break;
    }    
  }
  if ( PyDict_GetItemString(dict, "use_a6c_fits") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_a6c_fits"));
    for(EOBPars->use_a6c_fits=0; EOBPars->use_a6c_fits<=a6c_fits_NOPT; EOBPars->use_a6c_fits++){
      if (EOBPars->use_a6c_fits == a6c_fits_NOPT) EOBPars->use_a6c_fits = a6c_fits_HM;
      if (STREQUAL(val,use_a6c_fits_opt[EOBPars->use_a6c_fits])) break;
    }    
  }
  if ( PyDict_GetItemString(dict, "use_cN3LO_fits") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_cN3LO_fits"));
    for(EOBPars->use_cN3LO_fits=0; EOBPars->use_cN3LO_fits<=cN3LO_fits_NOPT; EOBPars->use_cN3LO_fits++){
      if (EOBPars->use_cN3LO_fits == cN3LO_fits_NOPT) EOBPars->use_cN3LO_fits = cN3LO_fits_HM_2023_432;
      if (STREQUAL(val,use_cN3LO_fits_opt[EOBPars->use_cN3LO_fits])) break;
    }    
  }
  if ( PyDict_GetItemString(dict, "use_tidal_fmode_model") != NULL ) { 
    EOBPars->use_tidal_fmode_model =  YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_tidal_fmode_model")));
  }
  if ( PyDict_GetItemString(dict, "use_speedytail") != NULL ) { 
    EOBPars->use_speedytail =  YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_speedytail")));
  }
  if ( PyDict_GetItemString(dict, "interp_uniform_grid") != NULL ) { 
    EOBPars->interp_uniform_grid =  YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "interp_uniform_grid")));
  }
  if ( PyDict_GetItemString(dict, "pGSF_tidal") != NULL ) { 
    EOBPars->pGSF_tidal = PyFloat_AsDouble(PyDict_GetItemString(dict, "pGSF_tidal"));
  }
  if ( PyDict_GetItemString(dict, "dt_merger_interp") != NULL ) { 
    EOBPars->dt_merger_interp = PyFloat_AsDouble(PyDict_GetItemString(dict, "dt_merger_interp"));
  }
  if ( PyDict_GetItemString(dict, "dt_interp") != NULL ) { 
    EOBPars->dt_interp = PyFloat_AsDouble(PyDict_GetItemString(dict, "dt_interp"));
  }
  if ( PyDict_GetItemString(dict, "dt") != NULL ) { 
    EOBPars->dt = PyFloat_AsDouble(PyDict_GetItemString(dict, "dt"));
  }
  if ( PyDict_GetItemString(dict, "srate_interp") != NULL ) { 
    EOBPars->srate_interp = PyFloat_AsDouble(PyDict_GetItemString(dict, "srate_interp"));
  }
  if ( PyDict_GetItemString(dict, "srate") != NULL ) { 
    EOBPars->srate = PyFloat_AsDouble(PyDict_GetItemString(dict, "srate"));
  }
  if ( PyDict_GetItemString(dict, "size") != NULL ) { 
    EOBPars->size = (int) PyLong_AsLong(PyDict_GetItemString(dict, "size"));
  }
  if ( PyDict_GetItemString(dict, "ringdown_extend_array") != NULL ) { 
    EOBPars->ringdown_extend_array = (int) PyLong_AsLong(PyDict_GetItemString(dict, "ringdown_extend_array"));
  }

  /* Adiabatic tidal ell>2 parameters */
  if ( PyDict_GetItemString(dict, "LambdaAl3") != NULL )
    EOBPars->LambdaAl3 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl3"));
  if ( PyDict_GetItemString(dict, "LambdaBl3") != NULL )
    EOBPars->LambdaBl3 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl3"));
  if ( PyDict_GetItemString(dict, "LambdaAl4") != NULL )
    EOBPars->LambdaAl4 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl4"));
  if ( PyDict_GetItemString(dict, "LambdaBl4") != NULL )
    EOBPars->LambdaBl4 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl4"));
  if ( PyDict_GetItemString(dict, "LambdaAl5") != NULL )
    EOBPars->LambdaAl5 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl5"));
  if ( PyDict_GetItemString(dict, "LambdaBl5") != NULL )
    EOBPars->LambdaBl5 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl5"));
  if ( PyDict_GetItemString(dict, "LambdaAl6") != NULL )
    EOBPars->LambdaAl6 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl6"));
  if ( PyDict_GetItemString(dict, "LambdaBl6") != NULL )
    EOBPars->LambdaBl6 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl6"));
  if ( PyDict_GetItemString(dict, "LambdaAl7") != NULL )
    EOBPars->LambdaAl7 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl7"));
  if ( PyDict_GetItemString(dict, "LambdaBl7") != NULL )
    EOBPars->LambdaBl7 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl7"));
  if ( PyDict_GetItemString(dict, "LambdaAl8") != NULL )
    EOBPars->LambdaAl8 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl8"));
  if ( PyDict_GetItemString(dict, "LambdaBl8") != NULL )
    EOBPars->LambdaBl8 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl8"));

  /* Extrinsic */
  if ( PyDict_GetItemString(dict, "distance") != NULL ) {
    EOBPars->distance = PyFloat_AsDouble(PyDict_GetItemString(dict, "distance"));
  }
  if ( PyDict_GetItemString(dict, "inclination") != NULL ) {
    EOBPars->inclination = PyFloat_AsDouble(PyDict_GetItemString(dict, "inclination"));
  }
  if ( PyDict_GetItemString(dict, "polarization") != NULL ) {
    EOBPars->polarization = PyFloat_AsDouble(PyDict_GetItemString(dict, "polarization"));
  }
  if ( PyDict_GetItemString(dict, "coalescence_angle") != NULL ) {
    EOBPars->coalescence_angle = PyFloat_AsDouble(PyDict_GetItemString(dict, "coalescence_angle"));
  }

  /* Modes */
  if ( PyDict_GetItemString(dict, "use_mode_lm") != NULL ) {
    if (EOBPars->use_mode_lm) free(EOBPars->use_mode_lm);
    PyListObject *tmp = PyDict_GetItemString(dict, "use_mode_lm");
    EOBPars->use_mode_lm_size = PyObject_Length(tmp);
    EOBPars->use_mode_lm = malloc ( EOBPars->use_mode_lm_size * sizeof(int) );
    for (int i = 0; i < EOBPars->use_mode_lm_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->use_mode_lm[i] = (int) PyLong_AsLong(item); 
    }
  }

  if (PyDict_GetItemString(dict, "use_mode_lm_inertial") != NULL) {
    if (EOBPars->use_mode_lm_inertial) free(EOBPars->use_mode_lm_inertial);
    PyListObject *tmp = PyDict_GetItemString(dict, "use_mode_lm_inertial");
    EOBPars->use_mode_lm_inertial_size = PyObject_Length(tmp);
    EOBPars->use_mode_lm_inertial = malloc ( EOBPars->use_mode_lm_inertial_size * sizeof(int) );
    for (int i = 0; i < EOBPars->use_mode_lm_inertial_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->use_mode_lm_inertial[i] = (int) PyLong_AsLong(item); 
    }     
  } else if (PyDict_GetItemString(dict, "use_mode_lm") != NULL) {
    /* use_mode_lm_inertial is not specified, but use_mode_lm is: use the latter
       to set the former as a default behavior
    */
    if (EOBPars->use_mode_lm_inertial) free(EOBPars->use_mode_lm_inertial);
    PyListObject *tmp = PyDict_GetItemString(dict, "use_mode_lm");
    EOBPars->use_mode_lm_inertial_size = PyObject_Length(tmp);
    EOBPars->use_mode_lm_inertial = malloc ( EOBPars->use_mode_lm_inertial_size * sizeof(int) );
    for (int i = 0; i < EOBPars->use_mode_lm_inertial_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->use_mode_lm_inertial[i] = (int) PyLong_AsLong(item); 
    }       
  }

  if ( PyDict_GetItemString(dict, "output_lm") != NULL ) {
    if (EOBPars->output_lm) free(EOBPars->output_lm);
    PyListObject *tmp = PyDict_GetItemString(dict, "output_lm");
    EOBPars->output_lm_size = PyObject_Length(tmp);
    EOBPars->output_lm = malloc ( EOBPars->output_lm_size * sizeof(int) );
    for (int i = 0; i < EOBPars->output_lm_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->output_lm[i] = (int) PyLong_AsLong(item); 
    }
  }

  /* k postpeak */
  if ( PyDict_GetItemString(dict, "kpostpeak") != NULL ) {
    if (EOBPars->kpostpeak) free(EOBPars->kpostpeak);
    PyListObject *tmp = PyDict_GetItemString(dict, "kpostpeak");
    EOBPars->kpostpeak_size = PyObject_Length(tmp);
    EOBPars->kpostpeak = malloc ( EOBPars->kpostpeak_size * sizeof(int) );
    for (int i = 0; i < EOBPars->kpostpeak_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->kpostpeak[i] = (int) PyLong_AsLong(item); 
    }
  }
  
  /* k nqc peak22 */
  if ( PyDict_GetItemString(dict, "knqcpeak22") != NULL ) {
    if (EOBPars->knqcpeak22) free(EOBPars->knqcpeak22);
    PyListObject *tmp = PyDict_GetItemString(dict, "knqcpeak22");
    EOBPars->knqcpeak22_size = PyObject_Length(tmp);
    EOBPars->knqcpeak22 = malloc ( EOBPars->knqcpeak22_size * sizeof(int) );
    for (int i = 0; i < EOBPars->knqcpeak22_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->knqcpeak22[i] = (int) PyLong_AsLong(item); 
    }
  }

  /* Post Adiabatic Dynamics */
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics") != NULL ) { 
    EOBPars->postadiabatic_dynamics =  YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "postadiabatic_dynamics")));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_N") != NULL ) { 
    EOBPars->postadiabatic_dynamics_N = (int) PyLong_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics_N"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_size") != NULL ) { 
    EOBPars->postadiabatic_dynamics_size = (int) PyLong_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics_size"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_stop") != NULL ) { 
    EOBPars->postadiabatic_dynamics_stop =  YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "postadiabatic_dynamics_stop")));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_rmin") != NULL ) {
    EOBPars->postadiabatic_dynamics_rmin = PyFloat_AsDouble(PyDict_GetItemString(dict, "postadiabatic_dynamics_rmin"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_dr") != NULL ) {
    EOBPars->postadiabatic_dynamics_dr = PyFloat_AsDouble(PyDict_GetItemString(dict, "postadiabatic_dynamics_dr"));
  }
  
  /* rc and flm */
  if ( PyDict_GetItemString(dict, "centrifugal_radius") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "centrifugal_radius"));
    for(EOBPars->centrifugal_radius=0; EOBPars->centrifugal_radius<=CENTRAD_NOPT; EOBPars->centrifugal_radius++){
      if (EOBPars->centrifugal_radius == CENTRAD_NOPT) EOBPars->centrifugal_radius = CENTRAD_NLO;
      if (STREQUAL(val,centrifugal_radius_opt[EOBPars->centrifugal_radius])) break;
    }     
  }
  if ( PyDict_GetItemString(dict, "use_flm") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "use_flm"));
    for(EOBPars->use_flm=0; EOBPars->use_flm<=USEFLM_NOPT; EOBPars->use_flm++){
      if (EOBPars->use_flm == USEFLM_NOPT) EOBPars->use_flm = USEFLM_HM;
      if (STREQUAL(val,use_flm_opt[EOBPars->use_flm])) break;
    }     
  }

  /* NQC */
  if ( PyDict_GetItemString(dict, "nqc") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "nqc"));
    for(EOBPars->nqc=0; EOBPars->nqc<=NQC_NOPT; EOBPars->nqc++){
      if (EOBPars->nqc == NQC_NOPT) EOBPars->nqc = NQC_AUTO;
      if (STREQUAL(val,nqc_opt[EOBPars->nqc])) break;
    }     
  }
  if ( PyDict_GetItemString(dict, "nqc_coefs_flx") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "nqc_coefs_flx"));
    for(EOBPars->nqc_coefs_flx=0; EOBPars->nqc_coefs_flx<=NQC_FLX_NOPT; EOBPars->nqc_coefs_flx++){
      if (EOBPars->nqc_coefs_flx == NQC_FLX_NOPT) EOBPars->nqc_coefs_flx = NQC_FLX_NONE;
      if (STREQUAL(val,nqc_flx_opt[EOBPars->nqc_coefs_flx])) break;
    }     
  }
  if ( PyDict_GetItemString(dict, "nqc_coefs_hlm") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "nqc_coefs_hlm"));
    for(EOBPars->nqc_coefs_hlm=0; EOBPars->nqc_coefs_hlm<=NQC_HLM_NOPT; EOBPars->nqc_coefs_hlm++){
      if (EOBPars->nqc_coefs_hlm == NQC_HLM_NOPT) EOBPars->nqc_coefs_hlm = NQC_HLM_NONE;
      if (STREQUAL(val,nqc_hlm_opt[EOBPars->nqc_coefs_hlm])) break;
    }     
  }

  /* LR and LSO */
  if ( PyDict_GetItemString(dict, "compute_LR") != NULL ) { 
    EOBPars->compute_LR = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "compute_LR")));
  }
  if ( PyDict_GetItemString(dict, "compute_LR_guess") != NULL ) { 
    EOBPars->compute_LR_guess = PyFloat_AsDouble(PyDict_GetItemString(dict, "compute_LR_guess"));
  }
  if ( PyDict_GetItemString(dict, "compute_LSO") != NULL ) { 
    EOBPars->compute_LSO = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "compute_LSO")));
  }
  if ( PyDict_GetItemString(dict, "compute_LSO_guess") != NULL ) { 
    EOBPars->compute_LSO_guess = PyFloat_AsDouble(PyDict_GetItemString(dict, "compute_LSO_guess"));
  }

  /* Spin Dynamics */
  if ( PyDict_GetItemString(dict, "project_spins") != NULL ) { 
    EOBPars->project_spins = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "project_spins")));
  }
  if ( PyDict_GetItemString(dict, "spin_interp_domain") != NULL ) { 
    EOBPars->spin_interp_domain = (int) PyLong_AsLong(PyDict_GetItemString(dict, "spin_interp_domain"));
  }
  if ( PyDict_GetItemString(dict, "spin_flx") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "spin_flx"));
    for(EOBPars->spin_flx=0; EOBPars->spin_flx<=SPIN_FLX_NOPT; EOBPars->spin_flx++){
      if (EOBPars->spin_flx == SPIN_FLX_NOPT) EOBPars->spin_flx = SPIN_FLX_PN;
      if (STREQUAL(val,spin_flx_opt[EOBPars->spin_flx])) break;
    }     
  }
  if ( PyDict_GetItemString(dict, "ringdown_eulerangles") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "ringdown_eulerangles"));
    for(EOBPars->ringdown_eulerangles=0; EOBPars->ringdown_eulerangles<=RD_EULERANGLES_NOPT; EOBPars->ringdown_eulerangles++){
      if (EOBPars->ringdown_eulerangles == RD_EULERANGLES_NOPT) EOBPars->ringdown_eulerangles = RD_EULERANGLES_QNMs;
      if (STREQUAL(val,ringdown_eulerangles_opt[EOBPars->ringdown_eulerangles])) break;
    }     
  }

  /* Output */
  if ( PyDict_GetItemString(dict, "output_hpc") != NULL ) { 
    EOBPars->output_hpc = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_hpc")));
  }
  if ( PyDict_GetItemString(dict, "output_multipoles") != NULL ) { 
    EOBPars->output_multipoles = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_multipoles")));
  }
  if ( PyDict_GetItemString(dict, "output_dynamics") != NULL ) { 
    EOBPars->output_dynamics = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_dynamics")));
  }
  if ( PyDict_GetItemString(dict, "output_nqc") != NULL ) { 
    EOBPars->output_nqc = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_nqc")));
  }
  if ( PyDict_GetItemString(dict, "output_nqc_coefs") != NULL ) { 
    EOBPars->output_nqc_coefs = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_nqc_coefs")));
  }
  if ( PyDict_GetItemString(dict, "output_ringdown") != NULL ) { 
    EOBPars->output_ringdown = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_ringdown")));
  }
  if ( PyDict_GetItemString(dict, "output_dir") != NULL ) {
    strcpy(EOBPars->output_dir, PyUnicode_AsUTF8(PyDict_GetItemString(dict, "output_dir")));
  }

  /* ODE */

  if ( PyDict_GetItemString(dict, "ode_timestep") != NULL ) { 
    char* val;
    val = PyUnicode_AsUTF8(PyDict_GetItemString(dict, "ringdown_eulerangles"));
    for(EOBPars->ode_timestep=0; EOBPars->ode_timestep<=ODE_TSTEP_NOPT; EOBPars->ode_timestep++){
      if (EOBPars->ode_timestep == ODE_TSTEP_NOPT) EOBPars->ode_timestep = ODE_TSTEP_ADAPTIVE;
      if (STREQUAL(val,ode_tstep_opt[EOBPars->ode_timestep])) break;
    }     
  }
  if ( PyDict_GetItemString(dict, "ode_abstol") != NULL ) { 
    EOBPars->ode_abstol = PyFloat_AsDouble(PyDict_GetItemString(dict, "ode_abstol"));
  }
  if ( PyDict_GetItemString(dict, "ode_reltol") != NULL ) { 
    EOBPars->ode_reltol = PyFloat_AsDouble(PyDict_GetItemString(dict, "ode_reltol"));
  }
  if ( PyDict_GetItemString(dict, "ode_tmax") != NULL ) { 
    EOBPars->ode_tmax = PyFloat_AsDouble(PyDict_GetItemString(dict, "ode_tmax"));
  }
  if ( PyDict_GetItemString(dict, "ode_stop_radius") != NULL ) { 
    EOBPars->ode_stop_radius = PyFloat_AsDouble(PyDict_GetItemString(dict, "ode_stop_radius"));
  }
  if ( PyDict_GetItemString(dict, "ode_stop_afterNdt") != NULL ) { 
    EOBPars->ode_stop_afterNdt = (int) PyLong_AsLong(PyDict_GetItemString(dict, "ode_stop_afterNdt"));
  }
  if ( PyDict_GetItemString(dict, "ode_stop_after_peak") != NULL ) { 
    EOBPars->ode_stop_after_peak = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "ode_stop_after_peak")));
  }

  /* FD */
  if ( PyDict_GetItemString(dict, "domain") != NULL ) { 
    EOBPars->domain = (int) PyLong_AsLong(PyDict_GetItemString(dict, "domain"));
  }
  if ( PyDict_GetItemString(dict, "df") != NULL ) { 
    EOBPars->df = PyFloat_AsDouble(PyDict_GetItemString(dict, "df"));
  }
  if ( PyDict_GetItemString(dict,"time_shift_FD") != NULL ) { 
    EOBPars->time_shift_FD = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "time_shift_FD")));
  }

  if ( PyDict_GetItemString(dict,"interp_freqs") != NULL ) {
    EOBPars->interp_freqs = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "interp_freqs")));
  }
  if ( PyDict_GetItemString(dict, "freqs") != NULL ) {
    if (EOBPars->freqs) free(EOBPars->freqs);
    PyListObject *tmp = PyDict_GetItemString(dict, "freqs");
    EOBPars->freqs_size = PyObject_Length(tmp);
    EOBPars->freqs = malloc ( EOBPars->freqs_size * sizeof(double) );
    for (int i = 0; i < EOBPars->freqs_size; i++){
      PyObject *item = PyList_GetItem(tmp, i);
      EOBPars->freqs[i] = PyFloat_AsDouble(item);
    }
  }

  return OK;
}

/*
 * Wrapper
 */

/* Wrapped functions */
static PyObject* EOBRunPy(PyObject* self, PyObject* args)
{
  PyObject* dict;

  /* Parse the input: from python float to c double or from dictionary to C pointer
     https://docs.python.org/3/c-api/arg.html 
     https://docs.python.org/2/c-api/dict.html
  */
  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) 
    return NULL;

  /* alloc output, set some defaults */

  Waveform *hpc = NULL; /* TD wvf */
  Waveform_lm *hmodes = NULL; /* modes */
  Waveform_lm *hTmodes = NULL; /* twisted modes */
  Waveform_lm *hTmmodes = NULL; /*twisted modes, m<0 */
  Waveform_lm *hT0modes = NULL; /*twisted modes, m<0 */

  WaveformFD *hfpc = NULL; /* FD wvf */
  WaveformFD_lm *hfmodes = NULL; /* modes */

  Dynamics *dynf = NULL;
  WaveformFD_lm *hfTmodes = NULL; /* twisted modes */

  int fc = 1;
  int default_choice = BINARY_BBH; /* default_choice, set to BBH */

  /* alloc EOBPars and set defaults based on Lambdas */
  EOBParameters_alloc ( &EOBPars ); 

  if ( PyDict_GetItemString(dict, "LambdaAl2") != NULL )
    EOBPars->LambdaAl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl2"));
  if ( PyDict_GetItemString(dict, "LambdaBl2") != NULL )
    EOBPars->LambdaBl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl2"));
  
  /* Add a warning for users */
  if ( PyDict_GetItemString(dict, "Lambda1") != NULL )
    errorexit("'Lambda1', 'Lambda2' are deprecated. Use LambdaAl2, LambdaBl2 instead.");
  if ( PyDict_GetItemString(dict, "Lambda2") != NULL )
    errorexit("'Lambda1', 'Lambda2' are deprecated. Use LambdaAl2, LambdaBl2 instead.");
  

  if(EOBPars->LambdaAl2 > 1. && EOBPars->LambdaBl2 > 1.) default_choice = BINARY_BNS;
  if(EOBPars->LambdaAl2 == 0. && EOBPars->LambdaBl2 > 1.) default_choice = BINARY_BHNS;
  EOBParameters_defaults (default_choice, EOBPars);  

  /* Read the dictionary in EOBPars */
  /* RG: there has to be a faster way...*/

  EOBPars->M = PyFloat_AsDouble(PyDict_GetItemString(dict, "M"));
  EOBPars->q = PyFloat_AsDouble(PyDict_GetItemString(dict, "q"));

  /* Tides*/
  if ( PyDict_GetItemString(dict, "LambdaAl2") != NULL )
    EOBPars->LambdaAl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaAl2"));
  if ( PyDict_GetItemString(dict, "LambdaBl2") != NULL )
    EOBPars->LambdaBl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "LambdaBl2"));
  
  /* Spins*/
  if (PyDict_GetItemString(dict, "chi1x") != NULL )
      EOBPars->chi1x = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi1x"));
  if (PyDict_GetItemString(dict, "chi1y") != NULL )
      EOBPars->chi1y = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi1y"));
  if (PyDict_GetItemString(dict, "chi1z") != NULL )
      EOBPars->chi1z = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi1z"));
  if (PyDict_GetItemString(dict, "chi1") != NULL )
      EOBPars->chi1 = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi1"));

  if (PyDict_GetItemString(dict, "chi2x") != NULL )
      EOBPars->chi2x = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi2x"));
  if (PyDict_GetItemString(dict, "chi2y") != NULL )
      EOBPars->chi2y = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi2y"));
  if (PyDict_GetItemString(dict, "chi2z") != NULL )
      EOBPars->chi2z = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi2z"));
  if (PyDict_GetItemString(dict, "chi2") != NULL )
     EOBPars->chi2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi2"));


  for (int k=0; k < NFIRSTCALL; k++){ 
    EOBPars->firstcall[k] = 1;
  }

  /* Optional arguments for the dictionary */

  /* Options */
  SetOptionalVariables(dict);

  /* output */
  const int output = EOBPars->output_dynamics
    + EOBPars->output_multipoles
    + EOBPars->output_hpc
    + EOBPars->output_nqc; 
  if (output) {
      if (system_mkdir(EOBPars->output_dir)) {
        printf("ERROR(TEOBResumS): %s\n",eob_error_msg[ERROR_MKDIR]);
        return ERROR_MKDIR;
      }
    }

  if (output){
    char outpar[STRLEN];
    strcpy(outpar,EOBPars->output_dir);
    EOBParameters_tofile(EOBPars,strcat(outpar,"/params.txt"));
  }

  eob_set_params(default_choice, fc);

  /* Overwrite spin-spin parameters, if required */
  if ( PyDict_GetItemString(dict, "C_Q1") != NULL ) { 
    EOBPars->C_Q1 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Q1"));
  }  
  if ( PyDict_GetItemString(dict, "C_Q2") != NULL ) { 
    EOBPars->C_Q2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Q2"));
  }
  if ( PyDict_GetItemString(dict, "C_Oct1") != NULL ) { 
    EOBPars->C_Oct1 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Oct1"));
  }  
  if ( PyDict_GetItemString(dict, "C_Oct2") != NULL ) { 
    EOBPars->C_Oct2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Oct2"));
  }
  if ( PyDict_GetItemString(dict, "C_Hex1") != NULL ) { 
    EOBPars->C_Hex1 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Hex1"));
  }  
  if ( PyDict_GetItemString(dict, "C_Hex2") != NULL ) { 
    EOBPars->C_Hex2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "C_Hex2"));
  }

  /* Overwrite a6c and cN3LO, if required */
  if ( PyDict_GetItemString(dict, "a6c") != NULL ) { 
    EOBPars->a6c = PyFloat_AsDouble(PyDict_GetItemString(dict, "a6c"));
  }   
  if ( PyDict_GetItemString(dict, "cN3LO") != NULL ) { 
    EOBPars->cN3LO = PyFloat_AsDouble(PyDict_GetItemString(dict, "cN3LO"));
  }   

  /* Run */
  int status = EOBRun(&hpc,    &hfpc, 
                      &hmodes, &hfmodes, &dynf,
                      &hTmodes, &hTmmodes, &hT0modes,
                      &hfTmodes,
                      default_choice, fc);

  if (status) printf("ERROR(TEOBResumS): %s\n",eob_error_msg[status]);  
  /*  Construct the output arrays */

  /* return modes? */
  int arg_out = 0; 
  if ( PyDict_GetItemString(dict, "arg_out") != NULL ) { 
    arg_out = YESNO2INT(PyUnicode_AsUTF8(PyDict_GetItemString(dict, "arg_out")));
  }
  
  /* build the dynamics dictionary */
  PyObject* dyndict  = PyDict_New(); /*dyn dictionary */
  npy_intp dims_dyn[1];
  dims_dyn[0] = dynf->size;

  /*time */
  double *pdt;
  PyArrayObject *pdto = (PyArrayObject *) PyArray_SimpleNew(1,dims_dyn,NPY_DOUBLE);
  pdt = pyvector_to_Carrayptrs(pdto);
  memcpy(pdt, dynf->time, dynf->size *sizeof(double));
  PyDict_SetItemString(dyndict, "t", pdto); 
  Py_DECREF(pdto);

  /* other variables */
  for(int v=0; v <EOB_DYNAMICS_NVARS; v++){
    double *pv;
    PyArrayObject *pvo = (PyArrayObject *) PyArray_SimpleNew(1,dims_dyn,NPY_DOUBLE);
    pv = pyvector_to_Carrayptrs(pvo);
    memcpy(pv, dynf->data[v], dynf->size *sizeof(double));
    PyDict_SetItemString(dyndict, eob_var[v], pvo); 
    Py_DECREF(pvo);
  }

  if(EOBPars->domain==DOMAIN_TD){
    
    double *pt, *php, *phc; /*t, h+ and hx */
    PyObject* hlmdict  = PyDict_New(); /*hlm dictionary */

    npy_intp dims[1];
    dims[0] = hpc->size;
    PyArrayObject *pto;
    PyArrayObject *phpo;
    PyArrayObject *phco;

    pto  = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phpo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phco = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);

    /* Cast py *arrays into C *arrays   */
    pt = pyvector_to_Carrayptrs(pto);
    php = pyvector_to_Carrayptrs(phpo);
    phc = pyvector_to_Carrayptrs(phco);

    /* Copy */
    memcpy(pt,  hpc->time, hpc->size * sizeof(double)); //t
    memcpy(php, hpc->real, hpc->size * sizeof(double)); //h+
    memcpy(phc, hpc->imag, hpc->size * sizeof(double)); //hx

    /*build hlm dictionary */
    if (EOBPars->use_spins ==  MODE_SPINS_ALIGNED){
      for(int k=0; k<KMAX; k++){
        if(hmodes->kmask[k]){
          double *pAhlm, *pphlm;
          PyArrayObject *pAhlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          PyArrayObject *pphlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          pAhlm = pyvector_to_Carrayptrs(pAhlmo);
          pphlm = pyvector_to_Carrayptrs(pphlmo);
          memcpy(pAhlm, hmodes->ampli[k], hmodes->size * sizeof(double));
          memcpy(pphlm, hmodes->phase[k], hmodes->size * sizeof(double));
          /* build dictionary entry*/
          PyObject *obj = Py_BuildValue("O:O", pAhlmo, pphlmo);
          char kst[12];
          sprintf(kst, "%i", k);
          /* populate the dictionary */
          PyDict_SetItemString(hlmdict, kst, obj);
          /* free */
          Py_DECREF(pAhlmo);
          Py_DECREF(pphlmo);
          Py_DECREF(obj);
        }
      }
    } else {

      int addzero = 1;
      for(int k=0; k<KMAX; k++){
        /* Add m > 0 */
        if(hTmodes->kmask[k]){
          double *pAhlm, *pphlm;
          PyArrayObject *pAhlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          PyArrayObject *pphlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          pAhlm = pyvector_to_Carrayptrs(pAhlmo);
          pphlm = pyvector_to_Carrayptrs(pphlmo);
          memcpy(pAhlm, hTmodes->ampli[k], hTmodes->size * sizeof(double));
          memcpy(pphlm, hTmodes->phase[k], hTmodes->size * sizeof(double));
          /* build dictionary entry*/
          PyObject *obj = Py_BuildValue("O:O", pAhlmo, pphlmo);
          char kst[12];
          sprintf(kst, "%i", k);
          /* populate the dictionary */
          PyDict_SetItemString(hlmdict, kst, obj); 
          /* free */
          Py_DECREF(pAhlmo);
          Py_DECREF(pphlmo);
          Py_DECREF(obj);
        }
        /* Add also the m<0 modes to the dictionary */
        if(hTmmodes->kmask[k]){
          double *pAhlm, *pphlm;
          PyArrayObject *pAhlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          PyArrayObject *pphlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          pAhlm = pyvector_to_Carrayptrs(pAhlmo);
          pphlm = pyvector_to_Carrayptrs(pphlmo);
          memcpy(pAhlm, hTmmodes->ampli[k], hTmmodes->size * sizeof(double));
          memcpy(pphlm, hTmmodes->phase[k], hTmmodes->size * sizeof(double));
          /* build dictionary entry*/
          PyObject *obj = Py_BuildValue("O:O", pAhlmo, pphlmo);
          char kst[12];
          sprintf(kst, "-%i", k);
          /* populate the dictionary */
          PyDict_SetItemString(hlmdict, kst, obj); 
          /* free */
          Py_DECREF(pAhlmo);
          Py_DECREF(pphlmo);
          Py_DECREF(obj);
        } 
        /* Finally, add the m=0 ones */  
        if(hT0modes->kmask[k] && addzero){
          /* do the thing */
          double *pAhlm, *pphlm;
          PyArrayObject *pAhlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          PyArrayObject *pphlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
          pAhlm = pyvector_to_Carrayptrs(pAhlmo);
          pphlm = pyvector_to_Carrayptrs(pphlmo);
          memcpy(pAhlm, hT0modes->ampli[k], hT0modes->size * sizeof(double));
          memcpy(pphlm, hT0modes->phase[k], hT0modes->size * sizeof(double));
          /* build dictionary entry*/
          PyObject *obj = Py_BuildValue("O:O", pAhlmo, pphlmo);
          char kst[12];
          sprintf(kst, "%i0", LINDEX[k]);
          /* populate the dictionary */
          PyDict_SetItemString(hlmdict, kst, obj); 
          /* free */
          Py_DECREF(pAhlmo);
          Py_DECREF(pphlmo);
          Py_DECREF(obj);
          addzero = 0;
        }
        if(LINDEX[k+1]-LINDEX[k]) addzero = 1;
      }
    }
    /* build the final object */
    PyObject *ret;
    if (arg_out == 0){
      ret = Py_BuildValue("OOO", pto, phpo, phco);
    } else if (arg_out == 1){
      ret = Py_BuildValue("OOOOO", pto, phpo, phco, hlmdict, dyndict);
    } else {
      printf("ERROR: arg_out has to be equal to 'yes' or 'no' ");
      ret = NULL;
    }
    /* Free C memory */
    Waveform_free (hpc);          
    WaveformFD_free (hfpc);       
    Waveform_lm_free(hmodes);    
    Waveform_lm_free (hTmodes);
    Waveform_lm_free (hTmmodes);
    Waveform_lm_free (hT0modes);
    WaveformFD_lm_free(hfmodes);  
    WaveformFD_lm_free (hfTmodes);
    EOBParameters_free (EOBPars);
    Dynamics_free (dynf);

    /* "Free" Python objects */
    Py_DECREF(pto);               
    Py_DECREF(phpo);
    Py_DECREF(phco);
    Py_DECREF(hlmdict);
    Py_DECREF(dyndict);

    return ret;  

  } else {

    double *pf, *phpr, *phpi, *phcr, *phci; /*FD: f,  Re and Im of h+, hx */
    PyObject* hflmdict = PyDict_New();      /* hlm FD dictionary */
    PyObject* htlmdict = PyDict_New();      /* hlm TD dictionary */

    npy_intp dims[1];
    dims[0] = hfpc->size;
    PyArrayObject *pfo;                 /* f  */
    PyArrayObject *phprealo, *phpimago; /* h+ */
    PyArrayObject *phcrealo, *phcimago; /* hx */

    pfo      = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phprealo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phpimago = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phcrealo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    phcimago = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    
    /* Cast py *arrays into C *arrays   */
    pf   = pyvector_to_Carrayptrs(pfo);
    phpr = pyvector_to_Carrayptrs(phprealo);
    phpi = pyvector_to_Carrayptrs(phpimago);
    phcr = pyvector_to_Carrayptrs(phcrealo);
    phci = pyvector_to_Carrayptrs(phcimago);

    memcpy(pf,   hfpc->freq,  hfpc->size * sizeof(double)); //f
    memcpy(phpr, hfpc->preal, hfpc->size * sizeof(double)); //Re h+
    memcpy(phpi, hfpc->pimag, hfpc->size * sizeof(double)); //Im h+
    memcpy(phcr, hfpc->creal, hfpc->size * sizeof(double)); //Re hx
    memcpy(phci, hfpc->cimag, hfpc->size * sizeof(double)); //Im hx
    
    /*build hflm dictionary */ 
    for(int k=0; k<KMAX; k++){
      if(hfmodes->kmask[k]){
        double *pAhflm, *pphflm;
        PyArrayObject  *pAhflmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
        PyArrayObject  *pphflmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
        pAhflm = pyvector_to_Carrayptrs(pAhflmo);
        pphflm = pyvector_to_Carrayptrs(pphflmo);
        memcpy(pAhflm, hfmodes->ampli[k], hfmodes->size * sizeof(double));
        memcpy(pphflm, hfmodes->phase[k], hfmodes->size * sizeof(double));
        /*build dictionary entry*/
        PyObject *obj = Py_BuildValue("O:O", pAhflmo, pphflmo);
        char kst[12];
        sprintf(kst, "%i", k);
        /*populate the dictionary*/
        PyDict_SetItemString(hflmdict, kst, obj);
        /*free*/
        Py_DECREF(pAhflmo);
        Py_DECREF(pphflmo);
        Py_DECREF(obj);
      }
    }

    /*build htlm dictionary */
    dims[0] = hmodes->size;
    for(int k=0; k<KMAX; k++){
      if(hmodes->kmask[k]){
        double *pAhlm, *pphlm;
        PyArrayObject *pAhlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
        PyArrayObject *pphlmo = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
        pAhlm = pyvector_to_Carrayptrs(pAhlmo);
        pphlm = pyvector_to_Carrayptrs(pphlmo);
        memcpy(pAhlm, hmodes->ampli[k], hmodes->size * sizeof(double));
        memcpy(pphlm, hmodes->phase[k], hmodes->size * sizeof(double));
        /* build dictionary entry*/
        PyObject *obj = Py_BuildValue("O:O", pAhlmo, pphlmo);
        char kst[12];
        sprintf(kst, "%i", k);
        /* populate the dictionary */
        PyDict_SetItemString(htlmdict, kst, obj); 
        /* free */
        Py_DECREF(pAhlmo);
        Py_DECREF(pphlmo);
        Py_DECREF(obj);
      }
    }
    /* add time to htlm */
    double *pt;
    PyArrayObject *pto = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    pt = pyvector_to_Carrayptrs(pto);
    memcpy(pt, hmodes->time, hmodes->size * sizeof(double));
    PyObject *obj = Py_BuildValue("O", pto);
    PyDict_SetItemString(htlmdict,"t", obj);
    Py_DECREF(obj);
    Py_DECREF(pto);

    /* build the final object */
    PyObject *ret;
    if (arg_out == 0){
      ret = Py_BuildValue("OOOOO", pfo, phprealo, phpimago, phcrealo, phcimago);
    } else if (arg_out == 1){
      ret = Py_BuildValue("OOOOOOOO", pfo, phprealo, phpimago, phcrealo, phcimago, hflmdict, htlmdict, dyndict);
    } else {
      printf("ERROR: arg_out has to be equal to 'yes' or 'no'");
      ret = NULL;
    }

    /* Free C memory */
    Waveform_free (hpc);          
    WaveformFD_free (hfpc); 
    Waveform_lm_free(hmodes); 
    Waveform_lm_free (hTmodes);
    Waveform_lm_free (hTmmodes);
    Waveform_lm_free (hT0modes);
    WaveformFD_lm_free(hfmodes);
    WaveformFD_lm_free (hfTmodes);
    EOBParameters_free (EOBPars);
    Dynamics_free (dynf);

    /* "Free" Python objects */
    Py_DECREF(pfo);             
    Py_DECREF(phprealo);
    Py_DECREF(phpimago);
    Py_DECREF(phcrealo);
    Py_DECREF(phcimago);
    Py_DECREF(hflmdict);
    Py_DECREF(htlmdict);
    Py_DECREF(dyndict);

    return ret;
  }
}

/* Define functions in module */
static PyMethodDef EOBRunMethods[] = {
  {"EOBRunPy", EOBRunPy, METH_VARARGS, "Generate a time or frequency domain TEOBResumS waveform"},
  {NULL} 
};

#if PY_MAJOR_VERSION >= 3

/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT,
  "EOBRun_module", "Some documentation",
  -1,
  EOBRunMethods
};

PyMODINIT_FUNC
PyInit_EOBRun_module(void)
{
  PyObject *module;
  module = PyModule_Create(&cModPyDem);
  if(module==NULL) return NULL;
  import_array();  /* IMPORTANT: this must be called */
  if (PyErr_Occurred()) return NULL;
  return module;
}

#else

/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initEOBRun_module(void)
{
  //(void) Py_InitModule("EOBRunTD_module", EOBRunTDMethods);
  PyObject *module;
  module = Py_InitModule("EOBRun_module", EOBRunMethods);
  if(module==NULL) return;
  import_array();  /* IMPORTANT: this must be called */
  return;
}

#endif
