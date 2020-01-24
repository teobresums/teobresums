/* Example of wrapping with the Python-C-API. 
 * See:
 * https://scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
 * And C_arraytest.c - SciPy Cookbook
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include "src/TEOBResumS.h"

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

/*
 * Wrapper
 */

/* Wrapped function */
static PyObject* EOBRunTD_func(PyObject* self, PyObject* args)
{
  PyObject* dict;
  double *pt, *php, *phc;

  /* Parse the input: from python float to c double or from dictionary to C pointer
     https://docs.python.org/3/c-api/arg.html 
     https://docs.python.org/2/c-api/dict.html
  */
  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) 
    return NULL;

  /* alloc output, set some defaults */
  Waveform *hpc;   
  int fc = 1;
  int default_choice = 0;

  //alloc EOBPars and set defaults based on Lambdas
  EOBParameters_alloc ( &EOBPars ); 

  EOBPars->LambdaAl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "Lambda1"));
  EOBPars->LambdaBl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "Lambda2"));

  if(EOBPars->LambdaAl2 > 1. && EOBPars->LambdaBl2 > 1.) default_choice = 1;
  EOBParameters_defaults (default_choice, EOBPars);  


  /* Read the dictionary in EOBPars */
  /* RG: there has to be a faster way...*/

  EOBPars->LambdaAl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "Lambda1"));
  EOBPars->LambdaBl2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "Lambda2"));
  EOBPars->M = PyFloat_AsDouble(PyDict_GetItemString(dict, "M"));
  EOBPars->q = PyFloat_AsDouble(PyDict_GetItemString(dict, "q"));
  EOBPars->chi1 = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi1"));
  EOBPars->chi2 = PyFloat_AsDouble(PyDict_GetItemString(dict, "chi2"));

  /* Optional arguments for the dictionary */
  // FIXME: reduce number of calls to PyDict_GetItemString
  // for integers, do the same but with int val = (int) PyInt_AsLong(pyObj_val);

  /* Options */
  if ( PyDict_GetItemString(dict, "r0") != NULL ) {
    EOBPars->r0 = PyFloat_AsDouble(PyDict_GetItemString(dict, "r0"));
  }
  if ( PyDict_GetItemString(dict, "initial_frequency") != NULL ) { 
    EOBPars->initial_frequency = PyFloat_AsDouble(PyDict_GetItemString(dict, "initial_frequency"));
  }
  if ( PyDict_GetItemString(dict, "use_geometric_units") != NULL ) { 
    EOBPars->use_geometric_units = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_geometric_units"));
  }
  if ( PyDict_GetItemString(dict, "use_spins") != NULL ) { 
    EOBPars->use_spins = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_spins"));
  }
  if ( PyDict_GetItemString(dict, "use_Yagi_fits") != NULL ) { 
    EOBPars->use_Yagi_fits = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_Yagi_fits"));
  }
  if ( PyDict_GetItemString(dict, "use_speedytail") != NULL ) { 
    EOBPars->use_speedytail = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_speedytail"));
  }
  if ( PyDict_GetItemString(dict, "interp_uniform_grid") != NULL ) { 
    EOBPars->interp_uniform_grid = (int) PyInt_AsLong(PyDict_GetItemString(dict, "interp_uniform_grid"));
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
  if ( PyDict_GetItemString(dict, "srate_interp") != NULL ) { 
    EOBPars->srate_interp = PyFloat_AsDouble(PyDict_GetItemString(dict, "srate_interp"));
  }

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

  /* Post Adiabatic Dynamics */
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics") != NULL ) { 
    EOBPars->postadiabatic_dynamics = (int) PyInt_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_N") != NULL ) { 
    EOBPars->postadiabatic_dynamics_N = (int) PyInt_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics_N"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_size") != NULL ) { 
    EOBPars->postadiabatic_dynamics_size = (int) PyInt_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics_size"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_stop") != NULL ) { 
    EOBPars->postadiabatic_dynamics_stop = (int) PyInt_AsLong(PyDict_GetItemString(dict, "postadiabatic_dynamics_stop"));
  }
  if ( PyDict_GetItemString(dict, "postadiabatic_dynamics_rmin") != NULL ) {
    EOBPars->postadiabatic_dynamics_rmin = PyFloat_AsDouble(PyDict_GetItemString(dict, "postadiabatic_dynamics_rmin"));
  }

  /* options yet to add

  eobp->centrifugal_radius=CENTRAD_LO; // {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES}
  eobp->use_flm=USEFLM_SSLO; // "SSLO", "SSNLO", "HM"
  
  eobp->compute_LR=0; // calculate LR ?
  eobp->compute_LSO=0; // calculate LSO ?
  eobp->compute_LR_guess=3.;
  eobp->compute_LSO_guess=6.;

  eobp->nqc=NQC_AUTO; // {"no", "auto", "manual"}
  eobp->nqc_coefs_flx=NQC_FLX_NONE; // {"none", "nrfit_nospin20160209", "fromfile"}
  eobp->nqc_coefs_hlm=NQC_HLM_NONE; // {"compute", "none", "nrfit_nospin20160209", "fromfile"}
  eobp->output_hpc= 0; // output h+,hx
  eobp->output_multipoles= 0; // output multipoles
  eobp->output_dynamics=0; // output dynamics
  eobp->output_nqc=0; // output NQC waveform
  eobp->output_nqc_coefs=0; // output multipolar NQC coefs (if determined)
  eobp->output_ringdown=0; // output ringdown waveform
  eobp->srate=4096.; // sampling rate, used if input is given in physical unit, reset based on tstep otherwise
  eobp->dt=0.5; // timestep, used if input is given in geometric unit, reset based on srate otherwise
  eobp->size=1; // size of the arrays (chunks, dynamically extended)
  eobp->ringdown_extend_array=500; // grid points to extend arrays for ringdown attachment
  eobp->ode_timestep=ODE_TSTEP_UNIFORM; // specify ODE solver timestep "uniform","adaptive","adaptive+uniform_after_LSO","undefined"
  eobp->ode_abstol=1e-13; // ODE solver absolute accuracy
  eobp->ode_reltol=1e-11; //  ODE solver relative accuracy
  eobp->ode_tmax=2e8; // max integration time
  eobp->ode_stop_radius=2.; // stop ODE integration at this radius (if > 0)
  eobp->ode_stop_afterNdt=4; // stop ODE N iters after the Omega peak

  */

  eob_set_params_EOBRun(default_choice, fc); 
  EOBRunTD(&hpc, default_choice, fc);
  
  /*  Construct the output arrays */
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
  memcpy(pt,  hpc->time, hpc->size * sizeof(double));
  memcpy(php, hpc->real, hpc->size * sizeof(double));
  memcpy(phc, hpc->imag, hpc->size * sizeof(double));

  Waveform_free (hpc); /* Free C memory */
  EOBParameters_free (EOBPars);

  return Py_BuildValue("OOO", pto, phpo, phco);  /* This also works, maybe better for multiple outputs? */
}

/*
 * Define module
 */

/* Define functions in module */
static PyMethodDef EOBRunTDMethods[] = {
  {"EOBRunTD", EOBRunTD_func, METH_VARARGS, "Generate a TD TEOBResumS waveform"},
  /* SB: Not understood following line, but uncommented version
  prevent a segfault after runtime ... */
  {NULL, NULL}  /* {NULL, NULL, 0, NULL} */ 
};

#if PY_MAJOR_VERSION >= 3

/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT,
  "EOBRunTD_module", "Some documentation",
  -1,
  EOBRunTDMethods
};

PyMODINIT_FUNC
PyInit_EOBRunTD_module(void)
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
initEOBRunTD_module(void)
{
  //(void) Py_InitModule("EOBRunTD_module", EOBRunTDMethods);
  PyObject *module;
  module = Py_InitModule("EOBRunTD_module", EOBRunTDMethods);
  if(module==NULL) return;
  import_array();  /* IMPORTANT: this must be called */
  return;
}

#endif
