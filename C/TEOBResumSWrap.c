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
  if ( PyDict_GetItemString(dict, "use_tidal") != NULL ) { 
    EOBPars->use_tidal = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_tidal"));
  }
  if ( PyDict_GetItemString(dict, "use_tidal_gravitomagnetic") != NULL ) { 
    EOBPars->use_tidal_gravitomagnetic = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_tidal_gravitomagnetic"));
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
    EOBPars->size = (int) PyInt_AsLong(PyDict_GetItemString(dict, "size"));
  }
  if ( PyDict_GetItemString(dict, "ringdown_extend_array") != NULL ) { 
    EOBPars->ringdown_extend_array = (int) PyInt_AsLong(PyDict_GetItemString(dict, "ringdown_extend_array"));
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

  /* rc and flm */
  if ( PyDict_GetItemString(dict, "centrifugal_radius") != NULL ) { 
    EOBPars->centrifugal_radius = (int) PyInt_AsLong(PyDict_GetItemString(dict, "centrifugal_radius"));
  }
  if ( PyDict_GetItemString(dict, "use_flm") != NULL ) { 
    EOBPars->use_flm = (int) PyInt_AsLong(PyDict_GetItemString(dict, "use_flm"));
  }

  /* NQC */
  if ( PyDict_GetItemString(dict, "nqc") != NULL ) { 
    EOBPars->nqc = (int) PyInt_AsLong(PyDict_GetItemString(dict, "nqc"));
  }
  if ( PyDict_GetItemString(dict, "nqc_coefs_flx") != NULL ) { 
    EOBPars->nqc_coefs_flx = (int) PyInt_AsLong(PyDict_GetItemString(dict, "nqc_coefs_flx"));
  }
  if ( PyDict_GetItemString(dict, "nqc_coefs_hlm") != NULL ) { 
    EOBPars->nqc_coefs_hlm = (int) PyInt_AsLong(PyDict_GetItemString(dict, "nqc_coefs_hlm"));
  }


  /* LR and LSO */
  if ( PyDict_GetItemString(dict, "compute_LR") != NULL ) { 
    EOBPars->compute_LR = (int) PyInt_AsLong(PyDict_GetItemString(dict, "compute_LR"));
  }
  if ( PyDict_GetItemString(dict, "compute_LR_guess") != NULL ) { 
    EOBPars->compute_LR_guess = PyFloat_AsDouble(PyDict_GetItemString(dict, "compute_LR_guess"));
  }
  if ( PyDict_GetItemString(dict, "compute_LSO") != NULL ) { 
    EOBPars->compute_LSO = (int) PyInt_AsLong(PyDict_GetItemString(dict, "compute_LSO"));
  }
  if ( PyDict_GetItemString(dict, "compute_LSO_guess") != NULL ) { 
    EOBPars->compute_LSO_guess = PyFloat_AsDouble(PyDict_GetItemString(dict, "compute_LSO_guess"));
  }
  
  /* Output */
  if ( PyDict_GetItemString(dict, "output_hpc") != NULL ) { 
    EOBPars->output_hpc = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_hpc"));
  }
  if ( PyDict_GetItemString(dict, "output_multipoles") != NULL ) { 
    EOBPars->output_multipoles = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_multipoles"));
  }
  if ( PyDict_GetItemString(dict, "output_dynamics") != NULL ) { 
    EOBPars->output_dynamics = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_dynamics"));
  }
  if ( PyDict_GetItemString(dict, "output_nqc") != NULL ) { 
    EOBPars->output_dynamics = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_dynamics"));
  }
  if ( PyDict_GetItemString(dict, "output_nqc") != NULL ) { 
    EOBPars->output_nqc = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_nqc"));
  }
  if ( PyDict_GetItemString(dict, "output_nqc_coefs") != NULL ) { 
    EOBPars->output_nqc_coefs = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_nqc_coefs"));
  }
  if ( PyDict_GetItemString(dict, "output_ringdown") != NULL ) { 
    EOBPars->output_ringdown = (int) PyInt_AsLong(PyDict_GetItemString(dict, "output_ringdown"));
  }

  /* ODE */

  if ( PyDict_GetItemString(dict, "ode_timestep") != NULL ) { 
    EOBPars->ode_timestep = (int) PyInt_AsLong(PyDict_GetItemString(dict, "ode_timestep"));
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
    EOBPars->ode_stop_afterNdt = (int) PyInt_AsLong(PyDict_GetItemString(dict, "ode_stop_afterNdt"));
  }

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
