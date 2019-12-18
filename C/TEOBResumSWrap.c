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
  double M, q, chi1, chi2, LA, LB;
  double *pt, *php, *phc;

  /* Parse the input, from python float to c double 
     https://docs.python.org/3/c-api/arg.html 
  */
  if (!PyArg_ParseTuple(args, "dddddd", &M, &q, &chi1, &chi2, &LA, &LB))
    return NULL;
  
  /* Call the C function */
  Waveform *hpc;   
  int fc = 1;
  int default_choice = 1;
  //alloc
  EOBParameters_alloc( &EOBPars );
  EOBPars->M = M;
  EOBPars->q = q;
  EOBPars->chi1 = chi1;
  EOBPars->chi2 = chi2;
  EOBPars->LambdaAl2 = LA;
  EOBPars->LambdaBl2 = LB;
  EOBRunTD(&hpc, default_choice, fc);
  //if (EOBRunTD(size, a, &p)) return NULL;
  
  /*  Construct the output arrays */
  int dims[2];
  dims[0] = hpc->size;
  PyArrayObject *pto;
  PyArrayObject *phpo;
  PyArrayObject *phco;

  pto  = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
  phpo = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
  phco = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
  
  /* Cast py *arrays into C *arrays   */
  pt = pyvector_to_Carrayptrs(pto);
  php = pyvector_to_Carrayptrs(phpo);
  phc = pyvector_to_Carrayptrs(phco);

  /* Copy */
  memcpy(pt,  hpc->time, hpc->size * sizeof(double));
  memcpy(php, hpc->real, hpc->size * sizeof(double));
  memcpy(phc, hpc->imag, hpc->size * sizeof(double));

  Waveform_free (hpc); /* Free C memory */

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
