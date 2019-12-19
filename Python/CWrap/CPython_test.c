/* Example of wrapping with the Python-C-API. 
 * See:
 * https://scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
 * And C_arraytest.c - SciPy Cookbook
 */

#include <Python.h>
#include <numpy/arrayobject.h>

/* 
 * My C function/program 
 */
#define OK (0);
#define ERROR (1);
int foo(int size, double a, double **p)
{
  *p = (double*)malloc(size * sizeof(double));
  if(*p == NULL) return ERROR;
  for (int i=0; i<size; ++i) (*p)[i] = a;
  return OK;
}

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
static PyObject* foo_func(PyObject* self, PyObject* args)
{
  int size;
  double a;
  double *p, *pa;

  /* Parse the input, from python float to c double 
     https://docs.python.org/3/c-api/arg.html 
  */
  if (!PyArg_ParseTuple(args, "id", &size, &a))
    return NULL;
  
  /* Call the C function */
  if (foo(size, a, &p)) return NULL;
  
  /*  Construct the output array */
  npy_intp dims[1];
  dims[0] = size;
  PyArrayObject *pao;
  pao = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
  
  /* Cast py *arrays into C *arrays   */
  pa =(double*)PyArray_DATA(pao);
  
  /* Copy */
  memcpy(pa, p, size * sizeof(double));
  free(p); /* Free C memory */

  return PyArray_Return(pao);
  /* return Py_BuildValue("O", pao); */ /* This also works, maybe better for multiple outputs? */
}

/*
 * Define module
 */

/* Define functions in module */
static PyMethodDef FooMethods[] = {
  {"foo", foo_func, METH_VARARGS, "Init and return a double array"},
  /* SB: Not understood following line, but uncommented version
  prevent a segfault after runtime ... */
  {NULL, NULL, 0, NULL} 
  //{NULL, NULL, 0, NULL} //*/ 
};

#if PY_MAJOR_VERSION >= 3

/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem = {
  PyModuleDef_HEAD_INIT,
  "foo_module", "Some documentation",
  -1,
  FooMethods
};

PyMODINIT_FUNC
PyInit_foo_module(void)
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
initfoo_module(void)
{
  //(void) Py_InitModule("foo_module", FooMethods);
  PyObject *module;
  module = Py_InitModule("foo_module", FooMethods);
  if(module==NULL) return;
  import_array();  /* IMPORTANT: this must be called */
  return;
}

#endif
