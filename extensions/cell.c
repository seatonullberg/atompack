#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "vector.h"

static PyObject *py_cell_contains(PyObject *self, PyObject *args);
static PyObject *py_cell_enforce(PyObject *self, PyObject *args);

// Returns 1 if the position is within the cell else 0
int cell_contains(double cell[3][3], double position[3], double tolerance) {
  double mag;
  int contains = 1; // true
  for (int i = 0; i < 3; i++) {
    mag = norm(cell[i]);
    if (position[i] > mag + tolerance) {
      contains = 0; // false
      break;
    }
    if (position[i] < -tolerance) {
      contains = 0; // false
      break;
    }
  }
  return contains;
}

// Enforces that a position is within a cell.
// The out parameter contains the initial position 
// transformed such that it is contained by the cell.
void cell_enforce(double cell[3][3], double position[3], double tolerance, double out[3]) {
  double mag;
  for (int i = 0; i < 3; i++) {
    mag = norm(cell[i]);
    out[i] = position[i];
    if (out[i] > mag - tolerance) {
      while (out[i] > mag - tolerance) {
        out[i] -= mag;
      }
    }
    if (out[i] < -tolerance) {
      while(out[i] < -tolerance) {
        out[i] += mag;
      }
    }
  }
}

/******************************
 *  Python Wrapper Functions  *
 ******************************/

static PyObject *py_cell_contains(PyObject *self, PyObject *args) {
  PyArrayObject *py_cell;
  double(*cell)[3];
  PyArrayObject *py_position;
  double(*position);
  double tolerance;

  if (!PyArg_ParseTuple(args, "OOd", &py_cell, &py_position, &tolerance)) {
    return NULL;
  }

  cell = (double(*)[3])PyArray_DATA(py_cell);
  position = (double(*))PyArray_DATA(py_position);

  if (cell_contains(cell, position, tolerance)) {
    Py_RETURN_TRUE;
  }
  Py_RETURN_FALSE;
}

static PyObject *py_cell_enforce(PyObject *self, PyObject *args) {
  PyArrayObject *py_cell;
  double(*cell)[3];
  PyArrayObject *py_position;
  double(*position);
  double tolerance;
  double out[3];
  PyArrayObject *res;

  if (!PyArg_ParseTuple(args, "OOd", &py_cell, &py_position, &tolerance)) {
    return NULL;
  }

  cell = (double(*)[3])PyArray_DATA(py_cell);
  position = (double(*))PyArray_DATA(py_position);
  cell_enforce(cell, position, tolerance, out);
  
  double *ptr;
  for (int i = 0; i < 3; i++) {
    ptr = PyArray_GETPTR1(py_position, i);
    *ptr = out[i];
  }
  Py_RETURN_NONE;

}


/*******************************
 *  Python Module Description  *
 *******************************/

static PyMethodDef method_def[] = {
    {"cell_contains", py_cell_contains, METH_VARARGS,
     "Returns True if point is within cell."},
    {"cell_enforce", py_cell_enforce, METH_VARARGS,
      "Enforces that a position is within a cell."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_cell",
    "Module `_cell` provides optimized C functions which operate on a 3x3 "
    "matrix representation of a prallelepiped cell.",
    -1,
    method_def,
};

PyMODINIT_FUNC PyInit__cell(void) {
  PyObject *module = PyModule_Create(&module_def);
  if (module == NULL) {
    return NULL;
  }
  return module;
}