#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "libatompack.h"

static PyObject *py_cell_contains(PyObject *self, PyObject *args);
static PyObject *py_nearest_neighbor(PyObject *self, PyObject *args);

static PyMethodDef libatompack_methods[] = {
    {"cell_contains", py_cell_contains, METH_VARARGS, "Checks if point is within cell."},
    {"nearest_neighbor", py_nearest_neighbor, METH_VARARGS, "Returns the index of and distance from the nearest neighbor."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef libatompack = {
    PyModuleDef_HEAD_INIT,
    "libatompack",
    "A pure C library for atomic structure generation.",
    -1,
    libatompack_methods,
};

PyMODINIT_FUNC Py_Init_libatompack(void)
{
    PyObject *module = PyModule_Create(&libatompack);
    if (module == NULL)
    {
        return NULL;
    }
    return module;
}