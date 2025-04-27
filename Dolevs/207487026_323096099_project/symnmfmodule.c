#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "./symnmf.h"




matrix parse_pylist_to_matrix(PyObject* py_matrix) {
    matrix mat;
    int i, j, rows, cols;

    if (!PyList_Check(py_matrix)) {
        return empty_matrix();
    }
    rows = PyList_Size(py_matrix);
    if (rows == 0) {
        return empty_matrix();
    }

    PyObject* first_row = PyList_GetItem(py_matrix, 0);
    if (!PyList_Check(first_row)) {
        return empty_matrix();
    }
    cols = PyList_Size(first_row);

    mat = init_zero_struct_matrix(rows, cols);
    if (!is_valid_matrix(&mat)) {
        return empty_matrix();
    }

    for (i = 0; i < rows; i++) {
        PyObject* row = PyList_GetItem(py_matrix, i);
        if (!PyList_Check(row)) {
            free_matrix_ptr(&mat);
            return empty_matrix();
        }
        for (j = 0; j < cols; j++) {
            PyObject* item = PyList_GetItem(row, j);
            if (!item || !PyFloat_Check(item)) {
                free_matrix_ptr(&mat);
                return empty_matrix();
            }
            mat.mat[i][j] = PyFloat_AsDouble(item);
        }
    }
    return mat;
}

PyObject* matrix_to_pylist(matrix* mat_ptr) {
    int i, j;
    PyObject* py_mat = PyList_New(mat_ptr->rows);
    if (!py_mat) return NULL;

    for (i = 0; i < mat_ptr->rows; i++) {
        PyObject* row = PyList_New(mat_ptr->cols);
        if (!row) {
            for (int k = 0; k < i; k++) {
                Py_DECREF(PyList_GetItem(py_mat, k));
            }
            Py_DECREF(py_mat);
            return NULL;
        }
        for (j = 0; j < mat_ptr->cols; j++) {
            PyObject* val = PyFloat_FromDouble(mat_ptr->mat[i][j]);
            if (!val) {
                for (int k = 0; k < j; k++) {
                    Py_DECREF(PyList_GetItem(row, k));
                }
                Py_DECREF(row);
                for (int k = 0; k < i; k++) {
                    Py_DECREF(PyList_GetItem(py_mat, k));
                }
                Py_DECREF(py_mat);
                return NULL;
            }
            PyList_SetItem(row, j, val);  
        }
        PyList_SetItem(py_mat, i, row);  
    }
    return py_mat;
}

static PyObject* py_sym(PyObject *self, PyObject *args) {
    PyObject *py_data;
    if (!PyArg_ParseTuple(args, "O", &py_data)) {
        PyErr_SetString(PyExc_TypeError, "Expected one argument: a list of lists (matrix).");
        return NULL;
    }

    matrix input = parse_pylist_to_matrix(py_data);
    if (!is_valid_matrix(&input)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid input matrix.");
        return NULL;
    }

    matrix result = sym(input);
    free_matrix_ptr(&input);
    PyObject* py_result = matrix_to_pylist(&result);
    free_matrix_ptr(&result);
    return py_result;
}

static PyObject* py_ddg(PyObject *self, PyObject *args) {
    PyObject *py_data;
    if (!PyArg_ParseTuple(args, "O", &py_data)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid input matrix.");
        return NULL;
    }

    matrix input = parse_pylist_to_matrix(py_data);
    if (!is_valid_matrix(&input)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid input matrix.");
        return NULL;
    }

    matrix result = ddg(input);
    free_matrix_ptr(&input);
        if (!is_valid_matrix(&result)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to compute DDG matrix.");
        return NULL;
    }
    PyObject* py_result = matrix_to_pylist(&result);
    free_matrix_ptr(&result);
    return py_result;
}

static PyObject* py_norm(PyObject *self, PyObject *args) {
    PyObject *py_data;
    if (!PyArg_ParseTuple(args, "O", &py_data)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to compute result matrix.");
        return NULL;
    }

    matrix input = parse_pylist_to_matrix(py_data);
    if (!is_valid_matrix(&input)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to compute result matrix.");
        return NULL;
    }

    matrix result = norm(input);
    free_matrix_ptr(&input);
    
    if (!is_valid_matrix(&result)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to compute normalized Laplacian.");
        return NULL;
    }
    PyObject* py_result = matrix_to_pylist(&result);
    free_matrix_ptr(&result);
    return py_result;
}

static PyObject* py_symnmf(PyObject *self, PyObject *args) {
    PyObject *py_H, *py_W;
    double beta, eps;
    int max_iter;

    if (!PyArg_ParseTuple(args, "OOddi", &py_H, &py_W, &beta, &eps, &max_iter)) { 
        PyErr_SetString(PyExc_TypeError, "Expected two matrices, a float, a float, and an int."); 
        return NULL;
    }

    matrix H = parse_pylist_to_matrix(py_H);
    matrix W = parse_pylist_to_matrix(py_W);
    if (!is_valid_matrix(&H) || !is_valid_matrix(&W)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid input matrices to SymNMF.");
        free_matrix_ptr(&H);
        free_matrix_ptr(&W);
        return NULL;
    }

    matrix result = convergence(H, W, beta, eps, max_iter);
    free_matrix_ptr(&H);
    free_matrix_ptr(&W);
    if (!is_valid_matrix(&result)) {
    PyErr_SetString(PyExc_RuntimeError, "SymNMF convergence failed.");
    return NULL;
    }
    PyObject* py_result = matrix_to_pylist(&result);
    free_matrix_ptr(&result);
    return py_result;
}

static PyMethodDef SymNmfMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Compute similarity matrix"},
    {"ddg", py_ddg, METH_VARARGS, "Compute degree matrix"},
    {"norm", py_norm, METH_VARARGS, "Compute normalized Laplacian"},
    {"symnmf", py_symnmf, METH_VARARGS, "Run SymNMF algorithm"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf_capi",
    NULL,
    -1,
    SymNmfMethods
};

PyMODINIT_FUNC PyInit_symnmf_capi(void) {
    return PyModule_Create(&symnmfmodule);
}
