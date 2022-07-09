#define NPY_NO_DEPRECATED_API NPY_1_22_API_VERSION
#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"

#include "rk23.h"
#include "verner65e.h"
#include "verner98e.h"
#include "lsoda.h"
#include "radau5.h"

#define PyArray_FContiguousFromObject(op,type,min_depth,max_depth) PyArray_FromAny(op, PyArray_DescrFromType(type), min_depth, max_depth, NPY_ARRAY_FARRAY | NPY_ARRAY_ENSUREARRAY, NULL)

#if defined(DEBUG)

#define LINE() std::cerr << "[INFO]: " << __LINE__ << '\n'

#else

#define LINE()

#endif

// declare lapack functions that will be retrieved from scipy
void (*sgetrf)(int *, int *, float *, int *, int *, int *);
void (*dgetrf)(int *, int *, double *, int *, int *, int *);
void (*cgetrf)(int *, int *, std::complex<float> *, int *, int *, int *);
void (*zgetrf)(int *, int *, std::complex<double> *, int *, int *, int *);

void (*sgetrs)(char *, int *, int *, const float *, int *, const int *, float *, int *, int *);
void (*dgetrs)(char *, int *, int *, const double *, int *, const int *, double *, int *, int *);
void (*cgetrs)(char *, int *, int *, const std::complex<float> *, int *, const int *, std::complex<float> *, int *, int *);
void (*zgetrs)(char *, int *, int *, const std::complex<double> *, int *, const int *, std::complex<double> *, int *, int *);


void *getPointer(PyObject *obj)
{
    if (PyCapsule_CheckExact(obj)) {
        return PyCapsule_GetPointer(obj, NULL);
    }
    PyErr_SetString(PyExc_ValueError, "Not an object containing a void ptr");
    return NULL;
}

bool importLapackRoutines()
{
    // scipy.linalg.blas.dgemm._cpointer
    PyObject *scipyModule = nullptr;
    PyObject *d = nullptr;
    PyObject *obj = nullptr;
    PyObject *c = nullptr;

    auto getPointerFromName = [&](const char *name, void **ptr) {
        obj = PyDict_GetItemString(d, name);
        c = PyObject_GetAttrString(obj, "_cpointer");
        *ptr = getPointer(c);
        Py_XDECREF(obj);
        Py_XDECREF(c);
    };

    scipyModule = PyImport_ImportModule("scipy.linalg.lapack");
    if (!scipyModule) {
        return false;
    }
    d = PyModule_GetDict(scipyModule);
    if (!d) {
        Py_XDECREF(scipyModule);
        return false;
    }

    getPointerFromName("sgetrf", (void **)(&sgetrf));
    getPointerFromName("dgetrf", (void **)(&dgetrf));
    getPointerFromName("cgetrf", (void **)(&cgetrf));
    getPointerFromName("zgetrf", (void **)(&zgetrf));

    getPointerFromName("sgetrs", (void **)(&sgetrs));
    getPointerFromName("dgetrs", (void **)(&dgetrs));
    getPointerFromName("cgetrs", (void **)(&cgetrs));
    getPointerFromName("zgetrs", (void **)(&zgetrs));

    return true;
}

template <typename funcT>
inline ODEResult callMethod(int method, funcT f, double *y0, double *t, double *y, double min_step, double max_step, double h0, double rtol, double atol, size_t m, double *sphere, PyObject *pyargs, size_t n)
{
    ODEResult result;
    LINE();
    switch (method) {
        case 0:
            // LSODA
            result = lsoda::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case 1:
            // RKV65E
            result = verner65e::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case 2:
            // RK23
            result = bogakiShampine::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case 3:
            // RKV98E
            result = verner98e::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case 5:
            // radau5IIA
            result = radau5iia::integrate<double, decltype(f),
                radau5iia::factorLU<double, dgetrf>,
                radau5iia::factorLU<std::complex<double>, zgetrf>,
                radau5iia::solveLU<double, dgetrs>,
                radau5iia::solveLU<std::complex<double>, zgetrs>>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Not a valid method");
            break;
    }
    return result;
}

PyObject *PySolveIVP(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *kwlist[] = { "f", "y0", "teval", "args", "out", "isvoid",
                             "atol", "rtol", "min_step", "max_step",
                             "window", "method", "full_output", "h0", NULL };

    PyObject *pyf = nullptr;
    PyObject *pyT = nullptr;
    PyObject *pyY = nullptr;
    PyObject *pyargs = nullptr;
    PyObject *pyW = nullptr;
    PyObject *pyY0 = nullptr;
    PyObject *npyY0 = nullptr;
    PyObject *npyT = nullptr;
    PyObject *npyY = nullptr;
    PyObject *npyW = nullptr;

    npy_intp shape[2];
    double *t = nullptr, *y = nullptr, *y0 = nullptr, *w = nullptr;
    double atol = 1.49012e-8, rtol = 1.49012e-8, min_step = 1e-6, h0 = 0.0, max_step = std::numeric_limits<double>::infinity();
    int isvoid = 0;
    size_t n, m;
    ODEResult result;
    int method = 0, full_output = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|OOpddddOipd", const_cast<char **>(kwlist), &pyf, &pyY0, &pyT, &pyargs, &pyY, &isvoid, &atol, &rtol, &min_step, &max_step, &pyW, &method, &full_output, &h0)) {
        if (PyErr_Occurred()) {
            return NULL;
        }
        PyErr_SetString(PyExc_RuntimeError, "Unable to parse args");
        return NULL;
    }

    // create required numpy arrays
    npyY0 = PyArray_ContiguousFromObject(pyY0, NPY_DOUBLE, 0, 0);
    if (!npyY0) {
        PyErr_SetString(PyExc_MemoryError, "continguousFromObject returned NULL");
        goto error;
    }
    npyT = PyArray_ContiguousFromObject(pyT, NPY_DOUBLE, 0, 0);
    if (!npyT) {
        PyErr_SetString(PyExc_MemoryError, "continguousFromObject returned NULL");
        goto error;
    }

    shape[0] = PyArray_SIZE((PyArrayObject *)npyT);
    shape[1] = PyArray_SIZE((PyArrayObject *)npyY0);
    if (!pyY) {
        npyY = PyArray_SimpleNew(2, shape, NPY_DOUBLE);
    }
    else {
        npyY = PyArray_ContiguousFromObject(pyY, NPY_DOUBLE, 0, 0);
        if (!npyY) {
            PyErr_SetString(PyExc_MemoryError, "continguousFromObject returned NULL");
            goto error;
        }
    }

    if (!pyW) {
        shape[1] += 1;
        npyW = PyArray_SimpleNew(1, &shape[1], NPY_DOUBLE);
        std::memset(PyArray_DATA((PyArrayObject *)npyW), 0, sizeof(double) * shape[1]);
        shape[1] -= 1;
        *((double *)PyArray_DATA((PyArrayObject *)npyW)) = std::numeric_limits<double>::infinity();
    }
    else {
        npyW = PyArray_ContiguousFromObject(pyW, NPY_DOUBLE, 0, 0);
        if (!npyW) {
            PyErr_SetString(PyExc_MemoryError, "continguousFromObject returned NULL");
            goto error;
        }
    }

    // check for shape errors
    if (PyArray_DIM((PyArrayObject *)npyY, 0) != PyArray_DIM((PyArrayObject *)npyT, 0)) {
        PyErr_SetString(PyExc_ValueError, "shape mismatch output array and time array");
        goto error;
    }
    else if ((PyArray_DIM((PyArrayObject *)npyY0, 0) + 1) != PyArray_DIM((PyArrayObject *)npyW, 0)) {
        PyErr_SetString(PyExc_ValueError, "window must have shape (n + 1, )");
        goto error;
    }
    else if (PyArray_DIM((PyArrayObject *)npyY0, 0) != PyArray_DIM((PyArrayObject *)npyY, 1)) {
        PyErr_SetString(PyExc_ValueError, "shape mismatch output array and y0 array");
        goto error;
    }

    // get pointers
    t = (double *)PyArray_DATA((PyArrayObject *)npyT);
    y = (double *)PyArray_DATA((PyArrayObject *)npyY);
    y0 = (double *)PyArray_DATA((PyArrayObject *)npyY0);
    w = (double *)PyArray_DATA((PyArrayObject *)npyW);
    m = shape[0];
    n = shape[1];

    if (isvoid) {
        auto f = [=](const double t, const RP(double) z, RP(double) o, const RP(void) a) {
            PyObject *pyz;
            PyObject *pyo;
            npy_intp dims[1] = { (npy_intp)n };
            pyz = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)z);
            pyo = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)o);
            PyArray_CLEARFLAGS((PyArrayObject *)pyz, NPY_ARRAY_OWNDATA);
            PyArray_CLEARFLAGS((PyArrayObject *)pyo, NPY_ARRAY_OWNDATA);

            if (a != nullptr) {
                Py_XDECREF(PyObject_CallFunction(pyf, "dOOO", t, pyz, pyo, (PyObject *)a));
            }
            else {
                Py_XDECREF(PyObject_CallFunction(pyf, "dOO", t, pyz, pyo));
            }

            Py_XDECREF(pyz);
            Py_XDECREF(pyo);
        };

        result = callMethod(method, f, y0, t, y, min_step, max_step, h0, rtol, atol, m, w, pyargs, n);
    }
    else {
        auto f = [=](const double t, const RP(double) z, RP(double) o, const RP(void) a) {
            PyObject *pyz = nullptr;
            PyObject *pyo = nullptr;
            PyArrayObject *pya = nullptr;
            npy_intp dims[1] = { (npy_intp)n };
            pyz = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)z);
            PyArray_CLEARFLAGS((PyArrayObject *)pyz, NPY_ARRAY_OWNDATA);

            if (a != nullptr) {
                pyo = PyObject_CallFunction(pyf, "dOO", t, pyz, (PyObject *)a);
            }
            else {
                pyo = PyObject_CallFunction(pyf, "dO", t, pyz);
            }

            PyObject *err = PyErr_Occurred();
            if (pyo && (err == nullptr)) {
                pya = (PyArrayObject *)PyArray_ContiguousFromObject(pyo, NPY_DOUBLE, 0, 0);
                std::memcpy(o, PyArray_DATA(pya), PyArray_SIZE(pya) * sizeof(double));
            }
            else {
                for (size_t i = 0; i < n; i++) {
                    o[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }

            Py_XDECREF(pyo);
            Py_XDECREF(pya);
            Py_XDECREF(pyz);
        };

        result = callMethod(method, f, y0, t, y, min_step, max_step, h0, rtol, atol, m, w, pyargs, n);
    }

error:
    Py_XDECREF(npyT);
    Py_XDECREF(npyY0);
    Py_XDECREF(npyW);

    if (PyErr_Occurred()) {
        Py_XDECREF(npyY);
        return NULL;
    }

    if (!full_output) {
        return Py_BuildValue("Ol", npyY, (long)result.k);
    }
    else {
        return Py_BuildValue("Ol{s:l,s:l,s:l,s:l}", npyY, (long)result.k,
            "feval", (long)result.feval,
            "jeval", (long)result.jeval,
            "fails", (long)result.fails,
            "steps", (long)result.steps);
    }
}


PyDoc_STRVAR(
    PySolveIVP_doc,
    "solveIVP(f, y0, teval, args=None, out=None, isvoid=False, atol=1e-6, rtol=1e-3, min_step=1e-6, max_step=1, window=None, method=1, full_output=False)\n"
    "--\n"
    "\n"
    "Integrate IVP and evaluate at user supplied points.\n"
    "Arguments: (timeout, flags=None)\n"
    "Doc blahblah doc doc doc.");


static PyMethodDef moduleMethods[] = {
    {"solveIVP", (PyCFunction)PySolveIVP, METH_VARARGS | METH_KEYWORDS, PySolveIVP_doc},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pyodepack",
    NULL,
    -1,
    moduleMethods };

extern "C" PyObject * PyInit__pyodepack(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    import_array();

    if (m == NULL) {
        return NULL;
    }

    if (!importLapackRoutines()) {
        Py_XDECREF(m);
        return NULL;
    }

    return m;
}