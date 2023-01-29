/*
 * A python wrapper for various methods for numerically integrating ODEs.
 *
 * Vincent Lovero
 *
 * [1] E. Hairer, S. P. Norsett G. Wanner, “Solving Ordinary Differential Equations I: Nonstiff Problems”
 * [2] E. Hairer, G. Wanner, “Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems”
 * [3] Li, Heng http://lh3lh3.users.sourceforge.net/download/lsoda.c
 * [4] Verner, J H (2013). Explicit Runge Kutta pairs with lower stage-order. Numerical Algorithms 65(3): 555–577
 * [5] P. Bogacki, L.F. Shampine, “A 3(2) Pair of Runge-Kutta Formulas”, Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
 *
 * TODO: 1 - add scalar version to module
 *       2 - clean up std::cout's in code
 *       3 - make steppers for LSODA and Radau
 */


#define NPY_NO_DEPRECATED_API NPY_1_22_API_VERSION
#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"

#include "rk12.h"
#include "rk23.h"
#include "verner65e.h"
#include "verner98e.h"
#include "lsoda.h"
#include "radau5.h"

 // Create Fortran contiguous numpy array from python object (increase ref count if already farray)
 // will be used later when support for user jacobian is added
#define PyArray_FContiguousFromObject(op,type,min_depth,max_depth) PyArray_FromAny(op, PyArray_DescrFromType(type), min_depth, max_depth, NPY_ARRAY_FARRAY | NPY_ARRAY_ENSUREARRAY, NULL)


 // declare lapack functions that will be retrieved from scipy
void (*sgetrf)(int *, int *, float *, int *, int *, int *);
void (*dgetrf)(int *, int *, double *, int *, int *, int *);
void (*cgetrf)(int *, int *, std::complex<float> *, int *, int *, int *);
void (*zgetrf)(int *, int *, std::complex<double> *, int *, int *, int *);

void (*sgetrs)(char *, int *, int *, const float *, int *, const int *, float *, int *, int *);
void (*dgetrs)(char *, int *, int *, const double *, int *, const int *, double *, int *, int *);
void (*cgetrs)(char *, int *, int *, const std::complex<float> *, int *, const int *, std::complex<float> *, int *, int *);
void (*zgetrs)(char *, int *, int *, const std::complex<double> *, int *, const int *, std::complex<double> *, int *, int *);


enum ODEMethods
{
    LSODA = 0,
    RKV65E,
    RK23,
    RKV98E,
    RADAU5,
    RK12
};


void destruct(PyObject *obj)
{
    void *stepper = PyCapsule_GetPointer(obj, NULL);
    if (stepper) {
        delete stepper;
    }
}

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
        *ptr = (c && obj) ? getPointer(c) : nullptr;
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

    if (PyErr_Occurred()) {
        return false;
    }

    return true;
}

template <typename funcT>
inline ODEResult<double> callMethod(int method, funcT f, double *y0, double *t, double *y, double min_step, double max_step, double h0, double rtol, double atol, size_t m, double *sphere, PyObject *pyargs, size_t n)
{
    ODEResult<double> result;

    switch (method) {
        case ODEMethods::LSODA:
            result = lsoda::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case ODEMethods::RKV65E:
            result = verner65e::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case ODEMethods::RK23:
            result = bogackiShampine::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case ODEMethods::RKV98E:
            result = verner98e::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case ODEMethods::RADAU5:
            result = radau5iia::integrate<double, decltype(f),
                radau5iia::factorLU<double, dgetrf>,
                radau5iia::factorLU<std::complex<double>, zgetrf>,
                radau5iia::solveLU<double, dgetrs>,
                radau5iia::solveLU<std::complex<double>, zgetrs>>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        case ODEMethods::RK12:
            result = Fehlberg12::integrate<double, decltype(f)>(f, y0, t, y, min_step, max_step, h0, rtol, atol, n, m, sphere, (void *)pyargs);
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Not a valid method");
            break;
    }
    return result;
}

template <typename funcT>
inline void callStepper(void **stepper, int method, funcT f, double *y0, double tn, double tf, double *y, double min_step, double max_step, double h0, double rtol, double atol, bool firstStep, double *sphere, double *work, PyObject *pyargs, size_t n)
{

    switch (method) {
        case ODEMethods::LSODA:
            PyErr_SetString(PyExc_NotImplementedError, "LSODA stepper is not implemented");
            break;
        case ODEMethods::RKV65E:
        {
            if (*stepper == nullptr) {
                *stepper = new verner65e::stepper64_t<funcT>;
            }
            ((verner65e::stepper64_t<funcT> *)(*stepper))->step(f, tn, tf, y0, y, atol, rtol, min_step, max_step, h0, n, firstStep, sphere, work, (void *)pyargs);
            break;
        }
        case ODEMethods::RK23:
        {
            if (*stepper == nullptr) {
                *stepper = new bogackiShampine::stepper64_t<funcT>;
            }
            ((bogackiShampine::stepper64_t<funcT> *)(*stepper))->step(f, tn, tf, y0, y, atol, rtol, min_step, max_step, h0, n, firstStep, sphere, work, (void *)pyargs);
            break;
        }
        case ODEMethods::RKV98E:
        {
            if (*stepper == nullptr) {
                *stepper = new verner98e::stepper64_t<funcT>;
            }
            ((verner98e::stepper64_t<funcT> *)(*stepper))->step(f, tn, tf, y0, y, atol, rtol, min_step, max_step, h0, n, firstStep, sphere, work, (void *)pyargs);
            break;
        }
        case ODEMethods::RADAU5:
            PyErr_SetString(PyExc_NotImplementedError, "RADAU5 stepper is not implemented");
            break;
        case ODEMethods::RK12:
        {
            if (*stepper == nullptr) {
                *stepper = new Fehlberg12::stepper64_t<funcT>;
            }
            ((Fehlberg12::stepper64_t<funcT> *)(*stepper))->step(f, tn, tf, y0, y, atol, rtol, min_step, max_step, h0, n, firstStep, sphere, work, (void *)pyargs);
            break;
        }
        default:
            PyErr_SetString(PyExc_ValueError, "Not a valid method");
            break;
    }
}

inline size_t getWorkArraySize(int method, size_t node)
{
    switch (method) {
        case ODEMethods::LSODA:
            PyErr_SetString(PyExc_NotImplementedError, "LSODA stepper is not implemented");
            break;
        case ODEMethods::RKV65E:
            return verner65e::getWorkArraySize(node);
        case ODEMethods::RK23:
            return bogackiShampine::getWorkArraySize(node);
        case ODEMethods::RKV98E:
            return verner98e::getWorkArraySize(node);
        case ODEMethods::RADAU5:
            PyErr_SetString(PyExc_NotImplementedError, "RADAU5 stepper is not implemented");
            break;
        case ODEMethods::RK12:
            return Fehlberg12::getWorkArraySize(node);
        default:
            PyErr_SetString(PyExc_ValueError, "Not a valid method");
            break;
    }
    return 0;
}

PyObject *PyParseWork(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *kwlist[] = { "work", "method", "node", NULL };

    PyObject *pyWork;
    int method;
    npy_intp node;

    npy_intp nwork;
    size_t s_star;
    PyArrayObject *npyWork = nullptr;
    double *work;
    double *yn;
    double *yn1;
    double *dy;
    double *temp;
    double *K;
    double *Q;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OiL", const_cast<char **>(kwlist), &pyWork, &method, &node)) {
        if (PyErr_Occurred()) {
            return NULL;
        }
        PyErr_SetString(PyExc_RuntimeError, "Unable to parse args");
        return NULL;
    }

    nwork = getWorkArraySize(method, node);
    if (!PyArray_Check(pyWork)) {
        PyErr_Format(PyExc_TypeError, "argument `work` is of type '%s' not np.ndarray", Py_TYPE(pyWork)->tp_name);
        return NULL;
    }

    npyWork = (PyArrayObject *)pyWork;
    if (nwork != PyArray_SIZE(npyWork)) {
        PyErr_Format(PyExc_ValueError, "`work` array size does not match required size for given method (%ld != %ld)", PyArray_SIZE(npyWork), nwork);
        return NULL;
    }

    switch (method) {
        case ODEMethods::LSODA:
            PyErr_SetString(PyExc_NotImplementedError, "LSODA stepper is not implemented");
            return NULL;
        case ODEMethods::RKV65E:
            s_star = 10;
            break;
        case ODEMethods::RK23:
            s_star = 4;
            break;
        case ODEMethods::RKV98E:
            s_star = 21;
            break;
        case ODEMethods::RADAU5:
            PyErr_SetString(PyExc_NotImplementedError, "RADAU5 stepper is not implemented");
            return NULL;
        case ODEMethods::RK12:
            s_star = 3;
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Not a valid method");
            return NULL;
    }

    work = (double *)PyArray_DATA(npyWork);
    yn = work;
    yn1 = work + node;
    dy = yn1 + node;
    temp = dy + node; // temp should contain actual yn on return
    K = temp + node;
    Q = K + (node * s_star);
    npy_intp shape[2] = { (npy_intp)s_star, node };

    return Py_BuildValue("{s:O,s:O,s:O}",
        "yn", PyArray_SimpleNewFromData(1, &node, NPY_DOUBLE, temp),
        "yn1", PyArray_SimpleNewFromData(1, &node, NPY_DOUBLE, yn1),
        "K", PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, K));
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
    ODEResult<double> result;
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
        PyErr_Format(PyExc_ValueError, "shape mismatch output array and time array (%li != %li)", PyArray_DIM((PyArrayObject *)npyY, 0), PyArray_DIM((PyArrayObject *)npyT, 0));
        goto error;
    }
    else if ((PyArray_DIM((PyArrayObject *)npyY0, 0) + 1) != PyArray_DIM((PyArrayObject *)npyW, 0)) {
        PyErr_Format(PyExc_ValueError, "window must have shape (n + 1, ), (%li != %li)", (PyArray_DIM((PyArrayObject *)npyY0, 0) + 1), PyArray_DIM((PyArrayObject *)npyW, 0));
        goto error;
    }
    else if (PyArray_DIM((PyArrayObject *)npyY0, 0) != PyArray_DIM((PyArrayObject *)npyY, 1)) {
        PyErr_Format(PyExc_ValueError, "shape mismatch output array and y0 array (%li != %li)", PyArray_DIM((PyArrayObject *)npyY0, 0), PyArray_DIM((PyArrayObject *)npyY, 1));
        goto error;
    }

    // get pointers
    t = (double *)PyArray_DATA((PyArrayObject *)npyT);
    y = (double *)PyArray_DATA((PyArrayObject *)npyY);
    y0 = (double *)PyArray_DATA((PyArrayObject *)npyY0);
    w = (double *)PyArray_DATA((PyArrayObject *)npyW);
    m = shape[0];
    n = shape[1];

    // is there a cleaner way of writing this??
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

            PyObject *err = PyErr_Occurred();
            if (err) {
                // fill with nan to stop integrator from proceeding
                for (size_t i = 0; i < n; i++) {
                    o[i] = std::numeric_limits<double>::quiet_NaN();
                }
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
        std::string message;
        switch (result.status) {
            case ODEExitCodes::success:
                message = "success";
                break;
            case ODEExitCodes::outOfBounds:
                message = "solution left computation window during between intervals";
                break;
            case ODEExitCodes::failedInterpolation:
                message = "interpolant left computation window";
                break;
            default:
                PyErr_Format(PyExc_ValueError, "Unknown exit code %lu", (unsigned long)(result.status));
                Py_XDECREF(npyY);
                return NULL;
        }
        return Py_BuildValue("Ol{s:l,s:l,s:l,s:l,s:l,s:s,s:d}", npyY, (long)result.k,
            "feval", (long)result.feval,
            "jeval", (long)result.jeval,
            "fails", (long)result.fails,
            "steps", (long)result.steps,
            "numlu", (long)result.numLU,
            "message", message.c_str(),
            "hout", result.stepOnExit);
    }
}

PyObject *PyODEStep(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *kwlist[] = { "f", "y0", "t0", "tf", "args", "out", "isvoid",
                             "atol", "rtol", "min_step", "max_step",
                             "window", "method", "full_output", "h0",
                             "first_step", "work", "stepper", NULL };

    PyObject *pyf = nullptr;
    PyObject *pyY = nullptr;
    PyObject *pyargs = nullptr;
    PyObject *pyW = nullptr;
    PyObject *pyWork = nullptr;
    PyObject *pyY0 = nullptr;
    PyObject *npyY0 = nullptr;
    PyObject *npyY = nullptr;
    PyObject *npyW = nullptr;
    PyObject *npyWork = nullptr;
    PyObject *pyStepper = nullptr;

    npy_intp node, nwork;
    double t0, tf;
    double *y = nullptr, *y0 = nullptr, *w = nullptr, *work = nullptr;
    double atol = 1.49012e-8, rtol = 1.49012e-8, min_step = 1e-6, h0 = 0.0, max_step = std::numeric_limits<double>::infinity();
    int isvoid = 0;
    ODEResult<double> result;
    void *stepper = nullptr;
    int method = ODEMethods::RKV65E, full_output = 0, first_step = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOdd|OOpddddOipdpOO", const_cast<char **>(kwlist),
        &pyf, &pyY0, &t0, &tf, &pyargs, &pyY, &isvoid, &atol, &rtol, &min_step, &max_step, &pyW, &method, &full_output, &h0, &first_step, &pyWork, &pyStepper)) {
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

    node = PyArray_SIZE((PyArrayObject *)npyY0);
    if (!pyY || (pyY == Py_None)) {
        npyY = PyArray_SimpleNew(1, &node, NPY_DOUBLE);
    }
    else {
        npyY = PyArray_ContiguousFromObject(pyY, NPY_DOUBLE, 0, 0);
        if (!npyY) {
            PyErr_SetString(PyExc_MemoryError, "continguousFromObject returned NULL");
            goto error;
        }
    }

    if (!pyW) {
        node += 1;
        npyW = PyArray_SimpleNew(1, &node, NPY_DOUBLE);
        std::memset(PyArray_DATA((PyArrayObject *)npyW), 0, sizeof(double) * node);
        node -= 1;
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
    if ((PyArray_DIM((PyArrayObject *)npyY0, 0) + 1) != PyArray_DIM((PyArrayObject *)npyW, 0)) {
        PyErr_Format(PyExc_ValueError, "window must have shape (n + 1, ), (%li != %li)", (PyArray_DIM((PyArrayObject *)npyY0, 0) + 1), PyArray_DIM((PyArrayObject *)npyW, 0));
        goto error;
    }
    else if (PyArray_DIM((PyArrayObject *)npyY0, 0) != PyArray_DIM((PyArrayObject *)npyY, 0)) {
        PyErr_Format(PyExc_ValueError, "shape mismatch output array and y0 array (%li != %li)", PyArray_DIM((PyArrayObject *)npyY0, 0), PyArray_DIM((PyArrayObject *)npyY, 1));
        goto error;
    }

    // get work array
    nwork = getWorkArraySize(method, node);
    if (!pyWork) {
        npyWork = PyArray_SimpleNew(1, &nwork, NPY_DOUBLE);
    }
    else {
        npyWork = PyArray_ContiguousFromObject(pyWork, NPY_DOUBLE, 0, 0);
        if (PyArray_SIZE((PyArrayObject *)npyWork) < nwork) {
            Py_XDECREF(npyWork);
            npyWork = PyArray_SimpleNew(1, &nwork, NPY_DOUBLE);
        }
    }

    // get stepper from PyCapsule
    if ((pyStepper != nullptr) && (pyStepper != Py_None)) {
        if (!PyCapsule_CheckExact(pyStepper)) {
            PyErr_Format(PyExc_TypeError, "supplied argument must be of type PyCapsule (%s != PyCapsule)", Py_TYPE(pyStepper)->tp_name);
            goto error;
        }
        else {
            stepper = PyCapsule_GetPointer(pyStepper, "__PyODEStep_local_stepper__");
            if (!stepper) {
                goto error;
            }
        }
    }
    if (first_step) {
        // don't use stepper when first_step is set to True
        stepper = nullptr;
    }

    // get pointers
    y = (double *)PyArray_DATA((PyArrayObject *)npyY);
    y0 = (double *)PyArray_DATA((PyArrayObject *)npyY0);
    w = (double *)PyArray_DATA((PyArrayObject *)npyW);
    work = (double *)PyArray_DATA((PyArrayObject *)npyWork);

    // is there a cleaner way of writing this??
    if (isvoid) {
        auto f = [=](const double t, const RP(double) z, RP(double) o, const RP(void) a) {
            PyObject *pyz;
            PyObject *pyo;
            pyz = PyArray_SimpleNewFromData(1, &node, NPY_DOUBLE, (void *)z);
            pyo = PyArray_SimpleNewFromData(1, &node, NPY_DOUBLE, (void *)o);
            PyArray_CLEARFLAGS((PyArrayObject *)pyz, NPY_ARRAY_OWNDATA);
            PyArray_CLEARFLAGS((PyArrayObject *)pyo, NPY_ARRAY_OWNDATA);

            if (a != nullptr) {
                Py_XDECREF(PyObject_CallFunction(pyf, "dOOO", t, pyz, pyo, (PyObject *)a));
            }
            else {
                Py_XDECREF(PyObject_CallFunction(pyf, "dOO", t, pyz, pyo));
            }

            PyObject *err = PyErr_Occurred();
            if (err) {
                // fill with nan to stop integrator from proceeding
                for (size_t i = 0; i < node; i++) {
                    o[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }

            Py_XDECREF(pyz);
            Py_XDECREF(pyo);
        };

        callStepper(&stepper, method, f, y0, t0, tf, y, min_step, max_step, h0, rtol, atol, first_step, w, work, pyargs, node);
        if (PyErr_Occurred()) {
            goto error;
        }
        result = ((verner65e::stepper64_t<decltype(f)> *)stepper)->toODEResult();
    }
    else {
        auto f = [=](const double t, const RP(double) z, RP(double) o, const RP(void) a) {
            PyObject *pyz = nullptr;
            PyObject *pyo = nullptr;
            PyArrayObject *pya = nullptr;
            pyz = PyArray_SimpleNewFromData(1, &node, NPY_DOUBLE, (void *)z);
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
                for (size_t i = 0; i < node; i++) {
                    o[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }

            Py_XDECREF(pyo);
            Py_XDECREF(pya);
            Py_XDECREF(pyz);
        };

        callStepper(&stepper, method, f, y0, t0, tf, y, min_step, max_step, h0, rtol, atol, first_step, w, work, pyargs, node);
        if (PyErr_Occurred()) {
            goto error;
        }
        result = ((verner65e::stepper64_t<decltype(f)> *)stepper)->toODEResult();
    }

error:
    Py_XDECREF(npyY0);
    Py_XDECREF(npyW);

    if (PyErr_Occurred()) {
        Py_XDECREF(npyY);
        Py_XDECREF(npyWork);
        return NULL;
    }

    if (!full_output) {
        delete stepper;
        Py_XDECREF(npyWork);
        return npyY;
    }
    else {

        std::string message;
        switch (result.status) {
            case ODEExitCodes::success:
                message = "success";
                break;
            case ODEExitCodes::outOfBounds:
                message = "solution left com putation window during between intervals";
                break;
            case ODEExitCodes::failedInterpolation:
                message = "interpolant left computation window";
                break;
            default:
                PyErr_Format(PyExc_ValueError, "Unknown exit code %lu", (unsigned long)(result.status));
                Py_XDECREF(npyY);
                return NULL;
        }
        return Py_BuildValue("O{s:l,s:l,s:l,s:l,s:l,s:s,s:d,s:O,s:O}", npyY,
            "feval", (long)result.feval,
            "jeval", (long)result.jeval,
            "fails", (long)result.fails,
            "steps", (long)result.steps,
            "numlu", (long)result.numLU,
            "message", message.c_str(),
            "hout", result.stepOnExit,
            "work", npyWork,
            "stepper", PyCapsule_New(stepper, "__PyODEStep_local_stepper__", destruct));
    }
}

PyObject *parseStepper(PyObject *self, PyObject *arg)
{
    // any namespace and function type should work for the cast
    using stepper_t = verner65e::stepper64_t<void (*)(double, double *, double *, const void *)>;

    stepper_t *stepper;
    if (!PyCapsule_CheckExact(arg)) {
        PyErr_Format(PyExc_TypeError, "supplied argument must be of type PyCapsule (%s != PyCapsule)", Py_TYPE(arg)->tp_name);
        return nullptr;
    }

    stepper = (stepper_t *)PyCapsule_GetPointer(arg, "__PyODEStep_local_stepper__");
    if (stepper) {
        return Py_BuildValue("{s:d,s:d,s:d,s:k,s:k,s:k,s:k,s:O,s:k}",
            "t", stepper->get_t(),
            "h", stepper->get_h(),
            "hprev", stepper->get_hprev(),
            "fails", stepper->get_fails(),
            "feval", stepper->get_feval(),
            "steps", stepper->get_steps(),
            "jeval", stepper->get_jeval(),
            "interpolate", (stepper->get_interpolate() ? Py_True : Py_False),
            "status", (size_t)(stepper->get_status()));
    }
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError, "PyCapsule_GetPointer returned NULL");
    }
    return nullptr;
}


PyDoc_STRVAR(
    PySolveIVP_doc,
    "solveIVP(f, y0, teval, args=None, out=None, isvoid=False, atol=1.49012e-8, rtol=1.49012e-8, min_step=1e-6, max_step=None, window=None, method=0, full_output=False)\n"
    "--\n"
    "\n"
    "Integrate a system of ODEs.\n"
    "\n"
    "    Parameters\n"
    "    ----------\n"
    "    f : Callable\n"
    "        RHS of system\n"
    "    y0 : ArrayLike\n"
    "        _Initial conditions\n"
    "    teval : ArrayLike\n"
    "        points to evaluate solution at (must include initial time point t0)\n"
    "    args : tuple, optional\n"
    "        extra arguments for f(t, y), by default None\n"
    "    out : ndarray, optional\n"
    "        preallocated output array, by default None\n"
    "    isvoid : bool, optional\n"
    "        f does not return anything (performance booster), by default False\n"
    "    atol : float, optional\n"
    "        absolute error tolerance, by default 1.49012e-8\n"
    "    rtol : float, optional\n"
    "        relative error tolerance, by default 1.49012e-8\n"
    "    min_step : float, optional\n"
    "        minimum allowable step, by default 1e-6\n"
    "    max_step : float, optional\n"
    "        maximum allowable step, by default inf\n"
    "    window : ArrayLike, optional\n"
    "        computational window solution must stay inside, by default None\n"
    "    method : int, optional\n"
    "        which integration method to use, by default methods.LSODA ( == 0)\n"
    "    full_output : bool, optional\n"
    "        return dict with integration stats, by default False\n"
    "\n"
    "    Returns\n"
    "    -------\n"
    "    Tuple[ndarray, int]\n"
    "        ndarray with shape (len(teval), len(y0)) containing solution\n"
    "        and number of successful time steps taken.\n"
    "    Tuple[ndarray, int, Dict[str, Any]]\n"
    "        If full_output == True, also return dictionary with integtation stats.\n");


static PyMethodDef moduleMethods[] = {
    {"solveIVP", (PyCFunction)PySolveIVP, METH_VARARGS | METH_KEYWORDS, PySolveIVP_doc},
    {"ODEStep", (PyCFunction)PyODEStep, METH_VARARGS | METH_KEYWORDS, PySolveIVP_doc},
    {"parse_stepper", parseStepper, METH_O, "parse stepper from capsule"},
    {"parse_work", (PyCFunction)PyParseWork, METH_VARARGS | METH_KEYWORDS, "parse work array into dict"},
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