#ifndef _ODE_COMMON_H
#define _ODE_COMMON_H

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <iostream>
#include <iomanip>
#include <complex>
#include <array>


// error print macros
#define PRINT_ERROR_NON_FINITE() std::cerr << "[ERROR]: NON-FINITE VALUES ENCOUNTERED\n"
#define PRINT_ERROR_MESSAGE(s) std::cerr << "[ERROR]: " << s << '\n';

#define RP(T) T * __restrict


enum ODEExitCodes
{
    success = 0,
    outOfBounds,

};

namespace std
{
    template <typename T>
    T sign(const T x)
    {
        return (x < 0) ? -1 : 1;
    }
}

struct ODEResult
{
    size_t k;
    size_t feval;
    size_t jeval;
    size_t fails;
    size_t steps;
    size_t status;
};


// RK error norm
// Using weighted root mean square where the weight is given by
// eps_abs + eps_rel * max(|yn|, |yn1|)
template <const size_t n, typename T>
__attribute__((always_inline)) T err_norm(const T delta_y[n], const T yn[n], const T yn1[n], const T atol, const T rtol)
{
    T _yn, _yn1, val, s = 0.0;
    for (size_t i = 0; i < n; i++) {
        _yn1 = std::abs(yn1[i]);
        _yn = std::abs(yn[i]);
        val = delta_y[i] / (atol + rtol * std::max(_yn, _yn1));
        s += val * val;
    }
    return std::sqrt(s / n);
}


// RK error norm
template <typename T>
__attribute__((always_inline)) T gerr_norm(const RP(T) delta_y, const RP(T) yn, const RP(T) yn1, const T atol, const T rtol, const size_t n)
{
    T _yn, _yn1, val, s = 0.0;
    for (size_t i = 0; i < n; i++) {
        _yn1 = std::abs(yn1[i]);
        _yn = std::abs(yn[i]);
        val = delta_y[i] / (atol + rtol * std::max(_yn, _yn1));
        s += val * val;
    }
    return std::sqrt(s / n);
}

// compute eulclidean norm
template <typename T>
__attribute__((always_inline)) T enorm(const RP(T) x, const size_t n)
{
    T s = 0.0;
    for (size_t i = 0; i < n; i++) {
        s = std::fma(x[i], x[i], s);
    }
    return std::sqrt(s);
}

template <typename T>
__attribute__((always_inline)) T rmsnorm(const RP(T) x, const RP(T) y, const T atol, const T rtol, const size_t n)
{
    T s = 0.0;
    for (size_t i = 0; i < n; i++) {
        s += x[i] * x[i] / (atol + rtol * std::abs(y[i]));
    }
    return std::sqrt(s / n);
}

template <typename T>
__attribute__((always_inline)) T rmsnorm(const RP(T) x, const RP(T) scale, const size_t n)
{
    T s = 0.0, val;
    for (size_t i = 0; i < n; i++) {
        val = x[i] / scale[i];
        s += val * val;
    }
    return std::sqrt(s / n);
}

// compute weighted max norm
template <typename T>
__attribute__((always_inline)) T wmnorm(const RP(T) x, const RP(T) w, const size_t n)
{
    T m = 0.0;
    for (size_t i = 0; i < n; i++) {
        m = std::max(m, std::abs(x[i] * w[i]));
    }
    return m;
}

// choose initial step-size
template <typename funcT, typename T>
T initStepSize(funcT f, const T t0, const RP(T) y0, const RP(T) f0, RP(T) y1, RP(T) f1, const T rtol, const T atol, const T direction, const T order, const size_t n)
{
    size_t i;
    T d0 = 0, d1 = 0, scale, h0, d2 = 0, h1;

    for (i = 0; i < n; i++) {
        scale = atol + rtol * std::abs(y0[i]);
        d0 += y0[i] * y0[i] / (scale * scale);
        d1 += f0[i] * f0[i] / (scale * scale);
    }
    d0 = std::sqrt(d0 / n);
    d1 = std::sqrt(d1 / n);


    if ((d0 < 1e-5) || (d1 < 1e-5)) {
        h0 = 1e-6;
    }
    else {
        h0 = 0.01 * d0 / d1;
    }

    for (i = 0; i < n; i++) {
        y1[i] = y0[i] + direction * h0 * f0[i];
    }
    f(t0 + h0, y1, f1);

    for (i = 0; i < n; i++) {
        scale = atol + rtol * std::abs(y0[i]);
        scale = (f1[i] - f0[i]) / scale;
        scale *= scale;
        d2 += scale;
    }
    d2 = std::sqrt(d2 / n) / h0;

    if (std::max(d1, d2) <= 1e-15) {
        h1 = std::max<T>(1e-6, h0 * 1e-3);
    }
    else {
        h1 = std::pow(0.01 / std::max(d1, d2), 1.0 / (order + 1));
    }

    return std::min<T>(100 * h0, h1);
}


/*
 * compute jacobian using forward difference, store in column major order
 * assumes a square jacobian
 */
template <typename funcT, typename T>
__attribute__((always_inline)) void njac1(funcT f, const T t, RP(T) x, const RP(T) fx, RP(T) J, const size_t n)
{
    size_t i, j;
    T temp, h;
    const T eps = std::sqrt(std::numeric_limits<T>::epsilon());

    for (i = 0; i < n; i++) {
        temp = x[i];
        h = (temp != 0) ? std::abs(temp) * eps : eps;
        x[i] = temp + h;
        f(t, x, &J[i * n]);
        x[i] = temp;
        for (j = 0; j < n; j++) {
            J[i * n + j] = (J[i * n + j] - fx[j]) / h;
        }
    }
}

#endif