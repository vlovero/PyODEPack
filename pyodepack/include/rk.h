/*
 * Generalized explicit Runge-Kutta solver that uses interpolation
 * to compute points between two successive steps. Algorithm inspired
 * by [1]
 *
 * The constexpr template arguments allows the compiler to unroll
 * and vectorize many of the loops which means performance will
 * be the same as if each method was hand coded.
 *
 * Vincent Lovero
 *
 * [1] E. Hairer, S. P. Norsett G. Wanner, “Solving Ordinary Differential Equations I: Nonstiff Problems”
 */

#ifndef ODEPACK_GENERAL_RK_H
#define ODEPACK_GENERAL_RK_H

#include "odepack_common.h"

namespace erk
{
    template <typename T, typename funcT, size_t s, size_t s_star, size_t p,
        const std::array<T, s_star> &C,
        size_t ldA, const std::array<T, s_star *ldA> &A,
        const std::array<T, s_star *p> &P,
        const std::array<T, s> &B,
        const std::array<T, s> &D,
        bool scalar = false>
    class RKStepper
    {
        T m_t;
        T m_h;
        T m_hprev;
        size_t m_fails;
        size_t m_feval;
        size_t m_steps;
        size_t m_jeval;
        bool m_interpolate;
        ODEExitCodes m_status;

    public:
        RKStepper()
        {
            // std::cout << "calling cunt RKStepper\n";
            m_t = std::numeric_limits<T>::quiet_NaN();
            m_h = 0;
            m_hprev = 0;
            m_fails = 0;
            m_feval = 0;
            m_steps = 0;
            m_jeval = 0;
            m_interpolate = 0;
            m_status = ODEExitCodes::success;
        }

        ODEResult<T> toODEResult() const
        {
            return { 1, m_feval, m_jeval, m_fails, m_steps, 0, m_status, m_h };
        }

        auto get_t() const { return m_t; }
        auto get_h() const { return m_h; }
        auto get_hprev() const { return m_hprev; }
        auto get_fails() const { return m_fails; }
        auto get_feval() const { return m_feval; }
        auto get_steps() const { return m_steps; }
        auto get_jeval() const { return m_jeval; }
        auto get_interpolate() const { return m_interpolate; }
        auto get_status() const { return m_status; }

        void step(funcT f, T tn, const T tf, const RP(T) y0, RP(T) y, const T atol, const T rtol, const T minStep, const T maxStep, const T initStep, const size_t node, const bool firstStep, const RP(T) window, RP(T) work, const RP(void) args)
        {
            // std::cout << "tn = " << tn << '\n';
            // std::cout << "tf = " << tf << '\n';
            size_t i, j, k;
            bool needToFree;
            const T backwards = std::sign<T>(tf - tn);
            const T maxL2 = window[0];
            const T *centers = window + 1;

            T h, hlo, val, norm, factor, t0, l2;
            T *yn, *yn1, *dy, *temp, *K, *Q;

            needToFree = work == nullptr;
            if (needToFree) {
                size_t nwork = node * (4 + s_star + p);
                work = (T *)std::malloc(nwork * sizeof(T));
            }

            yn = work;
            yn1 = work + node;
            dy = yn1 + node;
            temp = dy + node;
            K = temp + node;
            Q = K + (node * s_star);

            if (!std::isfinite(m_t)) {
                // std::cout << "m_t is nan -> set m_t = tn\n";
                m_t = tn;
            }

            if (((backwards * (tf - m_t)) <= 0) && m_interpolate && (!firstStep)) {
                this->interpolate(false, y, temp, K, Q, m_t - m_hprev, tf, m_hprev, node);
                return;
            }

            std::memcpy(yn, y0, node * sizeof(T));
            tn = m_t;

            if (firstStep || needToFree) {
                f(tn, yn, K, args);
            }
            else {
                std::memcpy(K, &K[(s - 1) * node], node * sizeof(T));
            }

            if (m_h == 0) {
                if (initStep) {
                    m_h = initStep;
                }
                else {
                    // std::cout << "choosing step size\n";
                    m_h = initStepSize<decltype(f), T>(f, tn, yn, K, temp, K + node, rtol, atol, backwards, p - 1, node, args);
                }
            }
            h = backwards * m_h;
            // std::cout << "tn = " << tn << '\n';
            // std::cout << "tf = " << tf << '\n';
            // std::cout << "h = " << h << '\n';
            // std::cout << "(backwards * (tf - tn)) = " << (backwards * (tf - tn)) << '\n';

            while (0 < (backwards * (tf - tn))) {
                m_steps++;

                hlo = std::abs(std::nextafter(tn, backwards * std::numeric_limits<T>::infinity()) - tn);
                m_h = std::max(m_h, hlo);
                h = backwards * m_h;
                // std::cout << "h = " << h << '\n';

                // compute K_2 - K_{s - 1}
                for (i = 1; i < s - 1; i++) {
                    for (j = 0; j < node; j++) {
                        val = 0;
                        for (k = 0; k < i; k++) {
                            // val += A[i * s + k] * K[k * n + j];
                            val = std::fma(A[i * ldA + k], K[k * node + j], val);
                        }
                        // temp[j] = yn[j] + dt * val;
                        temp[j] = std::fma(h, val, yn[j]);
                    }
                    f(std::fma(C[i], h, tn), temp, &K[i * node], args);
                }

                // compute y_{n + 1} and K_s
                for (j = 0; j < node; j++) {
                    val = 0;
                    for (k = 0; k < i; k++) {
                        // val += B[k] * K[k * n + j];
                        val = std::fma(B[k], K[k * node + j], val);
                    }
                    // yn1[j] = yn[j] + dt * val;
                    yn1[j] = std::fma(h, val, yn[j]);
                }
                f(std::fma(C[i], h, tn), yn1, &K[i * node], args);
                i++;

                // compute dy
                for (j = 0; j < node; j++) {
                    val = 0;
                    for (k = 0; k < i; k++) {
                        // val += D[k] * K[k * n + j];
                        val = std::fma(D[k], K[k * node + j], val);
                    }
                    dy[j] = h * val;
                }

                // compute error norm and scaling factor for step-size
                norm = errorNorm<T>(dy, yn, yn1, atol, rtol, node);
                factor = norm ? 0.9 * std::pow(norm, -1.0 / (p - 1)) : 10.0;
                factor = std::min<T>(10, factor);
                m_feval += (s - 1);

                if ((1.0 < norm) and (minStep < m_h)) {
                    m_h *= factor;
                    m_h = std::max(std::max(m_h, minStep), hlo);
                    h = backwards * m_h;
                    m_fails += 1;
                    continue;
                }
                else {
                    t0 = tn;
                    tn += h;
                    m_t = tn;
                    m_interpolate = 0 <= (backwards * (tn - tf));

                    if (m_interpolate) {
                        // std::cout << "i = " << i << '\n';
                        for (; i < s_star; i++) {
                            for (j = 0; j < node; j++) {
                                val = 0;
                                for (k = 0; k < i; k++) {
                                    // val += A[i * s + k] * K[k * n + j];
                                    val = std::fma(A[i * ldA + k], K[k * node + j], val);
                                }
                                // temp[j] = yn[j] + dt * val;
                                temp[j] = std::fma(h, val, yn[j]);
                            }
                            f(std::fma(C[i], h, tn), temp, &K[i * node], args);
                        }
                        m_feval += (s_star - s);
                        this->interpolate(true, y, yn, K, Q, t0, tf, h, node);

                        l2 = 0;
                        for (j = 0; j < node; j++) {
                            val = y[j] - centers[j];
                            l2 = std::fma(val, val, l2);
                        }
                        l2 = std::sqrt(l2);
                        if ((maxL2 < l2) || (!std::isfinite(l2))) {
                            m_status = ODEExitCodes::failedInterpolation;
                            if (needToFree) {
                                std::free(work);
                            }
                            return;
                        }
                    }

                    std::memcpy(temp, yn, node * sizeof(T));
                    std::memcpy(yn, yn1, node * sizeof(T));

                    m_hprev = h;
                    m_h *= factor;
                    m_h = std::min(maxStep, std::max(m_h, minStep));
                    h = backwards * m_h;

                    if (!m_interpolate) {
                        std::memcpy(K, &K[(s - 1) * node], node * sizeof(T));

                        l2 = 0;
                        for (j = 0; j < node; j++) {
                            val = yn[j] - centers[j];
                            l2 = std::fma(val, val, l2);
                        }
                        l2 = std::sqrt(l2);
                        if ((maxL2 < l2) || (!std::isfinite(l2))) {
                            m_status = ODEExitCodes::outOfBounds;
                            if (needToFree) {
                                std::free(work);
                            }
                            return;
                        }
                    }
                }
            }
            if (needToFree) {
                // std::cout << "freeing memory\n";
                std::free(work);
            }
        }

    private:
        void interpolate(const bool computeQ, RP(T) y, const RP(T) y0, const RP(T) K, RP(T) Q, const T t0, const T t, const T h, const size_t n)
        {
            size_t i, j, k;
            T z, acc, uk; // changed s to acc because size_t is already decleared

            if (computeQ) {
                // Q = K^T P
                for (i = 0; i < n; i++) {
                    for (j = 0; j < p; j++) {
                        z = 0;
                        for (k = 0; k < s_star; k++) {
                            z = std::fma(K[k * n + i], P[k * p + j], z);
                        }
                        Q[i * p + j] = z;
                    }
                }
            }

            // z = (uk, uk^2, ..., uk^p)^T
            // y = h * Q z + y0
            z = (t - t0) / h;
            for (i = 0; i < n; i++) {
                uk = z;
                acc = 0;
                for (k = 0; k < p; k++) {
                    // acc += Q[i * p + k] * uk;
                    acc = std::fma(Q[i * p + k], uk, acc);
                    uk *= z;
                }
                y[i] = std::fma(h, acc, y0[i]);
            }
        }
    };

    // interpolate points between t_n and t_n+1 using method specific interpolant P
    template <typename T, size_t s_star, size_t p, const std::array<T, s_star *p> &P>
    bool interpolate(RP(T) y, const RP(T) t, const RP(T) y0, const RP(T) K, const T t0, const T h, const T Dxy, const size_t n, const size_t m, RP(T) Q)
    {
        size_t i, j, k;
        T z, s, uk;
        // T Q[n * p];

        // Q = K^T P
        for (i = 0; i < n; i++) {
            for (j = 0; j < p; j++) {
                z = 0;
                for (k = 0; k < s_star; k++) {
                    // z += K[k * n + i] * P[k * p + j];
                    z = std::fma(K[k * n + i], P[k * p + j], z);
                }
                Q[i * p + j] = z;
            }
        }

        // z = (uk, uk^2, ..., uk^p)^T
        // y = Q z
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                z = (t[j] - t0) / h;
                uk = z;
                s = 0;
                for (k = 0; k < p; k++) {
                    // s += Q[i * p + k] * uk;
                    s = std::fma(Q[i * p + k], uk, s);
                    uk *= z;
                }
                y[j * n + i] = s;
            }
        }

        // y = y0 + h * y
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                // y[i * n + j] = y0[j] + h * y[i * n + j];
                y[i * n + j] = std::fma(h, y[i * n + j], y0[j]);
            }
        }
        // check for inf/nan
        for (i = 0; i < m; i++) {
            z = 0;
            for (j = 0; j < n; j++) {
                // z += y[i * n + j] * y[i * n + j];
                z = std::fma(y[i * n + j], y[i * n + j], z);
            }
            if ((Dxy < z) || (!std::isfinite(z))) {
                return true;
            }
        }
        return false;
    }

    template <typename T, typename funcT, size_t s, size_t s_star, size_t p,
        const std::array<T, s_star> &C,
        size_t ldA, const std::array<T, s_star *ldA> &A,
        const std::array<T, s_star *p> &P,
        const std::array<T, s> &B,
        const std::array<T, s> &D,
        bool scalar = false>
    ODEResult<T> integrate(funcT fn, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        const T Dxy = max_norm[0];
        const T *centers = &max_norm[1];

        auto f = [=](const T t, const T *input, T *out) { fn(t, input, out, fargs); };

        T *work = (T *)std::malloc((4 + s_star + p) * n * sizeof(T));
        T *yn = &work[0];
        T *yn1 = yn + n;
        T *dy = yn1 + n;
        T *temp = dy + n;
        T *K = temp + n;
        T *Q = K + (s_star * n);

        // step-size stuff and misc vars
        constexpr T err_exp = -1.0 / (T)p;
        T dt, dtabs, norm, factor;

        // keep track of function evaluations
        size_t fails = 0, feval = 0, jeval = 0, steps = 0;

        // indexing / indexing for interpolation
        size_t i, j, k, index_t0 = 1, index_t1;

        // values for interpolation
        bool performInterpolate;
        bool need2break;
        T t0;
        T t1 = teval[1];

        // values for keep track of current state
        T tn = teval[0];
        T tf = teval[m - 1];
        T val, dst;
        const T backwards = std::sign(t1 - tn);

        // start initializing things
        std::memcpy(yn, y0, n * sizeof(T));
        std::memcpy(Y, y0, n * sizeof(T));
        f(tn, yn, K);
        dtabs = (h0 == (T)0) ? initStepSize<decltype(f), T>(f, tn, y0, K, yn1, K + n, rtol, atol, backwards, p - 1, n) : std::abs(h0);
        dt = backwards * dtabs;
        feval = 2;

        while (backwards * (tf - tn) > 0) {
            steps++;
            if (backwards * (tn + dt - tf) > 0) {
                dt = tf - tn;
            }

            // compute K_2 - K_{s - 1}
            for (i = 1; i < s - 1; i++) {
                for (j = 0; j < n; j++) {
                    val = 0;
                    for (k = 0; k < i; k++) {
                        // val += A[i * s + k] * K[k * n + j];
                        val = std::fma(A[i * ldA + k], K[k * n + j], val);
                    }
                    // temp[j] = yn[j] + dt * val;
                    temp[j] = std::fma(dt, val, yn[j]);
                }
                f(std::fma(C[i], dt, tn), temp, &K[i * n]);
            }
            // compute y_{n + 1} and K_s
            for (j = 0; j < n; j++) {
                val = 0;
                for (k = 0; k < i; k++) {
                    // val += B[k] * K[k * n + j];
                    val = std::fma(B[k], K[k * n + j], val);
                }
                // yn1[j] = yn[j] + dt * val;
                yn1[j] = std::fma(dt, val, yn[j]);
            }
            f(std::fma(C[i], dt, tn), yn1, &K[i * n]);
            i++;

            // compute dy
            for (j = 0; j < n; j++) {
                val = 0;
                for (k = 0; k < i; k++) {
                    // val += D[k] * K[k * n + j];
                    val = std::fma(D[k], K[k * n + j], val);
                }
                dy[j] = dt * val;
            }

            // compute error norm and scaling factor for step-size
            norm = errorNorm<T>(dy, yn, yn1, atol, rtol, n);
            factor = norm ? 0.9 * std::pow(norm, err_exp) : 10.0;
            factor = std::min<T>(10, factor);
            feval += (s - 1);

            if ((1 < norm) && (min_step < dtabs)) {
                // if norm larger than one and step not too small then
                // redo step with smaller step size.
                dtabs *= factor;
                dtabs = std::max(dtabs, min_step);
                dt = dtabs * backwards;
                fails++;
                continue;
            }
            else {
                t0 = tn;
                tn += dt;

                // if we've reached / passed checkpoint, interpolate between current and last state
                performInterpolate = (0 <= (backwards * (tn - t1)));
                if (performInterpolate) {
                    // find first index out out of interpolant bounds
                    for (index_t1 = index_t0 + 1; index_t1 < m; index_t1++) {
                        if (0 < (backwards * (teval[index_t1] - tn))) {
                            break;
                        }
                    }

                    for (; i < s_star; i++) {
                        for (j = 0; j < n; j++) {
                            val = 0;
                            for (k = 0; k < i; k++) {
                                // val += A[i * s + k] * K[k * n + j];
                                val = std::fma(A[i * ldA + k], K[k * n + j], val);
                            }
                            // temp[j] = yn[j] + dt * val;
                            temp[j] = std::fma(dt, val, yn[j]);
                        }
                        f(std::fma(C[i], dt, tn), temp, &K[i * n]);
                    }
                    feval += (s_star - s);

                    need2break = interpolate<T, s_star, p, P>(&Y[n * index_t0], &teval[index_t0], yn, K, t0, dt, Dxy, n, index_t1 - index_t0, Q);
                    if (need2break) {
                        // PRINT_ERROR_MESSAGE("interpolant not in computation window or not finite");
                        std::free(work);
                        return { index_t0, feval, jeval, fails, steps, 0, ODEExitCodes::failedInterpolation, dt };
                    }
                    index_t0 = index_t1;
                    t1 = teval[index_t1];
                }
                else {
                    // compute norm if norm wasn't computed during interpolation
                    dst = 0;
                    for (i = 0; i < n; i++) {
                        val = yn1[i] - centers[i];
                        // dst += val * val;
                        dst = std::fma(val, val, dst);
                    }
                    if ((dst > Dxy) || !std::isfinite(dst)) {
                        // PRINT_ERROR_MESSAGE("solution not in computation window or not finite");
                        std::free(work);
                        return { index_t0, feval, jeval, fails, steps, 0, ODEExitCodes::outOfBounds, dt };
                    }
                }

                // update yn
                std::memcpy(yn, yn1, n * sizeof(T));
                std::memcpy(K, &K[(s - 1) * n], n * sizeof(T));

                dtabs *= factor;
                dtabs = std::min(max_step, std::max(dtabs, min_step));
                dt = backwards * dtabs;
            }
        }
        std::free(work);
        return { m, feval, jeval, fails, steps, 0, ODEExitCodes::success, dt };
    }

    template <typename T, typename funcT, size_t s, size_t s_star, size_t p,
        const std::array<T, s_star> &C,
        size_t ldA, const std::array<T, s_star *ldA> &A,
        const std::array<T, s_star *p> &P,
        const std::array<T, s> &B,
        const std::array<T, s> &D,
        bool scalar = false>
    ODEResult<T> integrateScalar(funcT fn, const T y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        const T Dxy = max_norm[0];
        const T center = max_norm[1];

        auto f = [=](const T t, const T *input, T *out) { fn(t, input, out, fargs); };

        T yn;
        T yn1;
        T dy;
        T temp;
        T K[s_star];
        T Q[p];

        // step-size stuff and misc vars
        constexpr T err_exp = -1.0 / (T)p;
        T dt, dtabs, norm, factor;

        // keep track of function evaluations
        size_t fails = 0, feval = 0, jeval = 0, steps = 0;

        // indexing / indexing for interpolation
        size_t i, j, k, index_t0 = 1, index_t1;

        // values for interpolation
        bool performInterpolate;
        bool need2break;
        T t0;
        T t1 = teval[1];

        // values for keep track of current state
        T tn = teval[0];
        T tf = teval[m - 1];
        T val, dst;
        const T backwards = std::sign(t1 - tn);

        // start initializing things
        yn = y0;
        Y[0] = y0;
        K[0] = f(tn, yn);
        dtabs = (h0 == (T)0) ? initStepSizeScalar<decltype(f), T>(f, tn, y0, K[0], rtol, atol, backwards, p - 1) : std::abs(h0);
        dt = backwards * dtabs;
        feval = 2;

        while (backwards * (tf - tn) > 0) {
            steps++;
            if (backwards * (tn + dt - tf) > 0) {
                dt = tf - tn;
            }

            // compute K_2 - K_{s - 1}
            for (i = 1; i < s - 1; i++) {
                val = 0;
                for (k = 0; k < i; k++) {
                    // val += A[i * s + k] * K[k * n + j];
                    val = std::fma(A[i * ldA + k], K[k], val);
                }
                // temp[j] = yn[j] + dt * val;
                temp = std::fma(dt, val, yn[j]);
                K[i] = f(std::fma(C[i], dt, tn), temp);
            }
            // compute y_{n + 1} and K_s
            val = 0;
            for (k = 0; k < i; k++) {
                // val += B[k] * K[k * n + j];
                val = std::fma(B[k], K[k], val);
            }
            // yn1[j] = yn[j] + dt * val;
            yn1 = std::fma(dt, val, yn);
            K[i] = f(std::fma(C[i], dt, tn), yn1);
            i++;

            // compute dy
            val = 0;
            for (k = 0; k < i; k++) {
                // val += D[k] * K[k * n + j];
                val = std::fma(D[k], K[k], val);
            }
            dy = dt * val;

            // compute error norm and scaling factor for step-size
            norm = std::abs(dy / (atol + rtol * std::max(std::abs(yn), std::abs(yn1))));
            factor = norm ? 0.9 * std::pow(norm, err_exp) : 10.0;
            factor = std::min<T>(10, factor);
            feval += (s - 1);

            if ((1 < norm) && (min_step < dtabs)) {
                // if norm larger than one and step not too small then
                // redo step with smaller step size.
                dtabs *= factor;
                dtabs = std::max(dtabs, min_step);
                dt = dtabs * backwards;
                fails++;
                continue;
            }
            else {
                t0 = tn;
                tn += dt;

                // if we've reached / passed checkpoint, interpolate between current and last state
                performInterpolate = (0 <= (backwards * (tn - t1)));
                if (performInterpolate) {
                    // find first index out out of interpolant bounds
                    for (index_t1 = index_t0 + 1; index_t1 < m; index_t1++) {
                        if (0 < (backwards * (teval[index_t1] - tn))) {
                            break;
                        }
                    }

                    for (; i < s_star; i++) {
                        val = 0;
                        for (k = 0; k < i; k++) {
                            // val += A[i * s + k] * K[k * n + j];
                            val = std::fma(A[i * ldA + k], K[k], val);
                        }
                        // temp[j] = yn[j] + dt * val;
                        temp = std::fma(dt, val, yn);
                        K[i] = f(std::fma(C[i], dt, tn), temp);
                    }
                    feval += (s_star - s);

                    need2break = interpolate<T, s_star, p, P>(&Y[index_t0], &teval[index_t0], yn, K, t0, dt, Dxy, 1, index_t1 - index_t0, Q);
                    if (need2break) {
                        // PRINT_ERROR_MESSAGE("interpolant not in computation window or not finite");
                        return { index_t0, feval, jeval, fails, steps, 0, ODEExitCodes::failedInterpolation, dt };
                    }
                    index_t0 = index_t1;
                    t1 = teval[index_t1];
                }
                else {
                    // compute norm if norm wasn't computed during interpolation
                    dst = std::abs(yn1 - center);
                    if ((dst > Dxy) || !std::isfinite(dst)) {
                        // PRINT_ERROR_MESSAGE("solution not in computation window or not finite");
                        return { index_t0, feval, jeval, fails, steps, 0, ODEExitCodes::outOfBounds, dt };
                    }
                }

                // update yn
                yn = yn1;
                K[0] = K[s - 1];

                dtabs *= factor;
                dtabs = std::min(max_step, std::max(dtabs, min_step));
                dt = backwards * dtabs;
            }
        }
        return { m, feval, jeval, fails, steps, 0, ODEExitCodes::success, dt };
    }
}

#endif // ODEPACK_GENERAL_RK_H