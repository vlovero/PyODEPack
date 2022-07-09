#ifndef ODEPACK_RADAU5_H
#define ODEPACK_RADAU5_H

#include "odepack_common.h"
#include <type_traits>


#define NEWTON_MAXITER 10

template <typename T>
constexpr bool is_supported_data_type()
{
    constexpr bool b = std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, std::complex<float>> || std::is_same_v<T, std::complex<double>>;
    return b;
}

namespace radau5iia
{

    constexpr double s6 = 2.4494897427831780982;
    constexpr double c1 = 0.15505102572168219018;
    constexpr double c2 = 0.64494897427831780982;
    constexpr double c3 = 1;
    constexpr double mu_real = 3.6378342527444957322;
    constexpr std::complex<double> mu_complex = { 2.6810828736277521339, -3.0504301992474105694 };

    constexpr double t11 = 0.09443876248897524;
    constexpr double t12 = -0.14125529502095421;
    constexpr double t13 = 0.03002919410514742;
    constexpr double t21 = 0.25021312296533332;
    constexpr double t22 = 0.20412935229379994;
    constexpr double t23 = -0.38294211275726192;
    constexpr double t31 = 1;
    constexpr double t32 = 1;
    constexpr double t33 = 0;

    constexpr double ti11 = 4.17871859155190428;
    constexpr double ti12 = 0.32768282076106237;
    constexpr double ti13 = 0.52337644549944951;
    constexpr double ti21 = -4.17871859155190428;
    constexpr double ti22 = -0.32768282076106237;
    constexpr double ti23 = 0.47662355450055044;
    constexpr double ti31 = 0.50287263494578682;
    constexpr double ti32 = -2.57192694985560522;
    constexpr double ti33 = 0.59603920482822492;

    constexpr double p11 = 10.048809399827415266;
    constexpr double p12 = -25.629591447076639683;
    constexpr double p13 = 15.580782047249223972;
    constexpr double p21 = -1.3821427331607491919;
    constexpr double p22 = 10.296258113743305757;
    constexpr double p23 = -8.9141153805825570096;
    constexpr double p31 = 1.0 / 3.0;
    constexpr double p32 = -8.0 / 3.0;
    constexpr double p33 = 10.0 / 3.0;

    constexpr double e1 = -10.048809399827415562;
    constexpr double e2 = 1.3821427331607488958;
    constexpr double e3 = -0.33333333333333333333;


    template <typename T>
    constexpr std::array<T, 3> getNodes()
    {
        return { c1, c2, c3 };
    }

    template <typename T>
    constexpr std::array<T, 3> getTiR()
    {
        return { ti11, ti12, ti13 };
    }

    template <typename T>
    constexpr std::array<std::complex<T>, 3> getTiC()
    {
        return { std::complex<T>{ ti21, ti31 }, std::complex<T>{ ti22, ti32 }, std::complex<T>{ ti23, ti33 } };
    }

    template <typename T>
    constexpr std::array<T, 9> getT()
    {
        // row major order
        return { t11, t12, t13, t21, t22, t23, t31, t32, t33 };
        // column major order
        // return { t11, t21, t31, t12, t22, t32, t13, t23, t33 };
    }

    template <typename T>
    constexpr std::array<T, 9> getTi()
    {
        // row major
        return { ti11, ti12, ti13, ti21, ti22, ti23, ti31, ti32, ti33 };
        // column major order
        // return { ti11, ti21, ti31, ti12, ti22, ti32, ti13, ti23, ti33 };
    }

    template <typename T>
    constexpr std::array<T, 9> getP()
    {
        return { p11, p12, p13, p21, p22, p23, p31, p32, p33 };
    }


    template <typename T, const auto &lapack_getrs>
    void solveLU(const T *A, T *b, const int *ipiv, int n, int nrhs)
    {
        static_assert(is_supported_data_type<T>(), "[solveLU]: unsupported type");
        int info;
        char trans = 'N';
        //solve(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        lapack_getrs(&trans, &n, &nrhs, A, &n, ipiv, b, &n, &info);
    }

    template <typename T, const auto &lapack_getrf>
    void factorLU(T *A, int *ipiv, int n)
    {
        static_assert(is_supported_data_type<T>(), "[solveLU]: unsupported type");
        int info;
        // void (*dgetrf)(int *, int *, double *, int *, int *, int *);
        //solve(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        lapack_getrf(&n, &n, A, &n, ipiv, &info);
    }

    // interpolate points between t_n and t_n+1 using method specific interpolant P
    template <typename T, bool skipQ = false>
    bool interpolate(RP(T) y, const RP(T) t, const RP(T) y0, const RP(T) Z, const T t0, const T h, const T Dxy, const size_t n, const size_t m, RP(T) Q)
    {
        constexpr size_t p = 3;
        constexpr size_t s_star = 3;
        constexpr std::array<T, 9> P = getP<T>();

        size_t i, j, k;
        T z, s, uk;
        // T Q[n * p];

        if constexpr (!skipQ) {
            // Q = Z^T P
            for (i = 0; i < n; i++) {
                for (j = 0; j < p; j++) {
                    z = 0;
                    for (k = 0; k < s_star; k++) {
                        // z += K[k * n + i] * P[k * p + j];
                        z = std::fma(Z[k * n + i], P[k * p + j], z);
                    }
                    Q[i * p + j] = z;
                }
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

        // y = y0 + y
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                y[i * n + j] += y0[j];
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


    template <typename T, typename funcT, const auto &solveLUR, const auto &solveLUC>
    T solveCollocationSystem(funcT f, const T t, const RP(T) y, T h, size_t n, RP(T) Z, const RP(T) scale, const T tol, const RP(T) LUR, const RP(std::complex<T>) LUC, const int *ipiv, RP(T) work, size_t &niter, bool &converged)
    {
        size_t i, j, k, m;
        T s, normdW, val, rate = 0, prevNormdW = 1;

        const T MR = mu_real / h;
        const std::complex<T> MC = mu_complex / h;
        const int *ipivR = ipiv;
        const int *ipivC = ipiv + n;

        constexpr std::array<T, 3> C = getNodes<T>();

        constexpr std::array<T, 9> T_ = getT<T>();
        constexpr std::array<T, 9> Ti = getTi<T>();

        constexpr std::array<T, 3> TiR = getTiR<T>();
        constexpr std::array<std::complex<T>, 3> TiC = getTiC<T>();

        // since matricies of collocation matrices are known (and small)
        // doing naive matrix multiplcation should be plenty fast
        // the compiler will be able to unroll and vectorize it.
        auto naiveDot = [&](auto *B, const auto *A1, const auto *A2, size_t r1, size_t c1, size_t r2, size_t c2) {
            using RT = typename std::remove_reference<decltype(B[0])>::type;

            RT sum;
            for (i = 0; i < r1; i++) {
                for (j = 0; j < c2; j++) {
                    sum = 0;
                    for (k = 0; k < r2; k++) {
                        sum += A1[i * c1 + k] * A2[k * c2 + j];
                    }
                    B[i * c2 + j] = sum;
                }
            }
        };

        // get pointers from work array (need (16 * n * sizeof(T)) bytes in work)
        T *W = work;
        T *F = W + (3 * n);
        T *dW = F + (3 * n);
        T *temp = dW + (3 * n);
        T *fR = temp + n;
        std::complex<T> *fC = (std::complex<T> *)(fR + n);
        T *dWR;
        std::complex<T> *dWC;

        converged = false;
        // Z has shape (3, n) -> W has shape (3, n)
        naiveDot(W, Ti.data(), Z, 3, 3, 3, n);

        for (m = 0; m < NEWTON_MAXITER; m++) {

            // eval f at nodes
            for (i = 0; i < 3; i++) {
                for (j = 0; j < n; j++) {
                    temp[j] = y[j] + Z[i * n + j];
                }
                f(std::fma(h, C[i], t), temp, &F[i * n]);
            }

            // fR = (TiR^T * F)^T - MR * W[0]
            //        ^ col vec -> row vec
            naiveDot(fR, TiR.data(), F, 1, 3, 3, n);
            for (i = 0; i < n; i++) {
                fR[i] = std::fma(-MR, W[i], fR[i]);
            }

            // fC = (TiC^T * F)^T - MC * (W[1] + 1j * W[2])
            naiveDot(fC, TiC.data(), F, 1, 3, 3, n);
            for (i = 0; i < n; i++) {
                fC[i] -= MC * std::complex<T>(W[n + i], W[2 * n + i]);
            }

            // solve real system
            solveLUR(LUR, fR, ipivR, n, 1);
            dWR = fR; // compiler should optimize this away

            // solve complex system
            solveLUC(LUC, fC, ipivC, n, 1);
            dWC = fC; // compiler should optimize this away

            // optimize this?
            for (i = 0; i < n; i++) {
                dW[i] = dWR[i];
                dW[n + i] = ((T *)dWC)[2 * i];
                dW[2 * n + i] = ((T *)dWC)[2 * i + 1];
            }

            // compute norm of change
            normdW = 0.0;
            for (i = 0; i < 3; i++) {
                for (j = 0; j < n; j++) {
                    val = dW[i * n + j] / scale[j];
                    normdW += val * val;
                }
            }
            normdW = std::sqrt(normdW / (3 * n));

            // if there's previous data (and a rate) check for failed convergence
            if (m != 0) {
                rate = normdW / prevNormdW;
                val = std::pow(rate, (NEWTON_MAXITER - m)) / (1.0 - rate) * normdW;
                if ((1 <= rate) || (tol < val)) {
                    converged = false;
                    niter = m + 1;
                    return rate;
                }
            }

            // if not not converged apply change
            for (i = 0; i < (3 * n); i++) {
                W[i] += dW[i];
            }

            // update Z
            naiveDot(Z, T_.data(), W, 3, 3, 3, n);

            // if converged break
            if ((normdW == 0.0) || ((m != 0) && ((rate / (1.0 - rate) * normdW) < tol))) {
                converged = true;
                niter = m + 1;
                return rate;
            }

            // update prev data
            prevNormdW = normdW;
        }

        niter = m;
        return rate;
    }

    template <typename T, typename funcT, const auto &factorLUR, const auto &factorLUC, const auto &solveLUR, const auto &solveLUC>
    ODEResult integrate(funcT fun, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        const T Dxy = max_norm[0];
        const T *centers = &max_norm[1];

        const T eta = std::numeric_limits<T>::epsilon();
        const T sqeta = std::sqrt(eta);

        auto f = [=](const T t, const T *input, T *out) { fun(t, input, out, fargs); };

        const size_t nrwork = (26 + 4 * n) * n;
        const size_t niwork = 2 * n;
        void *allwork = std::malloc(nrwork * sizeof(T) + niwork * sizeof(int));
        memset(allwork, 0, nrwork * sizeof(T) + niwork * sizeof(int));

        T *work = (T *)allwork;
        T *yn = work; // n
        T *yn1 = yn + n; // n
        T *dy = yn1 + n; // n
        T *Z = dy + n; // n
        T *fn = Z + (3 * n); // n
        T *LUR = fn + n; // n^2
        std::complex<T> *LUC = (std::complex<T> *)(LUR + (n * n)); // 2 n^2
        T *ewt = (T *)(LUC + (n * n));
        T *J = ewt + n;
        T *Z0 = J + (n * n);
        T *otherWork = Z0 + (3 * n);
        T *Q = otherWork;
        T *ZError = otherWork + (3 * n);
        T *temp = ZError + n;
        int *ipiv = (int *)(work + nrwork);

        // step-size stuff and misc vars
        constexpr T err_exp = -0.25;
        T dt, dtabs, norm, factor, dtprev, normprev = std::numeric_limits<T>::quiet_NaN(), rate;
        T fac, r0, r, yi, safety = 1;
        const T tol = std::max<T>(10 * eta / rtol, std::min<T>(0.03, std::sqrt(rtol)));

        // keep track of function evaluations
        size_t fails = 0, feval = 0, jeval = 0, steps = 0, nlu = 0;

        // indexing / indexing for interpolation
        size_t i, j, k, index_t0 = 1, index_t1, niter;

        // values for interpolation
        bool performInterpolate;
        bool need2break;
        bool converged;
        bool currentJac;
        bool rejected = false;
        bool recomputeJ;
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
        f(tn, yn, fn);
        dtabs = (h0 == (T)0) ? initStepSize<decltype(f), T>(f, tn, y0, fn, yn1, Z + n, rtol, atol, backwards, 3, n) : std::abs(h0);
        dt = backwards * dtabs;
        feval = 2;

        // init ewt
        for (i = 0; i < n; i++) {
            ewt[i] = std::fma(rtol, std::abs(yn[i]), atol);
        }

        // routine for creating jacobian
        auto genJ = [&]() __attribute__((always_inline))
        {
            fac = wmnorm<T>(ewt, fn, n);
            r0 = 1e3 * dtabs * eta * n * fac;
            if (r0 == 0) {
                r0 = 1;
            }
            for (i = 0; i < n; i++) {
                yi = yn[i];
                r = std::max(sqeta * std::abs(yi), r0 / ewt[i]);
                yn[i] += r;
                // pre-add negative sign to jacobian
                fac = -1.0 / r;
                f(tn, yn, Q);
                for (j = 0; j < n; j++) {
                    J[i * n + j] = (Q[j] - fn[j]) * fac;
                }
                yn[i] = yi;
            }
            jeval++;
            feval += n;
            currentJac = true;
        };

        // routine for real and complex lu factorizations
        auto updateLU = [&]() __attribute__((always_inline))
        {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    LUR[i * n + j] = (i != j) ? J[i * n + j] : (J[i * n + j] + mu_real / dt);
                    LUC[i * n + j] = (i != j) ? J[i * n + j] : (J[i * n + j] + mu_complex / dt);
                }
            }
            factorLUR(LUR, ipiv, n);
            factorLUC(LUC, ipiv + n, n);
            nlu += 2;
        };

        auto updateFactor = [&]() __attribute__((always_inline))
        {
            // no prev info on first step
            if ((normprev == std::numeric_limits<T>::quiet_NaN()) || (norm == 0)) {
                factor = 1;
            }
            else {
                factor = (dtabs / dtprev) * std::pow(normprev / norm, 0.25);
            }
            factor = (norm) ? (safety * std::min<T>(1, factor) * std::pow(norm, err_exp)) : 10;
            factor = std::max<T>(0.2, factor);
            factor = std::min<T>(10, factor);
        };

        // jacobian and factor real and complex
        genJ();
        updateLU();


        std::memset(Z0, 0, 3 * n * sizeof(T));

        while (backwards * (tf - tn) > 0) {
            steps++;
            if (backwards * (tn + dt - tf) > 0) {
                dt = tf - tn;
            }

            // init Z for solver
            std::memcpy(Z, Z0, 3 * n * sizeof(T));

            // solve for Z
            rate = solveCollocationSystem<T, decltype(f), solveLUR, solveLUC>(f, tn, yn, dt, n, Z, ewt, tol, LUR, LUC, ipiv, otherWork, niter, converged);
            feval += 3 * niter;

            if ((!converged) && (!currentJac)) {
                dtabs *= 0.5;
                dt *= 0.5;
                genJ();
                updateLU();
                fails++;
                continue;
            }

            // get yn1 from Z and vectors for error calculation
            for (i = 0; i < n; i++) {
                yn1[i] = yn[i] + Z[2 * n + i];
                ZError[i] = (e1 * Z[i] + e2 * Z[n + i] + e3 * Z[2 * n + i]) / dt;
                temp[i] = fn[i] + ZError[i];
            }

            // compute error norm and scaling factor for step-size
            solveLUR(LUR, temp, ipiv, n, 1); // error vector stored in temp
            norm = gerr_norm<T>(temp, yn, yn1, atol, rtol, n);
            safety = 0.9 * (T)(2 * NEWTON_MAXITER + 1) / (T)(2 * NEWTON_MAXITER + niter);

            if (rejected && ((T)1 < norm) && (min_step < dtabs)) {
                for (i = 0; i < n; i++) {
                    temp[i] += yn[i];
                }
                f(tn, temp, Q);
                for (i = 0; i < n; i++) {
                    temp[i] = ZError[i] + Q[i];
                }
                solveLUR(LUR, temp, ipiv, n, 1);
                norm = gerr_norm<T>(temp, yn, yn1, atol, rtol, n);
            }

            if (((T)1 < norm) && (min_step < dtabs)) {
                // if norm larger than one and step not too small then
                // redo step with smaller step size.
                fails++;
                rejected = true;
                updateFactor();
                dtabs *= factor;
                dtabs = std::max(dtabs, min_step);
                dt = dtabs * backwards;
                updateLU();
                continue;
            }
            else {
                rejected = false;
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

                    need2break = interpolate<T>(&Y[n * index_t0], &teval[index_t0], yn, Z, t0, dt, Dxy, n, index_t1 - index_t0, Q);
                    if (need2break) {
                        PRINT_ERROR_MESSAGE("interpolant not in computation window or not finite");
                        std::free(allwork);
                        return { index_t0, feval, jeval, fails, steps };
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
                        PRINT_ERROR_MESSAGE("solution not in computation window or not finite");
                        std::free(allwork);
                        return { index_t0, feval, jeval, fails, steps };
                    }
                }

                // recompute_jac = jac is not None and n_iter > 2 and rate > 1e-3
                recomputeJ = (2 < niter) && (1e-3 < rate);
                updateFactor();
                factor = std::min<T>(10, factor);

                if ((!recomputeJ) && (factor < 1.2)) {
                    factor = 1;
                }

                // update Z0;
                temp[0] = tn + c1 * dt;
                temp[1] = tn + c2 * dt;
                temp[2] = tn + c3 * dt;
                interpolate<T>(Z0, temp, yn, Z, t0, dt, Dxy, n, 3, Q);
                for (i = 0; i < 3; i++) {
                    for (j = 0; j < n; j++) {
                        Z0[i * n + j] -= yn1[j];
                    }
                }

                // update yn and fn
                std::memcpy(yn, yn1, n * sizeof(T));
                f(tn, yn, fn);

                // update ewt
                for (i = 0; i < n; i++) {
                    ewt[i] = std::fma(rtol, std::abs(yn[i]), atol);
                }

                // update values
                normprev = norm;
                dtprev = dtabs;
                dtabs *= factor;
                dtabs = std::min(max_step, std::max(dtabs, min_step));
                dt = backwards * dtabs;

                // update if necessary
                if (recomputeJ) {
                    genJ();
                    updateLU();
                }
                else {
                    currentJac = false;
                    if (factor != 1) {
                        updateLU();
                    }
                }
            }
        }
        std::free(allwork);
        return { m, feval, jeval, fails, steps };
    }
}

#endif // ODEPACK_RADAU5_H