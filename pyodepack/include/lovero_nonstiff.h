#ifndef ODEPACK_LOVERO_NS2_H
#define ODEPACK_LOVERO_NS2_H

#include "odepack_common.h"

#if 1
#define LOVERO_INLINE inline
#else
#define LOVERO_INLINE [[gnu::always_inline]]
#endif

#if 0
#define QUAD_INLINE
#else
#define QUAD_INLINE [[gnu::always_inline]]
#endif

#define TEST_FASTER_QUAD 1
#define TEST_SWITCHES_AT_TOP 1

namespace loveroNS
{
    // number of points is 13 -> only need up to 7 point gauss
    constexpr std::array<double, 0> gaussXK1 = {  };
    constexpr std::array<double, 1> gaussWK1 = { 2.000000000000000000000000000000 };

    constexpr std::array<double, 1> gaussXK2 = { 0.577350269189625731058868041146 };
    constexpr std::array<double, 1> gaussWK2 = { 1.000000000000000000000000000000 };

    constexpr std::array<double, 1> gaussXK3 = { 0.774596669241483404277914814884 };
    constexpr std::array<double, 2> gaussWK3 = { 0.888888888888888839545643349993, 0.555555555555555691249480787519 };

    constexpr std::array<double, 2> gaussXK4 = { 0.339981043584856257311344052141, 0.861136311594052572537805190223 };
    constexpr std::array<double, 2> gaussWK4 = { 0.652145154862546205798423670785, 0.347854845137453683179273866699 };

    constexpr std::array<double, 2> gaussXK5 = { 0.538469310105683107714469315397, 0.906179845938663963700321346550 };
    constexpr std::array<double, 3> gaussWK5 = { 0.568888888888888999417758896016, 0.478628670499366193347157150129, 0.236926885056189417966265864379 };

    constexpr std::array<double, 3> gaussXK6 = { 0.238619186083196932468553086437, 0.661209386466264481541088571248, 0.932469514203152050058065469784 };
    constexpr std::array<double, 3> gaussWK6 = { 0.467913934572691370128438848042, 0.360761573048138939334705810325, 0.171324492379169746048006572892 };

    constexpr std::array<double, 3> gaussXK7 = { 0.405845151377397184155881859624, 0.741531185599394460083999547351, 0.949107912342758486268223805382 };
    constexpr std::array<double, 4> gaussWK7 = { 0.417959183673468959163699310011, 0.381830050505118312464958307828, 0.279705391489276589123136318449, 0.129484966168870646585631334347 };

    /*
        Gaussian Quadrature for vector functions
        signature (scalar, pointer) -> void
    */
    template <const auto &nodes, const auto &weights, typename funcT, typename T>
    QUAD_INLINE void quadG(funcT &f, const T a, const T b, RP(T) result, const size_t n)
    {
        constexpr size_t nnodes = nodes.size();
        constexpr size_t nweights = weights.size();

        static_assert((1 <= nweights) && (nweights <= 4));
        constexpr bool even = nnodes == nweights;
        constexpr size_t shift = !even;

        const T center = (a + b) * 0.5;
        const T h2 = (b - a) * 0.5;
        size_t i, k;
        T fv1[n];
        T fv2[n];
        T fsum;
        T val;

        if constexpr (even) {
            for (i = 0; i < n; i++) {
                result[i] = 0;
            }
        }
        else {
            f(center, result);

            for (i = 0; i < n; i++) {
                result[i] *= weights[0];
            }
        }
        for (i = 0; i < nnodes; i++) {
            val = h2 * nodes[i];
            f(center - val, fv1);
            f(center + val, fv2);
            for (k = 0; k < n; k++) {
                fsum = fv1[k] + fv2[k];
                result[k] += weights[i + shift] * fsum;
            }
        }

        for (i = 0; i < n; i++) {
            result[i] *= h2;
        }
    }

    /*
        Gaussian Quadrature for scalar functions
        signature (scalar) -> scalar
    */
    template <const auto &nodes, const auto &weights, typename funcT, typename T>
    QUAD_INLINE T quadG(funcT &f, const T a, const T b)
    {
        constexpr size_t nnodes = nodes.size();
        constexpr size_t nweights = weights.size();

        static_assert((1 <= nweights) && (nweights <= 4));
        constexpr bool even = nnodes == nweights;
        constexpr size_t shift = !even;

        const T center = (a + b) * 0.5;
        const T h2 = (b - a) * 0.5;
        size_t i;
        T fv1;
        T fv2;
        T fsum;
        T val;
        T result;

        if constexpr (even) {
            result = 0;
        }
        else {
            result = weights[0] * f(center);
        }
        for (i = 0; i < nnodes; i++) {
            val = h2 * nodes[i];
            fv1 = f(center - val);
            fv2 = f(center + val);
            fsum = fv1 + fv2;
            result += weights[i + shift] * fsum;
        }

        result *= h2;
        return result;
    }

    /*
        Linearly weighted gaussian quadrature for scalar functions
        returns \int_a^b f(t)dt and wres will contain \int_a^b f(t) (t - x)dt
        This is used to make prediction step go from O(s^2 n) work to
        O(sn + s^2) work. So for larger systems, it will be much more efficient.
    */
    template <const auto &nodes, const auto &weights, typename funcT, typename T>
    QUAD_INLINE T quadGLW(funcT &f, const T a, const T b, T &wres, const T x)
    {
        constexpr size_t nnodes = nodes.size();
        constexpr size_t nweights = weights.size();

        static_assert((1 <= nweights) && (nweights <= 4));
        constexpr bool even = nnodes == nweights;
        constexpr size_t shift = !even;

        const T center = (a + b) * 0.5;
        const T h2 = (b - a) * 0.5;
        size_t i;
        T fv1;
        T fv2;
        T fsum;
        T val;
        T result;
        T weightedResult;

        if constexpr (even) {
            result = 0;
            weightedResult = 0;
        }
        else {
            result = weights[0] * f(center);
            weightedResult = result * (center - x);
        }
        for (i = 0; i < nnodes; i++) {
            val = h2 * nodes[i];
            fv1 = f(center - val);
            fv2 = f(center + val);
            fsum = fv1 + fv2;
            result += weights[i + shift] * fsum;
            // compute linearly weighted values
            fv1 *= (center - val - x);
            fv2 *= (center + val - x);
            fsum = fv1 + fv2;
            weightedResult += weights[i + shift] * fsum;
        }

        result *= h2;
        weightedResult *= h2;
        wres = weightedResult;
        return result;
    }


    template <typename T, typename funcT>
    struct LoveroNS
    {
        // typedef double T;
        // typedef void(*funcT)(double, double *, double *, void *);

        funcT f;
        size_t node;
        T *work, *yn, *dy, *yn1, *yn1_c, *ftmp, px[13], *A[13];
        T hn, tn, atol, rtol;
        int feval, fails, totalSteps, steps;
        T xhistory[13], *yhistory;  // might change to double*[13]
        size_t nhistory;
        size_t pnx, qnx;
        void *args;

    private:
        /*
            Newton divided difference polynomial
            no "solving" for coefficients
            only appending nodes and removing first node will be available
            adding Nodes will take O(ns) flops n is number of dimensions, s number of points
            removing first node will be O(1) time (just a few move instructions)
            same operator()
            new correction method
            no need for y vals since they will be in A
        */
        struct NewtonPolynomial
        {
            T *x, **A;
            size_t &nx, ny;

        public:
            NewtonPolynomial(size_t &n, size_t m, RP(T) xk, T **A)
                : x(xk), A(A), nx(n), ny(m)
            {
            }

            // quickly make a "copy" polynomial of degree n - k
            LOVERO_INLINE void referenceWithoutFirst(NewtonPolynomial &other, size_t k) const
            {
                other.x = x + k;
                other.A = A + k;
                other.nx = nx - k;
                other.ny = ny;
            }

            LOVERO_INLINE void addNode(const T xn, const T *yn)
            {
                ptrdiff_t i, j, k;
                x[nx] = xn;
                memcpy(A[nx], yn, ny * sizeof(T));

                j = 0;
                for (i = nx - 1; i > -1; i--) {
                    for (k = 0; k < (ptrdiff_t)ny; k++) {
                        A[i][(j + 1) * ny + k] = (A[i + 1][j * ny + k] - A[i][j * ny + k]) / (xn - x[i]);
                    }
                    j++;
                }
                nx++;
            }

            LOVERO_INLINE void removeLastNode()
            {
                nx--;
            }

            // only to be used at end of step (not for lower order polys)
            LOVERO_INLINE void removeFirstNode()
            {
                // since this class is only to be used in context
                // of solver, assume 13 points are always allocated
                size_t i;
                T *first = *A;
                nx--;
                for (i = 0; i < 12; i++) {
                    A[i] = A[i + 1];
                    x[i] = x[i + 1];
                }
                A[12] = first; // don't want memory leak and UB
            }

            QUAD_INLINE void operator() (T t, RP(T) res) const
            {
                size_t k, i;
                memset(res, 0, ny * sizeof(T));
                const T *a = *A;

                for (k = nx; k-- > 1; ) {
                    for (i = 0; i < ny; i++) {
                        res[i] += a[k * ny + i];
                        res[i] *= t - x[k - 1];
                    }
                }
                for (i = 0; i < ny; i++) {
                    res[i] += a[i];
                }
            }

            LOVERO_INLINE void integrate(T a, T b, RP(T) res) const
            {
                ptrdiff_t i, s = 0, m, n;
                bool r;
                T Q1, Q2;
                const T *c = *A;

                n = nx & (-2);
                r = (ptrdiff_t)nx != n;
                memset(res, 0, ny * sizeof(T));

                auto p = [&](T t) __attribute__((always_inline))
                {
                    ptrdiff_t j;
                    T tmp = 1;
                    for (j = 0; j < (s - 1); j++) {
                        tmp *= (t - x[j]);
                    }
                    return tmp;
                };

                for (s = 0; s < n; s += 2) {
                    // integrate in pairs since they require same number of nodes
                    s++;
                    m = (s >> 1) + (s & 1);
                    switch (m) {
                        case 1:
                            Q1 = quadGLW<gaussXK1, gaussWK1>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 2:
                            Q1 = quadGLW<gaussXK2, gaussWK2>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 3:
                            Q1 = quadGLW<gaussXK3, gaussWK3>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 4:
                            Q1 = quadGLW<gaussXK4, gaussWK4>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 5:
                            Q1 = quadGLW<gaussXK5, gaussWK5>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 6:
                            Q1 = quadGLW<gaussXK6, gaussWK6>(p, a, b, Q2, x[s - 1]);
                            break;
                        case 7:
                            Q1 = quadGLW<gaussXK7, gaussWK7>(p, a, b, Q2, x[s - 1]);
                            break;
                    }
                    s--;
                    for (i = 0; i < (ptrdiff_t)ny; i++) {
                        res[i] += (c[s * ny + i] * Q1) + (c[(s + 1) * ny + i] * Q2);
                    }
                }
                if (r) {
                    // if there was an odd number, integrate and add to result
                    s = nx;
                    m = (s >> 1) + (s & 1);
                    switch (m) {
                        case 1:
                            Q1 = quadG<gaussXK1, gaussWK1>(p, a, b);
                            break;
                        case 2:
                            Q1 = quadG<gaussXK2, gaussWK2>(p, a, b);
                            break;
                        case 3:
                            Q1 = quadG<gaussXK3, gaussWK3>(p, a, b);
                            break;
                        case 4:
                            Q1 = quadG<gaussXK4, gaussWK4>(p, a, b);
                            break;
                        case 5:
                            Q1 = quadG<gaussXK5, gaussWK5>(p, a, b);
                            break;
                        case 6:
                            Q1 = quadG<gaussXK6, gaussWK6>(p, a, b);
                            break;
                        case 7:
                            Q1 = quadG<gaussXK7, gaussWK7>(p, a, b);
                            break;
                    }
                    s--;
                    for (i = 0; i < (ptrdiff_t)ny; i++) {
                        res[i] += c[s * ny + i] * Q1;
                    }
                }
            }

            /*
                Corrector is of form p_{n+1} = p_{n} + a_{n+1}prod(x - x_k) = p_{n} + a_{n+1} q
                only need to integrate q to get Q = int q(t)
                instead of O(s^2 n) work, this will be O(s^2 + n) work, a HUGE savings for large
                systems.
            */
            LOVERO_INLINE T integrateCorrector(T a, T b, RP(T) dy, const RP(T) yn1, RP(T) yn1_c) const
            {
                const T *c = *A;
                const size_t m = (nx >> 1) + (nx & 1);

                T Q;
                auto p = [&, this](T t) __attribute__((always_inline))
                {
                    size_t j;
                    T tmp = (t - x[0]);
                    for (j = 1; j < (nx - 1); j++) {
                        tmp *= (t - x[j]);
                    }
                    return tmp;
                };

                switch (m) {
                    case 1:
                        Q = quadG<gaussXK1, gaussWK1>(p, a, b);
                        break;
                    case 2:
                        Q = quadG<gaussXK2, gaussWK2>(p, a, b);
                        break;
                    case 3:
                        Q = quadG<gaussXK3, gaussWK3>(p, a, b);
                        break;
                    case 4:
                        Q = quadG<gaussXK4, gaussWK4>(p, a, b);
                        break;
                    case 5:
                        Q = quadG<gaussXK5, gaussWK5>(p, a, b);
                        break;
                    case 6:
                        Q = quadG<gaussXK6, gaussWK6>(p, a, b);
                        break;
                    case 7:
                        Q = quadG<gaussXK7, gaussWK7>(p, a, b);
                        break;
                }

                for (size_t j = 0; j < ny; j++) {
                    dy[j] = c[(nx - 1) * ny + j] * Q;
                    yn1_c[j] = yn1[j] + dy[j];
                }
                return Q;
            }
        };

    public:
        LoveroNS(funcT f, size_t node, const RP(T) y0, T t0, T h0 = 0e-4, T atol = 1.49011612e-08, T rtol = 1.49011612e-08, void *args = nullptr)
            : f(f), node(node), hn(h0), tn(t0), atol(atol), rtol(rtol), feval(1), fails(0),
            totalSteps(0), steps(0), pnx(0), qnx(0), args(args)
        {
            // std::cout << __PRETTY_FUNCTION__ << '\n';
            size_t nwork = sizeof(T) * (node * (0
                + 1 // yn
                + 1 // dy
                + 1 // yn1
                + 1 // yn1_c
                + 13 // yhistory
                + 1 // ftmp
                + 13 * 13 // A
                ));

            // allocate and assign
            work = (T *)malloc(nwork); // need pointer to head of work array for free
            yn = work;  // yn is head of work -> free yn
            dy = yn + node;
            yn1 = dy + node;
            yn1_c = yn1 + node;
            yhistory = yn1_c + node;
            ftmp = yhistory + 13 * node;

            // this is overallocated but drastically simplifies the code
            T *head = ftmp + node;
            for (size_t i = 0; i < 13; i++) {
                A[i] = head;
                head += 13 * node;
            }

            // init
            *xhistory = tn;
            nhistory = 1;

            memcpy(yn, y0, node * sizeof(T));
            f(tn, yn, ftmp, args);
            NewtonPolynomial p(pnx, node, px, A);
            p.addNode(tn, ftmp);

            if (hn == 0) {
                hn = ::initStepSize(f, tn, yn, ftmp, yn1, yn1_c, rtol, atol, 1.0, 1.0, node, args);
            }
        }

        ~LoveroNS()
        {
            free(work);
        }

        inline void errorNorm(T &norm, T &hn1, const RP(T) dy, const RP(T) yn, const RP(T) yn1, size_t order) const
        {
            size_t i;
            T fac, tmp, scale, err;

            err = 0.0;
            for (i = 0; i < node; i++) {
                scale = atol + rtol * std::max(std::abs(yn[i]), std::abs(yn1[i]));
                tmp = dy[i] / scale;
                err += tmp * tmp;
            }
            norm = err = std::sqrt(err / node);
            fac = err ? 0.9 * std::pow(err, -1.0 / (order + 1)) : 10.0;
            fac = std::max<T>(0.2, std::min<T>(10.0, fac));
            hn1 = hn * fac;
        }

        // helper routine for adding this->yn to a vector x
        inline void addYnTo(RP(T) x) const
        {
            for (size_t i = 0; i < node; i++) {
                x[i] += yn[i];
            }
        }

        // helper routine for comuting the difference of two vectors
        inline void difference(RP(T) c, const RP(T) a, const RP(T) b) const
        {
            for (size_t i = 0; i < node; i++) {
                c[i] = a[i] - b[i];
            }
        }

        inline T nonstiffStep(T &norm0, T &ratio, bool &failed, const T direction)
        {
            T tn1, hn1_0, Q;
            NewtonPolynomial p(pnx, node, px, A);

            failed = false;

            while (1) {
                tn1 = tn + direction * hn;
                // O(sn + s^2) work
                p.integrate(tn, tn1, yn1);
                addYnTo(yn1);

                f(tn1, yn1, ftmp, args);
                feval++;
                // O(sn) work
                p.addNode(tn1, ftmp);

                // keeping track of Q leads to big reduction in repeated work 
                // O(s^2 + n) work
                Q = p.integrateCorrector(tn, tn1, dy, yn1, yn1_c);

                errorNorm(norm0, hn1_0, dy, yn, yn1_c, pnx - 1);
                ratio = hn1_0 / hn;
                hn = hn1_0;

                if (norm0 < 1) {
                    p.removeLastNode();
                    f(tn1, yn1_c, ftmp, args);
                    feval++;
                    p.addNode(tn1, ftmp);
                    break;
                }
                else {
                    T *a = &A[0][(pnx - 1) * node];
                    for (int i = 0; i < (failed ? 2 : 1); i++) {
                        std::swap(yn1, yn1_c); // just swap pointers instead of memcpy

                        if (i == 1) {
                            memcpy(yn1_c, a, node * sizeof(T)); // store previous coeff vals
                        }

                        p.removeLastNode();
                        f(tn1, yn1, ftmp, args);
                        feval++;
                        p.addNode(tn1, ftmp);

                        if (i == 0) {
                            for (size_t j = 0; j < node; j++) {
                                yn1_c[j] = yn1[j] + a[j] * Q - dy[j];
                                dy[j] = yn1_c[j] - yn1[j];
                            }
                        }
                        else {
                            for (size_t j = 0; j < node; j++) {
                                // yn1_c will be holding previous iteration's a values
                                yn1_c[j] = yn1[j] + (a[j] - yn1_c[j]) * Q;
                                dy[j] = yn1_c[j] - yn1[j];
                            }
                        }

                        errorNorm(norm0, hn1_0, dy, yn, yn1_c, pnx - 1);
                    }
                    if (norm0 < 1) {
                        break;
                    }
                }
                fails++;
                failed = true;
                p.removeLastNode();
            }

            // printf("%3d : %.15f | %s\n", totalSteps, norm0, failed ? "failed" : "passed");
            return tn1;
        }

        void orderSwitchNonstiff(T tn1, T norm0, T ratio)
        {
            NewtonPolynomial p(pnx, node, px, A);

            int order = p.nx - 1;
            int orderm1 = order - 1;
            int orderm2 = order - 2;
            T norm_m1, norm_m2, hn1_m1, hn1_m2;

            if (0 < orderm1) {
                NewtonPolynomial q(qnx, node, nullptr, nullptr);
                p.referenceWithoutFirst(q, 1);

                // yn1 holds yn1 order - 1
                q.integrate(tn, tn1, yn1);
                addYnTo(yn1);
                difference(dy, yn1_c, yn1);
                errorNorm(norm_m1, hn1_m1, dy, yn, yn1, orderm1);

                if (1 <= norm_m1) {
                    norm_m1 = std::numeric_limits<T>::infinity();
                    hn1_m1 = 0.0;
                }

                if (0 < orderm2) {
                    p.referenceWithoutFirst(q, 2);

                    // yn1 holds yn1 order - 2
                    q.integrate(tn, tn1, yn1);
                    addYnTo(yn1);
                    difference(dy, yn1_c, yn1);
                    errorNorm(norm_m2, hn1_m2, dy, yn, yn1, orderm2);

                    if (1 <= norm_m2) {
                        norm_m2 = std::numeric_limits<T>::infinity();
                        hn1_m2 = 0.0;
                    }
                }
                else {
                    norm_m2 = std::numeric_limits<T>::infinity();
                    hn1_m2 = 0.0;
                }
            }
            else {
                norm_m1 = std::numeric_limits<T>::infinity();
                norm_m2 = std::numeric_limits<T>::infinity();
                hn1_m1 = 0.0;
                hn1_m2 = 0.0;
            }

            if (std::max(norm_m1, norm_m2) <= (0.5 * norm0)) {
                T tmp = std::max(hn1_m1, hn1_m2);
                T hprev = hn / ratio;
                ratio = tmp / hprev;
                hn = tmp;

                if ((2 * order + 2) <= steps) {
                    p.removeFirstNode();
                    p.removeFirstNode();
                    steps = 0;
                }
                else {
                    p.removeFirstNode();
                }
            }
            else if (11 < order) {
                p.removeFirstNode();
            }
            else if (steps < order) {
                p.removeFirstNode();
            }
            else if (ratio < 0.9) {
                p.removeFirstNode();
            }
            else if (!norm0 && (1 < order)) {
                p.removeFirstNode();
                p.removeFirstNode();
                steps = 0;
            }
            else {
                steps = 0;
            }
        }

        bool step(const T direction = 1.0)
        {
            T norm0, ratio, tn1;
            bool failed;

            steps++;
            totalSteps++;

            constexpr bool stiff = false;

            if (stiff) {
                // stiff step
            }
            else {
                tn1 = nonstiffStep(norm0, ratio, failed, direction);
            }

            // method switch

            if (stiff) {
                // order switch stiff
            }
            else {
                orderSwitchNonstiff(tn1, norm0, ratio);
            }

            if (nhistory == 13) {
                memcpy(xhistory, xhistory + 1, (nhistory - 1) * sizeof(T));
                memcpy(yhistory, yhistory + node, node * (nhistory - 1) * sizeof(T));
                xhistory[nhistory - 1] = tn1;
                memcpy(yhistory + (nhistory - 1) * node, yn1_c, node * sizeof(T));
            }
            else {
                xhistory[nhistory] = tn1;
                memcpy(yhistory + nhistory * node, yn1_c, node * sizeof(T));
                nhistory++;
            }

            tn = tn1;
            std::swap(yn, yn1_c); // no memcpy
            return failed;
        }

        ODEFailCodes preStepChecks() const
        {
            size_t i;

            if (hn < (std::nextafter(tn, std::numeric_limits<T>::infinity()) - tn)) {
                return ODEFailCodes::stepTooSmall;
            }

            for (i = 0; i < node; i++) {
                if (!std::isfinite(yn[i]) || !std::isfinite(ftmp[i])) {
                    return ODEFailCodes::nonFiniteState;
                }
            }

            return ODEFailCodes::noFailures;
        }

        ptrdiff_t stepWithCheckpoints(const RP(T) teval, RP(T) yeval, const ptrdiff_t neval)
        {
            T norm0, ratio, tn1, ttarget, direction;
            ptrdiff_t k;
            bool failed;

            if (neval == 0) {
                return 0;
            }
            if (ODEFailCodes code = preStepChecks()) {
                return -code;
            }

            ttarget = *teval;
            direction = (ttarget < tn) ? -1.0 : 1.0;
            k = 0;

            steps++;
            totalSteps++;

            constexpr bool stiff = false;

            if (stiff) {
                // stiff step
            }
            else {
                tn1 = nonstiffStep(norm0, ratio, failed, direction);
            }

            // evaluate at checkpoints
            while ((k < neval) && (0 <= direction * (tn1 - ttarget))) {
                NewtonPolynomial(pnx, node, px, A).integrate(tn, ttarget, &yeval[node * k]);
                addYnTo(&yeval[node * k]);
                k++;
                if (k < neval) {
                    ttarget = teval[k];
                }
            }

            // method switch

            if (stiff) {
                // order switch stiff
            }
            else {
                orderSwitchNonstiff(tn1, norm0, ratio);
            }

            if (nhistory == 13) {
                memcpy(xhistory, xhistory + 1, (nhistory - 1) * sizeof(T));
                memcpy(yhistory, yhistory + node, node * (nhistory - 1) * sizeof(T));
                xhistory[nhistory - 1] = tn1;
                memcpy(yhistory + (nhistory - 1) * node, yn1_c, node * sizeof(T));
            }
            else {
                xhistory[nhistory] = tn1;
                memcpy(yhistory + nhistory * node, yn1_c, node * sizeof(T));
                nhistory++;
            }

            tn = tn1;
            std::swap(yn, yn1_c); // no memcpy
            return k;
        }
    };

    template <typename T, typename funcT>
    ODEResult<T> integrate(funcT fn, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, RP(void) fargs)
    {
        const T Dxy = max_norm[0];
        const T *centers = &max_norm[1];

        ptrdiff_t diff, k, i, j;
        ptrdiff_t neval = m, offset = 1;
        T val, dst, tf, direction;
        T *yj;

        memcpy(Y, y0, n * sizeof(T));
        diff = neval - offset;
        tf = teval[m - 1];
        direction = (tf < *teval) ? -1.0 : 1.0;

        LoveroNS<T, funcT> stepper(fn, n, y0, *teval, h0, atol, rtol, fargs);

        while (0 < (direction * (tf - stepper.tn))) {
            k = stepper.stepWithCheckpoints(&teval[offset], &Y[offset * n], diff);
            if (k < 0) {
                // pre step checks failed
                k = -k;
                switch (k) {
                    case ODEFailCodes::stepTooSmall:
                        return { (size_t)(offset), (size_t)stepper.feval, 0, (size_t)stepper.fails, (size_t)stepper.totalSteps, 0, ODEExitCodes::stepSizeTooSmall, stepper.hn, stepper.tn };
                    case ODEFailCodes::nonFiniteState:
                        return { (size_t)(offset), (size_t)stepper.feval, 0, (size_t)stepper.fails, (size_t)stepper.totalSteps, 0, ODEExitCodes::outOfBounds, stepper.hn, stepper.tn };
                    default:
                        assert((0 < k) && (k < 3));
                }
            }

            if (std::isfinite(Dxy)) {
                // do bounds check if Dxy is bounded
                for (j = 0; j < k; j++) {
                    yj = &Y[(offset + j)];
                    dst = 0;
                    for (i = 0; i < (ptrdiff_t)n; i++) {
                        val = yj[i] - centers[i];
                        dst += val * val;
                    }
                    if ((Dxy < dst) || !std::isfinite(dst)) {
                        return { (size_t)(offset + j), (size_t)stepper.feval, 0, (size_t)stepper.fails, (size_t)stepper.totalSteps, 0, ODEExitCodes::outOfBounds, stepper.hn, stepper.tn };
                    }
                }
            }
            offset += k;
            diff -= k;
        }
        return { m, (size_t)stepper.feval, 0, (size_t)stepper.fails, (size_t)stepper.totalSteps, 0, ODEExitCodes::success, stepper.hn, stepper.tn };
    }

}

#endif // ODEPACK_LOVERO_NS2_H