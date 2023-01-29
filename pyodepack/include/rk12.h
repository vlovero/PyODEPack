#ifndef ODEPACK_RK12_H
#define ODEPACK_RK12_H

#include "odepack_common.h"
#include "rk.h"

namespace Fehlberg12
{
    template <typename T>
    constexpr std::array<T, 9> _getInterpolant()
    {
        return { 1.0, -1.0, 1.0 / 512,
                 0.0, 1.0, -1.0 / 256,
                 0.0, 0.0, 1.0 / 512 };
        // return { 5.0 / 6.0, -1.0 / 2.0,
        //          1.0 / 3.0, 0.0,
        //          -1.0 / 6.0, 1.0 / 2.0 };
    }

    template <typename T>
    constexpr std::array<T, 3> _getNodes()
    {
        return { 0.0, 0.5, 1.0 };
    }

    template <typename T>
    constexpr std::array<T, 9> _getCoeffs()
    {
        return { 0.0, 0.0, 0.0,
                 0.5, 0.0, 0.0,
                 1.0 / 256, 255.0 / 256, 0.0 };
    }

    template <typename T>
    constexpr std::array<T, 3> _getWeights()
    {
        return { 1.0 / 512, 255.0 / 256, 1.0 / 512 };
    }

    template <typename T>
    constexpr std::array<T, 3> _getDeltas()
    {
        return { 1.0 / 512, 0.0, -1.0 / 512 };
    }

    constexpr std::array<double, 9> interpolant_f64 = _getInterpolant<double>();
    constexpr std::array<float, 9> interpolant_f32 = _getInterpolant<float>();

    constexpr std::array<double, 9> coefficients_f64 = _getCoeffs<double>();
    constexpr std::array<float, 9> coefficients_f32 = _getCoeffs<float>();

    constexpr std::array<double, 3> weights_f64 = _getWeights<double>();
    constexpr std::array<float, 3> weights_f32 = _getWeights<float>();

    constexpr std::array<double, 3> deltas_f64 = _getDeltas<double>();
    constexpr std::array<float, 3> deltas_f32 = _getDeltas<float>();

    constexpr std::array<double, 3> nodes_f64 = _getNodes<double>();
    constexpr std::array<float, 3> nodes_f32 = _getNodes<float>();

    template <typename T, typename funcT, const bool scalar = false>
    inline ODEResult<T> integrate(funcT fn, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value, "Only types float and double are currently supported");

        if constexpr (std::is_same<T, double>::value) {
            return erk::integrate<T, funcT, 3, 3, 3,
                nodes_f64,
                3, coefficients_f64,
                interpolant_f64,
                weights_f64,
                deltas_f64, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
        else if constexpr (std::is_same<T, float>::value) {
            return erk::integrate<T, funcT, 3, 3, 3,
                nodes_f32,
                3, coefficients_f32,
                interpolant_f32,
                weights_f32,
                deltas_f32, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
    }

    size_t getWorkArraySize(size_t node)
    {
        return (3 + 3 + 3) * node;
    }

    template <typename funcT, bool scalar = false>
    using stepper64_t = erk::RKStepper<double, funcT, 3, 3, 3,
        nodes_f64, 3, coefficients_f64,
        interpolant_f64,
        weights_f64,
        deltas_f64, scalar>;

    template <typename funcT, bool scalar = false>
    using stepper32_t = erk::RKStepper<float, funcT, 3, 3, 3,
        nodes_f32, 3, coefficients_f32,
        interpolant_f32,
        weights_f32,
        deltas_f32, scalar>;

    template <typename T, typename funcT, const bool scalar = false>
    inline auto step(funcT f, T tn, const T tf, const RP(T) y0, RP(T) y, const T atol, const T rtol,
        const T minStep, const T maxStep, const T initStep, const size_t node, const bool firstStep, const RP(T) window,
        RP(T) work, const RP(void) args)
    {
        static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value, "Only types float and double are currently supported");

        if constexpr (std::is_same<T, double>::value) {
            using stepper_t = erk::RKStepper<T, funcT, 3, 3, 3,
                nodes_f64, 3, coefficients_f64,
                interpolant_f64,
                weights_f64,
                deltas_f64, scalar>;

            stepper_t stepper;
            stepper.step(f, tn, tf, y0, y, atol, rtol, minStep, maxStep, initStep, node, firstStep, window, work, args);
            return stepper;
        }
        else if constexpr (std::is_same<T, float>::value) {
            using stepper_t = erk::RKStepper<T, funcT, 3, 3, 3,
                nodes_f32, 3, coefficients_f32,
                interpolant_f32,
                weights_f32,
                deltas_f32, scalar>;

            stepper_t stepper;
            stepper.step(f, tn, tf, y0, y, atol, rtol, minStep, maxStep, initStep, node, firstStep, window, work, args);
            return stepper;
        }
    }
}

#endif // ODEPACK_RK12_H