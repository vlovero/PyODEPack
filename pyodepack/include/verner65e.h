/*
 * Implementation of Verner's efficient 6(5) RK pair described in [1].
 *
 * Vincent Lovero
 *
 * [1] Verner, J H (2013). Explicit Runge Kutta pairs with lower stage-order. Numerical Algorithms 65(3): 555–577
 */

#ifndef ODEPACK_VERNER65E_H
#define ODEPACK_VERNER65E_H

#include "odepack_common.h"
#include "rk.h"

namespace verner65e
{
    // nodes
    constexpr double c1 = 0.0;
    constexpr double c2 = 0.6e-1;
    constexpr double c3 = 0.9593333333333333333333333333333333333333e-1;
    constexpr double c4 = 0.1439;
    constexpr double c5 = 0.4973;
    constexpr double c6 = 0.9725;
    constexpr double c7 = 0.9995;
    constexpr double c8 = 1.0;
    constexpr double c9 = 1.0;
    constexpr double c10 = 0.5;   // for 5th order interpolant
    constexpr double c11 = 0.828; // for 6th order interpolant
    constexpr double c12 = 0.28;  // for 6th order interpolant

    // high order weights
    constexpr double b1 = 0.3438957868357036009278820124728322386520e-1;
    constexpr double b2 = 0.0;
    constexpr double b3 = 0.0;
    constexpr double b4 = 0.2582624555633503404659558098586120858767;
    constexpr double b5 = 0.4209371189673537150642551514069801967032;
    constexpr double b6 = 4.405396469669310170148836816197095664891;
    constexpr double b7 = -176.4831190242986576151740942499002125029;
    constexpr double b8 = 172.3641334014150730294022582711902413315;
    constexpr double b9 = 0.0;

    // (low - high) order weights
    constexpr double d1 = 0.01471009780025453721628034803242903449959;
    constexpr double d2 = 0;
    constexpr double d3 = 0;
    constexpr double d4 = -0.03315123261169792512581627780462455292820;
    constexpr double d5 = 0.04853110633560248887893970109775854455210;
    constexpr double d6 = -3.598817244670423399385420635297573866747;
    constexpr double d7 = 176.4831190242986576151740942499002125029;
    constexpr double d8 = -172.9712528905928690091695534177158630437;
    constexpr double d9 = 0.05686113944047569241147603178766138153594;

    // coupling coefficients
    constexpr double a2_1 = 0.6e-1;
    constexpr double a3_1 = 0.1923996296296296296296296296296296296296e-1;
    constexpr double a3_2 = 0.7669337037037037037037037037037037037037e-1;
    constexpr double a4_1 = 0.35975e-1;
    constexpr double a4_2 = 0.0;
    constexpr double a4_3 = 0.107925;
    constexpr double a5_1 = 1.318683415233148260919747276431735612861;
    constexpr double a5_2 = 0.0;
    constexpr double a5_3 = -5.042058063628562225427761634715637693344;
    constexpr double a5_4 = 4.220674648395413964508014358283902080483;
    constexpr double a6_1 = -41.87259166432751461803757780644346812905;
    constexpr double a6_2 = 0.0;
    constexpr double a6_3 = 159.4325621631374917700365669070346830453;
    constexpr double a6_4 = -122.1192135650100309202516203389242140663;
    constexpr double a6_5 = 5.531743066200053768252631238332999150076;
    constexpr double a7_1 = -54.43015693531650433250642051294142461271;
    constexpr double a7_2 = 0.0;
    constexpr double a7_3 = 207.0672513650184644273657173866509835987;
    constexpr double a7_4 = -158.6108137845899991828742424365058599469;
    constexpr double a7_5 = 6.991816585950242321992597280791793907096;
    constexpr double a7_6 = -0.1859723106220323397765171799549294623692e-1;
    constexpr double a8_1 = -54.66374178728197680241215648050386959351;
    constexpr double a8_2 = 0.0;
    constexpr double a8_3 = 207.9528062553893734515824816699834244238;
    constexpr double a8_4 = -159.2889574744995071508959805871426654216;
    constexpr double a8_5 = 7.018743740796944434698170760964252490817;
    constexpr double a8_6 = -0.1833878590504572306472782005141738268361e-1;
    constexpr double a8_7 = -0.5119484997882099077875432497245168395840e-3;
    constexpr double a9_1 = 0.3438957868357036009278820124728322386520e-1;
    constexpr double a9_2 = 0.0;
    constexpr double a9_3 = 0.0;
    constexpr double a9_4 = 0.2582624555633503404659558098586120858767;
    constexpr double a9_5 = 0.4209371189673537150642551514069801967032;
    constexpr double a9_6 = 4.405396469669310170148836816197095664891;
    constexpr double a9_7 = -176.4831190242986576151740942499002125029;
    constexpr double a9_8 = 172.3641334014150730294022582711902413315;
    constexpr double a10_1 = 0.1652415901357280684383619367363197852645e-1;
    constexpr double a10_2 = 0.0;
    constexpr double a10_3 = 0.0;
    constexpr double a10_4 = 0.3053128187514178931377105638345583032476;
    constexpr double a10_5 = 0.2071200938201978848991082158995582390341;
    constexpr double a10_6 = -1.293879140655123187129665774327355723229;
    constexpr double a10_7 = 57.11988411588149149650780257779402737914;
    constexpr double a10_8 = -55.87979207510932290773937033203265749155;
    constexpr double a10_9 = 0.2483002829776601348057855515823731483430e-1;
    constexpr double a11_1 = 0.3815008181862774566390122945290392930051e-1;
    constexpr double a11_2 = 0.0;
    constexpr double a11_3 = 0.0;
    constexpr double a11_4 = 0.2502358252513705267411717435114864276070;
    constexpr double a11_5 = 0.3249441447817607965203943927614998837766;
    constexpr double a11_6 = 1.822460665832796218309507694895413940004;
    constexpr double a11_7 = -67.71372332692620404266416384783466979895;
    constexpr double a11_8 = 66.03587911808127345206414399075146898997;
    constexpr double a11_9 = -0.3638810874951269663495520353810337170236e-1;
    constexpr double a11_10 = 0.1064415999098880000000000000000000000000;
    constexpr double a12_1 = 0.1117816803966601180374383193782622585818;
    constexpr double a12_2 = 0.0;
    constexpr double a12_3 = 0.0;
    constexpr double a12_4 = 0.2575750510934521211743051533509279738017e-1;
    constexpr double a12_5 = 3.785140856363645858969575735188195099293;
    constexpr double a12_6 = 92.34088993695727177430761424880372341560;
    constexpr double a12_7 = -3819.461508432344125782385301764384773807;
    constexpr double a12_8 = 3732.492711530704227040059052662831877264;
    constexpr double a12_9 = -1.075694020996303327721922083840048543815;
    constexpr double a12_10 = -3.231539970732086111815287724591251178493;
    constexpr double a12_11 = -4.707539085458634781568599908721077305888;

    // 
    constexpr double p1_1 = 1.0;
    constexpr double p1_2 = -5.308169607103576297743491917539437544903;
    constexpr double p1_3 = 10.18168044895868030520877351032733768603;
    constexpr double p1_4 = -7.520036991611714828300683961994073691563;
    constexpr double p1_5 = 0.9340485368631160925057442706475838478288;
    constexpr double p1_6 = 0.7468671915770650884224462998058729264688;
    constexpr double p2_1 = 0.0;
    constexpr double p2_2 = 0.0;
    constexpr double p2_3 = 0.0;
    constexpr double p2_4 = 0.0;
    constexpr double p2_5 = 0.0;
    constexpr double p2_6 = 0.0;
    constexpr double p3_1 = 0.0;
    constexpr double p3_2 = 0.0;
    constexpr double p3_3 = 0.0;
    constexpr double p3_4 = 0.0;
    constexpr double p3_5 = 0.0;
    constexpr double p3_6 = 0.0;
    constexpr double p4_1 = 0.0;
    constexpr double p4_2 = 6.272050253212501244827865529084399503479;
    constexpr double p4_3 = -16.02618147467745958442607061022576892601;
    constexpr double p4_4 = 12.84435632451961742214954703737612797249;
    constexpr double p4_5 = -1.148794504476759027536609501260874665600;
    constexpr double p4_6 = -1.683168143014549714548776645115271798480;
    constexpr double p5_1 = 0.0;
    constexpr double p5_2 = 6.876491702846304590450466371720363234704;
    constexpr double p5_3 = -24.63576726084633318864583120149461053641;
    constexpr double p5_4 = 33.21078648379717088772133447477731248517;
    constexpr double p5_5 = -17.49461528263643828092150992351036511970;
    constexpr double p5_6 = 2.464041475806649706459795429914280132942;
    constexpr double p6_1 = 0.0;
    constexpr double p6_2 = -35.54445171059960218765875699270358093032;
    constexpr double p6_3 = 165.7016170190242105108446269288474846144;
    constexpr double p6_4 = -385.4635395491142731464726480659809841649;
    constexpr double p6_5 = 442.4324137015701845319394642134164121973;
    constexpr double p6_6 = -182.7206429912112095385038492673822360516;
    constexpr double p7_1 = 0.0;
    constexpr double p7_2 = 1918.654856698011449175045220651610014945;
    constexpr double p7_3 = -9268.121508966042500195164712930044037430;
    constexpr double p7_4 = 20858.33702877255011893787944928058522511;
    constexpr double p7_5 = -22645.82767158481047968149020787687967272;
    constexpr double p7_6 = 8960.474176055992754148556156624828257597;
    constexpr double p8_1 = 0.0;
    constexpr double p8_2 = -1883.069802132718312960582779305006646556;
    constexpr double p8_3 = 9101.025187200633795903395040749471730528;
    constexpr double p8_4 = -20473.18855195953591834830509979123557878;
    constexpr double p8_5 = 22209.76555125653413900516974418122400018;
    constexpr double p8_6 = -8782.168250963498630570274647563263264047;
    constexpr double p9_1 = 0.0;
    constexpr double p9_2 = 0.1190247963512364356614756628348873476063;
    constexpr double p9_3 = -0.1250269670503937512118264468821359362429;
    constexpr double p9_4 = 1.779956919394999075328101026471971070697;
    constexpr double p9_5 = -4.660932123043762639666625363637083723091;
    constexpr double p9_6 = 2.886977374347920879888875121212361241030;
    constexpr double p10_1 = 0.0;
    constexpr double p10_2 = -8.0;
    constexpr double p10_3 = 32.0;
    constexpr double p10_4 = -40.0;
    constexpr double p10_5 = 16.0;
    constexpr double p10_6 = 0.0;


    template <typename T>
    constexpr std::array<T, 60> _getOrder5Interpolant()
    {
        return { p1_1, p1_2, p1_3, p1_4, p1_5, p1_6,
                 p2_1, p2_2, p2_3, p2_4, p2_5, p2_6,
                 p3_1, p3_2, p3_3, p3_4, p3_5, p3_6,
                 p4_1, p4_2, p4_3, p4_4, p4_5, p4_6,
                 p5_1, p5_2, p5_3, p5_4, p5_5, p5_6,
                 p6_1, p6_2, p6_3, p6_4, p6_5, p6_6,
                 p7_1, p7_2, p7_3, p7_4, p7_5, p7_6,
                 p8_1, p8_2, p8_3, p8_4, p8_5, p8_6,
                 p9_1, p9_2, p9_3, p9_4, p9_5, p9_6,
                 p10_1, p10_2, p10_3, p10_4, p10_5, p10_6 };
    }

    template <typename T>
    constexpr std::array<T, 90> _getCoeffsForOrder5()
    {
        return { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,
                 a2_1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,
                 a3_1, a3_2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,
                 a4_1, a4_2, a4_3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,
                 a5_1, a5_2, a5_3, a5_4, 0.0, 0.0, 0.0, 0.0, 0.0 ,
                 a6_1, a6_2, a6_3, a6_4, a6_5, 0.0, 0.0, 0.0, 0.0 ,
                 a7_1, a7_2, a7_3, a7_4, a7_5, a7_6, 0.0, 0.0, 0.0 ,
                 a8_1, a8_2, a8_3, a8_4, a8_5, a8_6, a8_7, 0.0, 0.0 ,
                 a9_1, a9_2, a9_3, a9_4, a9_5, a9_6, a9_7, a9_8, 0.0 ,
                 a10_1, a10_2, a10_3, a10_4, a10_5, a10_6, a10_7, a10_8, a10_9 };
    }

    template <typename T>
    constexpr std::array<T, 10> _getNodesForOrder5()
    {
        return { c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 };
    }

    template <typename T>
    constexpr std::array<T, 9> _getWeights()
    {
        return { b1, b2, b3, b4, b5, b6, b7, b8, b9 };
    }

    template <typename T>
    constexpr std::array<T, 9> _getDeltas()
    {
        return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
    }

    constexpr std::array<double, 60> interpolant_f64 = _getOrder5Interpolant<double>();
    constexpr std::array<float, 60> interpolant_f32 = _getOrder5Interpolant<float>();

    constexpr std::array<double, 90> coefficients_f64 = _getCoeffsForOrder5<double>();
    constexpr std::array<float, 90> coefficients_f32 = _getCoeffsForOrder5<float>();

    constexpr std::array<double, 9> weights_f64 = _getWeights<double>();
    constexpr std::array<float, 9> weights_f32 = _getWeights<float>();

    constexpr std::array<double, 9> deltas_f64 = _getDeltas<double>();
    constexpr std::array<float, 9> deltas_f32 = _getDeltas<float>();

    constexpr std::array<double, 10> nodes_f64 = _getNodesForOrder5<double>();
    constexpr std::array<float, 10> nodes_f32 = _getNodesForOrder5<float>();


    template <typename T, typename funcT, const bool scalar = false>
    inline ODEResult integrate(funcT fn, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value, "Only types float and double are currently supported");

        if constexpr (std::is_same<T, double>::value) {
            return erk::integrate<T, funcT, 9, 10, 6,
                nodes_f64,
                9, coefficients_f64,
                interpolant_f64,
                weights_f64,
                deltas_f64, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
        else if constexpr (std::is_same<T, float>::value) {
            return erk::integrate<T, funcT, 9, 10, 6,
                nodes_f32,
                9, coefficients_f32,
                interpolant_f32,
                weights_f32,
                deltas_f32, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
    }
}

#endif // ODEPACK_VERNER65E_H