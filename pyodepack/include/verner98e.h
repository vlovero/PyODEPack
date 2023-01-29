/*
 * Implementation of Verner's efficient 9(8) RK pair described in [1].
 *
 * Vincent Lovero
 *
 * [1] Verner, J H (2013). Explicit Runge Kutta pairs with lower stage-order. Numerical Algorithms 65(3): 555–577
 */

#ifndef ODEPACK_VERNER98E_H
#define ODEPACK_VERNER98E_H

#include "odepack_common.h"
#include "rk.h"

namespace verner98e
{
    // nodes
    constexpr double c1 = 0.0;
    constexpr double c2 = 0.0346200000000000;
    constexpr double c3 = 0.09702435063878044594828361677100617517633;
    constexpr double c4 = 0.1455365259581706689224254251565092627645;
    constexpr double c5 = 0.561000000000000;
    constexpr double c6 = 0.2290079115904850126662751771814700052182;
    constexpr double c7 = 0.5449920884095149873337248228185299947818;
    constexpr double c8 = 0.645000000000000;
    constexpr double c9 = 0.4837500000000000000000000000000000000000;
    constexpr double c10 = 0.0675700000000000;
    constexpr double c11 = 0.250000000000000;
    constexpr double c12 = 0.6590650618730998549405331618649220295334;
    constexpr double c13 = 0.820600000000000;
    constexpr double c14 = 0.901200000000000;
    constexpr double c15 = 1.00000000000000;
    constexpr double c16 = 1.00000000000000;
    constexpr double c17 = 1.00000000000000;
    constexpr double c18 = 0.7404185470631561083004100761798676215811;
    constexpr double c19 = 0.888000000000000;
    constexpr double c20 = 0.696000000000000;
    constexpr double c21 = 0.487000000000000;
    constexpr double c22 = 0.0250000000000000;
    constexpr double c23 = 0.150000000000000;
    constexpr double c24 = 0.320000000000000;
    constexpr double c25 = 0.780000000000000;
    constexpr double c26 = 0.960000000000000;

    // high order weights
    constexpr double b1 = 0.01461197685842315252051541915018784713459;
    constexpr double b2 = 0.0;
    constexpr double b3 = 0.0;
    constexpr double b4 = 0.0;
    constexpr double b5 = 0.0;
    constexpr double b6 = 0.0;
    constexpr double b7 = 0.0;
    constexpr double b8 = -0.3915211862331339089410228267288242030810;
    constexpr double b9 = 0.2310932500289506415909675644868993669908;
    constexpr double b10 = 0.1274766769992852382560589467488989175618;
    constexpr double b11 = 0.2246434176204157731566981937082069688984;
    constexpr double b12 = 0.5684352689748512932705226972873692126743;
    constexpr double b13 = 0.05825871557215827200814768021863420902155;
    constexpr double b14 = 0.1364317403482215641609022744494239843327;
    constexpr double b15 = 0.03057013983082797397721005067920369646664;
    constexpr double b16 = 0.0;

    // (low - high) order weights
    constexpr double d1 = 0.005357988290444578334669665030800840430050;
    constexpr double d2 = 0;
    constexpr double d3 = 0;
    constexpr double d4 = 0;
    constexpr double d5 = 0;
    constexpr double d6 = 0;
    constexpr double d7 = 0;
    constexpr double d8 = 2.583020491182463963471769926039661727945;
    constexpr double d9 = -0.1425225315468662612826303441690357783613;
    constexpr double d10 = -0.01342065351268867600756328623975571429440;
    constexpr double d11 = 0.02867296291409493339975957974875822883630;
    constexpr double d12 = -2.624999655215792304429522291883350513167;
    constexpr double d13 = 0.2825509643291537215078618092038201722614;
    constexpr double d14 = -0.1364317403482215641609022744494239843327;
    constexpr double d15 = -0.03057013983082797397721005067920369646664;
    constexpr double d16 = 0.04834231373823958314376726739772871714902;

    // coupling coefficients
    constexpr double a2_1 = 0.0346200000000000;
    constexpr double a3_1 = -0.0389335438857287327017042687229284478532;
    constexpr double a3_2 = 0.1359578945245091786499878854939346230295;
    constexpr double a4_1 = 0.03638413148954266723060635628912731569111;
    constexpr double a4_2 = 0.0;
    constexpr double a4_3 = 0.1091523944686280016918190688673819470733;
    constexpr double a5_1 = 2.025763914393969636805657604282571047511;
    constexpr double a5_2 = 0.0;
    constexpr double a5_3 = -7.638023836496292020387602153091964592952;
    constexpr double a5_4 = 6.173259922102322383581944548809393545442;
    constexpr double a6_1 = 0.05112275589406060872792270881648288397197;
    constexpr double a6_2 = 0.0;
    constexpr double a6_3 = 0.0;
    constexpr double a6_4 = 0.1770823794555021537929910813839068684087;
    constexpr double a6_5 = 0.00080277624092225014536138698108025283759;
    constexpr double a7_1 = 0.1316006357975216279279871693164256985334;
    constexpr double a7_2 = 0.0;
    constexpr double a7_3 = 0.0;
    constexpr double a7_4 = -0.2957276252669636417685183174672273730699;
    constexpr double a7_5 = 0.0878137803564295237421124704053886667082;
    constexpr double a7_6 = 0.6213052975225274774321435005639430026100;
    constexpr double a8_1 = 0.07166666666666666666666666666666666666667;
    constexpr double a8_2 = 0.0;
    constexpr double a8_3 = 0.0;
    constexpr double a8_4 = 0.0;
    constexpr double a8_5 = 0.0;
    constexpr double a8_6 = 0.3305533578915319409260346730051472207728;
    constexpr double a8_7 = 0.2427799754418013924072986603281861125606;
    constexpr double a9_1 = 0.07180664062500000000000000000000000000000;
    constexpr double a9_2 = 0.0;
    constexpr double a9_3 = 0.0;
    constexpr double a9_4 = 0.0;
    constexpr double a9_5 = 0.0;
    constexpr double a9_6 = 0.3294380283228177160744825466257672816401;
    constexpr double a9_7 = 0.1165190029271822839255174533742327183599;
    constexpr double a9_8 = -0.03401367187500000000000000000000000000000;
    constexpr double a10_1 = 0.04836757646340646986611287718844085773549;
    constexpr double a10_2 = 0.0;
    constexpr double a10_3 = 0.0;
    constexpr double a10_4 = 0.0;
    constexpr double a10_5 = 0.0;
    constexpr double a10_6 = 0.03928989925676163974333190042057047002852;
    constexpr double a10_7 = 0.1054740945890344608263649267140088017604;
    constexpr double a10_8 = -0.02143865284648312665982642293830533996214;
    constexpr double a10_9 = -0.1041229174627194437759832813847147895623;
    constexpr double a11_1 = -0.02664561487201478635337289243849737340534;
    constexpr double a11_2 = 0.0;
    constexpr double a11_3 = 0.0;
    constexpr double a11_4 = 0.0;
    constexpr double a11_5 = 0.0;
    constexpr double a11_6 = 0.03333333333333333333333333333333333333333;
    constexpr double a11_7 = -0.1631072244872467239162704487554706387141;
    constexpr double a11_8 = 0.03396081684127761199487954930015522928244;
    constexpr double a11_9 = 0.1572319413814626097110769806810024118077;
    constexpr double a11_10 = 0.2152267478031879552303534778794770376960;
    constexpr double a12_1 = 0.03689009248708622334786359863227633989718;
    constexpr double a12_2 = 0.0;
    constexpr double a12_3 = 0.0;
    constexpr double a12_4 = 0.0;
    constexpr double a12_5 = 0.0;
    constexpr double a12_6 = -0.1465181576725542928653609891758501156785;
    constexpr double a12_7 = 0.2242577768172024345345469822625833796001;
    constexpr double a12_8 = 0.02294405717066072637090897902753790803034;
    constexpr double a12_9 = -0.0035850052905728761357394424889330334334;
    constexpr double a12_10 = 0.08669223316444385506869203619044453906053;
    constexpr double a12_11 = 0.4383840651968337846196219974168630120572;
    constexpr double a13_1 = -0.4866012215113340846662212357570395295088;
    constexpr double a13_2 = 0.0;
    constexpr double a13_3 = 0.0;
    constexpr double a13_4 = 0.0;
    constexpr double a13_5 = 0.0;
    constexpr double a13_6 = -6.304602650282852990657772792012007122988;
    constexpr double a13_7 = -0.281245618289472564778284183790118418111;
    constexpr double a13_8 = -2.679019236219849057687906597489223155566;
    constexpr double a13_9 = 0.518815663924157511565311164615012522024;
    constexpr double a13_10 = 1.365353187603341710683633635235238678626;
    constexpr double a13_11 = 5.885091088503946585721274891680604830712;
    constexpr double a13_12 = 2.802808786272062889819965117517532194812;
    constexpr double a14_1 = 0.4185367457753471441471025246471931649633;
    constexpr double a14_2 = 0.0;
    constexpr double a14_3 = 0.0;
    constexpr double a14_4 = 0.0;
    constexpr double a14_5 = 0.0;
    constexpr double a14_6 = 6.724547581906459363100870806514855026676;
    constexpr double a14_7 = -0.425444280164611790606983409697113064616;
    constexpr double a14_8 = 3.343279153001265577811816947557982637749;
    constexpr double a14_9 = 0.617081663117537759528421117507709784737;
    constexpr double a14_10 = -0.929966123939932833937749523988800852013;
    constexpr double a14_11 = -6.099948804751010722472962837945508844846;
    constexpr double a14_12 = -3.002206187889399044804158084895173690015;
    constexpr double a14_13 = 0.2553202529443445472336424602988558373637;
    constexpr double a15_1 = -0.779374086122884664644623040843840506343;
    constexpr double a15_2 = 0.0;
    constexpr double a15_3 = 0.0;
    constexpr double a15_4 = 0.0;
    constexpr double a15_5 = 0.0;
    constexpr double a15_6 = -13.93734253810777678786523664804936051203;
    constexpr double a15_7 = 1.252048853379357320949735183924200895136;
    constexpr double a15_8 = -14.69150040801686878191527989293072091588;
    constexpr double a15_9 = -0.494705058533141685655191992136962873577;
    constexpr double a15_10 = 2.242974909146236657906984549543692874755;
    constexpr double a15_11 = 13.36789380382864375813864978592679139881;
    constexpr double a15_12 = 14.39665048665068644512236935340272139005;
    constexpr double a15_13 = -0.7975813331776800379127866056663258667437;
    constexpr double a15_14 = 0.4409353709534277758753793068298041158235;
    constexpr double a16_1 = 2.058051337466886442151242368989994043993;
    constexpr double a16_2 = 0.0;
    constexpr double a16_3 = 0.0;
    constexpr double a16_4 = 0.0;
    constexpr double a16_5 = 0.0;
    constexpr double a16_6 = 22.35793772796803295519317565842520212899;
    constexpr double a16_7 = 0.90949810997556332745009198137971890783;
    constexpr double a16_8 = 35.89110098240264104710550686568482456493;
    constexpr double a16_9 = -3.442515027624453437985000403608480262211;
    constexpr double a16_10 = -4.865481358036368826566013387928704014496;
    constexpr double a16_11 = -18.90980381354342625688427480879773032857;
    constexpr double a16_12 = -34.26354448030451782929251177395134170515;
    constexpr double a16_13 = 1.264756521695642578827783499806516664686;
    constexpr double a16_14 = 0.0;
    constexpr double a16_15 = 0.0;
    constexpr double a17_1 = 0.01461197685842315252051541915018784713459;
    constexpr double a17_2 = 0.0;
    constexpr double a17_3 = 0.0;
    constexpr double a17_4 = 0.0;
    constexpr double a17_5 = 0.0;
    constexpr double a17_6 = 0.0;
    constexpr double a17_7 = 0.0;
    constexpr double a17_8 = -0.3915211862331339089410228267288242030810;
    constexpr double a17_9 = 0.2310932500289506415909675644868993669908;
    constexpr double a17_10 = 0.1274766769992852382560589467488989175618;
    constexpr double a17_11 = 0.2246434176204157731566981937082069688984;
    constexpr double a17_12 = 0.5684352689748512932705226972873692126743;
    constexpr double a17_13 = 0.05825871557215827200814768021863420902155;
    constexpr double a17_14 = 0.1364317403482215641609022744494239843327;
    constexpr double a17_15 = 0.03057013983082797397721005067920369646664;
    constexpr double a17_16 = 0.0;
    constexpr double a18_1 = 0.01549973668189559302279946863304789372788;
    constexpr double a18_2 = 0.0;
    constexpr double a18_3 = 0.0;
    constexpr double a18_4 = 0.0;
    constexpr double a18_5 = 0.0;
    constexpr double a18_6 = 0.0;
    constexpr double a18_7 = 0.0;
    constexpr double a18_8 = 0.3355153219059635054403439303177105512242;
    constexpr double a18_9 = 0.2003613944191860651552622660712101217322;
    constexpr double a18_10 = 0.1252060659283549312946162355194540994211;
    constexpr double a18_11 = 0.2298676393184206750544046308957155868736;
    constexpr double a18_12 = -0.2020250653476181447824906889122391003637;
    constexpr double a18_13 = 0.05917103230665456601422111997583025339897;
    constexpr double a18_14 = -0.02651834783047638681693835956996437528251;
    constexpr double a18_15 = -0.02384094602130971415278110567256446033405;
    constexpr double a18_16 = 0.0;
    constexpr double a18_17 = 0.02718171570208501807097257892166705118335;
    constexpr double a19_1 = 0.01302453943114338366054520296881099431474;
    constexpr double a19_2 = 0.0;
    constexpr double a19_3 = 0.0;
    constexpr double a19_4 = 0.0;
    constexpr double a19_5 = 0.0;
    constexpr double a19_6 = 0.0;
    constexpr double a19_7 = 0.0;
    constexpr double a19_8 = -0.7452850902413112085299330666038981625179;
    constexpr double a19_9 = 0.2643867896429300961465132150322749722129;
    constexpr double a19_10 = 0.1313961758372753932588328082078842388890;
    constexpr double a19_11 = 0.2167253815122927263092467187957410643315;
    constexpr double a19_12 = 0.8734117564076052559016338094938888451419;
    constexpr double a19_13 = 0.01185905643935776688228545787724340848142;
    constexpr double a19_14 = 0.05876002941689550612992712203494447529933;
    constexpr double a19_15 = 0.003266518630202087866399279690939423159022;
    constexpr double a19_16 = 0.0;
    constexpr double a19_17 = -0.008959308648417929824525368306101792182274;
    constexpr double a19_18 = 0.06941415157202692219907482080827253287034;
    constexpr double a20_1 = 0.01397089996925942721283716334050740168797;
    constexpr double a20_2 = 0.0;
    constexpr double a20_3 = 0.0;
    constexpr double a20_4 = 0.0;
    constexpr double a20_5 = 0.0;
    constexpr double a20_6 = 0.0;
    constexpr double a20_7 = 0.0;
    constexpr double a20_8 = -0.4665765335957674596054673402956853940520;
    constexpr double a20_9 = 0.2416372787216257077935214889875485248580;
    constexpr double a20_10 = 0.1290363341345674735721677437066933999929;
    constexpr double a20_11 = 0.2216700671735105311080225734522323922813;
    constexpr double a20_12 = 0.6257275123364644931771253383573999863003;
    constexpr double a20_13 = 0.04355312415679284117869124964829805160429;
    constexpr double a20_14 = 0.1011962491667290833450024852274278874501;
    constexpr double a20_15 = 0.01808582254679721049279369742685497400353;
    constexpr double a20_16 = 0.0;
    constexpr double a20_17 = -0.02079875587689169691156509689282083267654;
    constexpr double a20_18 = -0.09022232517086218976198252891464664868640;
    constexpr double a20_19 = -0.1212796735622254216011467740438097427634;
    constexpr double a21_1 = 0.01604638888318112738641232352800290501904;
    constexpr double a21_2 = 0.0;
    constexpr double a21_3 = 0.0;
    constexpr double a21_4 = 0.0;
    constexpr double a21_5 = 0.0;
    constexpr double a21_6 = 0.0;
    constexpr double a21_7 = 0.0;
    constexpr double a21_8 = 0.09517712399458336651642257453589397190702;
    constexpr double a21_9 = 0.1359187264655317806136927180199100622471;
    constexpr double a21_10 = 0.1237765280959854006935081364365637515893;
    constexpr double a21_11 = 0.2335656264102966047058755123098072346246;
    constexpr double a21_12 = -0.09051508172625873314662090873741762206189;
    constexpr double a21_13 = -0.02537574270006131028513276914038326155331;
    constexpr double a21_14 = -0.1359631696887162048002744757083947500478;
    constexpr double a21_15 = -0.04679214284145113075088049469061349990847;
    constexpr double a21_16 = 0.0;
    constexpr double a21_17 = 0.05177958859391748239949773879090325427473;
    constexpr double a21_18 = 0.09672595677476773313884172931875718705561;
    constexpr double a21_19 = 0.1477312690340742769720989417101989769314;
    constexpr double a21_20 = -0.1150750712958503934434410263732282100773;
    constexpr double a22_1 = 0.01802918623893620731908165792176564180038;
    constexpr double a22_2 = 0.0;
    constexpr double a22_3 = 0.0;
    constexpr double a22_4 = 0.0;
    constexpr double a22_5 = 0.0;
    constexpr double a22_6 = 0.0;
    constexpr double a22_7 = 0.0;
    constexpr double a22_8 = 0.06983601042028873702545973390560096201728;
    constexpr double a22_9 = -0.02541247660791663512384395986842781657182;
    constexpr double a22_10 = 0.008487827035463274491721441398893680307535;
    constexpr double a22_11 = -0.002427525516089801645451101966852425715128;
    constexpr double a22_12 = -0.1047839752893819879012607694745789515746;
    constexpr double a22_13 = -0.01473147795248041942353840372690095884761;
    constexpr double a22_14 = -0.03916338390816177165706892282751065537530;
    constexpr double a22_15 = -0.01005657343293959419073236542225421561652;
    constexpr double a22_16 = 0.0;
    constexpr double a22_17 = 0.01102510392204834322538452331445716455061;
    constexpr double a22_18 = 0.005092830749095398308703438556315975226108;
    constexpr double a22_19 = 0.04759715599420644505591133410826632557391;
    constexpr double a22_20 = 0.03386307003288382751110965442296681690349;
    constexpr double a22_21 = 0.02764422831404797700452373965825845732168;
    constexpr double a23_1 = 0.01677431640522778042988664067637191163626;
    constexpr double a23_2 = 0.0;
    constexpr double a23_3 = 0.0;
    constexpr double a23_4 = 0.0;
    constexpr double a23_5 = 0.0;
    constexpr double a23_6 = 0.0;
    constexpr double a23_7 = 0.0;
    constexpr double a23_8 = 0.6220437408820475326702539861577894278533;
    constexpr double a23_9 = -0.2060859809768841878234097076241307428139;
    constexpr double a23_10 = 0.1156394989766058889629372195583391792474;
    constexpr double a23_11 = 0.02664101793378358946544219293685167025971;
    constexpr double a23_12 = -0.9376810793418770527505892794460093668860;
    constexpr double a23_13 = -0.1367806466702160302637074581619101741312;
    constexpr double a23_14 = -0.3678480995268296672182605288991379118419;
    constexpr double a23_15 = -0.09547871314402478902820445838193201497337;
    constexpr double a23_16 = 0.0;
    constexpr double a23_17 = 0.1013492018422369748729008873270013785313;
    constexpr double a23_18 = -0.08911323084568593396468400926074881389560;
    constexpr double a23_19 = 0.4664140988974760478895528270623735057521;
    constexpr double a23_20 = 0.4502736292354579812232681662308722738519;
    constexpr double a23_21 = 0.1838522463326818655346135218242696774099;
    constexpr double a23_22 = 0.0;
    constexpr double a24_1 = 0.01071149731491444187554380927165768658192;
    constexpr double a24_2 = 0.0;
    constexpr double a24_3 = 0.0;
    constexpr double a24_4 = 0.0;
    constexpr double a24_5 = 0.0;
    constexpr double a24_6 = 0.0;
    constexpr double a24_7 = 0.0;
    constexpr double a24_8 = -0.07094336118221108191937165464264324417735;
    constexpr double a24_9 = 0.1002164900340091596740582334112699697590;
    constexpr double a24_10 = 0.1383453980468025108839271214703390659581;
    constexpr double a24_11 = 0.1796330633578163411338104055485109917477;
    constexpr double a24_12 = 0.09048246545576180974879274948815422276563;
    constexpr double a24_13 = -0.005460662294523338383345981122023862069115;
    constexpr double a24_14 = -0.03000457905119619782973021046143166498567;
    constexpr double a24_15 = -0.01145192026962799093665613252151017277867;
    constexpr double a24_16 = 0.0;
    constexpr double a24_17 = 0.01003394686109385076849515422360600302176;
    constexpr double a24_18 = -0.09506485282809046129031027932806241113157;
    constexpr double a24_19 = 0.04853358804093591445756711642658478691640;
    constexpr double a24_20 = 0.08013325919783924638483373011297347396327;
    constexpr double a24_21 = -0.1251643326835242045676140618774248455713;
    constexpr double a24_22 = 0.0;
    constexpr double a24_23 = 0.0;
    constexpr double a25_1 = 0.01410172088869221367153586187761994182069;
    constexpr double a25_2 = 0.0;
    constexpr double a25_3 = 0.0;
    constexpr double a25_4 = 0.0;
    constexpr double a25_5 = 0.0;
    constexpr double a25_6 = 0.0;
    constexpr double a25_7 = 0.0;
    constexpr double a25_8 = -0.3713379753704491105936205420001801316029;
    constexpr double a25_9 = 0.2231265548117180273161442520179150684520;
    constexpr double a25_10 = 0.1287005345918120122888629169443916280865;
    constexpr double a25_11 = 0.2224600659675494761192249831098918110654;
    constexpr double a25_12 = 0.5382853042550701952740528638168708946100;
    constexpr double a25_13 = 0.05417202616988763101781128062036252796548;
    constexpr double a25_14 = 0.1256968791308743925752109039299467082975;
    constexpr double a25_15 = 0.02784492789002054061504430663197543089132;
    constexpr double a25_16 = 0.0;
    constexpr double a25_17 = -0.03077409246205059733390460511525401688205;
    constexpr double a25_18 = 0.008569805293689777608077303071761466118035;
    constexpr double a25_19 = -0.1535174690587044615794997685221990516897;
    constexpr double a25_20 = -0.02179957030548196497189489878038029238243;
    constexpr double a25_21 = 0.01447128819737186799295514239727801525027;
    constexpr double a25_22 = 0.0;
    constexpr double a25_23 = 0.0;
    constexpr double a25_24 = 0.0;
    constexpr double a26_1 = 0.01424600411735646609296566581447532773183;
    constexpr double a26_2 = 0.0;
    constexpr double a26_3 = 0.0;
    constexpr double a26_4 = 0.0;
    constexpr double a26_5 = 0.0;
    constexpr double a26_6 = 0.0;
    constexpr double a26_7 = 0.0;
    constexpr double a26_8 = -0.3767107393295407091303982522049390741260;
    constexpr double a26_9 = 0.2252399780730421480874737297000189000070;
    constexpr double a26_10 = 0.1283603076292529988314451246143633426068;
    constexpr double a26_11 = 0.2230238705261692544876826347415151339678;
    constexpr double a26_12 = 0.5463127827750747224899202176094949607118;
    constexpr double a26_13 = 0.05526190791375779994553849469706124289752;
    constexpr double a26_14 = 0.1285613508749982456581494397108686240388;
    constexpr double a26_15 = 0.02857250681296406482698934635829147899039;
    constexpr double a26_16 = 0.0;
    constexpr double a26_17 = -0.02398761886357108720930416967644499057175;
    constexpr double a26_18 = 0.05556224458910509454379297181908734648749;
    constexpr double a26_19 = -0.01740675650762838674257930398070760254668;
    constexpr double a26_20 = -0.03815462365996979065575121886854199471011;
    constexpr double a26_21 = 0.01111878504898917877407531966545730451506;
    constexpr double a26_22 = 0.0;
    constexpr double a26_23 = 0.0;
    constexpr double a26_24 = 0.0;
    constexpr double a26_25 = 0.0;

    // interpolation coefficients
    constexpr double p1_1 = 1.;
    constexpr double p1_2 = -12.74966541771576102412726774271024205003;
    constexpr double p1_3 = 68.53080766672322360449931519438662010864;
    constexpr double p1_4 = -194.8119745354184115708662255742713867869;
    constexpr double p1_5 = 317.8426371352858465721400485549409045261;
    constexpr double p1_6 = -299.7155409593396451474764978718840789438;
    constexpr double p1_7 = 152.1119186420412302044289062135608459907;
    constexpr double p1_8 = -32.19357055471805948607776335487247499763;
    constexpr double p2_1 = 0.;
    constexpr double p2_2 = 0.;
    constexpr double p2_3 = 0.;
    constexpr double p2_4 = 0.;
    constexpr double p2_5 = 0.;
    constexpr double p2_6 = 0.;
    constexpr double p2_7 = 0.;
    constexpr double p2_8 = 0.;
    constexpr double p3_1 = 0.;
    constexpr double p3_2 = 0.;
    constexpr double p3_3 = 0.;
    constexpr double p3_4 = 0.;
    constexpr double p3_5 = 0.;
    constexpr double p3_6 = 0.;
    constexpr double p3_7 = 0.;
    constexpr double p3_8 = 0.;
    constexpr double p4_1 = 0.;
    constexpr double p4_2 = 0.;
    constexpr double p4_3 = 0.;
    constexpr double p4_4 = 0.;
    constexpr double p4_5 = 0.;
    constexpr double p4_6 = 0.;
    constexpr double p4_7 = 0.;
    constexpr double p4_8 = 0.;
    constexpr double p5_1 = 0.;
    constexpr double p5_2 = 0.;
    constexpr double p5_3 = 0.;
    constexpr double p5_4 = 0.;
    constexpr double p5_5 = 0.;
    constexpr double p5_6 = 0.;
    constexpr double p5_7 = 0.;
    constexpr double p5_8 = 0.;
    constexpr double p6_1 = 0.;
    constexpr double p6_2 = 0.;
    constexpr double p6_3 = 0.;
    constexpr double p6_4 = 0.;
    constexpr double p6_5 = 0.;
    constexpr double p6_6 = 0.;
    constexpr double p6_7 = 0.;
    constexpr double p6_8 = 0.;
    constexpr double p7_1 = 0.;
    constexpr double p7_2 = 0.;
    constexpr double p7_3 = 0.;
    constexpr double p7_4 = 0.;
    constexpr double p7_5 = 0.;
    constexpr double p7_6 = 0.;
    constexpr double p7_7 = 0.;
    constexpr double p7_8 = 0.;
    constexpr double p8_1 = 0.;
    constexpr double p8_2 = 141.0696092533712605710078918334708622849;
    constexpr double p8_3 = -1283.768450646593611469556198945305487267;
    constexpr double p8_4 = 4630.280061766681019456906394425884505766;
    constexpr double p8_5 = -8648.500976100316056391914391739346002386;
    constexpr double p8_6 = 8890.812161067019054072477301499300549282;
    constexpr double p8_7 = -4787.949212676938594146631544372229445471;
    constexpr double p8_8 = 1057.665286150543793998769524471496193587;
    constexpr double p9_1 = 0.;
    constexpr double p9_2 = -51.75101323451538200013913368910238945675;
    constexpr double p9_3 = 486.0412507312931363219502660627455126442;
    constexpr double p9_4 = -1777.475586368523869187924524440189207548;
    constexpr double p9_5 = 3345.498238649791484065923838383080369426;
    constexpr double p9_6 = -3455.762480007366110892744664331997981283;
    constexpr double p9_7 = 1867.081161290311461863227123612293652933;
    constexpr double p9_8 = -413.4004778109617695287019380323430573493;
    constexpr double p10_1 = 0.;
    constexpr double p10_2 = 16.32082008695896537661978157425676142665;
    constexpr double p10_3 = -118.7072740966746629137122793258036781419;
    constexpr double p10_4 = 379.8980653294655976666790006724113822620;
    constexpr double p10_5 = -659.1980179681815547557635890342282178072;
    constexpr double p10_6 = 645.0127094968870314571906430728856355521;
    constexpr double p10_7 = -335.3923630294779850989153429752631332406;
    constexpr double p10_8 = 72.19353685802189350615784496249014886653;
    constexpr double p11_1 = 0.;
    constexpr double p11_2 = -5.897787927512737700004074981094546266733;
    constexpr double p11_3 = 89.61427156602281202549839530979318529349;
    constexpr double p11_4 = -381.3887877052766946040448187450031021400;
    constexpr double p11_5 = 773.0964199867749581865429091250468401092;
    constexpr double p11_6 = -834.1212536283262120951621382995179788432;
    constexpr double p11_7 = 463.6209151933600203046608830911748528029;
    constexpr double p11_8 = -104.6991340674217303443344573066910439868;
    constexpr double p12_1 = 0.;
    constexpr double p12_2 = -211.5753578293942678386576895806496541353;
    constexpr double p12_3 = 1922.149747207956657441031349510100069933;
    constexpr double p12_4 = -6927.544647063403206822988172864704814066;
    constexpr double p12_5 = 12933.89131149183768334045388734576967421;
    constexpr double p12_6 = -13292.72252236109238830196845558874371599;
    constexpr double p12_7 = 7157.200591588665684047481512110694194478;
    constexpr double p12_8 = -1580.830687765595310572081908235178385210;
    constexpr double p13_1 = 0.;
    constexpr double p13_2 = -29.64315154973386448483853503830022916124;
    constexpr double p13_3 = 265.6161452668802722379486185155415065660;
    constexpr double p13_4 = -951.3146379152784513383019821793424939005;
    constexpr double p13_5 = 1769.876627954108227401999032086988914455;
    constexpr double p13_6 = -1814.926158732241586204906382926227921564;
    constexpr double p13_7 = 975.7252379518513874523768974027019911763;
    constexpr double p13_8 = -215.2758042600138267922695001811431333620;
    constexpr double p14_1 = 0.;
    constexpr double p14_2 = -78.71890822220322199073537821460790272065;
    constexpr double p14_3 = 702.2030963698086240509004430267820179684;
    constexpr double p14_4 = -2509.783618527791392508610051909561977746;
    constexpr double p14_5 = 4663.884687403353676982783161174879377253;
    constexpr double p14_6 = -4779.049763533834849959594190430001691071;
    constexpr double p14_7 = 2567.969360375736223208476377322945879725;
    constexpr double p14_8 = -566.3684221247208382190594586959862794242;
    constexpr double p15_1 = 0.;
    constexpr double p15_2 = -20.19281566680438095870219832836548040387;
    constexpr double p15_3 = 179.3637257949118586258743524266532193814;
    constexpr double p15_4 = -639.8126260732690700781771881696950701907;
    constexpr double p15_5 = 1187.623035597440146280793419574292438396;
    constexpr double p15_6 = -1216.083146478910882267036576696259876826;
    constexpr double p15_7 = 653.1305166034912224190607555907831343139;
    constexpr double p15_8 = -143.9981196370280660478353543467291609746;
    constexpr double p16_1 = 0.;
    constexpr double p16_2 = 0.;
    constexpr double p16_3 = 0.;
    constexpr double p16_4 = 0.;
    constexpr double p16_5 = 0.;
    constexpr double p16_6 = 0.;
    constexpr double p16_7 = 0.;
    constexpr double p16_8 = 0.;
    constexpr double p17_1 = 0.;
    constexpr double p17_2 = 22.25529172937004228561680254541914022993;
    constexpr double p17_3 = -202.2768882656563899475956661732714767089;
    constexpr double p17_4 = 741.4504035959390672891852702617386096470;
    constexpr double p17_5 = -1420.171375819646417700236265030810075574;
    constexpr double p17_6 = 1506.826759283909640808006726351811877831;
    constexpr double p17_7 = -842.0883145405746016477682230643054253648;
    constexpr double p17_8 = 194.0041240166586589127913551094173499391;
    constexpr double p18_1 = 0.;
    constexpr double p18_2 = 14.87763987360506280784259993181818867611;
    constexpr double p18_3 = -312.2840512785346745116969727542541022617;
    constexpr double p18_4 = 1842.506617216153294661664671292748129746;
    constexpr double p18_5 = -4868.899648863078786274115894868122667808;
    constexpr double p18_6 = 6508.346688203267682481628418765155233374;
    constexpr double p18_7 = -4307.866481530869189076138573916573603059;
    constexpr double p18_8 = 1123.319236379456609910815751549228821331;
    constexpr double p19_1 = 0.;
    constexpr double p19_2 = 94.98456650525377646383263451778568838647;
    constexpr double p19_3 = -820.5775654295434684704012573738365180124;
    constexpr double p19_4 = 2818.347256898649110692476937532304801453;
    constexpr double p19_5 = -4992.649725157522569698401937223054737557;
    constexpr double p19_6 = 4835.804647435804986125480300928971220935;
    constexpr double p19_7 = -2434.068718877444022356652060555528975267;
    constexpr double p19_8 = 498.1595386248021872436653821733585200624;
    constexpr double p20_1 = 0.;
    constexpr double p20_2 = 63.46082596146032232234944838983857033697;
    constexpr double p20_3 = -387.6244695782448432215454472613264482576;
    constexpr double p20_4 = 652.6008727039380999132826268449332985693;
    constexpr double p20_5 = 222.0020976000799700370477743563371275858;
    constexpr double p20_6 = -1696.555363345186459306426642315756771818;
    constexpr double p20_7 = 1674.058335196842891022210000150369785868;
    constexpr double p20_8 = -527.9422985388899807669177601643955622843;
    constexpr double p21_1 = 0.;
    constexpr double p21_2 = 57.55994643786018616993511878224123285359;
    constexpr double p21_3 = -588.2803453083489337731949182122044212454;
    constexpr double p21_4 = 2317.048600678134906430718062852747324934;
    constexpr double p21_5 = -4624.295311909926608047251992705773944827;
    constexpr double p21_6 = 5002.133263559409739230532157842265499369;
    constexpr double p21_7 = -2803.532946186995728195816710610623754886;
    constexpr double p21_8 = 639.3667927298664381850782820513480638022;


    template <typename T>
    constexpr std::array<T, 168> _getOrder5Interpolant()
    {
        return { p1_1, p1_2, p1_3, p1_4, p1_5, p1_6, p1_7, p1_8,
                p2_1, p2_2, p2_3, p2_4, p2_5, p2_6, p2_7, p2_8,
                p3_1, p3_2, p3_3, p3_4, p3_5, p3_6, p3_7, p3_8,
                p4_1, p4_2, p4_3, p4_4, p4_5, p4_6, p4_7, p4_8,
                p5_1, p5_2, p5_3, p5_4, p5_5, p5_6, p5_7, p5_8,
                p6_1, p6_2, p6_3, p6_4, p6_5, p6_6, p6_7, p6_8,
                p7_1, p7_2, p7_3, p7_4, p7_5, p7_6, p7_7, p7_8,
                p8_1, p8_2, p8_3, p8_4, p8_5, p8_6, p8_7, p8_8,
                p9_1, p9_2, p9_3, p9_4, p9_5, p9_6, p9_7, p9_8,
                p10_1, p10_2, p10_3, p10_4, p10_5, p10_6, p10_7, p10_8,
                p11_1, p11_2, p11_3, p11_4, p11_5, p11_6, p11_7, p11_8,
                p12_1, p12_2, p12_3, p12_4, p12_5, p12_6, p12_7, p12_8,
                p13_1, p13_2, p13_3, p13_4, p13_5, p13_6, p13_7, p13_8,
                p14_1, p14_2, p14_3, p14_4, p14_5, p14_6, p14_7, p14_8,
                p15_1, p15_2, p15_3, p15_4, p15_5, p15_6, p15_7, p15_8,
                p16_1, p16_2, p16_3, p16_4, p16_5, p16_6, p16_7, p16_8,
                p17_1, p17_2, p17_3, p17_4, p17_5, p17_6, p17_7, p17_8,
                p18_1, p18_2, p18_3, p18_4, p18_5, p18_6, p18_7, p18_8,
                p19_1, p19_2, p19_3, p19_4, p19_5, p19_6, p19_7, p19_8,
                p20_1, p20_2, p20_3, p20_4, p20_5, p20_6, p20_7, p20_8,
                p21_1, p21_2, p21_3, p21_4, p21_5, p21_6, p21_7, p21_8 };
    }

    template <typename T>
    constexpr std::array<T, 420> _getCoeffsForOrder5()
    {
        return { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a2_1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a3_1, a3_2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a4_1, a4_2, a4_3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a5_1, a5_2, a5_3, a5_4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a6_1, a6_2, a6_3, a6_4, a6_5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a7_1, a7_2, a7_3, a7_4, a7_5, a7_6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a8_1, a8_2, a8_3, a8_4, a8_5, a8_6, a8_7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a9_1, a9_2, a9_3, a9_4, a9_5, a9_6, a9_7, a9_8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a10_1, a10_2, a10_3, a10_4, a10_5, a10_6, a10_7, a10_8, a10_9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a11_1, a11_2, a11_3, a11_4, a11_5, a11_6, a11_7, a11_8, a11_9, a11_10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a12_1, a12_2, a12_3, a12_4, a12_5, a12_6, a12_7, a12_8, a12_9, a12_10, a12_11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a13_1, a13_2, a13_3, a13_4, a13_5, a13_6, a13_7, a13_8, a13_9, a13_10, a13_11, a13_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a14_1, a14_2, a14_3, a14_4, a14_5, a14_6, a14_7, a14_8, a14_9, a14_10, a14_11, a14_12, a14_13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a15_1, a15_2, a15_3, a15_4, a15_5, a15_6, a15_7, a15_8, a15_9, a15_10, a15_11, a15_12, a15_13, a15_14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                a16_1, a16_2, a16_3, a16_4, a16_5, a16_6, a16_7, a16_8, a16_9, a16_10, a16_11, a16_12, a16_13, a16_14, a16_15, 0.0, 0.0, 0.0, 0.0, 0.0,
                a17_1, a17_2, a17_3, a17_4, a17_5, a17_6, a17_7, a17_8, a17_9, a17_10, a17_11, a17_12, a17_13, a17_14, a17_15, a17_16, 0.0, 0.0, 0.0, 0.0,
                a18_1, a18_2, a18_3, a18_4, a18_5, a18_6, a18_7, a18_8, a18_9, a18_10, a18_11, a18_12, a18_13, a18_14, a18_15, a18_16, a18_17, 0.0, 0.0, 0.0,
                a19_1, a19_2, a19_3, a19_4, a19_5, a19_6, a19_7, a19_8, a19_9, a19_10, a19_11, a19_12, a19_13, a19_14, a19_15, a19_16, a19_17, a19_18, 0.0, 0.0,
                a20_1, a20_2, a20_3, a20_4, a20_5, a20_6, a20_7, a20_8, a20_9, a20_10, a20_11, a20_12, a20_13, a20_14, a20_15, a20_16, a20_17, a20_18, a20_19, 0.0,
                a21_1, a21_2, a21_3, a21_4, a21_5, a21_6, a21_7, a21_8, a21_9, a21_10, a21_11, a21_12, a21_13, a21_14, a21_15, a21_16, a21_17, a21_18, a21_19, a21_20 };
    }

    template <typename T>
    constexpr std::array<T, 21> _getNodesForOrder5()
    {
        return { c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21 };
    }

    template <typename T>
    constexpr std::array<T, 16> _getWeights()
    {
        return { b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16 };
    }

    template <typename T>
    constexpr std::array<T, 16> _getDeltas()
    {
        return { d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16 };
    }

    constexpr std::array<double, 168> interpolant_f64 = _getOrder5Interpolant<double>();
    constexpr std::array<float, 168> interpolant_f32 = _getOrder5Interpolant<float>();

    constexpr std::array<double, 420> coefficients_f64 = _getCoeffsForOrder5<double>();
    constexpr std::array<float, 420> coefficients_f32 = _getCoeffsForOrder5<float>();

    constexpr std::array<double, 16> weights_f64 = _getWeights<double>();
    constexpr std::array<float, 16> weights_f32 = _getWeights<float>();

    constexpr std::array<double, 16> deltas_f64 = _getDeltas<double>();
    constexpr std::array<float, 16> deltas_f32 = _getDeltas<float>();

    constexpr std::array<double, 21> nodes_f64 = _getNodesForOrder5<double>();
    constexpr std::array<float, 21> nodes_f32 = _getNodesForOrder5<float>();


    template <typename T, typename funcT, const bool scalar = false>
    inline ODEResult<T> integrate(funcT fn, const RP(T) y0, const RP(T) teval, RP(T) Y,
        const T min_step, const T max_step, const T h0, const T rtol, const T atol,
        const size_t n, const size_t m, const RP(T) max_norm, const RP(void) fargs)
    {
        static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value, "Only types float and double are currently supported");

        if constexpr (std::is_same<T, double>::value) {
            return erk::integrate<T, funcT, 16, 21, 8,
                nodes_f64,
                20, coefficients_f64,
                interpolant_f64,
                weights_f64,
                deltas_f64, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
        else if constexpr (std::is_same<T, float>::value) {
            return erk::integrate<T, funcT, 16, 21, 8,
                nodes_f32,
                20, coefficients_f32,
                interpolant_f32,
                weights_f32,
                deltas_f32, scalar>(fn, y0, teval, Y, min_step, max_step, h0, rtol, atol, n, m, max_norm, fargs);
        }
    }

    size_t getWorkArraySize(size_t node)
    {
        return (4 + 21 + 8) * node;
    }

    template <typename funcT, bool scalar = false>
    using stepper64_t = erk::RKStepper<double, funcT, 16, 21, 8,
        nodes_f64, 20, coefficients_f64,
        interpolant_f64,
        weights_f64,
        deltas_f64, scalar>;

    template <typename funcT, bool scalar = false>
    using stepper32_t = erk::RKStepper<float, funcT, 16, 21, 8,
        nodes_f32, 20, coefficients_f32,
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
            using stepper_t = erk::RKStepper<T, funcT, 16, 21, 8,
                nodes_f64, 20, coefficients_f64,
                interpolant_f64,
                weights_f64,
                deltas_f64, scalar>;

            stepper_t stepper;
            stepper.step(f, tn, tf, y0, y, atol, rtol, minStep, maxStep, initStep, node, firstStep, window, work, args);
            return stepper;
        }
        else if constexpr (std::is_same<T, float>::value) {
            using stepper_t = erk::RKStepper<T, funcT, 16, 21, 8,
                nodes_f32, 20, coefficients_f32,
                interpolant_f32,
                weights_f32,
                deltas_f32, scalar>;

            stepper_t stepper;
            stepper.step(f, tn, tf, y0, y, atol, rtol, minStep, maxStep, initStep, node, firstStep, window, work, args);
            return stepper;
        }
    }
}

#endif // ODEPACK_VERNER98E_H