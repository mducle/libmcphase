/*-----------------------------------------------------------------------------*\
| Matpack special functions - BesselExpI0(x)                          dbsi0e.cc |
|                                                                               |
| MatPack Library Release 1.0                                                   |
| Copyright (C) 1991,1995 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "../../include/mpspecfunp.h"

//-----------------------------------------------------------------------------//
//
// double BesselExpI0 (double x);
//
// BesselExpI0(x) calculates the double precision exponentially scaled 
// modified (hyperbolic) Bessel function of the first kind of order 
// zero for double precision argument x.  The result is the Bessel 
// function i0(X) multiplied by exp(-abs(x)). 
//
// This is a translation from the Fortran version of SLATEC, FNLIB,
// CATEGORY C10B1, REVISION 891214, originally written by Fullerton W.,(LANL)
// to C++.
//
// Series for BI0        on the interval  0.          to  9.00000E+00 
//                                        with weighted error   9.51E-34  
//                                         log weighted error  33.02 
//                               significant figures required  33.31 
//                                    decimal places required  33.65 
//
// Series for AI0        on the interval  1.25000E-01 to  3.33333E-01 
//                                        with weighted error   2.74E-32 
//                                         log weighted error  31.56 
//                               significant figures required  30.15 
//                                    decimal places required  32.39 
//
// Series for AI02       on the interval  0.          to  1.25000E-01 
//                                        with weighted error   1.97E-32 
//                                         log weighted error  31.71 
//                               significant figures required  30.15 
//                                    decimal places required  32.63 
//
//-----------------------------------------------------------------------------//

double BesselExpI0 (double x)
{
    static double bi0cs[18] = { 
	-0.07660547252839144951081894976243285,
	 1.927337953993808269952408750881196,
	 0.2282644586920301338937029292330415,
	 0.01304891466707290428079334210691888,
	 4.344270900816487451378682681026107e-4,
	 9.422657686001934663923171744118766e-6,
	 1.434006289510691079962091878179957e-7,
	 1.613849069661749069915419719994611e-9,
	 1.396650044535669699495092708142522e-11,
	 9.579451725505445344627523171893333e-14,
	 5.333981859862502131015107744e-16,
	 2.458716088437470774696785919999999e-18,
	 9.535680890248770026944341333333333e-21,
	 3.154382039721427336789333333333333e-23,
	 9.004564101094637431466666666666666e-26,
	 2.240647369123670016e-28,
	 4.903034603242837333333333333333333e-31,
	 9.508172606122666666666666666666666e-34 
    };

    static double ai0cs[46] = { 
	 0.07575994494023795942729872037438,
	 0.007591380810823345507292978733204,
	 4.153131338923750501863197491382e-4,
	 1.07007646343907307358242970217e-5,
	-7.90117997921289466075031948573e-6,
	-7.826143501438752269788989806909e-7,
	 2.783849942948870806381185389857e-7,
	 8.252472600612027191966829133198e-9,
	-1.204463945520199179054960891103e-8,
	 1.559648598506076443612287527928e-9,
	 2.292556367103316543477254802857e-10,
	-1.191622884279064603677774234478e-10,
	 1.757854916032409830218331247743e-11,
	 1.128224463218900517144411356824e-12,
	-1.146848625927298877729633876982e-12,
	 2.715592054803662872643651921606e-13,
	-2.415874666562687838442475720281e-14,
	-6.084469888255125064606099639224e-15,
	 3.145705077175477293708360267303e-15,
	-7.172212924871187717962175059176e-16,
	 7.874493403454103396083909603327e-17,
	 1.004802753009462402345244571839e-17,
	-7.56689536535053485342843588881e-18,
	 2.150380106876119887812051287845e-18,
	-3.754858341830874429151584452608e-19,
	 2.354065842226992576900757105322e-20,
	 1.11466761204792853022637335511e-20,
	-5.398891884396990378696779322709e-21,
	 1.439598792240752677042858404522e-21,
	-2.591916360111093406460818401962e-22,
	 2.23813318399858390743409229824e-23,
	 5.250672575364771172772216831999e-24,
	-3.249904138533230784173432285866e-24,
	 9.9242141032050379278572847104e-25,
	-2.164992254244669523146554299733e-25,
	 3.233609471943594083973332991999e-26,
	-1.184620207396742489824733866666e-27,
	-1.281671853950498650548338687999e-27,
	 5.827015182279390511605568853333e-28,
	-1.668222326026109719364501503999e-28,
	 3.6253095105415699757006848e-29,
	-5.733627999055713589945958399999e-30,
	 3.736796722063098229642581333333e-31,
	 1.602073983156851963365512533333e-31,
	-8.700424864057229884522495999999e-32,
	 2.741320937937481145603413333333e-32 
    };

    static double ai02cs[69] = {
	 0.0544904110141088316078960962268,
	 0.003369116478255694089897856629799,
	 6.889758346916823984262639143011e-5,
	 2.891370520834756482966924023232e-6,
	 2.048918589469063741827605340931e-7,
	 2.266668990498178064593277431361e-8,
	 3.396232025708386345150843969523e-9,
	 4.940602388224969589104824497835e-10,
	 1.188914710784643834240845251963e-11,
	-3.149916527963241364538648629619e-11,
	-1.321581184044771311875407399267e-11,
	-1.794178531506806117779435740269e-12,
	 7.180124451383666233671064293469e-13,
	 3.852778382742142701140898017776e-13,
	 1.540086217521409826913258233397e-14,
	-4.150569347287222086626899720156e-14,
	-9.554846698828307648702144943125e-15,
	 3.811680669352622420746055355118e-15,
	 1.772560133056526383604932666758e-15,
	-3.425485619677219134619247903282e-16,
	-2.827623980516583484942055937594e-16,
	 3.461222867697461093097062508134e-17,
	 4.465621420296759999010420542843e-17,
	-4.830504485944182071255254037954e-18,
	-7.233180487874753954562272409245e-18,
	 9.92147541217369859888046093981e-19,
	 1.193650890845982085504399499242e-18,
	-2.488709837150807235720544916602e-19,
	-1.938426454160905928984697811326e-19,
	 6.444656697373443868783019493949e-20,
	 2.886051596289224326481713830734e-20,
	-1.601954907174971807061671562007e-20,
	-3.270815010592314720891935674859e-21,
	 3.686932283826409181146007239393e-21,
	 1.268297648030950153013595297109e-23,
	-7.549825019377273907696366644101e-22,
	 1.502133571377835349637127890534e-22,
	 1.265195883509648534932087992483e-22,
	-6.100998370083680708629408916002e-23,
	-1.268809629260128264368720959242e-23,
	 1.661016099890741457840384874905e-23,
	-1.585194335765885579379705048814e-24,
	-3.302645405968217800953817667556e-24,
	 1.313580902839239781740396231174e-24,
	 3.689040246671156793314256372804e-25,
	-4.210141910461689149219782472499e-25,
	 4.79195459108286578063171401373e-26,
	 8.459470390221821795299717074124e-26,
	-4.03980094087283249314607937181e-26,
	-6.434714653650431347301008504695e-27,
	 1.225743398875665990344647369905e-26,
	-2.934391316025708923198798211754e-27,
	-1.961311309194982926203712057289e-27,
	 1.503520374822193424162299003098e-27,
	-9.588720515744826552033863882069e-29,
	-3.483339380817045486394411085114e-28,
	 1.690903610263043673062449607256e-28,
	 1.982866538735603043894001157188e-29,
	-5.317498081491816214575830025284e-29,
	 1.803306629888392946235014503901e-29,
	 6.213093341454893175884053112422e-30,
	-7.69218929277216186320072806673e-30,
	 1.858252826111702542625560165963e-30,
	 1.237585142281395724899271545541e-30,
	-1.102259120409223803217794787792e-30,
	 1.886287118039704490077874479431e-31,
	 2.16019687224365891314903141406e-31,
	-1.605454124919743200584465949655e-31,
	 1.965352984594290603938848073318e-32 
    };

    const double eta  = 0.5 * DBL_EPSILON * 0.1,
	         xsml = sqrt(0.5 * DBL_EPSILON * 4.5);

    double ret_val;

    static int nti0, ntai0, ntai02, first = 1;
    if (first) {
	nti0   = initds(bi0cs,  18, eta);
	ntai0  = initds(ai0cs,  46, eta);
	ntai02 = initds(ai02cs, 69, eta);
	first = 0;
    }

    double y = fabs(x);
    if (y > 3.0) goto L20;

    ret_val = 1.0 - x;
    if (y > xsml)
	ret_val = exp(-y) * (dcsevl(y * y / 4.5 - 1.0 , bi0cs, nti0) + 2.75);
    return ret_val;

  L20:
    if (y <= 8.0) 
	return (dcsevl((48.0 / y - 11.0)/5.0, ai0cs, ntai0) + 0.375) / sqrt(y);
    else // if (y > 8.0)
	return (dcsevl(16.0 / y - 1.0, ai02cs, ntai02) + 0.375) / sqrt(y);
}

//-----------------------------------------------------------------------------//