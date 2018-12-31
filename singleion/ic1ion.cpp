/* ic1ion.cpp
 * 
 * A class for calculating the crystal field Hamiltonian in intermediate coupling
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "ic1ion.hpp"

namespace libMcPhase {

#define Fk2F(F) F[1] *= 225.; F[2] *= 1089.; F[3] *= (184041./25); 

// --------------------------------------------------------------------------------------------------------------- //
// Setter/getter methods for cfpars class
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion::set(const Blm blm, double val) {
    cfpars::set(blm, val);
    m_ham_calc = false;
    m_ev_calc = false;
}

void ic1ion::set(int l, int m, double val) {
    cfpars::set(l, m, val);
    m_ham_calc = false;
    m_ev_calc = false;
}

void ic1ion::set_unit(cfpars::Units const newunit) {
    cfpars::set_unit(newunit);
    m_ham_calc = false;
    m_ev_calc = false;
}

void ic1ion::set_type(const cfpars::Type newtype) {
    cfpars::set_type(newtype);
    m_ham_calc = false;
    m_ev_calc = false;
}

void ic1ion::set_name(const std::string &ionname) {
    cfpars::set_name(ionname);
    getfromionname(ionname);
    m_ham_calc = false;
    m_ev_calc = false;
}

// --------------------------------------------------------------------------------------------------------------- //
// General methods for the ic1ion class
// --------------------------------------------------------------------------------------------------------------- //
RowMatrixXcd ic1ion::hamiltonian(bool upper) {
}

std::tuple<RowMatrixXcd, VectorXd> ic1ion::eigensystem() {
    if (!m_ev_calc) {
        SelfAdjointEigenSolver<RowMatrixXcd> es(hamiltonian(false));
        m_eigenvectors = es.eigenvectors();
        m_eigenvalues = es.eigenvalues();
        m_ev_calc = true;
    }
    return std::tuple<RowMatrixXcd, VectorXd>(m_eigenvectors, m_eigenvalues);
}

// --------------------------------------------------------------------------------------------------------------- //
// Looks up values of Slater and spin orbit radial integrals from published spectroscopic data
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion::getfromionname(const std::string &ionname)
{
    int n,i; 
    orbital l = (orbital)3;    // Defaults for f-electrons
    std::vector<double> F,a;
    double xi = 0.;
    double B=0.,C=0.; bool flg3d = false, flgBC = false; int S2;
    F.assign(4,0.); a.assign(3,0.);
    std::string ion = ionname;
    std::transform(ion.begin(), ion.end(), ion.begin(), [](unsigned char c) { return std::tolower(c); });
    ion.erase(ion.find("+")+1);
#define IONCMP ion.compare
    // Trivalent Lanthanides, from Carnall et al. J. Chem. Phys. v90, pp3443, 1989. All values in cm^{-1}, obtained from RE3+:LaCl_3
    //                               F^2            F^4            F^6        spin-orbit              alpha      beta        gamma
         if(IONCMP("ce3+")==0)                                              { xi = 647.3; n = 1; }
    else if(IONCMP("pr3+")==0) { F[1] = 68878.; F[2] = 50347.; F[3] = 32901.; xi = 751.7; n = 2;  a[0] = 16.23; a[1] = -566.6; a[2] = 1371.; }
    else if(IONCMP("nd3+")==0) { F[1] = 73018.; F[2] = 52789.; F[3] = 35757.; xi = 885.3; n = 3;  a[0] = 21.34; a[1] = -593.0; a[2] = 1445.; }
    else if(IONCMP("pm3+")==0) { F[1] = 76400.; F[2] = 54900.; F[3] = 37700.; xi = 1025.; n = 4;  a[0] = 20.50; a[1] = -560.;  a[2] = 1475.; }
    else if(IONCMP("sm3+")==0) { F[1] = 79805.; F[2] = 57175.; F[3] = 40250.; xi = 1176.; n = 5;  a[0] = 20.16; a[1] = -566.9; a[2] = 1500.; }
    else if(IONCMP("eu3+")==0) { F[1] = 83125.; F[2] = 59268.; F[3] = 42560.; xi = 1338.; n = 6;  a[0] = 20.16; a[1] = -566.9; a[2] = 1500.; }
    else if(IONCMP("gd3+")==0) { F[1] = 85669.; F[2] = 60825.; F[3] = 44776.; xi = 1508.; n = 7;  a[0] = 18.92; a[1] = -600.;  a[2] = 1575.; }
    else if(IONCMP("tb3+")==0) { F[1] = 88995.; F[2] = 62919.; F[3] = 47252.; xi = 1707.; n = 8;  a[0] = 18.40; a[1] = -590.9; a[2] = 1650.; }
    else if(IONCMP("dy3+")==0) { F[1] = 91903.; F[2] = 64372.; F[3] = 49386.; xi = 1913.; n = 9;  a[0] = 18.02; a[1] = -633.4; a[2] = 1790.; }
    else if(IONCMP("ho3+")==0) { F[1] = 94564.; F[2] = 66397.; F[3] = 52022.; xi = 2145.; n = 10; a[0] = 17.15; a[1] = -607.9; a[2] = 1800.; }
    else if(IONCMP("er3+")==0) { F[1] = 97483.; F[2] = 67904.; F[3] = 54010.; xi = 2376.; n = 11; a[0] = 17.79; a[1] = -582.1; a[2] = 1800.; }
    else if(IONCMP("tm3+")==0) { F[1] = 100134.;F[2] = 69613.; F[3] = 55975.; xi = 2636.; n = 12; a[0] = 17.26; a[1] = -624.5; a[2] = 1820.; }
    else if(IONCMP("yb3+")==0)                                              { xi = 2928.; n = 13;}
    // Trivalent Actinides, from Carnall et al. J. Chem. Phys. v90, pp3443, 1989. All values in cm^{-1}, from An3+:LaCl_3 and An3+:LaF_3
    else if(IONCMP("u3+")==0)  { F[1] = 39611.; F[2] = 32960.; F[3] = 23084.; xi = 1626.; n = 3;  a[0] = 29.26; a[1] = -824.6; a[2] = 1820.; }
    else if(IONCMP("np3+")==0) { F[1] = 45382.; F[2] = 37242.; F[3] = 25644.; xi = 1937.; n = 4;  a[0] = 31.78; a[1] = -728.0; a[2] = 1820.; }
    else if(IONCMP("pu3+")==0) { F[1] = 48679.; F[2] = 39333.; F[3] = 27647.; xi = 2242.; n = 5;  a[0] = 30.00; a[1] = -678.3; a[2] = 1820.; }
    else if(IONCMP("am3+")==0) { F[1] = 51900.; F[2] = 41600.; F[3] = 29400.; xi = 2564.; n = 6;  a[0] = 26.71; a[1] = -426.6; a[2] = 1820.; }
    else if(IONCMP("cm3+")==0) { F[1] = 55055.; F[2] = 43938.; F[3] = 32876.; xi = 2889.; n = 7;  a[0] = 29.42; a[1] = -362.9; a[2] = 1820.; }
    else if(IONCMP("bk3+")==0) { F[1] = 57697.; F[2] = 45969.; F[3] = 32876.; xi = 3210.; n = 8;  a[0] = 29.56; a[1] = -564.9; a[2] = 1820.; }
    else if(IONCMP("cf3+")==0) { F[1] = 60464.; F[2] = 48026.; F[3] = 34592.; xi = 3572.; n = 9;  a[0] = 27.36; a[1] = -587.5; a[2] = 1820.; }
    else if(IONCMP("es3+")==0) { F[1] = 63174.; F[2] = 50034.; F[3] = 36199.; xi = 3944.; n = 10; a[0] = 30.21; a[1] = -761.0; a[2] = 1820.; }
    else if(IONCMP("fm3+")==0) { F[1] = 65850.; F[2] = 52044.; F[3] = 37756.; xi = 4326.; n = 11; a[0] = 30.;   a[1] = -600.;  a[2] = 1820.; }
    else if(IONCMP("md3+")==0) { F[1] = 68454.; F[2] = 54048.; F[3] = 39283.; xi = 4715.; n = 12; a[0] = 30.;   a[1] = -600.;  a[2] = 1820.; }
    else if(IONCMP("no3+")==0)                                              { xi = 5144.; n = 13;}
    // Trivalent parameters from Sytsma et al., Phys. Rev. B, v52, pp12668, 1995, in cm^{-1}, obtained from An4+:LuPO_4
 //else if(IONCMP("cm3+")==0) { F[1] = 54669.1;F[2] = 44759.8; F[3]= 33021.4; xi= 2867.7;n = 7;  a[0] = 30.27; a[1] = -981.6; a[2] = 749.3; }
 //else if(IONCMP("gd3+")==0) { F[1] = 84075.; F[2] = 61410.8; F[3]= 44425.9; xi= 1494.; n = 7;  a[0] = 18.92; a[1] = -600.;  a[2] = 1575.; }
    // Tetravalent Actinides, from Conway, J. Chem. Phys. v41, pp904, 1964. All values in cm^{-1}, using hydrogenic wavefunctions.
    //                               F_2            F_4            F_6        spin-orbit
    else if(IONCMP("pa4+")==0)                                              {  xi= 1490.; n = 1; }
 //else if(IONCMP("u4+")==0)  { F[1] = 206.;   F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 1870.; n = 2; Fk2F(F); }
 //else if(IONCMP("np4+")==0) { F[1] = 223.8;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2193.; n = 3; Fk2F(F); }
 //else if(IONCMP("pu4+")==0) { F[1] = 242.9;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2429.; n = 3; Fk2F(F); }
    else if(IONCMP("am4+")==0) { F[1] = 282.1;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 2821.; n = 3; Fk2F(F); }
    else if(IONCMP("cm4+")==0) { F[1] = 307.0;  F[2]=F[1]*F425;F[3]=F[1]*F625; xi= 3042.; n = 3; Fk2F(F); }
    // Tetravalent Actinides parameters from Poirot et al., Phys. Rev. B, v39, pp6388, 1989, in cm^{-1}, obtained from An4+:ZrSiO_4
    //                               F^2            F^4            F^6        spin-orbit              alpha      beta        gamma
    else if(IONCMP("u4+")==0)  { F[1] = 44258.; F[2] = 40293.; F[3] = 31287.; xi = 1740.; n = 2;  a[0] = 23.; }
    else if(IONCMP("np4+")==0) { F[1] = 47949.; F[2] = 41455.; F[3] = 26528.; xi = 2088.; n = 3;  a[0] = 39.2;  a[1] = -610.;  a[2] = 1200.; }
    else if(IONCMP("pu4+")==0) { F[1] = 49394.; F[2] = 39495.; F[3] = 30684.; xi = 2366.; n = 3;  a[0] = 32.3;  a[1] = -783.;  a[2] = 1200.; }
    // U4+:CsCdBr_3 and other parameters, from Karbowiak et al., Chemical Physics, v308, p135, 2005 (Elsevier)
 //else if(IONCMP("u4+")==0)  { F[1] = 45601.; F[2] = 38622.; F[3] = 28423.; xi = 1718.; n = 2;  a[0] = 29.;   a[1] = -1440.; a[2] = 1532.; }
 //else if(IONCMP("u3+")==0)  { F[1] = 36918.; F[2] = 32942.; F[3] = 19906.; xi = 1601.; n = 3;  a[0] = 27.;   a[1] = -830.;  a[2] = 1093.; }
 //else if(IONCMP("cm3+")==0) { F[1] = 53309.; F[2] = 45993.; F[3] = 31047.; xi = 2800.; n = 7;  a[0] = 30.2;  a[1] = -947.;  a[2] = 910.;  }
 //else if(IONCMP("pr3+")==0) { F[1] = 67459.; F[2] = 49029.; F[3] = 32366.; xi = 741.07;n = 2;  a[0] = 22.98; a[1] = -682.98;a[2] = 1422.; }
  
    // All d-electron parameters from Appendix of AS Chakravarty, Introduction to Magnetic Properties of Solid, Wiley, 1980. Original work also cited
    // 3d ions parameters from JS Grifiths, The Theory of Transition Metal Ions, CUP, 1961 (A and B); lambda from TM Dunn, Trans. Faraday Soc. v57, 1441 (1961)
    //   where 0 is shown, parameter not known... (where 0. is shown parameter is zero!) where theory_lamdba shows "-" parameter is actually theoretical...
    //                                      lambda_experimental lambda_theory-[from M Blume and RE Watson, Proc. R. Soc. Lon. A v271, 565 (1963)]
    else if(IONCMP("sc2+")==0) { B = 0.;    C = 0.;     xi = 79.;  /*xi=86.;*/ n = 1; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("ti")==0)   { B = 560.;  C = 1840.;  xi = 0;    /*xi=   ;*/ n = 4; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("ti+")==0)  { B = 682.;  C = 2481.;  xi = 0;    /*xi=   ;*/ n = 3; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("ti2+")==0) { B = 718.;  C = 2629.;  xi = 60.;  /*xi=61.;*/ n = 2; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("ti3+")==0) { B = 0.;    C = 0.;     xi = 154.; /*xi=159;*/ n = 1; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("v")==0)    { B = 578.;  C = 2273.;  xi = 0;    /*xi=   ;*/ n = 5; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("v+")==0)   { B = 659.;  C = 2417.;  xi = 0;    /*xi=   ;*/ n = 4; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("v2+")==0)  { B = 766.;  C = 2855.;  xi = 55.;  /*xi=57.;*/ n = 3; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("v3+")==0)  { B = 861.;  C = 4165.;  xi = 106.; /*xi=104;*/ n = 2; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("v4+")==0)  { B = 0.;    C = 0.;     xi = 248.; /*xi=255;*/ n = 1; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("cr")==0)   { B = 790.;  C = 2520.;  xi = 0;    /*xi=   ;*/ n = 6; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("cr+")==0)  { B = 710.;  C = 2790.;  xi = 0;    /*xi=   ;*/ n = 5; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("cr2+")==0) { B = 830.;  C = 3430.;  xi = 58.;  /*xi=59.;*/ n = 4; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("cr3+")==0) { B = 1030.; C = 3850.;  xi = 91.;  /*xi=91.;*/ n = 3; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("cr4+")==0) { B = 1039.; C = 4238.;  xi = 164.; /*xi=163;*/ n = 2; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("mn")==0)   { B = 720.;  C = 3087.;  xi = 0;    /*xi=   ;*/ n = 7; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("mn+")==0)  { B = 873.;  C = 3130.;  xi = 64.;  /*xi=64.;*/ n = 6; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("mn2+")==0) { B = 960.;  C = 3325.;  xi = 68.57;/*xi=68.57*/n = 5; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("mn3+")==0) { B = 1140.; C = 3675.;  xi = 88.;  /*xi=87.;*/ n = 4; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("mn4+")==0) { F[1]=87044;F[2]=54316; xi = 134.; /*xi=135;*/ n = 3; l=D; flg3d=1;          } // (Expt.) Uylings et al., J. Phys. B. 17 (1984) 4103
    else if(IONCMP("fe")==0)   { B = 806.;  C = 3506.;  xi = 0;    /*xi=   ;*/ n = 8; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("fe+")==0)  { B = 869.;  C = 3638.;  xi = 119.; /*xi=115;*/ n = 7; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("fe2+")==0) { B = 1058.; C = 3091.;  xi = 103.; /*xi=114;*/ n = 6; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("fe3+")==0) { F[1]=97130;F[2]=60769; xi = 476.; /*xi=   ;*/ n = 5; l=D; flg3d=1;          } // Havercort Thesis 2p6 3d5
    else if(IONCMP("fe4+")==0) { B = 1144.; C = 4459.;  xi = 129.; /*xi=125;*/ n = 4; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("co")==0)   { B = 798.;  C = 4167.;  xi = 0;    /*xi=   ;*/ n = 9; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("co+")==0)  { B = 878.;  C = 3828.;  xi = 228.; /*xi=228;*/ n = 8; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("co2+")==0) { B = 1115.; C = 4366.;  xi = 178.; /*xi=189;*/ n = 7; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("co3+")==0) { B = 1065.; C = 5120.;  xi = 128.6;/*xi=145;*/ n = 6; l=D; flg3d=1; flgBC=1; } // Abragam Bleaney 1970 p 391 for B,C. xi from PRB 67 172401
// else if(IONCMP("ni")==0)   { B = 1025.; C = 4226.;  xi = 0;    /*xi=   ;*/ n =10; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("ni+")==0)  { B = 1037.; C = 4314.;  xi = 0;    /*xi=   ;*/ n = 9; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("ni2+")==0) { B = 1084.; C = 4831.;  xi = 324.; /*xi=343;*/ n = 8; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("ni3+")==0) { B = 1184.; C = 5105.;  xi = 272.; /*xi= - ;*/ n = 7; l=D; flg3d=1; flgBC=1; } // B,C, from fit to NIST data using Racah II, eqn 84
    else if(IONCMP("ni4+")==0) { F[1]=100185;F[2]=64787;xi = 197.; /*xi= - ;*/ n = 6; l=D; flg3d=1;          } // (Expt.) Uylings et al., J. Phys. B. 17 (1984) 4103
// else if(IONCMP("cu+")==0)  { B = 1216.; C = 4745.;  xi = 0;    /*xi=   ;*/ n =10; l=D; flg3d=1; flgBC=1; } //
    else if(IONCMP("cu2+")==0) { B = 1238.; C = 4659.;  xi = 830.; /*xi=830;*/ n = 9; l=D; flg3d=1; flgBC=1; }
    else if(IONCMP("cu3+")==0) { F[1]=111996;F[2]=69924;xi = 903.; /*xi= - ;*/ n = 8; l=D; flg3d=1;          } // Thesis Havercort Koeln Cu 2p6 3d8
    else if(IONCMP("zn3+")==0) { F[1]=116868;F[2]=72923;xi = 1097.;/*xi=   ;*/ n = 9; l=D; flg3d=1;          } // Havercort Thesis 2p6 3d9
     
    // 4d ions parameters from Richardson, Blackman and Ranschak, J. Chem. Phys. v58, 3010 (1973).
    //   xi from calculations of Blume, Freeman, Watson, Phys. Rev. v134, A320 (1964), or where not calculated from TM Dunn, Trans. Faraday Soc. v57, 1441 (1961)
    else if(IONCMP("y2+")==0)  { B = 0.;    C = 0.;     xi = 312.; /*xi=300;*/ n = 1; l=D; flgBC=1; }
    else if(IONCMP("zr2+")==0) { B = 333.;  C = 3.96*B; xi = 432.; /*xi=425;*/ n = 2; l=D; flgBC=1; }
    else if(IONCMP("zr3+")==0) { B = 0.;    C = 0.;     xi = 507.; /*xi=500;*/ n = 1; l=D; flgBC=1; }
    else if(IONCMP("nb+")==0)  { B = 324.;  C = 3.89*B; xi = 0;    /*xi=   ;*/ n = 4; l=D; flgBC=1; }
    else if(IONCMP("nb2+")==0) { B = 360.;  C = 3.97*B; xi = 560.; /*xi=555;*/ n = 3; l=D; flgBC=1; }
    else if(IONCMP("nb3+")==0) { B = 391.;  C = 4.03*B; xi = 644.; /*xi=670;*/ n = 2; l=D; flgBC=1; }
    else if(IONCMP("nb4+")==0) { B = 0.;    C = 0.;     xi = 750.; /*xi=   ;*/ n = 1; l=D; flgBC=1; }
    else if(IONCMP("mo+")==0)  { B = 355.;  C = 3.91*B; xi = 630.; /*xi=   ;*/ n = 5; l=D; flgBC=1; }
    else if(IONCMP("mo2+")==0) { B = 387.;  C = 3.98*B; xi = 717.; /*xi=695;*/ n = 4; l=D; flgBC=1; }
    else if(IONCMP("mo3+")==0) { B = 416.;  C = 4.03*B; xi = 812.; /*xi=800;*/ n = 3; l=D; flgBC=1; }
    else if(IONCMP("mo4+")==0) { B = 440.;  C = 4.08*B; xi = 950.; /*xi=   ;*/ n = 2; l=D; flgBC=1; }
    else if(IONCMP("mo5+")==0) { B = 0.;    C = 0.;     xi =1030.; /*xi=   ;*/ n = 1; l=D; flgBC=1; }
    else if(IONCMP("tc2+")==0) { B = 414.;  C = 3.99*B; xi = 850.; /*xi=   ;*/ n = 5; l=D; flgBC=1; }
    else if(IONCMP("tc3+")==0) { B = 440.;  C = 4.03*B; xi = 990.; /*xi=   ;*/ n = 4; l=D; flgBC=1; }
    else if(IONCMP("tc4+")==0) { B = 0;     C = 0;      xi =1150.; /*xi=   ;*/ n = 3; l=D; flgBC=1; } //
    else if(IONCMP("ru2+")==0) { B = 436.;  C = 3.99*B; xi =1077.; /*xi=1000*/ n = 6; l=D; flgBC=1; }
    else if(IONCMP("ru3+")==0) { B = 464.;  C = 4.04*B; xi =1197.; /*xi=1180*/ n = 5; l=D; flgBC=1; }
    else if(IONCMP("ru4+")==0) { B = 400;   C = 1270;   xi =1350.; /*xi=   ;*/ n = 4; l=D; flgBC=1; } //
    else if(IONCMP("rh+")==0)  { B = 427.;  C = 3.93*B; xi =0;     /*xi=   ;*/ n = 8; l=D; flgBC=1; } //
    else if(IONCMP("rh2+")==0) { B = 458.;  C = 3.98*B; xi =1664.; /*xi=1640*/ n = 7; l=D; flgBC=1; }
    else if(IONCMP("rh3+")==0) { B = 484.;  C = 4.03*B; xi =1416.; /*xi=1400*/ n = 6; l=D; flgBC=1; }
    else if(IONCMP("rh4+")==0) { B = 0;     C = 0;      xi =1570.; /*xi=   ;*/ n = 5; l=D; flgBC=1; } //
    else if(IONCMP("pd+")==0)  { B = 451.;  C = 3.94*B; xi =0;     /*xi=   ;*/ n = 9; l=D; flgBC=1; } //
    else if(IONCMP("pd2+")==0) { B = 480.;  C = 3.99*B; xi =1529.; /*xi=1600*/ n = 8; l=D; flgBC=1; }
    else if(IONCMP("pd3+")==0) { B = 506.;  C = 4.03*B; xi =1529.; /*xi=1600*/ n = 7; l=D; flgBC=1; }
    else if(IONCMP("ag2+")==0) { B = 502.;  C = 3.99*B; xi =1794.; /*xi=1840*/ n = 9; l=D; flgBC=1; }
    else if(IONCMP("ag3+")==0) { B = 528.;  C = 4.03*B; xi =1940.; /*xi=1930*/ n = 8; l=D; flgBC=1; }
    // 5d ions parameters from G Burns, J. Chem. Phys. v41, 1521 (1964) B,C only.
    else if(IONCMP("re2+")==0) { F[1]=30870;F[2]=22050; xi =0.;    /*xi=   ;*/ n = 5; l=D; }
    
    else { std::cerr << "getfromionname(): Name of ion " << ionname << " not recognised.\n"; return; }

    if(flg3d) // 3d ion - "xi" is actually lambda parameter -> |lambda| = xi/2S
    {
        S2 = (n<=(2*l+1)) ? n : ((4*l+2)-n);    // Finds 2S (maximise S from Hund's Rules)
        xi *= (double)S2;
    }
    if(flgBC) // Need to convert B and C parameters to F^2, F^4 slater integrals
    {
        // Eqn 77 of Racah II, A = F_0-49F_4 = F^0-F^4/9; B = F_2-5F_4 = (9F^2-5F^4)/441; C = 35F_4 = 5F^4/63;
        F[2] = (63./5)*C; F[1] = (441.*B+5*F[2])/9.;
    }

    m_n = n; 
    m_l = l;
    m_F = F;
    m_xi = xi;
    m_alpha = a;
    calc_stevfact();
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Stevens operator equivalent factors after Elliot, Judd and Runciman.
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion::calc_stevfact()
{
    if(m_l != F)  // <L||theta||L> from Abragam and Bleaney, EPR, Appendix B, Table 19, page 873
    {
        for(int i=0; i<3; i++) { m_stevfact[i]=0; }
        if(m_n==5) return;
        int S2 = (m_n<=5) ? m_n : (10-m_n);            // Maximise S.  (Hund's Rules with l=2 substituted in)
        int L = 0, ln = (m_n<=5) ? (2-m_n) : (7-m_n); 
        for(int ml=2; ml>ln; ml--) L+=ml;              // Maximise L.
        m_stevfact[0] = 2*(5-2*S2)/21./(2.*L-1);       // <L||alpha||L> = -/+ 2(2l+1-4S)/(2l-1)/(2l+3)/(2L-1)
        if(n<5) m_stevfact[0] = -m_stevfact[0];
        // <L||beta||L> = <L||alpha||L> * ( 3*(3(l-1)(l+2)-7(l-2S)(l+1-2S)) / 2(2l-3)(2l+5)(L-1)(2L-3) )
        m_stevfact[1] = m_stevfact[0]*3*(12-7*(2-S2)*(3-S2))/18./(L-1.)/(2*L-3.);
    }
    else {
        if(m_n==1) { 
            m_stevfact = {-2/35., 2/315., 0.};
        }
    
        m_stevfact = {-8.*sqrt(7./15.), 16.*sqrt(14./11.), -640.*sqrt(7./429.)};
    
        orbital L;
        int S2, J2;
        int ml, Li = 0, lm = abs(l), ln = (n<=(2*l+1)) ? (abs(l)-n) : (abs(l)-n+(2*l+1)); 
        fconf conf(n, l);
        fconf confp(n-1, l);
        std::vector<cfpls> cfps;
        double sumcfp, noncfpprod;
        int k, is, ic;
    
        // Calculates L,S,J from Hund's rules
        S2 = (n<=(2*l+1)) ? n : ((4*l+2)-n);                       // Maximise S.
        for(ml=lm; ml>ln; ml--) Li+=ml;                            // Maximise L.
        J2 = (n<=(2*l+1)) ? abs(2*Li-S2) : (2*Li+S2);              // J=|L-S| for <half-full, J=L+S for >half-full shell.
        L = (orbital)Li;
    
        // Finds the seniority and other quantum numbers from the list of states
        for(is=0; is<(int)conf.states.size(); is++) { if(conf.states[is].S2==S2 && conf.states[is].L==L) break; }
    
        noncfpprod = pow(-1.,-(double)l-Li) * (2.*Li+1.);
        if(l==D) cfps = racah_parents(n,conf.states[is].v,S2,L); else cfps = racah_parents(n,conf.states[is].v,conf.states[is].U,S2,L);
    
        // We use the formulae of Elliot, Judd and Runciman, Proc. R. Soc. Lon. A, v240, pp509, 1957 to calculate theta_k
        for(k=0; k<3; k++)
        {
            sumcfp = 0.;
            for(ic=0; ic<(int)cfps.size(); ic++) 
                sumcfp += m_racah.racahW(l*2,Li*2,l*2,Li*2,abs(confp.states[cfps[ic].ind].L)*2,(k+1)*4) * cfps[ic].cfp*cfps[ic].cfp 
    	            * pow(-1.,(double)abs(confp.states[cfps[ic].ind].L)+2.*(k+1.)) * noncfpprod;
            m_stevfact[k] *= sqrt(factorial(J2-2*(k+1))/factorial(J2+2*(k+1)+1)) * n * sumcfp * pow(-1.,2.*(k+1.)+Li+(S2+J2)/2.) 
                          * (J2+1) * sixj(Li*2,J2,S2,J2,Li*2,4*(k+1));
        }
    }
}

} // namespace libMcPhase
