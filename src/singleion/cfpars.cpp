/* cfpars.cpp
 * 
 * A class encapsulating crystal field parameters and their conversions
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cfpars.hpp"

namespace libMcPhase {

// --------------------------------------------------------------------------------------------------------------- //
// Reference tables (values taken from program cfield, by Peter Fabi, FZ Juelich, file theta.c)
// --------------------------------------------------------------------------------------------------------------- //

// Conversion factors from Stevens to Wybourn normalisation
static const std::array<double, 27> lambda = {
    {sqrt(6.)/2., sqrt(6.), 1./2., sqrt(6.), sqrt(6.)/2., // l=2
     sqrt(70.)/8., sqrt(35.)/2., sqrt(10.)/4., sqrt(5.)/2., 1./8., sqrt(5.)/2., sqrt(10.)/4., sqrt(35.)/2., sqrt(70.)/8., // l=4
     sqrt(231.)/16., 3*sqrt(77.)/8., 3*sqrt(14.)/16., sqrt(105.)/8., sqrt(105.)/16., sqrt(42.)/8., // l=6, m=-6 to m=-1
     1./16., sqrt(42.)/8., sqrt(105.)/16., sqrt(105.)/8., 3*sqrt(14.)/16., 3*sqrt(77.)/8., sqrt(231.)/16.} };

static const std::array<int, 27> id2l = { {0,0,0,0,0, 1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2,2,2,2,2,2,2} };

// Matching ion names to number of unfilled electrons in f-shell
static const std::unordered_map<std::string, int> ION_NUMBER_MAP = {
    {"ce3+",  1}, {"pr3+",  2}, {"nd3+",  3}, {"pm3+",  4}, {"sm3+",  5}, {"eu3+",  6}, {"gd3+",  7},
    {"tb3+",  8}, {"dy3+",  9}, {"ho3+", 10}, {"er3+", 11}, {"tm3+", 12}, {"yb3+", 13}, {"nd2+",  4},
    {"sm2+",  6}, {"eu2+",  7}, {"gd2+",  8}, {"tb2+",  9}, {"dy2+", 10}, {"ho2+", 11}, {"er2+", 12},
    {"tm2+", 13}, {"u4+",   2}, {"u3+",   3}, {"u2+",   4}, {"np4+",  3}, {"np3+",  4}, {"pu4+",  4},
    {"pu3+",  5} };

// Matching f-number to value of twice the total angular momentum J of the ground multiplet
static const std::array<int, 15> J2 = { { 0, 5, 8, 9, 8, 5, 0, 7, 12, 15, 16, 15, 12, 7, 0 } };

// Lande g-factor
static const std::array<double, 15> GJ = { { 0., 6./7, 0.8, 8./11, 0.6, 2./7, 0., 2., 1.5, 4./3, 1.25, 1.2, 7./6, 8./7 } };

static const std::array<double, 15> ALPHA_J = {{ // 2nd order Stevens Operator equivalent factor
    /* 4f_0        */    1.0 * 0                ,
    /* 4f_1 : Ce3+ */   -1.0 * 2/5/7            ,
    /* 4f_2 : Pr3+ */   -1.0 * 2*2*13/3/3/5/5/11,
    /* 4f_3 : Nd3+ */   -1.0 * 7/3/3/11/11      ,
    /* 4f_4 : Pm3+ */    1.0 * 2*7/3/5/11/11    ,
    /* 4f_5 : Sm3+ */    1.0 * 13/3/3/5/7       ,
    /* 4f_6 : Eu3+ */    1.0 * 0                ,
    /* 4f_7 : Gd3+ */    1.0 * 0                ,
    /* 4f_8 : Tb3+ */   -1.0 * 1/3/3/11         ,
    /* 4f_9 : Dy3+ */   -1.0 * 2/3/3/5/7        ,
    /* 4f_10: Ho3+ */   -1.0 * 1/2/3/3/5/5      ,
    /* 4f_11: Er3+ */    1.0 * 2*2/3/3/5/5/7    ,
    /* 4f_12: Tm3+ */    1.0 * 1/3/3/11         ,
    /* 4f_13: Yb3+ */    1.0 * 2/3/3/7          ,
    /* 4f_14       */    1.0 * 0
}};

static const std::array<double, 15> BETA_J = {{ // 4th order Stevens Operator equivalent factor
    /* 4f_0        */    1.0 * 0                              ,
    /* 4f_1 : Ce3+ */    1.0 * 2/3/3/5/7                      ,
    /* 4f_2 : Pr3+ */   -1.0 * 2*2/3/3/5/11/11                ,
    /* 4f_3 : Nd3+ */   -1.0 * 2*2*2*17/3/3/3/11/11/11/13     ,
    /* 4f_4 : Pm3+ */    1.0 * 2*2*2*7*17/3/3/3/5/11/11/11/13 ,
    /* 4f_5 : Sm3+ */    1.0 * 2*13/3/3/3/5/7/11              ,
    /* 4f_6 : Eu3+ */    1.0 * 0                              ,
    /* 4f_7 : Gd3+ */    1.0 * 0                              ,
    /* 4f_8 : Tb3+ */    1.0 * 2/3/3/3/5/11/11                ,
    /* 4f_9 : Dy3+ */   -1.0 * 2*2*2/3/3/3/5/7/11/13          ,
    /* 4f_10: Ho3+ */   -1.0 * 1/2/3/5/7/11/13                ,
    /* 4f_11: Er3+ */    1.0 * 2/3/3/5/7/11/13                ,
    /* 4f_12: Tm3+ */    1.0 * 2*2*2/3/3/3/3/5/11/11          ,
    /* 4f_13: Yb3+ */   -1.0 * 2/3/5/7/11                     ,
    /* 4f_14       */    1.0 * 0
}};

static const std::array<double, 15> GAMMA_J = {{ // 6th order Stevens Operator equivalent factor
    /* 4f_0        */    1.0 * 0                                 ,
    /* 4f_1 : Ce3+ */    1.0 * 0                                 ,
    /* 4f_2 : Pr3+ */    1.0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13   ,
    /* 4f_3 : Nd3+ */   -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13    ,
    /* 4f_4 : Pm3+ */    1.0 * 2*2*2*17*19/3/3/3/7/11/11/11/13/13,
    /* 4f_5 : Sm3+ */    1.0 * 0                                 ,
    /* 4f_6 : Eu3+ */    1.0 * 0                                 ,
    /* 4f_7 : Gd3+ */    1.0 * 0                                 ,
    /* 4f_8 : Tb3+ */   -1.0 * 1/3/3/3/3/7/11/11/13              ,
    /* 4f_9 : Dy3+ */    1.0 * 2*2/3/3/3/7/11/11/13/13           ,
    /* 4f_10: Ho3+ */   -1.0 * 5/3/3/3/7/11/11/13/13             ,
    /* 4f_11: Er3+ */    1.0 * 2*2*2/3/3/3/7/11/11/13/13         ,
    /* 4f_12: Tm3+ */   -1.0 * 5/3/3/3/3/7/11/11/13              ,
    /* 4f_13: Yb3+ */    1.0 * 2*2/3/3/3/7/11/13                 ,
    /* 4f_14       */    1.0 * 0
}};

using Map3 = std::unordered_map<std::string, std::array<double, 3>>;

// Expectation value of radial wavefunction <r^k>, from theta.c in McPhase cf1ion_module source code.
void createRkTable(Map3 &rk) {
    rk["ce3+"] = {1.309, 3.964, 23.31}; 
    rk["pr3+"] = {1.1963, 3.3335, 18.353}; /* U. Walter Diss.         */
    rk["nd3+"] = {1.114, 2.910, 15.03}; 
    rk["pm3+"] = {1.0353, 2.5390, 12.546}; /*          -"-            */
    rk["sm3+"] = {0.9743, 2.260, 10.55}; 
    rk["eu3+"] = {0.9175, 2.020, 9.039}; 
    rk["gd3+"] = {0.8671, 1.820, 7.831}; 
    rk["tb3+"] = {0.8220, 1.651, 6.852}; 
    rk["dy3+"] = {0.7814, 1.505, 6.048}; 
    rk["ho3+"] = {0.7446, 1.379, 5.379}; 
    rk["er3+"] = {0.7111, 1.270, 4.816}; 
    rk["tm3+"] = {0.6804, 1.174, 4.340}; 
    rk["yb3+"] = {0.6522, 1.089, 3.932}; 
    rk["nd2+"] = {1.392, 5.344, 45.450}; 
    rk["sm2+"] = {1.197, 3.861, 28.560}; 
    rk["eu2+"] = {1.098, 3.368, 23.580}; 
    rk["gd2+"] = {1.028, 2.975, 19.850}; 
    rk["tb2+"] = {0.968, 2.655, 16.980}; 
    rk["dy2+"] = {0.913, 2.391, 14.730}; 
    rk["ho2+"] = {0.866, 2.169, 12.920}; 
    rk["er2+"] = {0.824, 1.979, 11.450}; 
    rk["tm2+"] = {0.785, 1.819, 10.240}; 
    rk["u4+"] = {2.042, 7.632, 47.774};    /* Freeman et al. PRB 13 (1976) 1168 */
    rk["u3+"] = {2.346, 10.906, 90.544};   /* Freeman et al. PRB 13 (1976) 1168 */
    rk["u2+"] = {3.257, 26.82, 462.85};    /* Lewis et al. J. Chem Phys. 53 (1970) 809 */
    rk["np4+"] = {1.884, 6.504, 37.80};    /* Lewis et al. J. Chem Phys. 53 (1970) 809 */
    rk["np3+"] = {2.297, 11.00, 98.63};    /* Lewis et al. J. Chem Phys. 53 (1970) 809 */
    rk["pu4+"] = {1.838, 6.401, 38.77};    /* Lewis et al. J. Chem Phys. 53 (1970) 809 */
    rk["pu3+"] = {2.1025, 9.1775, 73.3};   /* Lewis et al. J. Chem Phys. 53 (1970) 809 */
}

const Map3 &RKTABLE() {
    static Map3 rk_table;
    if (rk_table.empty()) 
        createRkTable(rk_table);
    return rk_table;
}

// --------------------------------------------------------------------------------------------------------------- //
// General methods for cfpars class
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::getfromionname(const std::string &ionname) {
	std::string ion = ionname;
	std::transform(ion.begin(), ion.end(), ion.begin(), [](unsigned char c) { return std::tolower(c); });
    auto ion_number = ION_NUMBER_MAP.find(ion);
    if (ion_number == ION_NUMBER_MAP.end()) {
        throw std::runtime_error("Unknown ion");
    }
    m_ionname = ion;
    m_n = ion_number->second;
	m_J2 = J2[m_n];
    const Map3 &rktable = RKTABLE();
	auto rk_table = rktable.find(ion);
	m_rk = rk_table->second;
}

// --------------------------------------------------------------------------------------------------------------- //
// Setter/getter methods for cfpars class
// --------------------------------------------------------------------------------------------------------------- //
void cfpars::set(const Blm blm, double val) {
    int id = (int)blm;
    m_Bo[id] = val;
    m_Bi[id] = val / m_convfact[id] / m_econv;
}

void cfpars::set(int l, int m, double val) {
    int id;
    switch(l) {
        case 2: id = 2 + m; break;
        case 4: id = 9 + m; break;
        case 6: id = 20 + m; break;
        default:
            return;
    }
    m_Bo[id] = val;
    m_Bi[id] = val / m_convfact[id] / m_econv;
}

double cfpars::get(int l, int m) const {
    switch(l) {
        case 2: return m_Bo[2 + m];
        case 4: return m_Bo[9 + m];
        case 6: return m_Bo[20 + m];
        default:
            return 0.0;
    }
}

void cfpars::set_unit(cfpars::Units const newunit) {
    if (m_unit == newunit)
        return;
    m_unit = newunit;
    m_econv = ENERGYCONV[(int)m_unit];
    for (int id=0; id<27; id++) {
        m_Bo[id] = m_Bi[id] * m_econv * m_convfact[id];
    }
}

void cfpars::set_type(const cfpars::Type newtype) {
    if (!m_convertible) {
        throw std::runtime_error("Unknown ion, cannot set parameter type. Please set the ionname on construction.");
    }
    switch(newtype) {
        case cfpars::Type::Alm:
            for (int id=0; id<27; id++)
                if(m_stevfact[id2l[id]] != 0)
                    m_convfact[id] = lambda[id] / m_stevfact[id2l[id]] / m_rk[id2l[id]];
            break;
        case cfpars::Type::ARlm:
            for (int id=0; id<27; id++) 
                if(m_stevfact[id2l[id]] != 0)
                    m_convfact[id] = lambda[id] / m_stevfact[id2l[id]];
            break;
        case cfpars::Type::Blm:
            for (int id=0; id<27; id++) 
                    m_convfact[id] = lambda[id];
            break;
        case cfpars::Type::Llm:
            for (int id=0; id<27; id++)
                if(m_stevfact[id2l[id]] != 0)
                    m_convfact[id] = 1. / m_stevfact[id2l[id]];
            break;
    }
    for (int id=0; id<27; id++)
        m_Bo[id] = m_Bi[id] * m_convfact[id] * m_econv;
    m_type = newtype;
    if(newtype == cfpars::Type::Llm)
        m_norm = cfpars::Normalisation::Wybourne;
    else
        m_norm = cfpars::Normalisation::Stevens;
}

void cfpars::set_name(const std::string &ionname) {
    getfromionname(ionname);
    double alpha = ALPHA_J[m_n];
    double beta = BETA_J[m_n];
    double gamma = GAMMA_J[m_n];
    // Scales parameters by ratio of Stevens factor of old and new ion (if old ion defined).
    if (m_convertible) {
        std::array<double, 3> scale = { {alpha/m_stevfact[0], beta/m_stevfact[1], gamma/m_stevfact[2]} };
        for (int id=0; id<27; id++)
            m_Bi[id] *= scale[id2l[id]];
    } 
    else {
        for (int id=0; id<27; id++)
            m_Bi[id] = m_Bo[id] / m_econv;
    }
	m_stevfact = {alpha, beta, gamma};
	m_convertible = true;
    m_GJ = GJ[m_n];
    // Now reset the conversion table (from internal to external parameters)
    this->set_type(m_type);
}

void cfpars::set_J(const double J) {
    if (fmod(2 * J, 1.) == 0.) {
        m_J2 = (int)(2 * J);
    }
    else {
        throw std::runtime_error("Invalid value of J - must be integer or half-integer");
    }
    m_ionname = "";
    m_stevfact = {1., 1., 1.};
    m_rk = {0., 0., 0.};
    m_convertible = false;
    m_convfact = lambda;
    for (int id=0; id<27; id++) {
        m_Bi[id] = m_Bo[id] / m_econv;
    }
    m_norm = cfpars::Normalisation::Stevens;
    m_type = Type::Blm;
}

// --------------------------------------------------------------------------------------------------------------- //
// Constructor functions for cfpars class
// --------------------------------------------------------------------------------------------------------------- //
cfpars::cfpars() {
    m_convfact = lambda;
}

cfpars::cfpars(int J2) {
    m_J2 = J2;
    m_convfact = lambda;
}

cfpars::cfpars(const double J) {
    if (fmod(2 * J, 1.) == 0.) {
        m_J2 = (int)(2 * J);
    }
    else {
        throw std::runtime_error("Invalid value of J - must be integer or half-integer");
    }
    m_convfact = lambda;
}

cfpars::cfpars(const std::string &ionname) {
    getfromionname(ionname);
    double alpha = ALPHA_J[m_n];
    double beta = BETA_J[m_n];
    double gamma = GAMMA_J[m_n];
    m_convfact = lambda;
	m_stevfact = {alpha, beta, gamma};
	m_convertible = true;
}

} // namespace libMcPhase
