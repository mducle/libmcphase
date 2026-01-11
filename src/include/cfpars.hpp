/* cfpars.hpp
 * 
 * A class encapsulating crystal field parameters and their conversions
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#ifndef CFPARS_H
#define CFPARS_H

#include<array>
#include<algorithm>
#include<cctype>
#include<cmath>
#include<complex>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<limits>

#include "racah.hpp"
#include "eigen.hpp"

namespace libMcPhase {

// Conversion factors for different energy units[from][to], order: [meV, cm, K].
static const std::array<double, 3> ENERGYCONV = { {1., 8.065544005, 11.6045221} };

class cfpars {

    public:
    enum class Units {meV = 0, cm = 1, K = 2};
    enum class Normalisation {Stevens, Wybourne};
    enum class Type {Alm, Blm, Llm, ARlm, Nlm};
    enum class Blm {B22S = 0, B21S = 1, B20 = 2, B21 = 3, B22 = 4,
                    B44S = 5, B43S = 6, B42S = 7, B41S = 8, B40 = 9, B41 = 10, B42 = 11, B43 = 12, B44 = 13,
                    B66S = 14, B65S = 15, B64S = 16, B63S = 17, B62S = 18, B61S = 19, 
                    B60 = 20, B61 = 21, B62 = 22, B63 = 23, B64 = 24, B65 = 25, B66 = 26};
    enum class Sym {Ci=0, C1=1, C2=10, Cs=11, C2h=12, C2v=20, D2=21, D2h=22, C4=30, S4=31, C4h=32,
                    D4=40, C4v=41, D2d=42, D4h=43, C3=50, S6=51, D3=60, C3v=61, D3d=62, C6=71, C3h=72, C6h=73,
                    D6=80, C6v=81, D3h=82, D6h=83, T=90, Th=91, Td=92, O=93, Oh=94};

    protected:
        std::array<double, 27> m_Bi{};                        // Internal array of values (in Wybourne/theta_k in meV)
        std::array<double, 27> m_Bo{};                        // Output array of parameters.
        std::array<double, 3> m_rk{};                         // Radial integrals <r^k> (if constructed from ionname)
        std::array<double, 3> m_stevfact = {{1., 1., 1.}};    // Stevens factors \theta_k
        std::array<double, 27> m_convfact;                    // Conversion factor from internal to external parameters
        double m_econv = 1.;                                  // Conversion factor from internal to external energy
        Units m_unit = Units::meV;                            // Energy units of parameters
        Normalisation m_norm = Normalisation::Stevens;        // Normalisation (Stevens or Wybourne)
        Type m_type = Type::Blm;                              // Type of crystal field parameters
        std::string m_ionname;                                // Name of ion
        int m_J2 = 0;                                         // 2*J == twice the total angular momentum
        double m_GJ = -1.;                                    // Lande g-factor
        bool m_convertible = false;                           // True if can convert between types and normalisations
        racah m_racah{};                                      // Class to calc n-j symbols and cache factorials
        int m_n = 1;                                          // Number of open shell electrons in this configuration
        Sym m_sym = Sym::C1;                                  // Point symmetry of site
        
    public:
        // Methods
        virtual void getfromionname(const std::string &ionname);
        bool is_sym_allowed(const Blm blm);
        std::vector<Blm> get_sym_allowed_Blm();
        // Getters
        Units get_unit() const { return m_unit; }
        Normalisation get_normalisation() const { return m_norm; }
        Type get_type() const { return m_type; }
        std::string get_name() const { return m_ionname; }
        double get(const Blm blm) const { return m_Bo[(int)blm]; }
        double get(int l, int m) const;
        double get_GJ() const { return m_GJ; }
        double alpha() const { return m_stevfact[0]; }
        double beta() const { return m_stevfact[1]; }
        double gamma() const { return m_stevfact[2]; }
        double get_J() const { return (double)(m_J2 / 2.); }
        Sym get_sym() const { return m_sym; }
        // Setters
        virtual void set_unit(const Units newunit);
        virtual void set_type(const Type newtype);
        virtual void set_name(const std::string &ionname);
        virtual void set_J(const double J);
        virtual void set_GJ(const double GJ) { m_GJ = GJ; }
        virtual void set(const Blm blm, double val);
        virtual void set(int l, int m, double val);
        void set_sym(const Sym newsym) { m_sym = newsym; }
        // Constructors
        cfpars();
        cfpars(const int J2);
        cfpars(const double J);
        cfpars(const std::string &ionname);
    
}; // class cfpars

} // namespace libMcPhase

#endif
