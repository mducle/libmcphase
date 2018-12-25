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
#include<stdexcept>
#include<string>
#include<unordered_map>

namespace libMcPhase {

class cfpars {


    public:
    enum class Units {meV = 0, cm = 1, K = 2};
    enum class Normalisation {Stevens, Wybourne};
    enum class Type {Alm, Blm, Vlm, Wlm, Llm, ARlm};
    enum class Blm {B22S = 0, B21S = 1, B20 = 2, B21 = 3, B22 = 4,
                    B44S = 5, B43S = 6, B42S = 7, B41S = 8, B40 = 9, B41 = 10, B42 = 11, B43 = 12, B44 = 13,
                    B66S = 14, B65S = 15, B64S = 16, B63S = 17, B62S = 18, B61S = 19, 
                    B60 = 20, B61 = 21, B62 = 22, B63 = 23, B64 = 24, B65 = 25, B66 = 26};

    private:
        std::array<double, 27> m_Bi{};                        // Internal array of values (in Stevens in meV)
        std::array<double, 27> m_Bo{};                        // Output array of parameters.
        std::array<double, 3> m_rk{};                         // Radial integrals <r^k> (if constructed from ionname)
        std::array<double, 3> m_stevfact = {{1., 1., 1.}};    // Stevens factors \theta_k
        std::array<double, 3> m_invstevfact = {{1., 1., 1.}}; // Inverse Stevens factors
        std::array<double, 27> m_convfact;                    // Conversion factor from internal to external parameters
        Units m_unit = Units::meV;                            // Energy units of parameters
        Normalisation m_norm = Normalisation::Stevens;        // Normalisation (Stevens or Wybourne)
        Type m_type = Type::Blm;                              // Type of crystal field parameters
        std::string m_ionname;                                // Name of ion
        int m_J2 = 0;                                         // 2*J == twice the total angular momentum
        bool m_convertible = false;                           // True if can convert between types and normalisations
        
    public:
        // Getters
        const Units get_unit() const { return m_unit; }
        const Normalisation get_normalisation() const { return m_norm; }
        const Type get_type() const { return m_type; }
        const std::string get_name() const { return m_ionname; }
        const double get(const Blm blm) const { return m_Bo[(int)blm]; }
        const double get(int l, int m) const;
        const double alpha() { return m_stevfact[0]; }
        const double beta() { return m_stevfact[1]; }
        const double gamma() { return m_stevfact[2]; }
        const double get_J() { return (double)(m_J2 / 2.); }
        // Setters
        void set_unit(const Units newunit);
        void set_type(const Type newtype);
        void set_name(const std::string &ionname);
        void set_J(const double J);
        void set(const Blm blm, double val);
        void set(int l, int m, double val);
        // Constructors
        cfpars() { m_convfact.fill(1.); };
        cfpars(const int J2) : m_J2(J2) { m_convfact.fill(1.); };
        cfpars(const double J);
        cfpars(const std::string &ionname);
    
}; // class cfpars

} // namespace libMcPhase

#endif
