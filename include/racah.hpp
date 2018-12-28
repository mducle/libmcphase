/* racah.hpp
 * 
 * Class for calculating the n-j symbols of Wigner and Racah.
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#ifndef NJSYMS_H
#define NJSYMS_H

#include<algorithm>
#include<array>
#include<cmath>

namespace libMcPhase {

class racah {

    private:
        // Array to hold results of factorial calculations (170! is highest storable by a double)
        // The first row is i!, the second row is i!/1!, third is i!/2!, jth row is i!/j!
        std::array<std::array<double, 171>, 171> m_facts{};
        int m_highest = 0;                                                           // Highest factorial already calculated
        int m_max = 170;                                                             // Maximum value of factorial calculable
        void m_calc_f(int v);

    public:
        // Constructors
        racah() { m_facts[0][0] = 1; }; 
        // Methods
        double f(int v);                                                             // Returns v!
        double f_product(std::vector<int> v);                                        // Returns v[0]!v[1]!v[2]!... etc.
        double f_product_pz(std::vector<int> v, int z);                              // Returns (v[0]+z)!(v[1]+z)!... etc.
        double f_quotient(int x, int y);                                             // Returns x!/y!
        double tri(int a, int b, int c);                                             // Returns the triangular factor T(a,b,c)
        // Note the functions below use twice the value of any arguments, so as to 
        // represent half integer numbers
        double threej(int a, int b, int c, int d, int e, int f);                     // Calculates the 3j symbol (abd;def)
        double sixj(int a, int b, int c, int d, int e, int f);                       // Calculates the 6j symbol {abc;def}
        double ninej(int a, int b, int c, int d, int e, int f, int g, int h, int i); // Calculates the 9j symbol {abd;def;ghi}
        double wigner(int a, int b, int c, int d, int e, int f);                     // Calculates Wigner coeff. (abcd|ef)
        double wigner(int a, int b, int c, int d, int e, int f, int g, int h);       // Calculates Wigner coeff. (abcd|abef)
        double racahW(int a, int b, int c, int d, int e, int f);                     // Calculates Racah's W(abcd;ef) symbol

}; // class racah

} // namespace libMcPhase

#endif
