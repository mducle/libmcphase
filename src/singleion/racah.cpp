/* racah.cpp
 * 
 * Class for calculating the n-j symbols of Wigner and Racah.
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

// --------------------------------------------------------------------------------------------------------------- //
// NB. Note that all the functions in this file uses twice the value of any arguments - in order to represent
//     half-integral arguments as integers.
// --------------------------------------------------------------------------------------------------------------- //

#include "racah.hpp"

#define TRI(A,B,C) ( (C > (A+B)) || (C < abs(A-B)) )   // The triangular inequality

namespace libMcPhase {

// --------------------------------------------------------------------------------------------------------------- //
// Methods for racah class
// --------------------------------------------------------------------------------------------------------------- //
void racah::m_calc_f(int v) {
    // Calculates the factorial table up to v! (assuming that v > m_highest)
    for (int l=m_highest+1; l<=v; l++) {
        m_facts[0][l] = m_facts[0][l-1] * l;
        m_facts[l][0] = m_facts[l-1][0] / l;
        m_facts[l][l] = l;
        if(m_highest > 0) {
            for (int n=1; n<=m_highest; n++) {
                m_facts[n][l] = m_facts[n][l-1] * l;
                m_facts[l][n] = m_facts[l-1][n] / l;
            }
        }
        for (int m=l; m<=v; m++) {
            m_facts[l][m] = m_facts[l][m-1] * m;
            m_facts[m][l] = m_facts[m-1][l] / m;  // Lower triangle is inverse of upper triangle
        }
    }
    m_highest = v;
}

double racah::f(int v) {
    // Returns v! and populates the factorial table if needed.
    if (v < 0)
        return 0.;
    if (v > m_max) {
        throw std::runtime_error("Too high value of factorial.");
    }
    if (v > m_highest) {
        m_calc_f(v);
    }
    return m_facts[0][v];
}

double racah::f_product(std::vector<int> v) {
    // Returns product of factorials of elements of v
    double out = 1.;
    for (auto x: v) {
        out *= f(x);
    }
    return out;
}

double racah::f_product_pz(std::vector<int> v, int z) {
    // Returns product of factorials of elements of v plus a constant z
    double out = 1.;
    for (auto x: v) {
        out *= f(x+z);
    }
    return out;
}

double racah::f_quotient(int x, int y) {
    // Return x!/y!
    if (x < 0 || y < 0)
        return 0.;
    int mx = std::max(x, y);
    if (mx > m_max) {
        throw std::runtime_error("Too high value of factorial.");
    }
    if (mx > m_highest) {
        m_calc_f(mx);
    }
    return m_facts[y][x];
}

double racah::tri(int a, int b, int c) {
    // Return the triangular factor: (a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!
    // Note that in this implementation, a, b and c take twice their actual values
    return f((a+b-c)/2) * f((a-b+c)/2) * f_quotient((-a+b+c)/2, (a+b+c)/2+1);
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 3-j symbol.
// --------------------------------------------------------------------------------------------------------------- //
double racah::threej(int j1, int j2, int J, int m1, int m2, int M)
{
    // The analytical expression for the Wigner 3j symbol is given by the Racah formula
    // ( j1 j2 J )        (j1-j2-M)  
    // ( m1 m2 M ) = (-1)^          sqrt(T(j1,j2,J)) sqrt((j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(J+M)!(J-M)!)
    //                                               t
    //               * SUM _____________________(-1)^________________________________
    //                  t  t!(J-j2+t+m1)!(J-j1+t-m2)!(j1+j2-J-t)!(j1-t-m1)!(j2-t+m2)!

    // Selection rules:
    if (
        (m1 + m2 + M) != 0 ||                                        // 1. m1 + m2 + M = 0
        TRI(j1, j2, J) ||                                            // 2. |j1-j2| < J < |j1+j2|
        (m1>j1 || m1<-j1) || (m2>j2 || m2<-j2) || (-M>J || -M<-J) || // 3. -j1<=m1<=j1, -j2<=m2<=j2, -J<=M<=J
        ((j1 + j2 + J) % 2) != 0                                     // 4. Integer perimeter rule
       ) {
        return 0.;
    }

    std::vector<int> Fp = {(J-j2+m1)/2, (J-j1-m2)/2};
    std::vector<int> Fm = {(j1+j2-J)/2, (j1-m1)/2, (j2+m2)/2};

    // Sum over t runs over all positive integers of the arguments of the factorials.
    int mint = 0; 
    for(auto t: Fp)
        if(-t > mint) 
            mint = -t;

    int maxt = 171;
    for(auto t: Fm)
        if(t < maxt) 
            maxt = t;

    double sum_t = 0.;
    for (int t = mint; t <= maxt; t++) {
        sum_t += pow(-1., t) / ( f(t) * f_product_pz(Fp, t) * f_product_pz(Fm, -t) );
    }

    std::vector<int> P = {(j1+m1)/2, (j1-m1)/2, (j2+m2)/2, (j2-m2)/2, (J+M)/2, (J-M)/2};
    sum_t = pow(-1., (j1-j2-M)/2.) * sqrt(tri(j1, j2, J) * f_product(P)) * sum_t;
    return sum_t;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 6-j symbol.
// --------------------------------------------------------------------------------------------------------------- //
double racah::sixj(int a, int b, int c, int d, int e, int ff)
{
    // The 6-j symbols can be computed by the Racah formula:
    //
    //                                       1/2
    // { a b c } = [T(abc)T(aef)T(dbf)T(dec)]
    // { d e f }
    //                                                     z
    //              ---                                (-1)  * (z+1)!
    //            * >    ----------------------------------------------------------------------------
    //              ---z (z-a-b-c)!(z-a-e-f)!(z-d-b-f)!(z-d-e-c)!(a+b+d+e-z)!(b+c+e+f-z)!(a+c+d+f-z)!
 
    // Selection rules: The triads (a b c), (a e f), (d b f), (d e c) all:
    if (
        TRI(a,b,c) || TRI(a,e,ff) || TRI(d,b,ff) || TRI(d,e,c) ||              // 1. Satisfy |a-b| <= c <= |a+b|
        ((a+b+c)%2!=0) || ((a+e+ff)%2!=0) || ((d+b+ff)%2!=0) || ((d+e+c)%2!=0) // 2. Sum to integers
       ) {
        return 0.;
    }

    std::vector<int> Fp = {(-a-b-c)/2, (-a-e-ff)/2, (-d-b-ff)/2, (-d-e-c)/2};
    std::vector<int> Fm = {(a+b+d+e)/2, (b+c+e+ff)/2, (a+c+d+ff)/2};

    int minz = 0; 
    for(auto z: Fp)
        if(-z > minz) 
            minz = -z;

    int maxz = 171;
    for(auto z: Fm)
        if(z < maxz) 
            maxz = z;

    double sumz = 0.;
    for (int z = minz; z <= maxz; z++) {
        sumz += pow(-1., z) * f(z+1) / ( f_product_pz(Fp, z) * f_product_pz(Fm, -z) );
    }

    return sqrt(tri(a, b, c) * tri(a, e, ff) * tri(d, b, ff) * tri(d, e, c)) * sumz;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the 9j symbols
// --------------------------------------------------------------------------------------------------------------- //
double racah::ninej(int j1, int j2, int J12, int j3, int j4, int J34, int J13, int J24, int J)
{
    // The equation for the 9-j symbol is:
    //
    // { j1  j2  J12 }   ---     2g
    // { j3  j4  J34 } = >   (-1)    (2g+1)  { j1  j2 J12 } { j3 j4 J34 } { J13 J24 J  }
    // { J13 J24  J  }   ---                 { J34  J  g  } { j2 g  J24 } {  g  j1  j3 }
    //                    g
    //
    // where the two-row curly brackets indicate the 6-j symbols. The 6-j symbols are equal
    // to a permutation of their columns, so:
    //
    // { j3 j4 J34 } = { J34 j3 j4 }  and  { J13 J24 J  } = { J24 J  J13 }
    // { j2 g  J24 }   { J24 j2 g  }       {  g  j1  j3 }   { j1  j3  g  }
    //
    // Now, the 6-j symbols { j1 j2 j3 } have triads (j1 j2 j3), (j1 J2 J3), (J1 j2 J3) and (J1 J2 j3) which satisfy the triangular
    //                      { J1 J2 J3 }   inequalities. In particular, for those involving g [e.g. (j1 J2 J3)], |j1+J2| <= J3 <= j1+J2
    //
    // Thus for the 9j symbols we must have:   |J-j1|  <= g <=  J+j1
    //                                        |J34-j2| <= g <= J34+j2
    //                                        |J24-j3| <= g <= J24+j3

    int g, min_g, max_g;
    int G[] = {abs(J-j1), abs(J34-j2), abs(J24-j3), J+j1, J34+j2, J24+j3};
    double out = 0.;

    // Finds the allowed values of g
    min_g = G[0]; for(g=1; g<3; g++) if(G[g]<min_g) min_g = G[g];
    max_g = G[3]; for(g=4; g<6; g++) if(G[g]>max_g) max_g = G[g];
    for(g=min_g; g<=max_g; g++)  // g==2g, as all integers here represents twice their values (to accomodate half integral values).
        out += pow(-1.,g) * (g+1) * sixj(j1,j2,J12,J34,J,g) * sixj(j3,j4,J34,j2,g,J24) * sixj(J13,J24,J,g,j1,j3);

    return out;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the wigner coefficient (a b c d | a b e f) == (a  b  c  d | e  f) == (j1 j2 m1 m2 | j -m)
// --------------------------------------------------------------------------------------------------------------- //
double racah::wigner(int j1, int j2, int m1, int m2, int j, int m)
{
    // The analytical expression for the Wigner coefficient is given by: (cf. Racah II)
    //
    // (j j m m |j j jm)      ___________________
    //   1 2 1 2  1 2    =  \/ (2j+1) T(j1,j2,j) 
    //                                          ______________________________________________
    //                        ---     z       \/ (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j+m)!(j-m)! 
    //                      * >   (-1)   ----------------------------------------------------------
    //                        ---        z!(j1+j2-j-z)!(j1-m1-z)!(j2+m2-z)!(j-j2+m1+z)!(j-j1-m2+z)!
    //                         z

    // Selection rules:
    if (
        (m1 + m2 - m) != 0 ||                                        // 1. m1+m2-m = 0
        TRI(j1, j2, j) ||                                            // 2. |j1-j2| <= j <= |j1+j2|
        (m1>j1 || m1<-j1) || (m2>j2 || m2<-j2) || (-m>j || -m<-j) || // 3. -j1 <= m1 <= j1, -j2 <= m2 <= j2, -j <= -m <= j
        (j1 + j2 + j) % 2 != 0                                       // 4. Integer perimeter rule: j1+j2+j is integer
       ) {
        return 0.;
    }

    std::vector<int> Fm = {(j1+j2-j)/2, (j1-m1)/2, (j2+m2)/2};
    std::vector<int> Fp = {(j-j2+m1)/2, (j-j1-m2)/2};

    // Sum over z runs over all positive integers of the arguments of the factorials.
    int minz = 0; 
    for(auto z: Fp)
        if(-z > minz) 
            minz = -z;

    int maxz = 171;
    for(auto z: Fm)
        if(z < maxz) 
            maxz = z;

    double sumz = 0.;
    for (int z = minz; z <= maxz; z++) {
        sumz += pow(-1., z) * ( f(z) * f_product_pz(Fm, -z) * f_product_pz(Fp, z) );
    }

    std::vector<int> P = {(j1+m1)/2, (j1-m1)/2, (j2+m2)/2, (j2-m2)/2, (j+m)/2, (j-m)/2};
    return sqrt((j+1) * tri(j1, j2, j) * f_product(P)) / sumz;
}

double racah::wigner(int a, int b, int c, int d, int e, int f, int g, int h)
{  
    if(e==a && f==b) 
        return wigner(a,b,c,d,g,h); 
    else 
        return 0.; 
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the value of Racah's W symbols W(abcd;ef) which are related to the 6-j symbols
// --------------------------------------------------------------------------------------------------------------- //
double racah::racahW(int a, int b, int c, int d, int e, int f)
{
    // The formula for the Racah W function is:       where the triangle functions are:
    //                                           
    //                                        1/2                                        (a+b-c)!(a-b+c)!(-a+b+c)!
    // W(abcd;ef) = [T(abe)T(cde)T(acf)T(bdf)]                                  T(abc) = -------------------------
    //                                                                                       (a + b + c + 1)!
    //               ---                     (a + b + c + d + 1 - z)!
    //             * >    ------------------------------------------------------------------
    //               ---z z!(a+b-e-z)!(c+d-e-z)!(a+c-f-z)!(b+d-f-z)!(e+f-a-d+z)!(e+f-b-c+z)!
    //
    // The W functions are only defined for integer or half integer arguments that satisfy the selection rule that 
    // the four triads (abe), (cde), (acf), and (bdf) has an integral sum and satisfy the triangular inequality 
    // |a-b| <= c <= |a+b|. Otherwise W is zero.

    // Selection rules: The triads (a b e), (c d e), (a c f), (b d f)
    if( 
        // 1. All satisfy the triangular inequality: |a-b| <= c <= a+b 
        (TRI(a,b,e) || TRI(c,d,e) || TRI(a,c,f) || TRI(b,d,f)) ||
        // 2. Elements of each triad sum to an integer  [a,b,c,d,e,f are actually twice their nominal values]
        (((a+b+e)%2!=0) || ((c+d+e)%2!=0) || ((a+c+f)%2!=0) || ((b+d+f)%2!=0) )
        ) {
        return 0.;
    }
    
    std::vector<int> Fm = {(a+b-e)/2, (c+d-e)/2, (a+c-f)/2, (b+d-f)/2};
    std::vector<int> Fp = {(e+f-a-d)/2, (e+f-b-c)/2};

    // Sum over z runs over all positive integers of the arguments of the factorials.
    int minz = 0; 
    for(auto z: Fp)
        if(-z > minz) 
            minz = -z;

    int maxz = a + b + c + d + 1; 
    for(auto z: Fm)
        if(z < maxz) 
            maxz = z;

    double sumz = 0.;
    for (int z = minz; z <= maxz; z++) {
        sumz += pow(-1., z) * f_quotient((a + b + c + d)/2 + 1 - z, z) / f_product_pz(Fm, -z) / f_product_pz(Fp, z);
    }

    return sumz * sqrt(tri(a, b, e) * tri(c, d, e) * tri(a, c, f) * tri(b, d, f));
}

} // namespace libMcPhase

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ njsyms.cpp && ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*#include<iostream>
int main(int argc, char *argv[])
{
   int a,b,c,d,e,f,g,h,i;

   if(argc<7)
   {
      a = 3; b = 4; c = 1; d = 2; e = 1; f = 4;  // for racahW and sixj
   }
   else
   {
      a = atoi(argv[1])*2; b = atoi(argv[2])*2; c = atoi(argv[3])*2;
      d = atoi(argv[4])*2; e = atoi(argv[5])*2; f = atoi(argv[6])*2;
   }
   
   // Hmmm... problems with sizes of ints with factorials > 12 - had to modify factorial() to use long ints - again probs at >20!
   std::cout << "3! = " << factorial(3) << "; 0! = " << factorial(0) << "; 1! = " << factorial(1) << "\n";
   std::cout << "11! = " << factorial(11); std::cout << "; 12! = " << factorial(12) << "; 13! = " << factorial(13) << "\n";
   std::cout << "14! = " << factorial(14); std::cout << "; 15! = " << factorial(15) << "; 16! = " << factorial(16) << "\n";
   std::cout << "20! = " << factorial(20); std::cout << "; 21! = " << factorial(21) << "; 22! = " << factorial(22) << "\n";

   std::cout << "{" << a << "/2 " << b << "/2 " << c << "/2}\n" << "{" << d << "/2 " << e << "/2 " << f << "/2} = ";
   std::cout << sixj(a,b,c,d,e,f) << "\n\n";

   std::cout << "W(" << a << "/2 " << b << "/2 " << c << "/2 " << d << "; " << e << "/2 " << f << "/2) = " << racahW(a,b,c,d,e,f) << "\n";

   std::cout << "(" << a << "/2 " << b << "/2 " << c << "/2 " << d << "/2 | " << e << "/2 " << f << "/2) = " << wigner(a,b,c,d,e,f) << "\n\n";

   std::cout << "{" << a << "/2 " << b << "/2 " << c << "/2}\n" << "{" << e << "/2 " << d << "/2 " << c << "/2}\n" << "{" << f << "/2 " << f << "/2 " << " 0} = ";
   std::cout << ninej(a,b,c,e,d,c,f,f,0) << "\n\n";

   if(argc<7)
   {
      a = 1; b = 2; c = 1; d = 1; e = 0; f = -1;   // for threej
   }
   std::cout << "(" << a << "/2 " << b << "/2 " << c << "/2)\n" << "(" << d << "/2 " << e << "/2 " << f << "/2) = ";
   std::cout << threej(a,b,c,d,e,f) << "\n\n";

   return 0;
}*/
