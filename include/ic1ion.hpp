/* ic1ion.hpp
 * 
 * Header file for the ic1ionmodule
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#ifndef IC1ION_H
#define IC1ION_H

#include "cfpars.hpp"
#include "eigen.hpp"
#include "ic_states.hpp"

namespace libMcPhase {

class ic1ion : public cfpars {

    protected:
        bool m_ham_calc = false;                              // Flag to indicate if Hamiltonian calculated
        bool m_ev_calc = false;                               // Flag to indicate if eigenvectors/values calculated
        RowMatrixXcd m_hamiltonian;                           // Cached Hamiltonian
        RowMatrixXcd m_eigenvectors;                          // Cached eigenvectors
        VectorXd m_eigenvalues;                               // Cached eigenvalues
		// Declarations for functions in so_cf.cpp
		RowMatrixXd racah_so(int n, double xi, orbital e_l=F);                   // Calculates the spin-orbit matrix
		RowMatrixXd racah_Umat(int n, int k, orbital e_l=F);                     // Calculates the reduced matrix U^k
		RowMatrixXd racah_uJ(int n, int k, orbital e_l=F);                       // Calculates the U^k redmat in |LSJ>
		RowMatrixXd racah_ukq(int n, int k, int q, orbital e_l=F);               // Calculates the tensor operator U^k_q 
		RowMatrixXd fast_ukq(int n, int k, int q, orbital e_l=F);                // Calculates the tensor operator U^k_q 
		RowMatrixXd racah_mumat(int n, int q, orbital e_l=F);                    // Calculates the magnetic moment operator
		void racah_mumat(int n,int q,RowMatrixXd&L, RowMatrixXd&S, orbital l=F); // Calculates the magnetic moment operators
		void chanlam_mumat(int n,int q,RowMatrixXd&mu, orbital l=F);             // Calculates the magnetic moment operator
		// Declarations for functions in lovesey.cpp
		bool lovesey_aKK(RowMatrixXd &aKK, int K, int Kp, int n, orbital l);     // Calculates the matrix a(K,K')
		bool lovesey_cKK(RowMatrixXd &aKK, int K, int Kp, int n, orbital l);     // Calculates the matrix c(K,K')
		void lovesey_Qq(std::vector<RowMatrixXd >&Q, int q, int n, orbital l,    // Calculates the transition matrix Qq
		       std::vector<double>&);
		RowMatrixXd balcar_MSq(int q, int K, int Q, int n, orbital l);           // Calculates the coeff. of the spin density
		RowMatrixXd balcar_MLq(int q, int K, int Q, int n, orbital l);           // Calculates the coeff. of the orbital dens.
		RowMatrixXcd balcar_Mq(int xyz,int K,int Q,int n,orbital l);             // Driver for calculation of density coeff.

    public:
        // Setters
        virtual void set_unit(const Units newunit);
        virtual void set_type(const Type newtype);
        virtual void set_name(const std::string &ionname);
        virtual void set_J(const double J);
        virtual void set(const Blm blm, double val);
        virtual void set(int l, int m, double val);
        // Constructors
        ic1ion() : cfpars() {};
        ic1ion(const int J2) : cfpars(J2) {};
        ic1ion(const double J) : cfpars(J) {};
        ic1ion(const std::string &ionname) : cfpars(ionname) {};
        // Methods
        RowMatrixXcd hamiltonian(bool upper=true);
        std::tuple<RowMatrixXcd, VectorXd> eigensystem();

}; // class ic1ion

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in coulomb.cpp
// --------------------------------------------------------------------------------------------------------------- //
double racah_xwu(qR7 W, qG2 U, qG2 Up);                                   // Looks up values of x(W|UU')
double racah_chi(orbital L, orbital Lp, qG2 U, qG2 Up);                   // Looks up values of (U|chi(L)|U')
int racah_e2sign(int S2, int v);                                          // Looks up sign (+/-) of e2 (Racah 4, eq 73)
double racah_e2prod(qR7 W, qG2 U, qG2 Up, orbital L, orbital Lp);         // Calculates sum of x(W|UU')*(U|chi(L)|U')
RowMatrixXd racah_e2(int n);                                              // Calculates the matrix elements of e2
double racah_yfn(int n, int v, int S2, qG2 U, int vp, qG2 Up);            // Looks up values of y(f^n,vSU|v'S'U')
double racah_phi(qG2 U, qG2 Up, orbital Lp, orbital L);                   // Looks up values of (U|phi(L)|U')
double racah_g(qG2 U, bool R5flag=false);                                 // Calcs. eigenvalue of Casimir's op for G2 or R5
double racah_g(qR7 W);                                                    // Calcs. eigenvalue of Casimir's op for R7
RowMatrixXd racah_e3(int n);                                              // Calculates the matrix elements of e3
RowMatrixXd racah_emat(int n, double E0, double E1, double E2, double E3);//Calculates the Coulomb interaction matrix
RowMatrixXd racah_emat(int n, double F0, double F2, double F4);           // Calculates the Coulomb matrix for d-electrons
RowMatrixXd racah_emat(int n, double F0, double F2);                      // Calculates the Coulomb matrix for p-electrons
std::vector<double> racah_FtoE(std::vector<double> F);                    // Converts from F_k to E
std::vector<double> racah_EtoF(std::vector<double> E);                    // Converts from E to F_k
std::vector<double> racah_FtoF_k(std::vector<double> F);                  // Converts from F^k to F_k
std::vector<double> racah_F_ktoF(std::vector<double> F_k);                // Converts from F_k to F^k
RowMatrixXd racah_ci(int n, double alpha, double beta, double gamma);     // Calcs. conf. interaction matrix for f-elec.
RowMatrixXd racah_ci(int n, double alpha, double beta);                   // Calcs. conf. interaction matrix for d-elec.
RowMatrixXd racah_ci(int n, double alpha);                                // Calcs. conf. interaction matrix for p-elec.

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in cfp.cpp
// --------------------------------------------------------------------------------------------------------------- //
double racah_ulf(qG2 U, orbital L, qG2 Up, orbital Lp);                   // Calculates (UL|U'L'+f)
double racah_wupf(qR7 W, qG2 U, qR7 Wp, qG2 Up);                          // Calculates (WU|W'U'+f)
double racah_cfp(int n, qG2 U, int v, int S2, orbital L, qG2 Up, int vp, int S2p, orbital Lp);
double racah_cfp(int n, std::string child, std::string parent);           // Calculates the coeff. fractional parentage.
double racah_cfp(int n, int v, int S2, orbital L, int vp, int S2p, orbital Lp); // cfp for d-electrons
double racah_cfp(int n, int S2, orbital L, int S2p, orbital Lp);          // cfp for p-electrons
std::vector<cfpls> racah_parents(int n, int v, qG2 U, int S2, orbital L); // Calculates the cfp's of all the parents
std::vector<cfpls> racah_parents(int n, std::string state);               //    of a particular state.
std::vector<cfpls> racah_parents(int n, int v, int S2, orbital L);        // Calc. the cfp's of all parents for d-elec.
std::vector<cfpls> racah_parents(int n, int S2, orbital L);               // Calc. the cfp's of all parents for p-elec.
//  RowMatrixXd cfp_orthog_test(int n, const char* LS);                   // These two functions are for testing the
//  RowMatrixXd cfp_cowan_test(int n, const char* LSp);                   //    calculations of the cfp's.


/*
// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in ic_hmltn.cpp
// --------------------------------------------------------------------------------------------------------------- //
RowMatrixXd convH2H(RowMatrixXd H,int l,std::vector<std::vector<int> > c); // Converts from one basis to another
void convH2H(complexdouble *Hin, complexdouble *Hout, int lnIn, int lnOut, std::vector< std::vector<int> > cv);
bool ic_parseheader(const char *filename, icpars &pars);                  // Determines the paramters from saved file
RowMatrixXd ic_Hcso(icpars &pars);                                       // Calculates the Hamiltonian H=H_c+H_so
RowMatrixXd ic_hmltn(RowMatrixXd &H_cfi, icpars &pars);                 // Calculates the IC Hamiltonian matrix
std::vector<double> ic_mag(RowMatrixXd &Hic, RowMatrixXd &iHic,         // Calculates the magnetisation
       RowMatrixXd &Jmat, RowMatrixXd &iJmat, std::vector<double> &T, 
       double H_mag=1., int nev=30);

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in icpars.cpp
// --------------------------------------------------------------------------------------------------------------- //
int  getdim(int n, orbital l);                                            // Number of states = ^{4l+2}C_{n}
void strtolower(std::string &instring);                                   // Converts a string to lower case
void conv_e_units(icpars &flags, std::string &newunits);                  // Converts 1-ion pars to diff. energy units
std::vector<double> stev_thetak(int n, orbital l);                        // Calculates the Stevens factors.
std::vector<double> rk_int(std::string &ionname);                         // Looks up value of radial integrals

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions ic_diag.cpp
// --------------------------------------------------------------------------------------------------------------- //
int ic_diag(RowMatrixXd&Hic, RowMatrixXd&iH, complexdouble*V, double*E);// Diagonalises complex hermitian Hic+iH
int ic_diag(int n, complexdouble *zm, complexdouble *z, double *eigval);
int ic_diag(RowMatrixXd&Hic, double *V, double *E);                      // Diagonalises real symmetric Hic 
int ic_diag(double *mz, int lda, int n, double *m, double *eigval);
int ic_leig(RowMatrixXd&H,RowMatrixXd&i,complexdouble*V,double*E,int n);// Finds only the n lowest eigenval/vec
int ic_leig(int n, complexdouble *zm, complexdouble *z, double*E, int iu);
int ic_leig(RowMatrixXd&Hic, double *V, double *E, int n);               // Finds only the n lowest eigenval/vec
#ifndef NO_ARPACK
int ic_arpackeig(int n, complexdouble*m,complexdouble*z, double*E,int iu);// Finds only lowest eigval/vecs, use ARPACK
int ic_arpackeig(int n, double *m, double *z, double *E, int iu);         // Finds only lowest eigval/vecs, use ARPACK
#endif

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in spectre.cpp
// --------------------------------------------------------------------------------------------------------------- //
void spectre_vrot(icpars &pars, complexdouble *Vrot, double &elim);
void spectre_vrot(icpars &pars, double *Vrot, double &elim);
iceig spectre_eig(RowMatrixXd Hic, double *Vrot, int cb);
iceig spectre_eig(RowMatrixXd Hic, RowMatrixXd iHic, complexdouble *Vrot, int cb);

// --------------------------------------------------------------------------------------------------------------- //
// Declarations for functions in ic1ion.cpp
// --------------------------------------------------------------------------------------------------------------- //
void getfromionname(std::string &ion, icpars &flags);                     // Gets free ion parameters from tables
void ic_parsecfpars(std::string &n, std::string &v, icpars &p, int l=1);  // Parses CF parameter for k and q
void ic_parseinput(const char *file, icpars &flags);                      // Parses file for 1-ion pars & phys prop.
void ic_conv_basis(icpars &pars, iceig &VE, fconf &conf);                 // Converts eigenvectors to different basis
void ic_printheader(const char *filename, icpars &pars);                  // Prints header to file
void ic_showoutput(const char *file, icpars&pars, iceig&VE, int iconf=1,  // Prints calculated spectra to file
                   std::vector<int> ikeepJ = std::vector<int>());
void ic_cmag(const char *filename, icpars &pars, double elim=-DBL_MAX);   // Calcs. magnetisation using icmfmat::
*/

} // namespace libMcPhase

#endif
