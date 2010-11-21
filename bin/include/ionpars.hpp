#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>
#include <mpspecfunp.h>

// ionpars: class to 
//          - read single ion parameterfile for cfield module
//          - load and store matrices for internal module cfield
//          - diagonalize cfproblem and calculate moment and transition matrix elements

class ionpars  
{private: 
   // calculates scattering operator 
   void MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc);
  int calcmagdensity;  // 0 ... normal mode, 1,2,3 calc <J'i>=gJ/2 (<J1,2,3 * Ji>+<Ji*J1,2,3>) ... gives magnetisationdensity in a b c dir instead
                        // of chargedensiy in chrgplt,charges ...
  
 public:
   int so1ion; // switch to identify the coordinate system orientation with respect to the axes abc
              // 0 ...  cfield xyz||cab
              // 1 .... so1ion xyz||abc
   char * iontype; // description string
   int nof_electrons; // nof 4f electrons
   double J;// momentum quantum number
   double gJ; // Lande factor
   double alpha;
   double beta;  // stevens factors
   double gamma;
   double r2;
   double r4;  // radial wave function exp values
   double r6;
  
   Matrix Ja; Matrix Jb; Matrix Jc; Matrix Hcf;
   ComplexMatrix Jaa;
   ComplexMatrix Jbb;
   ComplexMatrix Jcc;
 

   Matrix **Olm; ComplexMatrix **OOlm; // array of matrices
 
   Vector Blm; // Cf parameters  
   Vector Llm; // Cf parameters  

   // functions needed to calculate thermal expectation value of moment  
   Vector & cfield (double & T,Vector & H, double & Z,double & U, ComplexMatrix & ests);
   void cfieldJJ(Vector & JJ,double & T, Vector & gjmbH, double & lnZs, double & U, ComplexMatrix & ests);
   void cfeigenstates (ComplexMatrix *est, Vector & H, double & T);
   // and transition matrix elements
   int  cfielddm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta,ComplexMatrix & ests);
   int cfielddn(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexMatrix & nat);
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est
   // calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
   // system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
   ComplexVector & MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est);


   void savBlm(FILE * file); // saving Blm to file 
   void savLlm(FILE * file); // saving Blm to file 
   void save(FILE * file); // save ion parameters to file 


   ionpars(int dimj);
   ionpars(FILE * cf_file);
   ionpars (char * iontype); // constructor from iontype (mind:no matrices filled with values !)
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif
