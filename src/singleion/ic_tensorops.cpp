/* ic_tensorops.cpp
 *
 * Contains methods to calculate the tensor operators for the ic1ion class.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 (icmf.cpp), 2020 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"

namespace libMcPhase {

// The tensor operators are saved within ic1ion as a vector of matrices. The following constants are indices into
// the rank k and order q of these matrices.
static const std::array<int, 51> k = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
static const std::array<int, 51> q = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the tensor operators up to 
// --------------------------------------------------------------------------------------------------------------- //
void ic1ion::calc_tensorops(int num_op) {
    if (num_op < m_tensorops.size() && num_op <= 0)
        return;
    // Calculates the L and S operator matrix for the each directions
    // Note operators up to 6 are the magnetic moment operators in x, y, z basis 
    // rather than tensor operators in the radial basis.
    // They are, in order: Sx, Lx, Sy, Ly, Sz, Lz
    RowMatrixXd Sp1, Sm1, Lp1, Lm1; 
    racah_mumat(1, Lp1, Sp1);
    racah_mumat(-1, Lm1, Sm1);
    // Note also that m_tensorops is a vector of _real_ matrices.
    // BUT negative order (-q) operators are purely imaginary, as are the Sy and Ly matrices.
    m_tensorops.push_back((Sm1 - Sp1) / sqrt(2.));  // Sx
    m_tensorops.push_back((Lm1 - Lp1) / sqrt(2.));  // Lx
    m_tensorops.push_back((Sm1 + Sp1) / sqrt(2.));  // Sy
    m_tensorops.push_back((Lm1 + Lp1) / sqrt(2.));  // Ly
    if(num_op > 4) {
        RowMatrixXd Sz, Lz;
        racah_mumat(0, Lz, Sz);
        m_tensorops.push_back(Sz);
        m_tensorops.push_back(Lz);
    }
    // For higher rank tensor operators use racah_ukq to calculate in radial basis
    for (int i=6; i<num_op; i++) {
        m_tensorops.push_back(racah_ukq(k[i], q[i]));
    }
}

RowMatrixXcd ic1ion::zeeman_hamiltonian(double H, std::vector<double> Hdir) {
    if (Hdir.size() != 3) {
        throw std::runtime_error("ic1ion::zeeman_hamiltonian(): Hdir must be a 3-vector");
    }
    // Normalise the field direction vector
    double Hnorm = sqrt(Hdir[0] * Hdir[0] + Hdir[1] * Hdir[1] + Hdir[2] * Hdir[2]);
    if (fabs(Hnorm) < 1.e-6) {
        throw std::runtime_error("ic1ion::magnetisation(): Direction vector cannot be zero");
    }
    std::vector<double> nHdir;
    std::transform(Hdir.begin(), Hdir.end(), std::back_inserter(nHdir), [Hnorm](double Hd){ return Hd / Hnorm; });
    // Calculates the L, S operators if not done already.
    calc_tensorops(6);
    RowMatrixXcd zeeman = RowMatrixXcd::Zero(m_tensorops[0].rows(), m_tensorops[0].cols());
    zeeman.real() = (GS * m_tensorops[0] + m_tensorops[1]) * nHdir[0]   // gSx + Lx
                   +(GS * m_tensorops[4] + m_tensorops[5]) * nHdir[2];  // gSz + Lz
    zeeman.imag() = (GS * m_tensorops[2] + m_tensorops[3]) * nHdir[1];  // gSy + Ly
    double H_in_energy = H * MU_Bc * m_econv; // MU_B in internal units (cm/T); want H in external energy units
    return (zeeman * H_in_energy);
}

std::vector<RowMatrixXcd> ic1ion::calculate_moments_matrix(RowMatrixXcd ev) {
    std::vector<RowMatrixXcd> moments;  // Vector of x, y, z moment (matrix elements) matrices; magnetisation
                                        //   is calculated from diagonal; susceptibility from all elements
    calc_tensorops(6);                  // Calculates the L, S operators if not done already.
    for (int ii=0; ii<6; ii+=2) {
        RowMatrixXcd Jmat = RowMatrixXcd::Zero(ev.rows(), ev.cols());
        if (ii == 2)
            Jmat.imag() = GS * m_tensorops[ii] + m_tensorops[ii+1];
        else
            Jmat.real() = GS * m_tensorops[ii] + m_tensorops[ii+1];
        moments.push_back((ev.adjoint()) * (Jmat * ev));
    }
    return moments;
}

/*
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the mean field matrix sum_i (H_i*J_i)
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::Jmat(sMat<double>&Jmat, sMat<double>&iJmat, std::vector<double>&gjmbH, bool save_matrices)
{
   int i; Jmat.zero(J[0].nr(),J[0].nc()); iJmat.zero(J[0].nr(),J[0].nc()); 
   if(_num_op<(int)gjmbH.size()) { iflag.resize(_num_op,0); _num_op = (int)gjmbH.size(); }
   for(i=0; i<(_num_op>6?6:_num_op); i++)
      if(fabs(gjmbH[i])>DBL_EPSILON*100) { if(iflag[i]==1) iJmat += J[i]*gjmbH[i]; else Jmat += J[i]*gjmbH[i]; }
   // Higher order than dipole operators needed
   if(_num_op>6)
   {
      char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
      if(save_matrices) {
      #ifndef _WINDOWS
      struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
         if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::Jmat(): Can't create mms dir, " << strerror(errno) << "\n";
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
         if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::Jmat(): Cannot create mms directory\n";
      #endif
      nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
      strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
      } else { strcpy(basename,"nodir/"); }
#define NSTR(K,Q) nstr[1] = K+48; nstr[2] = Q+48; nstr[3] = 0
#define MSTR(K,Q) nstr[1] = K+48; nstr[2] = 109;  nstr[3] = Q+48; nstr[4] = 0
      // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
      int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
      int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
    //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
      sMat<double> Upq,Umq; double redmat; int n = _n; //if(n>(2*_l+1)) n = 4*_l+2-n; 

      for(i=6; i<_num_op; i++)
      {
         if(q[i]<0) iflag[i]=1; 
         if (fabs(gjmbH[i])>DBL_EPSILON) 
         {
            redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[i],2*_l,0,0,0);// * wy2stev(i);
            if(k[i]%2==1) continue;   // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
            if(k[i]>4 && _l==D) continue;
            NSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[i],abs(q[i]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
            if(q[i]==0) { Jmat += Upq * (gjmbH[i]*redmat); continue; }
            MSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[i],-abs(q[i]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
            if(q[i]<0) { 
            // if((q[i]%2)==0) iJmat += (Upq - Umq) * (gjmbH[i]*redmat); else iJmat += (Upq + Umq) * (gjmbH[i]*redmat); } changed MR 15.12.09
               if((q[i]%2)==0) iJmat += (Umq - Upq) * (gjmbH[i]*redmat); else iJmat += (Umq + Upq) * (gjmbH[i]*redmat); }
            else {
               if((q[i]%2)==0)  Jmat += (Umq + Upq) * (gjmbH[i]*redmat); else  Jmat += (Umq - Upq) * (gjmbH[i]*redmat); } 
         } 
      }
   }
}

#define MAXNOFCHARINLINE 144

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the expectation values <V|J|V>exp(-beta*T) given a set of eigenstates
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> icmfmat::expJ(iceig &VE, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   double *vt=0, Z=0., U=0.; complexdouble *zt=0, zme;
   std::vector<double> E, ex((_num_op>6?_num_op:6)+2,0.), me, eb; matel.clear();
   int iJ, ind_j, Esz, Hsz=VE.Hsz(), incx=1; 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::expJ() - Hamiltonian matrix size not same as mean field operator!\n"; return E; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   // Checks that the eigenvalues are orthonormal
// char transa='C', transb='N'; double summm=0.;
// if(VE.iscomplex())
// {
//    complexdouble *zmm = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); 
//    complexdouble *vet = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); memcpy(vet,VE.zV(0),Hsz*Hsz*sizeof(complexdouble));
//    F77NAME(zgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &zalpha, vet, &Hsz, VE.zV(0), &Hsz, &zbeta, zmm, &Hsz);
//    for(int ii=0; ii<Hsz; ii++) { zmm[ii*Hsz+ii].r-=1.; summm += F77NAME(dzasum)(&Hsz, &zmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
//    std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
//    free(zmm); free(vet);
// }
// else
// {
//    double *dmm = (double*)malloc(Hsz*Hsz*sizeof(double)); 
//    double *vet = (double*)malloc(Hsz*Hsz*sizeof(double)); memcpy(vet,VE.V(0),Hsz*Hsz*sizeof(double));
//    F77NAME(dgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &alpha, vet, &Hsz, VE.V(0), &Hsz, &beta, dmm, &Hsz);
//    for(int ii=0; ii<Hsz; ii++) { dmm[ii*Hsz+ii]-=1.; summm += F77NAME(dasum)(&Hsz, &dmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
//    std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
//    free(dmm); free(vet);
// }

   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<J[0].nr(); Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0 || VE.E(Esz+1)==-DBL_MAX) break; }

   if (T<0){Esz=(int)(-T);printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Excitation Energy\n");
     for (ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } // MR 10.9.2010

   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) 
      if(fabs(VE.zV(ii,jj).r*VE.zV(ii,jj).r+VE.zV(ii,jj).i*VE.zV(ii,jj).i)<DBL_EPSILON*100000) 
      {
         VE.zV(ii,jj) = 0;
      }  
   // For first run calculate also the partition function and internal energy
   me.assign(Esz,0.); eb.assign(Esz,0.); Z=0.;
   if(!VE.iscomplex()) 
   {
      double *fJmat=J[0].f_array(); vt = (double*)malloc(Hsz*sizeof(double)); 
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
#ifdef _G77 
         F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
#else
         me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
#endif
//MR 10.9.2010
         if (T<0)
         {  char instr[MAXNOFCHARINLINE];
            printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
            if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
            eb[ind_j]=strtod(instr,NULL);
         }
        else
         { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 10.9.2010                                                                    !!!!    -----------------!!!!
      }
      free(fJmat); free(vt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }
   else
   {
      complexdouble *zJmat;
      zeroes.zero(J[0].nr(),J[0].nc()); if(iflag[0]==0) zJmat=zmat2f(J[0],zeroes); else zJmat = zmat2f(zeroes,J[0]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77 
         F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
         zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
         me[ind_j] = zme.r;
//MR 10.9.2010
         if (T<0)
         {  char instr[MAXNOFCHARINLINE];
            printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
            if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
            eb[ind_j]=strtod(instr,NULL);
         }
         else
         { eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 10.9.2010
      }
      free(zJmat); free(zt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }

   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::expJ(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::expJ(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   sMat<double> Upq,Umq; double redmat; int n = _n; //if(n>(2*_l+1)) n = 4*_l+2-n; 

   if(!_density.empty()) 
   { 
      std::cout << "Calculating the expectation of the moment density operator ";
      if(_density.find("l")!=std::string::npos || _density.find("s")!=std::string::npos) std::cout << _density << "\n"; 
      else
      {
         if(_density.find("1")!=std::string::npos) std::cout << "Sx\n";
         if(_density.find("2")!=std::string::npos) std::cout << "Lx\n";
         if(_density.find("3")!=std::string::npos) std::cout << "Sy\n";
         if(_density.find("4")!=std::string::npos) std::cout << "Ly\n";
         if(_density.find("5")!=std::string::npos) std::cout << "Sz\n";
         if(_density.find("6")!=std::string::npos) std::cout << "Lz\n";
      }
   }

   // Rest of the runs only calculate the new matrix elements
   for(iJ=1; iJ<(_num_op>6?_num_op:6); iJ++)
   {
      me.assign(Esz,0.);
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if(iJ>6 && (k[iJ]%2==1 || k[iJ]>_l*2)) { matel.push_back(me); continue; }
      if(!VE.iscomplex() && _density.empty())
      {
         if(iflag[iJ]==0) {
            double *fJmat; vt = (double*)malloc(Hsz*sizeof(double)); 
            if(iJ<6) 
               fJmat=J[iJ].f_array(); 
            else 
            {
               NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
               MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
               redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);// * wy2stev(iJ);
//             if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
               if(q[iJ]<0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
               #ifdef JIJCONV
               if(jijconv.size()>1) redmat*=jijconv[iJ+1];
               #endif
               Umq *= redmat; fJmat = Umq.f_array();
            }
            for(ind_j=0; ind_j<Esz; ind_j++)
            {  // Calculates the matrix elements <Vi|J.H|Vi>
               F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
#ifdef _G77 
               F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
#else
               me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
#endif
               ex[iJ]+=me[ind_j]*eb[ind_j];
            }
            free(fJmat); free(vt); matel.push_back(me); ex[iJ]/=Z; 
         } 
         else { me.assign(Esz,0.); matel.push_back(me); }
      }
      else
      {
         complexdouble *zJmat; zeroes.zero(J[0].nr(),J[0].nc());
         if(iJ<6) 
         {
            if(iflag[iJ]==0) zJmat=zmat2f(J[iJ],zeroes); else zJmat = zmat2f(zeroes,J[iJ]);
         }
         else 
         {
         // if(!_density.empty()) { zJmat = balcar_Mq(_density,k[iJ],q[iJ],_n,_l); } else {
            NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
            MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
            redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);// * wy2stev(iJ);
//          if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
            if(q[iJ]<0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
            #ifdef JIJCONV
            if(jijconv.size()>1) redmat*=jijconv[iJ+1];
            #endif
            Umq *= redmat; if(iflag[iJ]==0) zJmat=zmat2f(Umq,zeroes); else zJmat = zmat2f(zeroes,Umq); //}
         }
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J.H|Vi>
            F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
#ifdef _G77 
            F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
#else
            zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
#endif
            me[ind_j] = zme.r;
            ex[iJ]+=me[ind_j]*eb[ind_j];
         }
         free(zJmat); free(zt); matel.push_back(me); ex[iJ]/=Z;
      }
      if(fabs(ex[iJ])<DBL_EPSILON) ex[iJ]=0.; 
   }
   ex[iJ] = log(Z)-VE.E(0)/(KB*T); ex[iJ+1] = U;
   return ex;
}
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix M_ab=<i|Ja|j><j|Jb|i>{exp(-beta_i*T)-exp(-beta_j*T)} for some state i,j
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::u1(std::vector<double>&u, std::vector<double>&iu, iceig&VE, double T, int i, int j,int pr,float & delta, bool save_matrices)
{  double *vt=0, Z=0., therm; complexdouble *zt=0, zme; zme.r=0; zme.i=0.;
   int sz = (_num_op>6?_num_op:6);
   std::vector<double> mij(sz,0.);//, mji(6,0.);
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
// u.zero(sz); iu.zero(sz);
   int iJ, Hsz=VE.Hsz(), incx=1; 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::u1() - Hamiltonian matrix size not same as mean field operator!\n"; return; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   complexdouble *zJmat=0;
   char uplo = 'U';

   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::u1(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::u1(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
                 
   sMat<double> Upq,Umq; double redmat; int n = _n; if(n>(2*_l+1)) n = 4*_l+2-n; 

   // Calculates the matrix elements: <i|Ja|j> and <j|Ja|i> for each of the six Ja's
   for(iJ=0; iJ<sz; iJ++)
   {
//    if(k[iJ]%2==1) { if(VE.iscomplex()) { zij[iJ].r=0.; zij[iJ].i=0.; } else mij[iJ]=0.; continue; }
      if(iJ>=6)
      {
         NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
         MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
         #ifdef JIJCONV
         if(jijconv.size()>1) redmat*=jijconv[iJ];
         #endif
         redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);
//       if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
         if(q[iJ]<0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
         Umq *= redmat;
      }

      if(!VE.iscomplex() && im[iJ]==0)
      {
         vt = (double*)malloc(Hsz*sizeof(double)); 
         double *fJmat; if(iJ>=6) fJmat=Upq.f_array(); else fJmat=J[iJ].f_array();
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(j), &incx, &beta, vt, &incx);
         #ifdef _G77 
         F77NAME(ddot)(mij[iJ], &Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
         #else
         mij[iJ] = F77NAME(ddot)(&Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
         #endif
         free(fJmat); free(vt);
      } 
      else
      {
         zeroes.zero(J[0].nr(),J[0].nc());
         if(iJ>=6) { if(im[iJ]==0) zJmat=zmat2f(Umq,zeroes);   else zJmat = zmat2f(zeroes,Umq); }
         else      { if(im[iJ]==0) zJmat=zmat2f(J[iJ],zeroes); else zJmat = zmat2f(zeroes,J[iJ]); }
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(j), &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, VE.zV(i), &incx, zt, &incx);
         #else
         zij[iJ] = F77NAME(zdotc)(&Hsz, VE.zV(i), &incx, zt, &incx);
         #endif
//       int k;for(k=0;k<Hsz;++k)printf("%6.3f %+6.3f i  ",VE.zV(j)[k].r,VE.zV(j)[k].i);
         free(zJmat); free(zt);
      }
   }

   if(i==j&&T>0) {//subtract thermal expectation value from zij=zii
            std::vector< std::vector<double> > matel;
            std::vector<double> vJ = expJ(VE,T,matel,save_matrices);
            for(iJ=0; iJ<sz; iJ++)zij[iJ].r-=vJ[iJ];
            }
   if (T<0){T=-T;}

   // Calculates the matrix M_ab and iM_ab
   for(iJ=0; iJ<sz; iJ++)
      {  
         u[iJ+1] = zij[iJ].r;
         iu[iJ+1]= zij[iJ].i;
      }

   delta = VE.E(j)-VE.E(i);
   if(delta<-0.000001)
   {
      std::cerr << "ERROR module ic1ion - du1calc: energy gain delta gets negative\n"; 
      exit(EXIT_FAILURE);
   }
   if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // Calculates the partition function
   for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(VE.E(iJ)-VE.E(0))/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T)) - exp(-(VE.E(j)-VE.E(0))/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ia|%i>|^2=%6.3f\n |<%i|Ib|%i>|^2=%6.3f\n |<%i|Ic|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
         printf(" |<%i|Id|%i>|^2=%6.3f\n |<%i|Ie|%i>|^2=%6.3f\n |<%i|If|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ia-<Ia>|%i>|^2=%6.3f\n |<%i|Ib-<Ib>|%i>|^2=%6.3f\n |<%i|Ic-<Ic>|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
         printf(" |<%i|Id-<Id>|%i>|^2=%6.3f\n |<%i|Ie-<Ie>|%i>|^2=%6.3f\n |<%i|If-<If>|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
         printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
      }
   }

   // multiply matrix Mab by occupation factor
   for(iJ=0; iJ<sz; iJ++)
      { u[iJ+1] *= sqrt(therm/Z); iu[iJ+1] *= sqrt(therm/Z); }

}

//--------------------------------------------------------------------------------------------------------------
std::vector<double> icmfmat::orbmomdensity_expJ(iceig &VE,int xyz, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   return spindensity_expJ(VE,-xyz, T, matel, save_matrices);
}

//--------------------------------------------------------------------------------------------------------------
std::vector<double> icmfmat::spindensity_expJ(iceig &VE,int xyz, double T, std::vector< std::vector<double> > &matel, bool save_matrices)
{
   double *vt=0, Z=0., U=0.; complexdouble *zt=0, zme;
   std::vector<double> E, ex((_num_op>6?_num_op:6)+2,0.), me, eb; matel.clear();
   int iJ, ind_j, Esz, Hsz=VE.Hsz(), incx=1;
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::expJ() - Hamiltonian matrix size not same as mean field operator!\n"; return E; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   double alpha = 1, beta = 0; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   char uplo = 'U';
   // Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<J[0].nr(); Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0) break; }
   if (T<0){Esz=(int)(-T);printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                          printf ("Number   Excitation Energy\n");
   for (ind_j=0;ind_j<Esz;++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
   } // MR 17.9.2010
   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++)
      if(fabs(VE.zV(ii,jj).r*VE.zV(ii,jj).r+VE.zV(ii,jj).i*VE.zV(ii,jj).i)<DBL_EPSILON*100000)
      {
         VE.zV(ii,jj) = 0.;
      }

   // For first run calculate also the partition function and internal energy
   me.assign(Esz,0.); eb.assign(Esz,0.); Z=0.;
   if(!VE.iscomplex())
   {
      double *fJmat=J[0].f_array(); vt = (double*)malloc(Hsz*sizeof(double));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(ind_j), &incx, &beta, vt, &incx);
         #ifdef _G77
         F77NAME(ddot)(me[ind_j],&Hsz, VE.V(ind_j), &incx, vt, &incx);
         #else
         me[ind_j] = F77NAME(ddot)(&Hsz, VE.V(ind_j), &incx, vt, &incx);
         #endif
 //MR 17.9.2010
      if (T<0)
      {  char instr[MAXNOFCHARINLINE];
         printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
         if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
         eb[ind_j]=strtod(instr,NULL);
      }
       else
      {  eb[ind_j] = exp(-E[ind_j]/(KB*T));} ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 17.9.2010
      }
      free(fJmat); free(vt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }
   else
   {
      complexdouble *zJmat;
      zeroes.zero(J[0].nr(),J[0].nc()); if(iflag[0]==0) zJmat=zmat2f(J[0],zeroes); else zJmat = zmat2f(zeroes,J[0]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      for(ind_j=0; ind_j<Esz; ind_j++)
      {  // Calculates the matrix elements <Vi|J.H|Vi>
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
         #ifdef _G77
         F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
         #else
         zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
         #endif
         me[ind_j] = zme.r;
//       eb[ind_j] = exp(-E[ind_j]/(KB*T)); ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MR 17.9.2010
         if (T<0)
         {  char instr[MAXNOFCHARINLINE];
            printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
            if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
            eb[ind_j]=strtod(instr,NULL);
         }
         else
         {  eb[ind_j] = exp(-E[ind_j]/(KB*T)); } ex[0]+=me[ind_j]*eb[ind_j]; Z+=eb[ind_j]; U+=(E[ind_j]+VE.E(0))*eb[ind_j];
//MRend 17.9.2010                                                                   !!!!    -----------------!!!!
      }
      free(zJmat); free(zt); matel.push_back(me); ex[0]/=Z; U/=Z;
   }

   char nstr[6]; char basename[255]; strcpy(basename,"results/mms/");
   if(save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::expJ(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY))
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::expJ(): Cannot create mms directory\n";
   #endif
   nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }
   int k[] = {0,1, 1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   sMat<double> Upq,Umq;  //if(n>(2*_l+1)) n = 4*_l+2-n;
   char xyzstr[] = "xyz";
   if(xyz>0) { std::cout << "Calculating the expectation values of the spin density operator S" << xyzstr[xyz-1] << "\n"; }
   else      { std::cout << "Calculating the expectation values of the orbital moment density operator L" << xyzstr[-xyz-1] << "\n"; }

   // Rest of the runs only calculate the new matrix elements
   for(iJ=0; iJ<(_num_op>6?_num_op:6); iJ++)
   {  ex[iJ]=0;
      me.assign(Esz,0.);
      // Using the above reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
      if((k[iJ]%2==1) || (k[iJ]>4 && _l==D)) { matel.push_back(me); continue; }
         complexdouble *zJmat; zeroes.zero(J[0].nr(),J[0].nc());
         
         zJmat = balcar_Mq(xyz,k[iJ],q[iJ],_n,_l); // minus sign stands for orbital density coeff
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J.H|Vi>
            F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(ind_j), &incx, &zbeta, zt, &incx);
            #ifdef _G77
            F77NAME(zdotc)(&zme, &Hsz, VE.zV(ind_j), &incx, zt, &incx);
            #else
            zme = F77NAME(zdotc)(&Hsz, VE.zV(ind_j), &incx, zt, &incx);
            #endif
            me[ind_j] = zme.r;
            ex[iJ]+=me[ind_j]*eb[ind_j];
         }
         free(zJmat); free(zt); matel.push_back(me); ex[iJ]/=Z;
      
      if(fabs(ex[iJ])<DBL_EPSILON) ex[iJ]=0.;
   }
   //ex[iJ] = log(Z)-VE.E(0)/(KB*T); ex[iJ+1] = U;
   return ex;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix M_ab=<i|M(q)|j><j|:(q)|i>{exp(-beta_i*T)-exp(-beta_j*T)} for some state i,j
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::dod_u1(int xyz, std::vector<double>&u, std::vector<double>&iu, iceig&VE, double T, int i, int j,int pr,float & delta, bool save_matrices)
{  
// double *vt=0;
   double Z=0., therm; complexdouble *zt=0, zme; zme.r=0; zme.i=0.;
   int sz = (_num_op>6?_num_op:6);
   std::vector<double> mij(sz,0.);//, mji(6,0.);
   std::vector<complexdouble> zij(sz,zme);//, zji(6,zme);
// u.zero(sz); iu.zero(sz);
   int iJ, Hsz=VE.Hsz(), incx=1; 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::u1() - Hamiltonian matrix size not same as mean field operator!\n"; return; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
// double alpha = 1, beta = 0;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   complexdouble *zJmat=0;
   char uplo = 'U';

// char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
// if(save_matrices) {
// #ifndef _WINDOWS
// struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
//    if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::u1(): Can't create mms dir, " << strerror(errno) << "\n";
// #else
// DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
//    if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::u1(): Cannot create mms directory\n";
// #endif
// nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
// strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
// } else { strcpy(basename,"nodir/"); }
   // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
   int k[] = {0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
                 
// sMat<double> Upq,Umq; double redmat; int n = _n; if(n>(2*_l+1)) n = 4*_l+2-n; 

   // Calculates the matrix elements: <i|Ja|j> and <j|Ja|i> for each of the six Ja's
   for(iJ=0; iJ<sz; iJ++)
   {
//    if(k[iJ]%2==1) { if(VE.iscomplex()) { zij[iJ].r=0.; zij[iJ].i=0.; } else mij[iJ]=0.; continue; }
//    if(iJ>=6)
//    {
//       NSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
//       Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[iJ],abs(q[iJ]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
//       MSTR(k[iJ],abs(q[iJ])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
//       Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[iJ],-abs(q[iJ]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
//       #ifdef JIJCONV
//       if(jijconv.size()>1) redmat*=jijconv[iJ];
//       #endif
//       redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[iJ],2*_l,0,0,0);
////     if(q[iJ]<0) { if((q[iJ]%2)==0) Upq -= Umq; else Upq += Umq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Upq += Umq; else Upq -= Umq; } changed MR 15.12.09
//       if(q[iJ]<0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; } else if(q[iJ]>0) { if((q[iJ]%2)==0) Umq += Upq; else Umq -= Upq; }
//       Umq *= redmat;
//    }
//
//    if(!VE.iscomplex() && im[iJ]==0)
//    {
//       vt = (double*)malloc(Hsz*sizeof(double)); 
//       double *fJmat; if(iJ>=6) fJmat=Upq.f_array(); else fJmat=J[iJ].f_array();
//       F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, VE.V(j), &incx, &beta, vt, &incx);
//       #ifdef _G77 
//       F77NAME(ddot)(mij[iJ], &Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
//       #else
//       mij[iJ] = F77NAME(ddot)(&Hsz, VE.V(i), &incx, vt, &incx); zij[iJ].r = mij[iJ];
//       #endif
//       free(fJmat); free(vt);
//    } 
//    else
//    {
//       zeroes.zero(J[0].nr(),J[0].nc());
//       if(iJ>=6) { if(im[iJ]==0) zJmat=zmat2f(Umq,zeroes);   else zJmat = zmat2f(zeroes,Umq); }
//       else      { if(im[iJ]==0) zJmat=zmat2f(J[iJ],zeroes); else zJmat = zmat2f(zeroes,J[iJ]); }
         zJmat = balcar_Mq(xyz,k[iJ],q[iJ],_n,_l);
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, VE.zV(j), &incx, &zbeta, zt, &incx);
         #ifdef _G77 
         F77NAME(zdotc)(&zij[iJ], &Hsz, VE.zV(i), &incx, zt, &incx);
         #else
         zij[iJ] = F77NAME(zdotc)(&Hsz, VE.zV(i), &incx, zt, &incx);
         #endif
//       int k;for(k=0;k<Hsz;++k)printf("%6.3f %+6.3f i  ",VE.zV(j)[k].r,VE.zV(j)[k].i);
         free(zJmat); free(zt);
//    }
   }

   if(i==j&&T>0) {//subtract thermal expectation value from zij=zii
            std::vector< std::vector<double> > matel;
            std::vector<double> vJ = spindensity_expJ(VE,xyz,T,matel,save_matrices);
            for(iJ=0; iJ<sz; iJ++)zij[iJ].r-=vJ[iJ];
            }
   if (T<0){T=-T;}

   // Calculates the matrix M_ab and iM_ab
   for(iJ=0; iJ<sz; iJ++)
      {  
         u[iJ+1] = zij[iJ].r;
         iu[iJ+1]= zij[iJ].i;
      }

   delta = VE.E(j)-VE.E(i);
   if(delta<-0.000001)
   {
      std::cerr << "ERROR module ic1ion - du1calc: energy gain delta gets negative\n"; 
      exit(EXIT_FAILURE);
   }
   if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // Calculates the partition function
   for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(VE.E(iJ)-VE.E(0))/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T)) - exp(-(VE.E(j)-VE.E(0))/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ia|%i>|^2=%6.3f\n |<%i|Ib|%i>|^2=%6.3f\n |<%i|Ic|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
         printf(" |<%i|Id|%i>|^2=%6.3f\n |<%i|Ie|%i>|^2=%6.3f\n |<%i|If|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         printf(" |<%i|Ia-<Ia>|%i>|^2=%6.3f\n |<%i|Ib-<Ib>|%i>|^2=%6.3f\n |<%i|Ic-<Ic>|%i>|^2=%6.3f\n",i+1,j+1,u[1]*u[1]+iu[1]*iu[1],i+1,j+1,u[2]*u[2]+iu[2]*iu[2],i+1,j+1,u[3]*u[3]+iu[3]*iu[3]);
         printf(" |<%i|Id-<Id>|%i>|^2=%6.3f\n |<%i|Ie-<Ie>|%i>|^2=%6.3f\n |<%i|If-<If>|%i>|^2=%6.3f\n",i+1,j+1,u[4]*u[4]+iu[4]*iu[4],i+1,j+1,u[5]*u[5]+iu[5]*iu[5],i+1,j+1,u[6]*u[6]+iu[6]*iu[6]);
         printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
      }
   }

   // multiply matrix Mab by occupation factor
   for(iJ=0; iJ<sz; iJ++)
      { u[iJ+1] *= sqrt(therm/Z); iu[iJ+1] *= sqrt(therm/Z); }

}
*/

} // namespace libMcPhase
