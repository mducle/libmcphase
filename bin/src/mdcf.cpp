 // *************************************************************************
 // ****** md cf - single ion transition matrix elements storage ************
 // *************************************************************************
// methods for class mdcf 
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "mdcf.hpp"

#define MAXNOFSPINS  200

// returns sum of all components of intvector
int sum(IntVector & v)
 {int i,sum=0;
 for (i=v.Lo();i<=v.Hi();++i)
 {sum+=v(i);}
  return sum;
}

// returns md of cf [i=na,j=nb,k=nc] 
ComplexMatrix & mdcf::U(int na, int nb, int nc) const
{ return (*s[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::V(int na, int nb, int nc) const
{ return (*sb[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::M(int na, int nb, int nc)
{ return (*m[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::N(int na, int nb, int nc)
{ return (*mb[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::sqrt_gamma(int na, int nb, int nc) const
{ return (*l[in(na,nb,nc)]);
}
ComplexMatrix & mdcf::sqrt_Gamma(int na, int nb, int nc) const
{ return (*lb[in(na,nb,nc)]);
}
Vector & mdcf::delta(int na, int nb, int nc)
{ return (*d[in(na,nb,nc)]);
}
// the same but for cf number "i"
ComplexMatrix & mdcf::Ui(int i)
{ return (*s[i]);
}
ComplexMatrix & mdcf::Vi(int i)
{ return (*sb[i]);
}
ComplexMatrix & mdcf::Mi(int i)
{ return (*m[i]);
}
ComplexMatrix & mdcf::Ni(int i)
{ return (*mb[i]);
}
ComplexMatrix & mdcf::sqrt_gammai(int i)
{ return (*l[i]);
}
ComplexMatrix & mdcf::sqrt_Gammai(int i)
{ return (*lb[i]);
}
Vector  & mdcf::deltai(int i)
{ return (*d[i]);
}
// get index ijk=iv(1-3)  of cf configuration number in 
int * mdcf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of cf from indizes i,j,k
int mdcf::in(int i, int j, int k) const
{return ((i*mxb+j)*mxc+k);}

// get number of cf from indizes i,j,k,l
int mdcf::ind(int i, int j, int k, int l)
{int indd=(((i*mxb+j)*mxc+k)*nofatoms+l);
 if(indd<0||indd>mxa*mxb*mxc*(nofatoms+1)) {fprintf(stderr,"mdcf indexing error");exit(EXIT_FAILURE);}
 return indd;}

// return number of cfs
int mdcf::n()
{return (nofa*nofb*nofc);
}
int mdcf::na()
{return nofa;
}
int mdcf::nb()
{return nofb;
}
int mdcf::nc()
{return nofc;
}

/**************************************************************************/

//constructors
mdcf::mdcf (int n1,int n2,int n3,int n,int nc)
{  int i;
   nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   nofatoms=n;nofcomponents=nc;

//dimension arrays
  s = new ComplexMatrix * [mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (m == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (l == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  sb = new ComplexMatrix * [mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (sb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  mb = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (mb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb = new ComplexMatrix * [mxa*mxb*mxc+1];
  if (lb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector * [mxa*mxb*mxc+1]; //(1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  nt= new IntVector * [mxa*mxb*mxc+1];
  if (nt == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  for(i=0;i<=mxa*mxb*mxc;++i){nt[i]=new IntVector(1,nofatoms);}

  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*(nofatoms+1)+1];   
  if (eigenstates == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  

}

ComplexMatrix & mdcf::est(int i, int j, int k, int l)
{return (*eigenstates[ind(i,j,k,l)]);}


void mdcf::est_ini(int i, int j, int k, int l,ComplexMatrix & M) // initialize est
{eigenstates[ind(i,j,k,l)]=new ComplexMatrix(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
 (*eigenstates[ind(i,j,k,l)])=M;
}

// has to be called before mdcf object can be used for calculation
void mdcf::set_noftransitions(int i, int j, int k, IntVector & notr)
{      
      (*nt[in(i,j,k)])=notr;
     s[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     m[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     l[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     sb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     mb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     lb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
     d[in(i,j,k)]= new Vector(1,sum((*nt[in(i,j,k)])));
      
}

int mdcf::baseindex(int i, int j, int k, int l, int tn) const
{// the baseindex is used to number the rows an columns of the
 // matrices associated with the crystallographic unit number ijk
 // it starts at one and combines indices l (atom number) and t (number of
 // transition for this atom) into a single index bi, its maximum value is
 // baseindex_max
  int bi=0;
  int i1;
for(i1=1;i1<=l;++i1)
{bi+=(*nt[in(i,j,k)])(i1);}
bi-=(*nt[in(i,j,k)])(l);
bi+=tn;
return bi;
}

//void mdcf::baseindex2ltn(int baseindex, int i, int j, int k, int & l, int & tn)
//{// the inverse of function baseindex -> returns l and tn for given baseindex, ijk
//  int bi=0;
//  int i1;
//  for(i1=1;i1<=(*nt[in(i,j,k)]).Hi;++i1)
//  {if(bi+(*nt[in(i,j,k)])(i1);}
//   bi-=(*nt[in(i,j,k)])(l);
//   bi+=tn;
//?????
//
//}

int mdcf::baseindex_max(int i, int j, int k) const
{return sum((*nt[in(i,j,k)]));}

int mdcf::noft(int i, int j, int k,int l) const
{return (*nt[in(i,j,k)])(l);}

//destruktor
mdcf::~mdcf ()
{int i,j,k;
 for (i=1;i<=nofa;++i){
 for (j=1;j<=nofb;++j){
 for (k=1;k<=nofc;++k){
 delete s[in(i,j,k)];
 delete m[in(i,j,k)];
 delete l[in(i,j,k)];
 delete sb[in(i,j,k)];
 delete mb[in(i,j,k)];
 delete lb[in(i,j,k)];
 delete d[in(i,j,k)];
 }}}
 delete []s;delete []m;delete []d;delete []l;delete []nt;
 delete []sb;delete []mb;delete []lb;
}



//*
//kopier-konstruktor
mdcf::mdcf (const mdcf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  nofatoms=p.nofatoms;nofcomponents=p.nofcomponents;
//dimension arrays
  s = new ComplexMatrix*[mxa*mxb*mxc+1]; //(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix*[mxa*mxb*mxc+1];
  if (m == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix*[mxa*mxb*mxc+1];
  if (l == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  sb = new ComplexMatrix *[mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (sb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  mb = new ComplexMatrix *[mxa*mxb*mxc+1];
  if (mb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  lb = new ComplexMatrix *[mxa*mxb*mxc+1];
  if (lb == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector*[mxa*mxb*mxc+1];//(1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  nt = new IntVector*[mxa*mxb*mxc+1];//(1,nofatoms);
  if (nt == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  eigenstates= new ComplexMatrix * [mxa*mxb*mxc*nofatoms+1];   
  if (eigenstates == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 


 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {
      nt[in(i,j,k)]=p.nt[in(i,j,k)];

      s[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      m[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      l[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      sb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      mb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      lb[in(i,j,k)]= new ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      d[in(i,j,k)]= new Vector(1,sum((*nt[in(i,j,k)])),1,sum((*nt[in(i,j,k)])));

      *d[in(i,j,k)]=*p.d[in(i,j,k)];
      *s[in(i,j,k)]=*p.s[in(i,j,k)];
      *m[in(i,j,k)]=*p.m[in(i,j,k)];
      *l[in(i,j,k)]=*p.l[in(i,j,k)];
      *sb[in(i,j,k)]=*p.sb[in(i,j,k)];
      *mb[in(i,j,k)]=*p.mb[in(i,j,k)];
      *lb[in(i,j,k)]=*p.lb[in(i,j,k)];
     } 
    }
  }           

}
//*/
/*
//zuweisung
mdcf & mdcf::operator= (const mdcf & op2)
{int i,j,k;
 nofa=op2.nofa; nofb=op2.nofb; nofc=op2.nofc;
 mxa=op2.mxa; mxb=op2.mxb; mxc=op2.mxc;
 nofatoms=op2.nofatoms;
 nofcomponents=op2.nofcomponents;
  
  delete[]s;delete []m;delete []d;delete []l;delete[]nt;
//dimension arrays
  s = new ComplexMatrix[mxa*mxb*mxc+1];//(1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  m = new ComplexMatrix[mxa*mxb*mxc+1];
  if (m == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  l = new ComplexMatrix[mxa*mxb*mxc+1];
  if (l == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  d = new Vector[mxa*mxb*mxc+1];//(1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  nt = new IntVector[mxa*mxb*mxc+1](1,nofatoms);
  if (nt == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {
      (*nt[in(i,j,k)])=op2.(*nt[in(i,j,k)]);
      
      s[in(i,j,k)]= ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      m[in(i,j,k)]= ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      l[in(i,j,k)]= ComplexMatrix(1,nofcomponents*sum((*nt[in(i,j,k)])),1,nofcomponents*sum((*nt[in(i,j,k)])));
      d[in(i,j,k)]= Vector(1,sum((*nt[in(i,j,k)])),1,sum((*nt[in(i,j,k)])));
           
      d[in(i,j,k)]=op2.d[in(i,j,k)];

      s[in(i,j,k)]=op2.s[in(i,j,k)];
      m[in(i,j,k)]=op2.m[in(i,j,k)];
      l[in(i,j,k)]=op2.l[in(i,j,k)];
     } 
    }
  }           
  return *this;
}

*/