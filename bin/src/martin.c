// martins functions
// to get on better in c


#include<martin.h>
#include<myev.h>

#ifndef __linux__
#include<cfloat>
#endif


// ********************************************************************************************************
// extract parameter 'parameter'  from string instr (z.B. "blabla dmin=0.2 blabla") -
// output: var ... value of parameter
// returns 1 on error and 0 if successful
int extract(char * instr,const char * parameter,double & var)
{ //const char delimiters[] = " =:\n";
  char *token,*td,*te;
  
// check if line is comment line -> if yes return 1
if (instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!') return 1; //removed 26.5.02 in order to be able to place parameters in comment lines
                                 // inserted again 27.8.09 to be able to have real comment lines ignored
                                 // by mcphase - however "#!" will be treated as comment with variable to be read
 
  // strip /r (dos line feed) from line if necessary
 // while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
 
  if ((token = strstr (instr, parameter))==NULL) return 1; // parameter string not found
  if (token>instr&&1!=strspn(token-1," \t!?$%&*()[]{}\"��/><;@:=+-~|"))return 1; // no space etc before parameter
 td=instr;while ((te=strstr(td,"#!"))!=NULL){td=te+1;} // skip all "#!" signs and
 if ((td=strchr(td,'#'))!=NULL){if(td<token) return 1;} // check if comment sign "#" appears before parameter - if yes return 1
 
  //extract parameter  
  token+=strlen(parameter);
  if (strstr (token, "=")==NULL) return 1;  // no '=' found after parameter string
  while(strstr(token," ")==token||strstr(token,"\t")==token)++token;
  if (strstr(token,"=")!=token) return 1; // there are other characters than tab or spaces between parameter and =
  ++token;
  var = strtod (token, NULL);
  return 0;
}

// same for int and double 
// extract a variable named [parmeter] into var out of a string [instr]
  int extract(char * instr,const char * parameter,int & var)
      {double dd;if(0==extract(instr,parameter,dd)){var=(int)dd;return 0;}else{return 1;}}
  int extract(char * instr,const char * parameter,float & var)
      {double dd;if(0==extract(instr,parameter,dd)){var=(float)dd;return 0;}else{return 1;}}


//the same for a string ... maximal n characters are copied
int extract(char * instr,const char * parameter,char * var,size_t n)
{ const char delimiters[] = " \n";
  char *token,*td,*te;
  
// check if line is comment line -> if yes return 1
if (instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!') return 1; //removed 26.5.02 in order to be able to place parameters in comment lines
                                 // inserted again 27.8.09 to be able to have real comment lines ignored
                                 // by mcphase - however "#!" will be treated as comment with variable to be read
 
// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

  if ((token = strstr (instr, parameter))==NULL) return 1; // check if parameter is found - if not return 1
  if (token>instr&&1!=strspn(token-1," \t!?$%&*()[]{}\"��/><;@:=+-~|")){return 1;} // no space etc before parameter
 td=instr;while ((te=strstr(td,"#!"))!=NULL){td=te+1;} // skip all "#!" signs and
 if ((td=strchr(td,'#'))!=NULL){if(td<token) return 1;} // check if comment sign "#" appears before parameter - if yes return 1

  //extract parameter  
  token+=strlen(parameter);
  token+=strspn(token, " \t=");// look for '=' but continue if not present
                                // remove starting spaces
  strncpy (var,token, n);
  //remove from string var all characters after delimiters
  strtok(var,delimiters);
  return 0;
}

// same for variable with prefix
int extract_with_prefix(char * instr,char * prefix, const char * parameter,double & var)
{if(0==extract(instr,parameter,var)){return 0;}
 else{char par[MAXNOFCHARINLINE];sprintf(par,"%s%s",prefix,parameter);  if(0==extract(instr,par,var)){return 0;}}
 return 1;
}
// same for int variable with prefix
int extract_with_prefix(char * instr,char * prefix, const char * parameter,int & var)
{if(0==extract(instr,parameter,var)){return 0;}
 else{char par[MAXNOFCHARINLINE];sprintf(par,"%s%s",prefix,parameter);  if(0==extract(instr,par,var)){return 0;}}
 return 1;
}
// same for char variable with prefix
int extract_with_prefix(char * instr,char * prefix, const char * parameter,char * var,size_t n)
{if(0==extract(instr,parameter,var,n)){return 0;}
 else{char par[MAXNOFCHARINLINE];sprintf(par,"%s%s",prefix,parameter);  if(0==extract(instr,par,var,n)){return 0;}}
 return 1;
}

//****************************************************************************************************************
//open a file: similar fopen but with error check 
FILE * fopen_errchk (const char * filename,const char * mode)
{ FILE *file;
 
  errno = 0;
  file = fopen (filename,mode);
  if (file == NULL)
    { 
      fprintf (stderr, "Couldn't open file %s: %s\n",filename, strerror (errno));
      exit (EXIT_FAILURE);
     }
  return file;
} 

//input file from string with error check
char * fgets_errchk (char * instr,int size, FILE * file)
{char * s;
 errno=0;
 s=fgets(instr,size,file);
 if (s==NULL)
 {fprintf (stderr,"ERROR fgets_errchk: end of file\n"); 
  exit(EXIT_FAILURE);}
 return s;
}

// **************************************************************************
// input line and split into numbers
// on input: nn[0]... size of array nn[] (not changed in routine)
//           fin_coq..filepointer
// returns nof numbers read and
// 0 if no numbers have been read (e.g. empty line or line is comment line starting with # or EOF)

#define maxnofcharinline  7024


char * mystrtok (char * s, char * delimiters)
{char * pointer;
 pointer=s+strcspn(s,delimiters);
 pointer+=strspn(pointer,delimiters);
 if ((*pointer)=='\0')pointer=NULL;
return pointer;
}

// function to input a line of numbers separated by delimiters
// example:
// 3 23 542 23
// returns:0 .... it is a comment line (starting with #) or empty line
//         n .... number of numbers read
int inputline (FILE * fin_coq, float *nn)
{
  char instr[maxnofcharinline];
  char delimiters[] = " \n\t";
  char *token;
  int i;
  errno=0;

  if (fgets (instr, sizeof (instr), fin_coq) == NULL)
    { return 0;}

 
  if (feof(fin_coq)==0&&strchr(instr,'\n')==NULL)
    { fprintf (stderr, "Error in function inputline: input string too long");
      exit (EXIT_FAILURE);
     }

// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

// check if line is comment line -> if yes return 0
  if (instr[strspn(instr," \t")]=='#') return 0;
 

  //initialize token to first nonspace character
  token = instr+strspn(instr," \t");
  i=0;
  if (*token!='\n'){
  for (i = 1; token != NULL&&(*token!='#') ; ++i)
    {
if(i>=(int)nn[0])
        { fprintf (stderr, "Error in function inputline: maximum value of numbers in line exceeded,more numbers in line (>%i) to be read.\n",i);
           exit (EXIT_FAILURE);
        }
    
      nn[i] = strtod (token, NULL);
      token = mystrtok (token, delimiters);
 // printf("i=%i token=%g\n",i,nn[i]);
    }
  }        
//  printf("i=%i token=%g\n",i,nn[i]);

  if (i<1) {return 0;}
  return i-1;
}

// function to input a line of numbers separated by delimiters
// example:
// parname=3 23 542 23
// returns:0 .... it is a comment line (starting with #) or parameter parname not found
//         n .... number of numbers read
int inputparline (const char * parname, FILE * fin_coq, float *nn)
{
  char instr[maxnofcharinline];
  char delimiters[] = " \n";
  char *token;
  int i;
  errno=0;

  if (fgets (instr, sizeof (instr), fin_coq) == NULL)
    { return 0;}
  
  if (feof(fin_coq)==0&&strchr(instr,'\n')==NULL)
    { fprintf (stderr, "Error in function inputparline: input string too long");
      exit (EXIT_FAILURE);
     }
// check if line is comment line -> if yes return 0
  if (instr[strspn(instr," \t")]=='#') return 0;

// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
 
// check line contains parname -> if no return 0
  if ((token = strstr (instr, parname))==NULL) return 0;
  if (strstr(instr,"#")!=NULL){if(token>strstr(instr,"#")) return 0;}
  token+=strlen(parname);
// check if there is a "=" character -> if no return 0
  if((token = strstr (token, "="))==NULL) return 0;
  token +=1;
  while (*token==' '){++token;} // remove starting spaces
  //+ initialize token to first nonspace character
  i=0;
  if (*token!='\n'){
  for (i = 1; token != NULL ; ++i)
    {
if(i>=(int)nn[0])
        { fprintf (stderr, "Error in function inputparline: maximum value of numbers in line exceeded,more numbers in line (>%i) to be read.\n",i);
           exit (EXIT_FAILURE);
        }
    
      nn[i] = strtod (token, NULL);
      token = mystrtok (token, delimiters);
 // printf("i=%i token=%g\n",i,nn[i]);
    }
  }        
//  printf("i=%i token=%g\n",i,nn[i]);

  if (i<1) {return 0;}
  return i-1;
}

 // *************************************************************************

// return random number between 0 and z
float rnd(float z){return  z*rand()/RAND_MAX;}

// return integer of floating number (like basic integer function)
float integer (float s){  double result;modf(s,&result);  return result;}

//return rounded integer of floating number
int cint (float s){  double result;if(modf(s,&result)<0.5){return (int)result;}{return (int)result+1;}}

// round to 1e-11 precision
double myround(double s){return myround(1e-10,s);}
double myround(double prec,double s){
double result;result=prec*rint(s/prec);
if (result==0){result*=result;}
return result;
}

// factorial of an integer number
double factorial(double number) {
	double temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}
// Above recursive function too slow. For pure int return values, use a look
// up table for up to 12! (as much as 32-bit in can hold) inlined in <martin.h>

// return threej symbol 
float threej (float AJ1,float  AJ2,float  AJ3,float AM1,float AM2,float AM3)
{//  THIS IS A PROGRAM TO CALCULATE THE 3_J SYMBOL
//    (AJ1,AJ2,AJ3;AM1,AM2,AM3)
//    F(I)=I! IS THE FACTORIAL FUNCTION
//    I = 0 TO 30 IS USUSALLY ENOUGH
//    PLEASE NOTE THAT I STARTS FROM ZERO
//
//    "CRYSTAL FIELD HANDBOOK", EDITED BY D.J. NEWMAN AND BETTY NG
//    CAMBRIDGE UNIVERSITY PRESS 
//    ... and transfered to c and improved by M. Rotter 2008

        int I,J7;
	double R0,R4,R5,R6,R8,R9;
	double F[31];      

//  DEFINE THE FACTORIAL FUNCTION
        F[0] = 1;for(I=1;I<=30;++I){F[I]=F[I-1]*I;}

//   SELECTION RULES
     if (AM1 + AM2 + AM3 != 0) { return 0;}
     if (AJ1 + AJ2 - AJ3 < 0) { return 0;}
     if (AJ3 + AJ1 - AJ2 < 0) { return 0;}
     if (AJ3 + AJ2 - AJ1 < 0) { return 0;}
     if (AJ1 + AJ2 + AJ3 + 1 < 0) { return 0;}
     if (AJ1 - AM1 < 0) { return 0;}
     if (AJ2 - AM2 < 0) { return 0;}
     if (AJ3 - AM3 < 0) { return 0;}
     if (AJ1 + AM1 < 0) { return 0;}
     if (AJ2 + AM2 < 0) { return 0;}
     if (AJ3 + AM3 < 0) { return 0;}

   R4 = F[int(AJ1 + AJ2 - AJ3)] * F[int(AJ1 - AM1)] * F[int(AJ2 - AM2)] * F[int(AJ3 - AM3)] * F[int(AJ3 + AM3)];
   R5 = F[int(AJ1 + AJ2 + AJ3 + 1)] * F[int(AJ3 + AJ1 - AJ2)] * F[int(AJ3 + AJ2 - AJ1)] * F[int(AJ1 + AM1)] * F[int(AJ2 + AM2)];
   R6 = 0;
   for (J7 = 0;J7<=25;++J7)
     { if (AJ1 + AM1 + J7 >= 0 && AJ2 + AJ3 - AM1 - J7 >= 0 && AJ3 + AM3 - J7 >= 0 && AJ1 - AM1 - J7 >= 0
           && AJ2 - AJ3 + AM1 + J7 >= 0) 
       {R8 = F[int(AJ1 + AM1 + J7)] * F[int(AJ2 + AJ3 - AM1 - J7)];
        R8 *= odd(int(AJ1 - AM1 - J7)) ? -1 : 1; // .... equiv to  * (-1) ^ (int(AJ1 - AM1 - J7));
        R9 = F[J7] * F[int(AJ3 + AM3 - J7)] * F[int(AJ1 - AM1 - J7)] * F[int(AJ2 - AJ3 + AM1 + J7)];
        R6 +=  R8 / R9;
       }    
     }
     R0 = sqrt(R4 / R5) * R6 ;
     R0 *= odd(int(AJ1 - AJ2 - AM3)) ? -1 : 1; //.... equiv to * (-1) ^ int((AJ1 - AJ2 - AM3))
  return R0;
}



#ifndef __linux__
// /* not needed any more in mingw3.1
// return rounded value
//double rint(double value)
//{if (fabs(floor(value)-value)<fabs(ceil(value)-value))
//  {return floor(value);}
// else
//  {return ceil(value);}
//}




// MaCOS,linux etc check ... compare http://sourceforge.net/p/predef/wiki/OperatingSystems/
//sleep function
int sleep(int a)
{
#if defined CYGWIN || defined __CYGWIN__ || (__APPLE__ & __MACH__)
#include<unistd.h>
  return usleep(a*1000);
#else
  _sleep(a*1000);
  return true;
#endif
}

#endif

// some vector functions

/**********************************************************************/
void xproduct(Vector & result,Vector a, Vector b)
{
 result(1)=a(2)*b(3)-a(3)*b(2);
 result(2)=-a(1)*b(3)+a(3)*b(1);
 result(3)=a(1)*b(2)-a(2)*b(1);

 return ;
}

// calculate transition matrix element of eigenvectr <i|op|j> with eigenvectors
// given as column vectors i and j of Matrix (zr + i zi), the Hermitian operator op is described
// by   matrix op.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
double matelr (int k,int l,Matrix & zr, Matrix & zi, Matrix & op)
{double sumr=0;//,p1r,p1i,p2r,p2i;
  //xisyj,xjsyi,im(0,1);
  int dl,dh;
 if((dl=op.Rlo())!=op.Clo()||op.Rlo()!=zr.Rlo()||op.Rlo()!=zi.Rlo()||
    (dh=op.Rhi())!=op.Chi()||op.Rhi()!=zr.Rhi()||op.Rhi()!=zi.Rhi()){fprintf(stderr,"Error matel: operator not sqaure and/or eigenvector dimension do not match\n");exit(EXIT_FAILURE);}
 for(int i=dl;i<=dh;++i){//complex <double> xis(zr[i][k],-zi[i][k]);
                         //complex <double> yi(zr[i][l],zi[i][l]);
                       sumr+=(zr[i][k]*zr[i][l]+ zi[i][k]*zi[i][l])*op[i][i]; //xis*yi*op[i][i];
  for(int j=i+1;j<=dh;++j){//complex <double> yj(zr[j][l],zi[j][l]);
                           //complex <double> opc1(op[j][i],-op[i][j]);                          
                                       // sum+=opc1*xis*yj;
                           //complex <double> xjs(zr[j][k],-zi[j][k]);
                           //complex <double> opc2(op[j][i],op[i][j]);
                                        //sum+=opc2*xjs*yi;
                          if(op[i][j]!=0){
                             //p1i=zr[i][k]*zi[j][l]- zi[i][k]*zr[j][l];
                             //p2i=zr[j][k]*zi[i][l]- zi[j][k]*zr[i][l];
                           sumr+=(zr[i][k]*zi[j][l]- zi[i][k]*zr[j][l]-zr[j][k]*zi[i][l]+ zi[j][k]*zr[i][l])*op[i][j];
                                          }
                          if(op[j][i]!=0){
                            // p1r=zr[i][k]*zr[j][l]+ zi[i][k]*zi[j][l];// xis*yj
                            // p2r=zr[j][k]*zr[i][l]+ zi[j][k]*zi[i][l];// xjs*yi
                           sumr+=(zr[i][k]*zr[j][l]+ zi[i][k]*zi[j][l]+zr[j][k]*zr[i][l]+ zi[j][k]*zi[i][l])*op[j][i];
                                          }
                                                         
                         }
        }
 return sumr;
}
double mateli (int k,int l,Matrix & zr, Matrix & zi, Matrix & op)
{double sumi=0;//,p1r,p1i,p2r,p2i;
  //xisyj,xjsyi,im(0,1);
  int dl,dh;
 if((dl=op.Rlo())!=op.Clo()||op.Rlo()!=zr.Rlo()||op.Rlo()!=zi.Rlo()||
    (dh=op.Rhi())!=op.Chi()||op.Rhi()!=zr.Rhi()||op.Rhi()!=zi.Rhi()){fprintf(stderr,"Error matel: operator not sqaure and/or eigenvector dimension do not match\n");exit(EXIT_FAILURE);}
 for(int i=dl;i<=dh;++i){//complex <double> xis(zr[i][k],-zi[i][k]);
                         //complex <double> yi(zr[i][l],zi[i][l]);
                       sumi+=(zr[i][k]*zi[i][l]- zi[i][k]*zr[i][l])*op[i][i];
  for(int j=i+1;j<=dh;++j){//complex <double> yj(zr[j][l],zi[j][l]);
                           //complex <double> opc1(op[j][i],-op[i][j]);                          
                                       // sum+=opc1*xis*yj;
                           //complex <double> xjs(zr[j][k],-zi[j][k]);
                           //complex <double> opc2(op[j][i],op[i][j]);
                                        //sum+=opc2*xjs*yi;
                          if(op[i][j]!=0){
                             //p1r=zr[i][k]*zr[j][l]+ zi[i][k]*zi[j][l];// xis*yj
                             //p2r=zr[j][k]*zr[i][l]+ zi[j][k]*zi[i][l];// xjs*yi
                           sumi+=(zr[j][k]*zr[i][l]+ zi[j][k]*zi[i][l]-zr[i][k]*zr[j][l]- zi[i][k]*zi[j][l])*op[i][j];
                                          }
                          if(op[j][i]!=0){
                             //p1i=zr[i][k]*zi[j][l]- zi[i][k]*zr[j][l];
                             //p2i=zr[j][k]*zi[i][l]- zi[j][k]*zr[i][l];
                           sumi+=(zr[i][k]*zi[j][l]- zi[i][k]*zr[j][l]+zr[j][k]*zi[i][l]- zi[j][k]*zr[i][l])*op[j][i];                           
                                          }
                                                         
                         }
        }
 return sumi;
}

// calculate some column i of ComplexMatrix opZ which is defined as the 
// product of two complex matrices op * Z, however Z is given as
// as real matrices  Matrix zr + i Matrix zi, the Hermitian operator op is described
// by   Matrix op.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].
void opZcol (int i,ComplexMatrix & opZ, Matrix & op,Matrix & zr, Matrix & zi)
{int k,l,d=op.Rhi(); // no check of dimensions to improve speed !!
 for(k=1;k<=d;++k){opZ[k][i]=complex<double>(op[k][k]*zr[k][i],op[k][k]*zi[k][i]);
                   for(l=1;l<k;++l)    opZ[k][i]+=complex<double>(op[l][k]*zr[l][i]-op[k][l]*zi[l][i],op[l][k]*zi[l][i]+op[k][l]*zr[l][i]);
                   for(l=k+1;l<=d;++l) opZ[k][i]+=complex<double>(op[k][l]*zr[l][i]+op[l][k]*zi[l][i],op[k][l]*zi[l][i]-op[l][k]*zr[l][i]);
 
                  } 
}

void rezcalc(Vector r1,Vector  r2,Vector  r3,Vector  rez1,Vector  rez2,Vector  rez3)
{// calculate reciprocal lattice rezi from real lattice ri
 float vol;
 xproduct(rez1,r2,r3); vol=rez1*r1; rez1*=2.0*PI/vol;
 xproduct(rez2,r1,r3); vol=rez2*r2; rez2*=2.0*PI/vol;
 xproduct(rez3,r1,r2); vol=rez3*r3; rez3*=2.0*PI/vol;
return;}

void get_abc_in_ijk(Matrix & abc_in_ijk,Vector & abc)
{// get lattice vectors in terms of coordinate system ijk
 //defined by  j||b, k||(a x b) and i normal to k and j
double alpha=90; if(abc.Hi()==6) alpha=abc(4); // this is to ensure backward compatibility
double beta=90; if(abc.Hi()==6)  beta=abc(5);
double gamma=90; if(abc.Hi()==6) gamma=abc(6);

abc_in_ijk(1,1)=abc(1)*sin(gamma*PI/180);  // vector a in terms of ijk
abc_in_ijk(2,1)=abc(1)*cos(gamma*PI/180);
abc_in_ijk(3,1)=0;

abc_in_ijk(1,2)=0;   // vector b in terms of ijk
abc_in_ijk(2,2)=abc(2);
abc_in_ijk(3,2)=0;

abc_in_ijk(2,3)=abc(3)*cos(alpha*PI/180);  //vector c in terms of ijk
abc_in_ijk(1,3)=abc(3)*(cos(beta*PI/180)-cos(alpha*PI/180)*cos(gamma*PI/180))/sin(gamma*PI/180);
if (fabs(abc_in_ijk(1,3))>abc(3)){fprintf(stderr,"ERROR getabc_in_ijk: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
abc_in_ijk(3,3)=abc(3)*abc(3)-abc_in_ijk(2,3)*abc_in_ijk(2,3)-abc_in_ijk(1,3)*abc_in_ijk(1,3);
if (abc_in_ijk(3,3)<=0){fprintf(stderr,"ERROR getabc_in_ijk: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
abc_in_ijk(3,3)=sqrt(abc_in_ijk(3,3));
}

void dadbdc2ijk(Matrix & rijk,Matrix & rdadbdc, Vector & abc)
{// transforms primitive lattice vector matrix r given in terms of abc
 // to ijk coordinate system
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system
                           //defined by  j||b, k||(a x b) and i normal to k and j
get_abc_in_ijk(rtoijk,abc);

rijk=rtoijk*rdadbdc;
}

void dadbdc2ijk(Vector & rijk,Vector & dadbdc, Vector & abc)
{// transforms vector r given in terms of abc
 // to ijk coordinate system
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system
                           //defined by  j||b, k||(a x b) and i normal to k and j
get_abc_in_ijk(rtoijk,abc);

rijk=rtoijk*dadbdc;
}

void ijk2dadbdc(Vector & dadbdc,Vector & rijk, Vector & abc)
{// transforms vector r given in terms of ijk coordinates
// to da db dc
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system
                           //defined by  j||b, k||(a x b) and i normal to k and j
get_abc_in_ijk(rtoijk,abc);
dadbdc=rtoijk.Inverse()*rijk;
}

void hkl2ijk(Vector & qijk,Vector & hkl, Vector & abc)
{// transforms Miller indices (in terms of reciprocal lattice abc*)
// to Q vector in ijk coordinate system
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system
                           //defined by  j||b, k||(a x b) and i normal to k and j
get_abc_in_ijk(rtoijk,abc);
qijk=rtoijk.Inverse().Transpose()*hkl;
qijk*=2*PI;
}

void ijk2hkl(Vector & hkl, Vector & qijk,Vector & abc)
{// transforms Q vector in ijk coordinate system to Miller indices (in terms of reciprocal lattice abc*)
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system
                           //defined by  j||b, k||(a x b) and i normal to k and j
get_abc_in_ijk(rtoijk,abc);
hkl=rtoijk.Transpose()*qijk;
hkl/=2*PI;
}


void nlimits_calc(Vector & nmin, Vector & nmax, double radius, Matrix & a)
{// problem: we want to find all lattice vectors Rn=ni*ai which are within a
 // sphere of radius r from the origin (ai = column vectors of matrix a)
 // this routine returns the maximum and minimum values of ni i=1,2,3
 // by probing the corners of a cube
 Vector ddd(1,8),dd0(1,3),dd(1,3),minv(1,3),maxv(1,3);
  int i;
 minv(1)=-radius;
 minv(2)=-radius;
 minv(3)=-radius;
 maxv(1)=radius;
 maxv(2)=radius;
 maxv(3)=radius;
 // determine ijkmin ijkmax by calculating the 8 corners of the  quader
  // in terms of primitive lattice
  // i*p.Column(1)+j*p.Column(2)+k*p.Column(3)=cornerpointvector ... i,j,k =?
  // ijk=p^-1*cornerpointvector
  for (i=1;i<=3;++i)
  {dd0=minv;               dd=a.Inverse()*dd0;ddd(1)=dd(i);
   dd0=minv;dd0(1)=maxv(1);dd=a.Inverse()*dd0;ddd(2)=dd(i);
   dd0=minv;dd0(2)=maxv(2);dd=a.Inverse()*dd0;ddd(3)=dd(i);
   dd0=minv;dd0(3)=maxv(3);dd=a.Inverse()*dd0;ddd(4)=dd(i);
   dd0=maxv;               dd=a.Inverse()*dd0;ddd(5)=dd(i);
   dd0=maxv;dd0(1)=minv(1);dd=a.Inverse()*dd0;ddd(6)=dd(i);
   dd0=maxv;dd0(2)=minv(2);dd=a.Inverse()*dd0;ddd(7)=dd(i);
   dd0=maxv;dd0(3)=minv(3);dd=a.Inverse()*dd0;ddd(8)=dd(i);
   nmin(i)=Min(ddd)-1;nmax(i)=Max(ddd)+1;
  }
}

// some matrix functions for hermitian matrices in
// real notation: The real parts of the elements must be
 //  stored in the lower triangle of z,the imaginary parts (of the elements
 //  corresponding to the lower triangle) in the positions
 //  of the upper triangle of z[lo..hi,lo..hi].
Matrix herm_dirprod(Matrix  R, Matrix  T) // direct product
{if(R.Rhi()!=R.Chi()){fprintf(stderr,"Error martin.c herm_dirproduct: Matrix A not square\n");exit(EXIT_FAILURE);}
 if(T.Rhi()!=T.Chi()){fprintf(stderr,"Error martin.c herm_dirproduct: Matrix B not square\n");exit(EXIT_FAILURE);}
//myPrintMatrix(stdout,R);
//myPrintMatrix(stdout,T);

 Matrix P(1,R.Rhi()*T.Rhi(),1,R.Rhi()*T.Rhi());
 //P=0;
 for(int a=1;a<=R.Rhi();++a)
 for(int r=1;r<=T.Rhi();++r)
 {for(int b=1;b<=a;++b)
  {for(int s=1;s<r;++s)
   {if(a==b)
    {int i=(a-1)*T.Rhi()+r;
     int j=(a-1)*T.Rhi()+s;
     if(i>=j)P(i,j)=R(a,a)*T(r,s);
     else    P(i,j)=-R(a,a)*T(s,r);
         i=(a-1)*T.Rhi()+s;
         j=(a-1)*T.Rhi()+r;
     if(i>=j)P(i,j)=R(a,a)*T(r,s);
     else    P(i,j)=R(a,a)*T(s,r);
    }
    else
    {
    int i=(a-1)*T.Rhi()+r;
    int j=(b-1)*T.Rhi()+s;
    if(i>=j)P(i,j)=R(a,b)*T(r,s)-R(b,a)*T(s,r);
    else    P(i,j)=-R(a,b)*T(s,r)-R(b,a)*T(r,s);
        i=(b-1)*T.Rhi()+r;
        j=(a-1)*T.Rhi()+s;
    if(i>=j)P(i,j)=R(a,b)*T(r,s)+R(b,a)*T(s,r);
    else    P(i,j)=-R(a,b)*T(s,r)+R(b,a)*T(r,s);
        i=(a-1)*T.Rhi()+s;
        j=(b-1)*T.Rhi()+r;
    if(i>=j)P(i,j)=R(a,b)*T(r,s)+R(b,a)*T(s,r);
    else    P(i,j)=-R(b,a)*T(r,s)+R(a,b)*T(s,r);
        i=(b-1)*T.Rhi()+s;
        j=(a-1)*T.Rhi()+r;
    if(i>=j)P(i,j)=R(a,b)*T(r,s)-R(b,a)*T(s,r);
    else    P(i,j)=R(b,a)*T(r,s)+R(a,b)*T(s,r);
    }
   }
   //s=r
   int i=(a-1)*T.Rhi()+r;
   int j=(b-1)*T.Rhi()+r;
   if(i>=j)P(i,j)=R(a,b)*T(r,r);
   else    P(i,j)=-R(b,a)*T(r,r);
       i=(b-1)*T.Rhi()+r;
       j=(a-1)*T.Rhi()+r;
   if(i>=j)P(i,j)=R(a,b)*T(r,r);
   else    {P(i,j)=R(b,a)*T(r,r);
            if(a==b)P(i,j)=0;
           }
  }
 }

//myPrintMatrix(stdout,P);
return P;
}

double aMb_real(Matrix & M, Matrix & zr,Matrix & zc, int ia, int ib) // transition matrix element
{double real=0.0;                                                    // <a|M|b>  a,b are columns ia and ib
                                                                     // of zr+izc

  if(M.Rhi()!=M.Chi()){fprintf(stderr,"Error martin.c aMb_real: Matrix M not square\n");exit(EXIT_FAILURE);}
  for(int a=1;a<=M.Rhi();++a)
  {for(int b=1;b<a;++b)
   {real+=zr(a,ia)*M(a,b)*zr(b,ib)-zr(a,ia)*M(b,a)*zc(b,ib)+zc(a,ia)*M(b,a)*zr(b,ib)+zc(a,ia)*M(a,b)*zc(b,ib);
   }
    real+=zr(a,ia)*M(a,a)*zr(a,ib)+zc(a,ia)*M(a,a)*zc(a,ib);
   for(int b=a+1;b<=M.Rhi();++b)
   {real+=zr(a,ia)*M(b,a)*zr(b,ib)+zr(a,ia)*M(a,b)*zc(b,ib)+zc(a,ia)*M(b,a)*zc(b,ib)-zc(a,ia)*M(a,b)*zr(b,ib);
   }

  }

 return real;
}

double aMb_imag(Matrix & M, Matrix & zr,Matrix & zc, int ia, int ib)
{double imag=0.0;
 if(M.Rhi()!=M.Chi()){fprintf(stderr,"Error martin.c aMb_imag: Matrix M not square\n");exit(EXIT_FAILURE);}
  for(int a=1;a<=M.Rhi();++a)
  {for(int b=1;b<a;++b)
   {imag+=zc(a,ia)*M(b,a)*zc(b,ib)-zc(a,ia)*M(a,b)*zr(b,ib)+zr(a,ia)*M(b,a)*zr(b,ib)+zr(a,ia)*M(a,b)*zc(b,ib);
   }
    imag+=zr(a,ia)*M(a,a)*zc(a,ib)-zc(a,ia)*M(a,a)*zr(a,ib);
   for(int b=a+1;b<=M.Rhi();++b)
   {imag+=zr(a,ia)*M(b,a)*zc(b,ib)-zr(a,ia)*M(a,b)*zr(b,ib)-zc(a,ia)*M(b,a)*zr(b,ib)-zc(a,ia)*M(a,b)*zc(b,ib);
   }

  }

 return imag;
}


Matrix MatrixfromVectors(Vector & v1,Vector & v2,Vector & v3)
{static Matrix m(1,3,v1.Lo(),v1.Hi());
 for(int i=v1.Lo();i<=v1.Hi();++i){m(i,1)=v1(i);m(i,2)=v2(i);m(i,3)=v3(i);}
return m;
}

void set_zlm_constants(Matrix & cnst)
{// cnst is the Zlm constants - put them into the matrix ... (same code is reused in ionpars.cpp)
if(cnst.Rlo()!=0||cnst.Rhi()!=6||cnst.Clo()!=-6||cnst.Chi()!=6){fprintf(stderr,"Error matrix in Zlm constants wrong dimension");exit(EXIT_FAILURE);}
cnst(0,0)=0.28209479177387814347403972578039;

cnst(1,0)=0.48860251190291992158638462283835;
cnst(1,1)=0.48860251190291992158638462283835;

cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;

cnst(3,0)=0.373176332590115391414395913199;
cnst(3,1)=0.45704579946446573615802069691665;
cnst(3,2)=1.445305721320277027694690077199;
cnst(3,3)=0.5900435899266435103456102775415;


cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;

cnst(5,0)=0.11695032245342359643971519209028;
cnst(5,1)=0.45294665119569692129844165821715;
cnst(5,2)=2.3967683924866618775009505697816;
cnst(5,3)=0.48923829943525038768400871584296;
cnst(5,4)=2.0756623148810412789957985225952;
cnst(5,5)=0.65638205684017010281411876637614;


cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
int l,m;
for(l=1;l<=6;l+=1){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}
}
