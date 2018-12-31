/* pario.cpp
 *
 * Parameter input/output routines
 *
 * Functions:
 *   int  getdim(int n, orbital l);                                            // Number of states = ^{4l+2}C_{n}
 *   void getfromionname(std::string &ion, icpars &flags);                     // Gets free ion parameters from tables
 *   void ic_parsecfpars(std::string &n, std::string &v, icpars &p, int l=1);  // Parses CF parameter for k and q
 *   void ic_parseinput(const char *file, icpars &flags);                      // Parses file for 1-ion pars & phys prop.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2010 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <fstream>
#include <iomanip>
#include <ctime>

void ic_parsecfpars(std::string &varname, std::string &varval, icpars &pars, int length)
{
   std::string tmpstr = varname.substr(length,1);
   int k = atoi(tmpstr.c_str()), q;
   if((int)varname.length()==(length+3))
   {
      if(varname.compare(length+1,1,"m")==0 || varname.compare(length+1,1,"-")==0) tmpstr = varname.substr(length+2,1);
      else if(varname.compare(length+2,1,"s")==0) tmpstr = varname.substr(length+1,1);
      else { std::cerr << "ic_parseinput: CF parameter name " << varname << " not recognised.\n"; return; }
      q = -atoi(tmpstr.c_str());
   }
   else { tmpstr = varname.substr(length+1,1); q = atoi(tmpstr.c_str()); }
   tmpstr = varname.substr(0,length);
   pars.B.assign(tmpstr,k,q,atof(varval.c_str()));
}

// --------------------------------------------------------------------------------------------------------------- //
// Parses input file for IC parameters
// --------------------------------------------------------------------------------------------------------------- //
void ic_parseinput(const char *filename, icpars &pars)
{
   std::string strline,varname,varval,tmpstr,issstr;
   size_t equal_sym_pos,first_nonspace_pos,first_num_pos=std::string::npos,last_num_pos=std::string::npos;
   std::istringstream iss;
   int k;

   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) { std::cerr << "ic_parseinput(): Cannot open file " << filename << "\n"; return; }

   // If the parameters file is from cfield, assume parameters are in meV, in Stevens normalisation, with theta_k=1.
   getline(FILEIN,strline);
   if(strline.compare(0,8,"#!cfield")==0)
   {
      tmpstr.assign("meV"); conv_e_units(pars,tmpstr); tmpstr.assign("B"); pars.B.conv(tmpstr);
   }

   while(!FILEIN.eof())
   {
      strline.clear(); getline(FILEIN,strline);
      if(strline.find("#",0,1)!=std::string::npos)              // Line is a comment, ignores it
         continue; 

      first_nonspace_pos = strline.find_first_not_of(" \t");
      equal_sym_pos = strline.find("=");
    
      varname.clear(); varval.clear();
      if(first_nonspace_pos==std::string::npos || equal_sym_pos==std::string::npos) 
         varname.assign(strline);
      else
      {
         iss.clear();
         varname = strline.substr(first_nonspace_pos,equal_sym_pos-first_nonspace_pos);
         varname = varname.substr(0,varname.find_first_of(" \t"));
         strtolower(varname);
         varval = strline.substr(equal_sym_pos+1);
         size_t subfirst = varval.find_first_not_of(" \t"); 
         size_t sublast  = varval.find_last_not_of(" \t")+1; 
         if(subfirst!=std::string::npos)
         {
            varval = varval.substr(subfirst,sublast);
            first_num_pos = varval.find_first_of("-0123456789.EDed");
            last_num_pos = varval.find_first_not_of("-0123456789.EDed",first_num_pos);
            if(first_num_pos!=std::string::npos)
               iss.str(varval.substr(first_num_pos,last_num_pos));
            else
               iss.str(varval);
         }
         else
            iss.str(varval);
      }

      if(varname.find("iontype")!=std::string::npos)
         getfromionname(varval,pars);                           // Looks up values of Fk and xi from tables
      else if(varname.find("unit")!=std::string::npos)
         conv_e_units(pars,varval);
      else if(varname.find("conf")!=std::string::npos)
      {  tmpstr = varval.substr(0,1); pars.l = Lin(tmpstr); tmpstr = varval.substr(1,2); pars.n = atoi(tmpstr.c_str()); }
      else if(varname.compare("n")==0)
      {  iss >> pars.n; pars.B.calc_stevfact(pars.n,pars.l); }
      else if(varname.compare("l")==0)
      {  tmpstr = varval.substr(0,1); pars.l = Lin(tmpstr); pars.B.calc_stevfact(pars.n,pars.l); }
      else if(varname.compare("f")==0)                          // Slater integrals expressed as four floats
         for(k=0; k<4; k++)
         {
            iss >> pars.F[k];
            first_num_pos = varval.find_first_of("0123456789.-EeDd",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.-EeDd",first_num_pos);
            iss.str(varval.substr(first_num_pos,last_num_pos));
         }
      else if(varname.compare("f0")==0) { iss >> pars.F[0]; pars._F[0]=pars.F[0]*pars._econv; }
      else if(varname.compare("f2")==0) { iss >> pars.F[1]; pars._F[1]=pars.F[1]*pars._econv; }
      else if(varname.compare("f4")==0) { iss >> pars.F[2]; pars._F[2]=pars.F[2]*pars._econv; }
      else if(varname.compare("f6")==0) { iss >> pars.F[3]; pars._F[3]=pars.F[3]*pars._econv; }
      else if(varname.compare("xi")==0) { iss >> pars.xi; pars._xi = pars.xi*pars._econv; }
      else if(varname.compare("zeta")==0) { iss >> pars.xi; pars._xi = pars.xi*pars._econv; }
      else if(varname.compare("alpha")==0)
      {
         if(varval.find(",")!=std::string::npos)
            for(k=0; k<3; k++)
            {
               iss >> pars.alpha[k]; pars._alpha[k] = pars.alpha[k]*pars._econv;
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               iss.str(varval.substr(first_num_pos,last_num_pos));
            }
         else                              { iss >> pars.alpha[0]; pars._alpha[0]=pars.alpha[0]*pars._econv; }
      }
      else if(varname.compare("beta")==0)  { iss >> pars.alpha[1]; pars._alpha[1]=pars.alpha[1]*pars._econv; }
      else if(varname.compare("gamma")==0) { iss >> pars.alpha[2]; pars._alpha[2]=pars.alpha[2]*pars._econv; }

      else if(varname.find_first_of("awbvld")==0 && varname.find_first_of("0123456789")==1) 
         ic_parsecfpars(varname, varval, pars);
      else if(varname.compare(0,2,"ar")==0 && varname.find_first_of("0123456789")==2)
         ic_parsecfpars(varname, varval, pars, 2);
      else if(varname.find("density")!=std::string::npos)
         pars.density = varval;
//    else if(varname.find("observable")!=std::string::npos)
//       pars.observable = varval;
      else if(varname.find("eigenvectors")!=std::string::npos)
         iss >> pars.num_eigv;
      else if(varname.find("basis")!=std::string::npos)
         pars.basis = varval;
      else if(varname.find("calc")!=std::string::npos)
      {
         if(varname.find("mag")!=std::string::npos)  // Physical properties calculation flags.
            pars.calcphys += PHYSPROP_MAGBIT; 
         else if(varname.find("cp")!=std::string::npos || varname.find("heat")!=std::string::npos)
            pars.calcphys += PHYSPROP_CP_BIT; 
         else if(varname.find("susc")!=std::string::npos) {
            if(varname.find("inv")!=std::string::npos) 
               pars.calcphys += PHYSPROP_INVBIT;
            else 
               pars.calcphys += PHYSPROP_SUSBIT; }
      }
      else if(varname.find("partial_standalone")!=std::string::npos)
         pars.partial_standalone = true;
      else if(varname.find("partial")!=std::string::npos)
         pars.partial = true;
      #ifndef NO_ARPACK
      else if(varname.find("arnoldi_standalone")!=std::string::npos)
         pars.arnoldi_standalone = true;
      else if(varname.find("arnoldi")!=std::string::npos)
         pars.arnoldi = true;
      #endif
      else if(varname.find("save_matrices")!=std::string::npos)
         pars.save_matrices = true;
      else if(varname.find("truncate_matrix")!=std::string::npos)
      {  
         iss >> pars.truncate_level; if(pars.truncate_level<=0 || pars.truncate_level>=1) pars.truncate_level=0.5;
      }  
      else if(varname.find("use_J_operator_equivalent")!=std::string::npos)
      {  pars.B.op_equiv = Jt; pars.B.calc_stevfact(pars.n,pars.l); pars.B.convback(); }
      else if(varname.find("use_L_operator_equivalent")!=std::string::npos)
      {  pars.B.op_equiv = Lt; pars.B.calc_stevfact(pars.n,pars.l); pars.B.convback(); }
      else if(varname.find("dx2")!=std::string::npos) iss >> pars.Dx2;
      else if(varname.find("dy2")!=std::string::npos) iss >> pars.Dy2;
      else if(varname.find("dz2")!=std::string::npos) iss >> pars.Dz2;
      else if(varname.find("emu")!=std::string::npos)
         pars.mag_units = 1;
      else if(varname.find("simag")!=std::string::npos)
         pars.mag_units = 2;
      else if(varname.compare("xval")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xT = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
         if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHa = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
            if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHb = atof(tmpstr.c_str());
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
               {
                  tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xHc = atof(tmpstr.c_str());
               }
            }
         }
      }
      else if(varname.compare("xrange")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xMin = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xStep = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	    if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.xMax = atof(tmpstr.c_str());
            }
         }
      }
      else if(varname.compare("yval")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yT = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHa = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
            if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHb = atof(tmpstr.c_str());
               first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
               last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
               if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
               {
                  tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yHc = atof(tmpstr.c_str());
               }
            }
         }
      }
      else if(varname.compare("yrange")==0)
      {
         tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yMin = atof(tmpstr.c_str());
         first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
         last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	 if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
         {
            tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yStep = atof(tmpstr.c_str());
            first_num_pos = varval.find_first_of("0123456789.",last_num_pos);
            last_num_pos = varval.find_first_not_of("0123456789.",first_num_pos);
	    if(first_num_pos!=std::string::npos || last_num_pos!=std::string::npos)
            {
               tmpstr = varval.substr(first_num_pos,last_num_pos-first_num_pos); pars.yMax = atof(tmpstr.c_str());
            }
         }
      }
      else if(varname.compare("xt")==0)    { tmpstr=iss.str(); pars.xT    = atof(tmpstr.c_str()); }
      else if(varname.compare("xha")==0)   { tmpstr=iss.str(); pars.xHa   = atof(tmpstr.c_str()); }
      else if(varname.compare("xhb")==0)   { tmpstr=iss.str(); pars.xHb   = atof(tmpstr.c_str()); }
      else if(varname.compare("xhc")==0)   { tmpstr=iss.str(); pars.xHc   = atof(tmpstr.c_str()); }
      else if(varname.compare("yt")==0)    { tmpstr=iss.str(); pars.yT    = atof(tmpstr.c_str()); }
      else if(varname.compare("yha")==0)   { tmpstr=iss.str(); pars.yHa   = atof(tmpstr.c_str()); }
      else if(varname.compare("yhb")==0)   { tmpstr=iss.str(); pars.yHb   = atof(tmpstr.c_str()); }
      else if(varname.compare("yhc")==0)   { tmpstr=iss.str(); pars.yHc   = atof(tmpstr.c_str()); }
      else if(varname.compare("xmin")==0)  { tmpstr=iss.str(); pars.xMin  = atof(tmpstr.c_str()); }
      else if(varname.compare("xstep")==0) { tmpstr=iss.str(); pars.xStep = atof(tmpstr.c_str()); }
      else if(varname.compare("xmax")==0)  { tmpstr=iss.str(); pars.xMax  = atof(tmpstr.c_str()); }
      else if(varname.compare("ymin")==0)  { tmpstr=iss.str(); pars.yMin  = atof(tmpstr.c_str()); }
      else if(varname.compare("ystep")==0) { tmpstr=iss.str(); pars.yStep = atof(tmpstr.c_str()); }
      else if(varname.compare("ymax")==0)  { tmpstr=iss.str(); pars.yMax  = atof(tmpstr.c_str()); }
      else if(varname.compare("bx")==0)    { tmpstr=iss.str(); pars.Bx    = atof(tmpstr.c_str()); }
      else if(varname.compare("by")==0)    { tmpstr=iss.str(); pars.By    = atof(tmpstr.c_str()); }
      else if(varname.compare("bz")==0)    { tmpstr=iss.str(); pars.Bz    = atof(tmpstr.c_str()); }
   }  // while(!FILEIN.eof())
   FILEIN.close();
}

