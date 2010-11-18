/***************************************************
 * orbmomdensities - display orbmomdensities at given htpoint
 * Author: M. Rotter
 **************************************************/
#include "../../version"
#include "spincf.hpp"
#include "martin.h"
#include "myev.h"
#include<par.hpp>
#include<graphic_parameters.hpp>
#include <cryststruct.hpp>
#include "densities_func.c"

/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("# **********************************************************\n");
printf("# * orbmomdensities - display orbmomdensities at given H and T *\n");
printf("# * Reference: M. Rotter PRB 79 (2009) 140405R             *\n");
printf("# * %s                                     *\n",MCPHASVERSION);
printf("# **********************************************************\n");
// check command line
  if (argc < 6)
    { printf ("\n \
program orbmomdensities - display orbmomdensities at HT point\n\n \
use as: orbmomdensities threshhold T Ha Hb Hc [i j k] [file.mf]\n \
    or: orbmomdensities threshhold T Ha Hb Hc -div [file.mf]\n \
                        (default input file is results/mcphas.mf)\n\n \
This program outputs orbmomdensities on of magnetic ions in the magnetic\n \
unit cell the graphics output format can be fine tuned in .mf and\n \
and results/graphic_parameters.set\n\n \
jvx files can be viewed by: java javaview results/orbmomdensities.jvx \n \
                options: \n\
                -div triggers calculation of divergence of the vector field\n\
\n\n");
      exit (1);
    }
FILE * fin_coq, * fout;
   double T,ha,hb,hc,xx=0,yy=0,zz=0;
int doijk=0;// to be implemented : orbmomdensity-component along specific direction (doijk=1,2,3)
 graphic_parameters gp;
 gp.threshhold=strtod(argv[1],NULL);
 gp.scale_density_vectors=1;

   spincf savmf;
   int i,n;
   char outstr[MAXNOFCHARINLINE];
   cryststruct cs;
// read input file with mf configuration
if (argc>=9){
  xx=strtod(argv[6],NULL);
  yy=strtod(argv[7],NULL);
  zz=strtod(argv[8],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;
  doijk=3;
 }
if (argc>7&&strncmp(argv[6],"-div",4)==0)
{doijk=1;
}
 if (argc==7+doijk)
 { fin_coq = fopen_errchk (argv[6+doijk], "rb");}
 else
 { fin_coq = fopen_errchk ("./results/mcphas.mf", "rb");}

//  fout = fopen_errchk ("./results/orbmomdensities.out", "w");
   // input file header and mfconf------------------------------------------------------------------
   n=headerinput(fin_coq,stderr,gp,cs);

   // check for spinfconfiguration which is nearest to the T/H values chosen by user in command line
   check_for_best(fin_coq,strtod(argv[2],NULL),strtod(argv[3],NULL),strtod(argv[4],NULL),strtod(argv[5],NULL),savmf,T,ha,hb,hc,outstr);
  fclose (fin_coq);

// create plot of  -----------------------------------------------------------
  int ii,nt,k,j;
  //double lnz,u;
  

  par inputpars("./mcphas.j");

  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  int dim=49; if(doijk<3){dim=3*49;}
  spincf extendedspincf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,dim);

  // determine primitive magnetic unit cell
     Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
     savmf.calc_prim_mag_unitcell(p,cs.abc,cs.r);
  // .............................................................................

//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: mcalc has to be used to calculate all the <Olm>.
hh=0;for(ii=1;ii<=inputpars.nofatoms;++ii)
{(*inputpars.jjj[ii]).mcalc_parameter_storage_init(hh,T);} // initialize mcalc module parameter storage

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  extendedspincf.m(i,j,k)=0;
  for(ii=1;ii<=inputpars.nofatoms;++ii)
 {  Vector h(1,savmf.nofcomponents);
    Vector moms(1,savmf.nofcomponents);
    Vector moments(1,dim);
    Vector momentsx(1,49);
    Vector momentsy(1,49);
    Vector momentsz(1,49);
    h=0;
   for(nt=1;nt<=savmf.nofcomponents;++nt){h(nt)=hh(nt+savmf.nofcomponents*(ii-1));}
            if((*inputpars.jjj[ii]).module_type==0)
            {//(*inputpars.jjj[ii]).mcalc(moms,T,h,lnz,u,(*inputpars.jjj[ii]).mcalc_parstorage); // here we trigger single ion
                                                           // module to calculate all dim moments of orbmomdensity
if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_mcalc (momentsx,1, T, h, (*inputpars.jjj[ii]).mcalc_parstorage);
if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_mcalc (momentsy,2, T, h, (*inputpars.jjj[ii]).mcalc_parstorage);
if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_mcalc (momentsz,3, T, h, (*inputpars.jjj[ii]).mcalc_parstorage);
if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
else       {
            for(nt=1;nt<=49;++nt){moments(nt)=momentsx(nt);moments(nt+49)=momentsy(nt);moments(nt+2*49)=momentsz(nt);}
            }
            }
            else
            {moments=0;}

          // output atoms and moments in primitive unit cell to stdout
              Vector dd3(1,3);
              dd3=savmf.pos(i,j,k,ii, cs);
              dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
              printf("{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	              cs.cffilenames[ii],dd3(1)/cs.abc(1),dd3(2)/cs.abc(2),dd3(3)/cs.abc(3),dd0(1),dd0(2),dd0(3));

                     for(nt=1;nt<=dim;++nt)
		        {extendedspincf.m(i,j,k)(nt+dim*(ii-1))=moments(nt);
                    }

  }}}}
 gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
 gp.spins_scale_moment=0;
 gp.spins_show_static_moment_direction=0;
 gp.spins_wave_amplitude=0;
 gp.spins_show_oscillation=0;

//print out the long vector of moments 1-48
  printf("%s - spin configuration moments(i)\n",outstr);
  extendedspincf.print(stdout);

             Vector hkl(1,3);hkl=0;
             spincf savev_real(extendedspincf*0.0);
             spincf savev_imag(extendedspincf*0.0);
             gp.showprim=0;
  if(doijk==3) sprintf(gp.title,"projection of orbmomdensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of orbmomdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of orbmomdensity |ML(r)|");
  printf("%s\n",gp.title);

  fout = fopen_errchk ("./results/orbmomdensities.grid", "w");
     extendedspincf.cd(fout,cs,gp,savev_real,savev_imag,0.0,hkl,T,hh);
    fclose (fout);

  fout = fopen_errchk ("./results/orbmomdensities.jvx", "w");
    gp.showprim=0;
     extendedspincf.jvx_cd(fout,outstr,cs,gp,
                  0.0,savev_real,savev_imag,hkl,T,hh);
    fclose (fout);

  fout = fopen_errchk ("./results/orbmomdensities_prim.jvx", "w");
     gp.showprim=1;
    extendedspincf.jvx_cd(fout,outstr,cs,gp,
                  0.0,savev_real,savev_imag,hkl,T,hh);
    fclose (fout);



fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program orbmomdensities\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n\n results/orbmomdensities.grid created");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/orbmomdensities.jvx\n");
fprintf(stderr,"# * java javaview \"model=results/orbmomdensities.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n");
fprintf(stderr,"# * saved density mesh in results/obrmomdensities.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  for(i=1;i<=cs.nofatoms;++i){  delete cs.cffilenames[i];}

  return 0;



}