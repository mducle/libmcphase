/***************************************************
 * charges - display charges at given htpoint
 * Author: M. Rotter
 **************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#define MUB  5.788378E-02 // Bohrmagneton in meV/tesla

#include "../../version"
#include "spincf.hpp"
#include "martin.h"
#include "myev.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>
#include<par.hpp>




/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ spincf savmf;
 FILE * fin_coq, * fout;
 float delta,dd,ddT,ddHa,ddHb,ddHc,alpha,beta,gamma;
 int i,n=0,nofatoms=0,nofcomponents=3;
 long int pos=0,j;
 float numbers[13];numbers[9]=1;numbers[10]=3;
 numbers[0]=13;
 char instr[MAXNOFCHARINLINE];
 char outstr[MAXNOFCHARINLINE];
 char filename[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS],gJ[MAXNOFATOMS];
 char * cffilenames[MAXNOFATOMS];
// ComplexMatrix * eigenstates[MAXNOFATOMS];
  Matrix r(1,3,1,3);
  Vector abc(1,3);
  int max_ext_nof_components=51;
  int ext_nof_components[MAXNOFATOMS];
printf("#***************************************************\n");
printf("# * charges - display charges at given htpoint\n");
printf("# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("# * %s\n",MCPHASVERSION);
printf("# **************************************************\n");
// check command line
  if (argc < 5)
    { printf (" program charges - display charges at HT point\n\
                use as: charges T Ha Hb Hc [file.mf]\n \
                            (default input file is results/mcphas.mf)\n \
                    or: charges T Ha Hb Hc h k l E\n \
                        (reads from results/mcphas.mf and results/mcdisp.qev\n\n \
                This program outputs a magnetic/charge structure (and magnetic/orbital excitation)\n \
                graphic/movie in the output files results/charges*.jvx (javaview)\n\n \
                the graphics output format can be fine tuned in .mf and .qev input files\n\n \
                jvx files can be viewed by: java javaview results/charges.jvx \n \
                                            java javaview \"model=results/charges.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n\n");
      exit (1);
    }

if (argc==6) 
 { fin_coq = fopen_errchk (argv[5], "rb");}
 else
 { fin_coq = fopen_errchk ("./results/mcphas.mf", "rb");}
    
 fout = fopen_errchk ("./results/charges.out", "w");



double show_abc_unitcell=1.0,show_primitive_crystal_unitcell=1.0,show_magnetic_unitcell=1.0,show_atoms=1.0,scale_view_1=1.0,scale_view_2=1.0,scale_view_3=1.0;
double show_chargedensity=1.0,show_spindensity=1.0;
abc=0;char *token;
 // input file header ------------------------------------------------------------------
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
   if (pos==-1) 
       {fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq);
   // strip /r (dos line feed) from line if necessary
    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,"%s",instr);}
   if(abc[1]==0){extract(instr,"a",abc[1]);extract(instr,"b",abc[2]); extract(instr,"c",abc[3]); 
                 extract(instr,"alpha",alpha);  extract(instr,"beta",beta);extract(instr,"gamma",gamma); 
   }
   extract(instr,"show_abc_unitcell",show_abc_unitcell);
   extract(instr,"show_primitive_crystal_unitcell",show_primitive_crystal_unitcell);
   extract(instr,"show_magnetic_unitcell",show_magnetic_unitcell);
   extract(instr,"show_atoms",show_atoms);
   extract(instr,"show_spindensity",show_spindensity);
   extract(instr,"show_chargedensity",show_chargedensity);

   extract(instr,"scale_view_1",scale_view_1);
   extract(instr,"scale_view_2",scale_view_2);
   extract(instr,"scale_view_3",scale_view_3);

   extract(instr,"r1x",r[1][1]);extract(instr,"r2x",r[1][2]); extract(instr,"r3x",r[1][3]); 
   extract(instr,"r1y",r[2][1]); extract(instr,"r2y",r[2][2]); extract(instr,"r3y",r[2][3]);
   extract(instr,"r1z",r[3][1]); extract(instr,"r2z",r[3][2]); extract(instr,"r3z",r[3][3]);
   extract(instr,"r1a",r[1][1]);extract(instr,"r2a",r[1][2]); extract(instr,"r3a",r[1][3]); 
   extract(instr,"r1b",r[2][1]); extract(instr,"r2b",r[2][2]); extract(instr,"r3b",r[2][3]);
   extract(instr,"r1c",r[3][1]); extract(instr,"r2c",r[3][2]); extract(instr,"r3c",r[3][3]);
   extract(instr,"nofatoms",nofatoms);    extract(instr,"nofcomponents",nofcomponents); 
   if (nofatoms>0&&(extract(instr,"x",x[n+1])+
                   extract(instr,"y",y[n+1])+
  		       extract(instr,"z",z[n+1])==0)||
		       (extract(instr,"da",x[n+1])+
                   extract(instr,"db",y[n+1])+
		       extract(instr,"dc",z[n+1])==0))
		  {++n;if(n>nofatoms||nofatoms>MAXNOFATOMS) 
                    {fprintf(stderr,"ERROR charges.c reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"cffilename",cffilenames[n],(size_t)MAXNOFCHARINLINE);
                   extract(instr,"gJ",gJ[n]);
  ext_nof_components[n]=48;
  if (gJ[n]==0)
  {ext_nof_components[n]=51;  // here set for 3+48 components, module ic1ion
   fprintf(stderr,"WARNING program charges: gJ=0 for ion %i - intermediate coupling calculations not supported yet, will create only charges.out and no chargeplot charges.jvx\n",n);}     
//		   printf("%s\n",cffilenames[n]);
                  }
  }
  if (alpha!=90||beta!=90||gamma!=90)
  {fprintf(stderr,"ERROR: non orthogonal lattice not supported yet\n");exit(EXIT_FAILURE);}
   Vector gJJ(1,n); for (i=1;i<=n;++i){gJJ(i)=gJ[i];}
  
  
// load mfconfigurations and check which one is nearest -------------------------------   
  double T,ha,hb,hc;

   j=fseek(fin_coq,pos,SEEK_SET); 
    if (j!=0){fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
   
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      spins.load(fin_coq);
      ddT=strtod(argv[1],NULL)-numbers[3];ddT*=ddT;
      ddHa=strtod(argv[2],NULL)-numbers[5];ddHa*=ddHa;
      ddHb=strtod(argv[3],NULL)-numbers[6];ddHb*=ddHb;
      ddHc=strtod(argv[4],NULL)-numbers[7];ddHc*=ddHc;
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if (dd<delta)
       {delta=dd;
        sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",numbers[3],numbers[5],numbers[6],numbers[7],numbers[8],(int)numbers[9],(int)numbers[10]);
        savmf=spins;T=numbers[3];ha=numbers[5];hb=numbers[6];hc=numbers[7];
       }
    }
  fclose (fin_coq);
  
// create plot of spin+chargeconfiguration -----------------------------------------------------------
  int ii,nt,k,l,m;
  double lnz,u;
  float d;

  par inputpars("./mcphas.j");
  
  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  spincf extendedspincf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,max_ext_nof_components);

          // the following is for the printout of charges.out ...........................
           fprintf(fout,"#!T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,ha,hb,hc,savmf.na(),savmf.nb(),savmf.nc(),inputpars.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
            fprintf(fout,"#J=value {atom-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] <Ja> <Jb> <Jc> ...\n");
          fprintf(fout,"#{corresponding effective fields gjmbHeff [meV]- if passed to mcdiff only these are used for caculation (not the magnetic moments)}\n");

  // determine primitive magnetic unit cell
     Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
     savmf.calc_prim_mag_unitcell(p,abc,r);
  // .............................................................................                                
	       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: mcalc has to be used to calculate all the <Olm>.
hh=0;for(ii=1;ii<=inputpars.nofatoms;++ii)
{(*inputpars.jjj[ii]).eigenstates(hh,T);} // initialize eigenstate matrices

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  extendedspincf.m(i,j,k)=0;
  for(ii=1;ii<=inputpars.nofatoms;++ii)
 {  Vector h(1,ext_nof_components[ii]);
    Vector moments(1,ext_nof_components[ii]);
    h=0;
   for(nt=1;nt<=savmf.nofcomponents;++nt){h(nt)=hh(nt+savmf.nofcomponents*(ii-1));}
            if((*inputpars.jjj[ii]).module_type!=1&&(*inputpars.jjj[ii]).module_type!=3)
            {(*inputpars.jjj[ii]).mcalc(moments,T,h,lnz,u,(*inputpars.jjj[ii]).est); // here we trigger single ion
                                                           // module to calculate all 48 (ext_nof_components)
                                                           // higher order moments 
            }
            else
            {moments=0;}

          // output atoms and moments in primitive unit cell to stdout
              Vector dd3(1,3);
              dd3=savmf.pos(i,j,k,ii, abc, r,x,y,z);   
              dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
              fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              cffilenames[ii],dd3(1)/abc(1),dd3(2)/abc(2),dd3(3)/abc(3),dd0(1),dd0(2),dd0(3));
              printf("{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	              cffilenames[ii],dd3(1)/abc(1),dd3(2)/abc(2),dd3(3)/abc(3),dd0(1),dd0(2),dd0(3));
                     for(nt=1;nt<=3;++nt){if(gJ[ii]!=0){fprintf(fout," %4.4f",gJ[ii]*moments(nt));}
                                          else         {fprintf(fout," %4.4f",2*moments(2*nt-1)+moments(2*nt));}
                                         }
                     for(nt=1;nt<=ext_nof_components[ii];++nt)                                                               // this else is when gJ=0: means intermediate coupling
		        {extendedspincf.m(i,j,k)(nt+max_ext_nof_components*(ii-1))=moments(nt);
                         fprintf(fout," %4.4f",extendedspincf.m(i,j,k)(nt+max_ext_nof_components*(ii-1)));
                        }
                         fprintf(fout,"\n");
                      fprintf(fout,"                  corresponding effective fields gjmbHeff [meV]-->          ");
                      for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout meanfields
                        {fprintf(fout," %4.4f",h(nt));}
                         fprintf(fout,"\n");
                             
//	                 myPrintComplexMatrix(fout,(*inputpars.jjj[ii]).eigenstates(h));      
							   // ... and the eigenvalues + eigenvectors !

  }}}}

fclose(fout);

//print out the long vector of moments 1-48
  printf("%s - spin configuration <Olm>(i)\n",outstr);
  extendedspincf.print(stdout);

             Vector hkl(1,3);hkl=0;
             spincf savev_real(extendedspincf*0.0);
             spincf savev_imag(extendedspincf*0.0);
             
  fout = fopen_errchk ("./results/charges.grid", "w");
     extendedspincf.cd(fout,abc,r,x,y,z,cffilenames,0,50,50,50,scale_view_1,scale_view_2,scale_view_3,
                       savev_real,savev_imag,0.0,0.0,hkl);
    fclose (fout);

  fout = fopen_errchk ("./results/charges.jvx", "w");
     extendedspincf.jvx_cd(fout,outstr,abc,r,x,y,z,gJJ,show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell,show_atoms,scale_view_1,scale_view_2,scale_view_3,
                  0,0.0,savev_real,savev_imag,0.0,hkl,0.0,0.0,cffilenames,show_chargedensity,show_spindensity);
    fclose (fout);

  fout = fopen_errchk ("./results/charges_prim.jvx", "w");
     extendedspincf.jvx_cd(fout,outstr,abc,r,x,y,z,gJJ,show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell,show_atoms,scale_view_1,scale_view_2,scale_view_3,
                  1,0.0,savev_real,savev_imag,0.0,hkl,0.0,0.0,cffilenames,show_chargedensity,show_spindensity);
    fclose (fout);


if (argc>=9){// try a spinwave picture
             double h,k,l,E,ddh,ddk,ddl,ddE;
             int extended_eigenvector_dimension;
             double spins_wave_amplitude=1.0,spins_show_ellipses=1.0,spins_show_direction_of_static_moment=1.0; 
             fin_coq = fopen_errchk ("./results/mcdisp.qee", "rb");
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin_coq); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;} 
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",spins_show_ellipses);
               extract(instr,"spins_show_direction_of_static_moment",spins_show_direction_of_static_moment);
               extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
              }
//               if(extended_eigenvector_dimension!=48){fprintf(stderr,"Error program charges - extended_eigenvector_dimension in results/mcdisp.qee not equal to 48\n");exit(1);}
               j=fseek(fin_coq,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 spincf ev_real(extendedspincf.na(),extendedspincf.nb(),extendedspincf.nc(),extendedspincf.nofatoms,extended_eigenvector_dimension);
                 spincf ev_imag(extendedspincf.na(),extendedspincf.nb(),extendedspincf.nc(),extendedspincf.nofatoms,extended_eigenvector_dimension);
                 ev_real.load(fin_coq);ev_imag.load(fin_coq);
                 ddT=strtod(argv[1],NULL)-numbers[4];ddT*=ddT;
                 ddHa=strtod(argv[2],NULL)-numbers[1];ddHa*=ddHa;
                 ddHb=strtod(argv[3],NULL)-numbers[2];ddHb*=ddHb;
                 ddHc=strtod(argv[4],NULL)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];E=numbers[9];                                        
                  savev_real=ev_real;
                  savev_imag=ev_imag;                  
                 }
               }
              fclose (fin_coq);
              fprintf(stdout,"#%s - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              savev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              savev_imag.print(stdout);

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
              complex <double> im(0,1);
              for(i=0;i<16;++i)
              {phase=2*3.1415*i/15;
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(16)\n",i+1);
               printf("********************************************\n");
               
               sprintf(filename,"./results/charges.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");
                     extendedspincf.jvx_cd(fin_coq,outstr,abc,r,x,y,z,gJJ,show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell,show_atoms,scale_view_1,scale_view_2,scale_view_3,
                                  0,phase,savev_real,savev_imag,spins_wave_amplitude,hkl,spins_show_ellipses,spins_show_direction_of_static_moment,cffilenames,show_chargedensity,show_spindensity);
               fclose (fin_coq);
               sprintf(filename,"./results/charges_prim.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");
                     extendedspincf.jvx_cd(fin_coq,outstr,abc,r,x,y,z,gJJ,show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell,show_atoms,scale_view_1,scale_view_2,scale_view_3,
                                  1,phase,savev_real,savev_imag,spins_wave_amplitude,hkl,spins_show_ellipses,spins_show_direction_of_static_moment,cffilenames,show_chargedensity,show_spindensity);
               fclose (fin_coq);
              }
          printf("# %s\n",outstr);
          }

fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program charges\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/charges.jvx\n");
fprintf(stderr,"# * java javaview \"model=results/charges.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n");
fprintf(stderr,"# ************************************************************************\n");

  for(i=1;i<=nofatoms;++i){  delete cffilenames[i];}

  return 0;



}


