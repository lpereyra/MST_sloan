#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesloan.h"

int    ngrup;
struct SnapST snap;
struct particle_data *P;
struct grup_data *Gr;
struct cosmoparam cp;
struct galsloang  *gal;

// AUXILIARES TABLA //
double mt1min                    ;  
double mt1max                    ;
double mt2min                    ;  
double mt2max                    ;
double dmt1                      ;
double dmt2                      ;
double **rinttabla               ;
double *dis2red                  ;
double *rt                       ;

void readoutsloan()
{

  FILE *pf;
  char filename[200];
  int i,j,iden_01,iden_02;
  float disred, dlum;
  float mt1,mt2;
  double dismin,dismax;
  #ifdef GAL_LUM
  float mtest;
  ngrup = 0;
  #endif

  cp.omegam = 0.3                       ;  /* OMEGA MATERIA                              */
  cp.omegal = 0.7                       ;  /* OMEGA LAMBDA                               */
  cp.omegak = 1.0-cp.omegam-cp.omegal   ;  /* OMEGA CURVATURA                            */
  cp.h0     = 100.                      ;  /* ESTO DEJA TODO EN UNIDADES DE H^-1         */
  cp.dlummax = -1.E26                   ;
  mt1min     = 1.E26                    ;
  mt1max     = -1.E26                   ;
  mt2min     = 1.E26                    ;
  mt2max     = -1.E26                   ;
  dismin     = 1.E26                    ;
  dismax     = -1.E26                   ;

  RED("Read OUT Sloan...\n");

  sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");

  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&cp.npart, sizeof(int),1,pf);

  /* ALOCATACION Y LECTURA */
  gal = (struct galsloang *) malloc(cp.npart*sizeof(struct galsloang));
  P   = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  #ifdef GAL_LUM
  Gr = (struct grup_data *) malloc(cp.npart*sizeof(struct grup_data));
  #endif

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  iden_01 = iden_02 = 0;
  for(i=0;i<cp.npart;i++)
  {
    fread(&gal[i].gal.id,sizeof(long int),1,pf);
    fread(&gal[i].gal.targetid,sizeof(long int),1,pf);
    fread(&gal[i].gal.specid,sizeof(long int),1,pf);
    fread(&gal[i].gal.alfa,sizeof(float),1,pf);
    fread(&gal[i].gal.delta,sizeof(float),1,pf);
    fread(&gal[i].gal.red,sizeof(float),1,pf);
    fread(&gal[i].gal.rederr,sizeof(float),1,pf);
    fread(&gal[i].gal.m[0],sizeof(float),1,pf);
    fread(&gal[i].gal.m[1],sizeof(float),1,pf);
    fread(&gal[i].gal.m[2],sizeof(float),1,pf);
    fread(&gal[i].gal.m[3],sizeof(float),1,pf);
    fread(&gal[i].gal.m[4],sizeof(float),1,pf);
    fread(&gal[i].gal.merr[0],sizeof(float),1,pf);
    fread(&gal[i].gal.merr[1],sizeof(float),1,pf);
    fread(&gal[i].gal.merr[2],sizeof(float),1,pf);
    fread(&gal[i].gal.merr[3],sizeof(float),1,pf);
    fread(&gal[i].gal.merr[4],sizeof(float),1,pf);
    fread(&gal[i].gal.modelm[0],sizeof(float),1,pf);
    fread(&gal[i].gal.modelm[1],sizeof(float),1,pf);
    fread(&gal[i].gal.modelm[2],sizeof(float),1,pf);
    fread(&gal[i].gal.modelm[3],sizeof(float),1,pf);
    fread(&gal[i].gal.modelm[4],sizeof(float),1,pf);
    fread(&gal[i].gal.modelmerr[0],sizeof(float),1,pf);
    fread(&gal[i].gal.modelmerr[1],sizeof(float),1,pf);
    fread(&gal[i].gal.modelmerr[2],sizeof(float),1,pf);
    fread(&gal[i].gal.modelmerr[3],sizeof(float),1,pf);
    fread(&gal[i].gal.modelmerr[4],sizeof(float),1,pf);
    fread(&gal[i].gal.ext[0],sizeof(float),1,pf);
    fread(&gal[i].gal.ext[1],sizeof(float),1,pf);
    fread(&gal[i].gal.ext[2],sizeof(float),1,pf);
    fread(&gal[i].gal.ext[3],sizeof(float),1,pf);
    fread(&gal[i].gal.ext[4],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosg[0],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosg[1],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosg[2],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosgerr[0],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosgerr[1],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosgerr[2],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosr[0],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosr[1],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosr[2],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosrerr[0],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosrerr[1],sizeof(float),1,pf);
    fread(&gal[i].gal.petrosrerr[2],sizeof(float),1,pf);
    fread(&gal[i].gal.k[0],sizeof(float),1,pf);
    fread(&gal[i].gal.k[1],sizeof(float),1,pf);
    fread(&gal[i].gal.k[2],sizeof(float),1,pf);
    fread(&gal[i].gal.k[3],sizeof(float),1,pf);
    fread(&gal[i].gal.k[4],sizeof(float),1,pf);  
    fread(&gal[i].grupo[0],sizeof(int),1,pf);  
    fread(&gal[i].grupo[1],sizeof(int),1,pf);  
    fread(&gal[i].npgrupo[0],sizeof(int),1,pf);  
    fread(&gal[i].npgrupo[1],sizeof(int),1,pf);  

    gal[i].gal.alfa  *= PI180;
    gal[i].gal.delta *= PI180;

    #ifdef LIM_VOL
    assert(gal[i].gal.red<zcut);
    #endif

    disred = red2dis(gal[i].gal.red); //EN MPC
    dlum = disred*(1.0+gal[i].gal.red); 

    if(dlum>cp.dlummax) cp.dlummax = dlum; // DISTANCIA LUMINOSIDAD MAXIMA

    mt2  = rmaplim-25.0-5.0*log10(disred*(1.0+gal[i].gal.red)) ;
    if(mt2>mt2max) mt2max = mt2; //magnitud absoluta limite maxima
    if(mt2<mt2min) mt2min = mt2; //magnitud absoluta limite minima

    mt1  = rmapmin-25.0-5.0*log10(disred*(1.0+gal[i].gal.red)) ;
    if(mt1>mt1max) mt1max = mt1; //magnitud absoluta minima maxima
    if(mt1<mt1min) mt1min = mt1; //magnitud absoluta minima minima

    P[i].Pos[0] = disred*cos(gal[i].gal.delta)*cos(gal[i].gal.alfa) ;
    P[i].Pos[1] = disred*cos(gal[i].gal.delta)*sin(gal[i].gal.alfa) ;
    P[i].Pos[2] = disred*sin(gal[i].gal.delta)                      ;
    P[i].Dis    = disred;
    P[i].sub    = gal[i].grupo[0];

    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];

    if(P[i].Dis > dismax) dismax = P[i].Dis;
    if(P[i].Dis < dismin) dismin = P[i].Dis;

    if(gal[i].grupo[0] > iden_01) iden_01 = gal[i].grupo[0];
    if(gal[i].grupo[1] > iden_02) iden_02 = gal[i].grupo[1];

    #ifdef GAL_LUM

    if(P[i].sub==0) continue;

    mtest = (gal[i].gal.m[2]-gal[i].gal.k[2])-25.0-5.0*log10(disred*(1.0+gal[i].gal.red));
    if(mtest<mcut)
    {
      Gr[ngrup].save = P[i].sub;
      Gr[ngrup].id = i;
      memcpy(Gr[ngrup].Pos,P[i].Pos,3*sizeof(float));
      Gr[ngrup].mr = mtest;
      ngrup++;
    }

    #endif
  }

  fprintf(stdout,"Num Total %d\n",cp.npart);
  fprintf(stdout,"Num Total de grupos en la primera identificacion %d\n",iden_01);
  fprintf(stdout,"Num Total de grupos en la segunda identificacion %d\n",iden_02);
  fflush(stdout);

  cp.vol = (2.210)*(pow(dismax,3.)-pow(dismin,3.))/3.0;

  #ifdef GAL_LUM
  fprintf(stdout,"%d Num de glx mas lum que %f\n",ngrup,mcut);
  fflush(stdout);
  Gr = (struct grup_data *) realloc (Gr,ngrup*sizeof(struct grup_data));
  #endif

  RED("End Read OUT Sloan\n");

  fclose(pf);

  #ifdef LIM_VOL
  rmablim = rmaplim-25.0-5.0*log10(red2dis(zcut)*(1.0+zcut)); // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = MAGMENOSINF;  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  #else
  rmablim = rmaplim-25.0-5.0*log10(vfid);  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = rmapmin-25.0-5.0*log10(vfid);  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  #endif

  rintlim = intfl(rmabmin,rmablim);        // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
  fprintf(stdout, "INTEGRAL LIMITES\n");
  fprintf(stdout, "%f MABS correspondiente a %f MAPA con una velocidad fiducial %f\n",rmabmin,rmapmin,vfid);
  fprintf(stdout, "%f MABS correspondiente a %f MAPA con una velocidad fiducial %f\n",rmablim,rmaplim,vfid);
  fflush(stdout);

  mt2min *= (1.0-1.0e-2*mt2min/fabs(mt2min)); // agranda y achica los limites 
  mt2max *= (1.0+1.0e-2*mt2max/fabs(mt2max));
  mt1min *= (1.0-1.0e-2*mt1min/fabs(mt1min));
  mt1max *= (1.0+1.0e-2*mt1max/fabs(mt1max));

  // TABLA DE LA INTEGRAL DEL FACTOR DE ESCALA //
  dmt1 = (mt1max-mt1min)/(float)NTABLA;
  dmt2 = (mt2max-mt2min)/(float)NTABLA;
  dlum = (redmax-redmin)/(float)NTABLADIS;

  fprintf(stdout,"mt1min mt1max %f %f\n", mt1min, mt1max);
  fprintf(stdout,"mt2min mt2max %f %f\n", mt2min, mt2max);
  fprintf(stdout,"redmin redmax %f %f\n", redmin, redmax);
  fprintf(stdout,"Crea las tablas\n");

  rt = (double *) malloc(NTABLADIS*sizeof(double));
  dis2red = (double *) malloc(NTABLADIS*sizeof(double));
  rinttabla = (double **) malloc(NTABLA*sizeof(double *));
  for(i=0;i<NTABLA;i++)
    rinttabla[i] = (double *) malloc(NTABLA*sizeof(double));

  #ifdef NTHREADS
    #pragma omp parallel for default(none) \
    num_threads(NTHREADS) private(i,j,mt1,mt2) \
    shared(mt1min,dmt1,mt2min,dmt2,rinttabla)   
  #endif
  for(i=0;i<NTABLA;i++ )
  {
    mt1 = mt1min + (double)i*dmt1;
    for(j=0;j<NTABLA;j++ )
    {
      mt2 = mt2min + (double)j*dmt2;
      rinttabla[i][j] = intfl(mt1,mt2);
    }
  }

  #ifdef NTHREADS
    #pragma omp parallel for default(none) \
    num_threads(NTHREADS) private(i,disred) \
    shared(redmin,dlum,rt,dis2red)   
  #endif
  for(i=0;i<NTABLADIS;i++)
  {
    disred     = redmin + (double)i*dlum;
    rt[i]      = disred;
    dis2red[i] = red2dis(disred);
  }

  fprintf(stdout,"Termina las tablas\n");
  fflush(stdout);

}

#ifndef GAL_LUM

void read_grup_fof(double *fof)
{
  char filename[200];
  int  i;
  FILE *pfin;
 
  GREEN("Read Grups...\n");
  #ifdef LIM_VOL
    sprintf(filename,"%s%.2f_%.2f_%.2f_centros.bin",snap.root,zcut,fof[1],fof[0]);
  #else
    sprintf(filename,"%s%.2f_%.2f_centros.bin",snap.root,fof[1],fof[0]);
  #endif
  pfin = fopen(filename,"rb"); 

  fread(&ngrup,sizeof(int),1,pfin);

  fprintf(stdout,"Grupos %d\n",ngrup);
  fflush(stdout);

  Gr = (struct grup_data *) malloc(ngrup*sizeof(struct grup_data));

  for(i=0;i<ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(int),1,pfin);
    fread(&Gr[i].id,sizeof(int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[1],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[2],sizeof(float),1,pfin);
    fread(&Gr[i].NumPart,sizeof(int),1,pfin);
  }

  fclose(pfin);

  GREEN("End Read Grups\n");

}

#endif

void change_positions(int n)
{
  int ip, idim;

  RED("Inicio Change Positions\n");

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip = 0; ip < n; ip++)
    for(idim = 0; idim < 3; idim++)
      P[ip].Pos[idim] -= pmin[idim];

  for(ip = 0; ip < ngrup; ip++)
    for(idim = 0; idim < 3; idim++)
      Gr[ip].Pos[idim] -= pmin[idim];

  cp.lbox = pmax[0] - pmin[0];
  for(idim = 1; idim < 3; idim++)
    if(cp.lbox < (pmax[idim] - pmin[idim])) cp.lbox = (pmax[idim] - pmin[idim]);

  cp.lbox *= 1.001;
  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
}

double funlum(double x, void *p)
{
  struct paramfl *par=(struct paramfl *)p ;
  double ma=(par->ma)                     ;
  double alfa=(par->alfa)                 ;
  double t                                ;

  t=pow(10.0,0.4*(ma-x)) ;
  t=pow(t,1.0+alfa)*exp(-t) ;

  return(t) ;
}

double intfl(double x1, double x2)
{
  double resultado, error;
  size_t neval;
  struct paramfl pfl;
  gsl_function F; 

  pfl.ma=flma;
  pfl.alfa=flalfa;
  pfl.fia=flfia;

  F.function = &funlum;
  F.params = &pfl;

  if(x1<MAGMENOSINF)x1=MAGMENOSINF;
  gsl_integration_qng(&F,x1,x2,1.0e-7,1.0e-7,&resultado,&error,&neval);
  return(resultado);
}

double f(double z, void *p) /* FUNCION A INTEGRAR PARA LA DISTANCIA EN FUNCION DE Z*/
{
  struct paramcos *par=(struct paramcos *)p ;
  double om=(par->omegam);
  double ol=(par->omegal);
  double ok=(par->omegak);
  double q;
  q=pow(1.+z,3.)*om+pow(1.+z,2.)*ok+ol;
  return(1.0/sqrt(q));
}

double red2dis(double z)
{
  double resultado, error;
  size_t neval;
  struct paramcos pcos;
  gsl_function F;

  pcos.omegam=cp.omegam;
  pcos.omegal=cp.omegal;
  pcos.omegak=cp.omegak;

  F.function=&f;
  F.params=&pcos;

  gsl_integration_qng(&F,0.0,z,1.0e-7,1.0e-7,&resultado,&error,&neval);
 
  return(CVEL/cp.h0*resultado);
}
