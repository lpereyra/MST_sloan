#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesloan.h"
#include "colores.h"

int     nfrac               ;
char    message[200]        ;
double  rmaplim             ;  // MAGNITUD APARENTE LIMITE DEL CATALOGO    
double  rmapmin             ;  // MAGNITUD APARENTE MINIMA DEL CATALOGO    
double  redmax              ;  // REDSHIT MAXIMO                           
double  redmin              ;  // REDSHIT MINIMO                           
double  flfia               ;  // AMPLITUD DE LA FL                        
double  flma                ;  // MAGNITUD CARACTERISTICA DE LA FL         
double  flalfa              ;  // PENDIENTE EN EL EXTREMO DEBIL DE LA FL     
double  vfid                ;  // VELOCIDAD FIDUCIAL                   
double  v0                  ;  // EN MPC H^-1                          
double  rmablim             ;  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO
double  rmabmin             ;  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO
double  rintlim             ;  // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
float   pmin[3]             ;
float   pmax[3]             ;
double  *fof                ;
#ifdef LIM_VOL
double  zcut                ;  // REDSHIT CUT                              
#endif
#ifdef GAL_LUM
float   mcut                ;
#endif

void init_variables(int argc, char **argv){
  FILE *pfin;
  char filename[200];

  RED("Initializing variables...\n");

  nfrac = 2;
  fof = (double *) malloc(nfrac*sizeof(double));
 
  sprintf(filename,"%s",argv[1]);
  if(!(pfin=fopen(filename,"r")))
  {
    sprintf(message,"can't open file `%s` \n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%d  \n",&snap.nfiles))
  {
    sprintf(message,"can't read file `%s`\nneed # of snapshots\n",filename);RED(message);
    exit(0);
  }
    
  if(!fscanf(pfin,"%s  \n",&snap.root))
  {
    sprintf(message,"can't read file `%s`\nneed snapshots directory\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%s  \n",&snap.name))
  {
    sprintf(message,"can't read file `%s`\nneed snapname\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf  \n",&vfid))
  {
    sprintf(message,"can't read file `%s`\nvel fiducial\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf \n",&fof[0]))
  {
    sprintf(message,"can't read file `%s`\nlen fof[0]\n",filename);RED(message);
    exit(0);
  }

  if(!fscanf(pfin,"%lf \n",&fof[1]))
  {
    sprintf(message,"can't read file `%s`\nneed len fof[1]\n",filename);RED(message);
    exit(0);
  }

  #ifdef GAL_LUM
  if(!fscanf(pfin,"%f \n",&mcut))
  {
    sprintf(message,"can't read file `%s`\nneed mcut\n",filename);RED(message);
    exit(0);
  }

  if(mcut>=0.)
  {
    sprintf(message,"error mcut\n");RED(message);
    exit(0);
  }
  #endif

  #ifdef LIM_VOL
  if(!fscanf(pfin,"%lf \n",&zcut))
  {
    sprintf(message,"can't read file `%s`\n redshift cut\n",filename);RED(message);
    exit(0);
  }

  if(zcut<0.)
  {
    sprintf(message,"error zcut\n");RED(message);
    exit(0);
  }
  #endif

  fclose(pfin);  

  rmaplim  = 17.77    ;  // MAGNITUD APARENTE LIMITE DEL CATALOGO      
  rmapmin  = 14.5     ;  // MAGNITUD APARENTE MINIMA DEL CATALOGO      
  redmax   = 0.2      ;  // REDSHIT MAXIMO                             
  redmin   = 0.01     ;  // REDSHIT MINIMO                             
  flfia    = 0.0149   ;  // AMPLITUD DE LA FL                          
  flma     = -20.44   ;  // MAGNITUD CARACTERISTICA DE LA FL           
  flalfa   = -1.05    ;  // PENDIENTE EN EL EXTREMO DEBIL DE LA FL     
  v0       = 2.0      ;  // EN MPC H^-1   VELOCIDAD FIDUCIAL           

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
  sprintf(message,"Magnitud aparente limite %f\n",rmaplim);BLUE(message);
  sprintf(message,"Magnitud aparente minima %f\n",rmapmin);BLUE(message);
  #ifdef LIM_VOL
  sprintf(message,"redshift cut             %f\n",zcut);BLUE(message);
  #endif
  BLUE("************* Options ************\n");
  sprintf(message,"primera identificacion   %f\n",fof[0]);RED(message);
  sprintf(message,"segunda identificacion   %f\n",fof[1]);RED(message);
  #ifdef GAL_LUM
  sprintf(message,"glx mas luminosas que    %f\n",mcut);RED(message);
  #endif
  BLUE("************* LEN FOF ************\n");
  #ifdef LEN_FOF_MERCHAN
  GREEN("Usa LEN FOF MERCHAN\n");
  #else
  GREEN("Usa LEN FOF APROX VOL\n");
  #endif
  BLUE("**********************************\n");

  #ifdef LOCK
  RED("USING LOCKS\n");
  #endif


  GREEN("END\n");
}
