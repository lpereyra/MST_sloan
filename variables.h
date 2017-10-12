#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
  #define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
  #define VELFACTOR 1.0
#endif

#ifndef NPARTMIN
  #define NPARTMIN 4 
#endif

#ifndef NTABLA
  #define NTABLA 1000
#endif

#ifndef NTABLADIS
  #define NTABLADIS 10000
#endif

/* Posiciones, velocidades de las partículas */
struct particle_data 
{
  float          Pos[3];
  float          Dis;
  int            sub;
  int            gr;
};

struct grup_data
{
  int        save;
  int        id;
  float      Pos[3];
  #ifdef GAL_LUM
    float      mr; 
  #else
    int        NumPart;
  #endif
};

extern int     nfrac               ;
extern int     ngrup               ;
extern struct  particle_data *P    ;
extern struct  grup_data *Gr       ;
extern double  rmaplim             ;  // MAGNITUD APARENTE LIMITE DEL CATALOGO     
extern double  rmapmin             ;  // MAGNITUD APARENTE MINIMA DEL CATALOGO     
extern double  redmax              ;  // REDSHIT MAXIMO                            
extern double  redmin              ;  // REDSHIT MINIMO                            
extern double  flfia               ;  // AMPLITUD DE LA FL                         
extern double  flma                ;  // MAGNITUD CARACTERISTICA DE LA FL          
extern double  flalfa              ;  // PENDIENTE EN EL EXTREMO DEBIL DE LA FL    
extern double  vfid                ;  // VELOCIDAD FIDUCIAL                       
extern double  v0                  ;  // EN MPC H^-1                              
extern double  rmablim             ;  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
extern double  rmabmin             ;  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
extern double  rintlim             ;  // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
extern float   pmin[3]             ;
extern float   pmax[3]             ;
extern double  *fof                ;
#ifdef LIM_VOL
  extern double  zcut              ;  // REDSHIT CUT                               
#endif
#ifdef GAL_LUM
  extern float mcut;
  #ifdef BRANCH_SURVIVE
    extern float mr_survive;
  #endif
#else
  #ifdef BRANCH_SURVIVE
    extern int N_part_survive;
  #endif
#endif

// AUXILIARES TABLA //
extern double mt1min                    ;  
extern double mt1max                    ;
extern double mt2min                    ;  
extern double mt2max                    ;
extern double dmt1                      ;
extern double dmt2                      ;
extern double **rinttabla               ;
extern double *dis2red                  ;
extern double *rt                       ;

void init_variables(int argc, char **argv);

#endif
