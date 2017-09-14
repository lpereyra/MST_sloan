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

/* Posiciones, velocidades y energias de las partículas */
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
extern double  rmaplim             ;  /* MAGNITUD APARENTE LIMITE DEL CATALOGO      */
extern double  rmapmin             ;  /* MAGNITUD APARENTE MINIMA DEL CATALOGO      */
extern double  redmax              ;  /* REDSHIT MAXIMO                             */
extern double  redmin              ;  /* REDSHIT MINIMO                             */
extern double  zcut                ;  /* REDSHIT CUT                                */
extern double  flfia               ;  /* AMPLITUD DE LA FL                          */
extern double  flma                ;  /* MAGNITUD CARACTERISTICA DE LA FL           */
extern double  flalfa              ;  /* PENDIENTE EN EL EXTREMO DEBIL DE LA FL     */
extern float   pmin[3]             ;
extern float   pmax[3]             ;
extern double  *fof                ;
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

void init_variables(int argc, char **argv);

#endif
