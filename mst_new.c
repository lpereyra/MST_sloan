#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>

#include "variables.h"
#include "cosmoparam.h"
#include "leesloan.h"
#include "timer.h"
#include "colores.h"
#include "voronoi.h"
#include "mst_kruskal.h"

#include <algorithm>

#ifdef CALCULA_MEDIA
  #include "grid.h"
#endif

#ifdef GAL_LUM
  std::vector<std::pair<float,int> > orden;
  #ifdef BRANCH_SURVIVE
    float mr_survive;
  #endif
#else
  std::vector<std::pair<int,int> > orden;
  #ifdef BRANCH_SURVIVE
    int N_part_survive;
  #endif
#endif

void Write_Segments(int *Padre, int *Rank, double *fof);

int main(int argc, char **argv)
{
  int    i;
  int    *Padre, *Rank;
  double start,end;
  std::vector<std::vector<int> > adjacency_list;
  std::vector<std::pair<float,std::pair<int,int> > > edges;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  #ifdef BRANCH_SURVIVE
  BLUE("********** Important *************\n");
    #ifdef GAL_LUM
     mr_survive = atof(argv[2]);
     if(mr_survive>=-1) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Mr Mag %g\n",mr_survive);
    #else
      N_part_survive = atoi(argv[2]);
      if(N_part_survive==0) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Numpar %g\n",N_part_survive);
    #endif
  RED(message);
  BLUE("**********************************\n");
  fflush(stdout);
  #endif

  readoutsloan();
 
  #ifndef GAL_LUM
    read_grup_fof(fof);
  #endif 

  /// Importante sino los grupos quedan fuera del box ///
  change_positions(cp.npart);

  Voronoi_Grupos(fof[0],edges);

  fprintf(stdout,"%d NumEdges\n",(int)edges.size());
  fflush(stdout);

  Padre = (int *) malloc(ngrup*sizeof(int));
  Rank =  (int *) malloc(ngrup*sizeof(int));

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = i;
    Rank[i] = 0;
    adjacency_list.push_back(std::vector<int>());
  }

  Kruskal(Padre,Rank,edges,adjacency_list);

  Podado(4,adjacency_list);  

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = -1;
    Rank[i]  = 0;

    if(adjacency_list[i].size()>0)
    #ifdef GAL_LUM
      orden.push_back(std::make_pair(fabs(Gr[i].mr),i));
    #else   
      orden.push_back(std::make_pair(Gr[i].NumPart,i));
    #endif
  }

  sort(orden.begin(),orden.end());
  
  for(i=(int)orden.size()-1;i>=0;i--)
  {
    int k = orden[i].second;

    if(Padre[k]==-1)
      DLU(k,Padre[k],adjacency_list,Padre,Rank);
  }

  fprintf(stdout,"Escribe\n");
  fflush(stdout);

  Write_Segments(Padre,Rank,fof);

  #ifdef CALCULA_MEDIA
  free(P);
  #endif  
  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}

void Write_Segments(int *Padre, int *Rank, double *fof)
{
  char filename[200];
  FILE *pfout, *pfpropiedades;
  int i,j,k,id;
  float dx,dy,dz;
  float dux,duy,duz;
  float r,lenr,elong,rms;
  std::vector<int> aux;
  #ifdef CALCULA_MEDIA
  float rcil2 = 1000.;
  float rlong = 500.;
  #endif

  j = 0;
  #ifdef GAL_LUM
    #ifdef BRANCH_SURVIVE    
      sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
    #else
      sprintf(filename,"segmentos_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
    #endif
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);

    #ifdef BRANCH_SURVIVE    
      sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
    #else
      sprintf(filename,"propiedades_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
    #endif      
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);
  #else
    #ifdef BRANCH_SURVIVE    
      sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
    #else
      sprintf(filename,"segmentos_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
    #endif
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);

    #ifdef BRANCH_SURVIVE    
      sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
    #else
      sprintf(filename,"propiedades_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
    #endif
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);
  #endif

  #ifdef CALCULA_MEDIA
  
  fprintf(stdout,"Build grid\n");

  grid.step = 0; /// IMPORTANTE
  grid.nobj = cp.npart;
  r = sqrt(rcil2);
 
  if(rlong>r){
    grid.ngrid = (int)(cp.lbox/rlong);
  }else{ 
    grid.ngrid = (int)(cp.lbox/r);
  }

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid_init();
  grid_build();

  fprintf(stdout,"End build grid\n");

  #endif

  while(!orden.empty())
  {

    i = orden.back().second;
 
    if(Rank[i]==1 && Padre[i]>=0)
    {
      aux.push_back(i);
      id = Padre[i];
      Rank[i] *= -1;

      while(id>=0)
      {
        if(Rank[id]<0)
        {
          aux.push_back(id);
          break;
        }

        //if(Rank[id]>2)
        #ifdef GAL_LUM
        if(Gr[id].mr<mr_survive && Padre[id]!=-1)
        #else
        if(Gr[id].NumPart>N_part_survive && Padre[id]!=-1)
        #endif
        {
          aux.push_back(id);

          #ifdef SORT_DERECHA
            #ifdef GAL_LUM
            if(Gr[aux[0]].mr<Gr[aux.back()].mr) // las magnitudes son negativas recordar 
              reverse(aux.begin(),aux.end());
            #else
            if(Gr[aux[0]].NumPart>Gr[aux.back()].NumPart)
              reverse(aux.begin(),aux.end());
            #endif
          #endif
 
          k = 0;
          #ifdef BRANCH_SURVIVE    

            #ifdef GAL_LUM
              if(Gr[aux[0]].mr<mr_survive) k++;
              if(Gr[id].mr<mr_survive)     k++;
            #else
              if(Gr[aux[0]].NumPart>N_part_survive) k++;
              if(Gr[id].NumPart>N_part_survive)     k++;
            #endif

          #else

            if(abs(Rank[id])>2) k++;
            if(abs(Rank[aux[0]])>2) k++;

          #endif


          fwrite(&k,sizeof(int),1,pfpropiedades);

          k = (int)aux.size();
          fwrite(&k,sizeof(int),1,pfout);
          fwrite(&k,sizeof(int),1,pfpropiedades);

          dux = Gr[id].Pos[0] - Gr[aux[0]].Pos[0];
          duy = Gr[id].Pos[1] - Gr[aux[0]].Pos[1];
          duz = Gr[id].Pos[2] - Gr[aux[0]].Pos[2];

          #ifdef PERIODIC
          dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
          dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;

          duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
          duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;

          duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
          duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
          #endif

          elong = dux*dux+duy*duy+duz*duz;

          lenr = rms = 0.0;

          for(k=0;k<(int)aux.size();k++)
          {

          #ifdef GAL_LUM
            fwrite(&Gr[aux[k]].id,sizeof(int),1,pfout);
          #else
            fwrite(&aux[k],sizeof(int),1,pfout);
          #endif

            if(k==0) continue;

            dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
            dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
            dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

            #ifdef PERIODIC
            dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
            dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

            dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
            dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

            dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
            dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
            #endif
            
            r = sqrt(dx*dx+dy*dy+dz*dz);

            lenr += r;

            if(k==(int)aux.size()-1) continue;

            dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
            dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
            dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

            #ifdef PERIODIC
            dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
            dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

            dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
            dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

            dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
            dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;

            #endif
                
            r = pow(dy*duz-dz*duy,2);
            r += pow(dz*dux-dx*duz,2);
            r += pow(dx*duy-dy*dux,2);
            r /= elong;
            
            rms+=r;

          }

          r = sqrt(elong)/lenr;

          k = (int)aux.size();
          rms /= (float)k;
          rms = sqrt(rms);

          fwrite(&lenr,sizeof(float),1,pfpropiedades);
          fwrite(&r,sizeof(float),1,pfpropiedades);
          fwrite(&rms,sizeof(float),1,pfpropiedades);

          aux.clear();
          j++;
        }

        aux.push_back(id);
        Rank[id] *= -1;
        id = Padre[id];               
      }

      #ifdef SORT_DERECHA
        #ifdef GAL_LUM
          if(Gr[aux[0]].mr<Gr[aux.back()].mr) // las magnitudes son negativas recordar 
            reverse(aux.begin(),aux.end());
        #else
          if(Gr[aux[0]].NumPart>Gr[aux.back()].NumPart)
            reverse(aux.begin(),aux.end());
        #endif
      #endif
 
      id = aux.back();

      k = 0;
      #ifdef BRANCH_SURVIVE    
        #ifdef GAL_LUM
          if(Gr[aux[0]].mr<mr_survive) k++;
          if(Gr[id].mr<mr_survive)     k++;
        #else
          if(Gr[aux[0]].NumPart>N_part_survive) k++;
          if(Gr[id].NumPart>N_part_survive)     k++;
        #endif
      #else
        if(abs(Rank[id])>2) k++;
        if(abs(Rank[aux[0]])>2) k++;
      #endif

      fwrite(&k,sizeof(int),1,pfpropiedades);

      k = (int)aux.size();
      fwrite(&k,sizeof(int),1,pfout);
      fwrite(&k,sizeof(int),1,pfpropiedades);

      dux = Gr[id].Pos[0] - Gr[aux[0]].Pos[0];
      duy = Gr[id].Pos[1] - Gr[aux[0]].Pos[1];
      duz = Gr[id].Pos[2] - Gr[aux[0]].Pos[2];

      #ifdef PERIODIC
      dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
      dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;
    
      duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
      duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;
    
      duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
      duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
      #endif

      elong = dux*dux+duy*duy+duz*duz;

      lenr = rms = 0.0;

      for(k=0;k<(int)aux.size();k++)
      {

      #ifdef GAL_LUM
        fwrite(&Gr[aux[k]].id,sizeof(int),1,pfout);
      #else
        fwrite(&aux[k],sizeof(int),1,pfout);
      #endif

        if(k==0) continue;

        dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
        dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
        dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif
        
        r = sqrt(dx*dx+dy*dy+dz*dz);
        
        lenr += r;

        if(k==(int)aux.size()-1) continue;

        dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
        dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
        dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif
            
        r = pow(dy*duz-dz*duy,2);
        r += pow(dz*dux-dx*duz,2);
        r += pow(dx*duy-dy*dux,2);
        r /= elong;

        rms += r;

      }

      r = sqrt(elong)/lenr;

      k = (int)aux.size();
      rms /= (float)k;
      rms = sqrt(rms);

      fwrite(&lenr,sizeof(float),1,pfpropiedades);
      fwrite(&r,sizeof(float),1,pfpropiedades);
      fwrite(&rms,sizeof(float),1,pfpropiedades);

      aux.clear();
      j++;
    }

    orden.pop_back();
  }

  rewind(pfout);
  fwrite(&j,sizeof(int),1,pfout);
  fclose(pfout);

  rewind(pfpropiedades);
  fwrite(&j,sizeof(int),1,pfpropiedades);
  fclose(pfpropiedades);

  #ifdef CALCULA_MEDIA
  grid_free();
  #endif
}
