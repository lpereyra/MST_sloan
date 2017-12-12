#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.h"
#include "cosmoparam.h"
#include "leesloan.h"
#include "timer.h"
#include "voronoi.h"
#include "mst_kruskal.h"
#include "colores.h"

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
#ifdef CALCULA_MEDIA
 void calc_media(float xc, float yc, float zc, \
  float vx, float vy, float vz, int *numpart, float *red_media, \
  float rcil2, float rlong_2);
 bool point_inside(float dot, float rlong_2);
#endif

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
     if(mr_survive>=0) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Mr Mag %g\n",mr_survive);
    #else
      N_part_survive = atoi(argv[2]);
      if(N_part_survive==0) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Numpar %d\n",N_part_survive);
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

  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();
  #ifdef CALCULA_MEDIA
  free(P);
  #endif
  for(i=0;i<NTABLA;i++)
    free(rinttabla[i]);
  free(rinttabla);
  free(dis2red);
  free(rt);

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
  int *contador;
  float dx,dy,dz;
  float dux,duy,duz;
  float r,lenr,elong,rms;
  std::vector<int> aux;

  #ifdef CALCULA_MEDIA
  std::vector<std::vector<int> > segmentos;
  float xc,yc,zc;
  float rbin,racum;
  float redmed;
  int   npart;
  int   bin,lbin;
  float rcil  = 1.0;
  float rcil2 = rcil*rcil;
  float rlong = 0.5;
  #endif

  j = 0;

  #ifdef GAL_LUM

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"%.2f_segmentos_%2.2f_lum_%2.2f.bin",fabs(mr_survive),fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"segmentos_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"segmentos_%2.2f_lum_%2.2f.bin",fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"%.2f_propiedades_%2.2f_lum_%2.2f.bin",fabs(mr_survive),fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"propiedades_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"propiedades_%2.2f_lum_%2.2f.bin",fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);

  #else

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"%.2f_segmentos_%2.2f_%2.2f.bin",(float)N_part_survive,fof[0],fof[1]);
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"segmentos_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"segmentos_%2.2f_%2.2f.bin",fof[0],fof[1]);
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"%.2f_propiedades_%2.2f_%2.2f.bin",(float)N_part_survive,fof[0],fof[1]);
      #endif
    #else
      #ifdef LIM_VOL
        sprintf(filename,"propiedades_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"propiedades_%2.2f_%2.2f.bin",fof[0],fof[1]);
      #endif
    #endif   // close BRANCH_SURVIVE
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);

  #endif // close GAL_LUM

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

  contador = (int *) calloc(3,sizeof(int));

  while(!orden.empty())
  {

    i = orden.back().second;
 
    //if(Rank[i]==1 && Padre[i]>=0)
    if(Padre[i]>=0)
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
          
          #ifdef CALCULA_MEDIA
          segmentos.push_back(aux);
          #endif
 
          k = 0;
          #ifdef BRANCH_SURVIVE    

            #ifdef GAL_LUM
              if(Gr[aux[0]].mr<mr_survive) k++;
              if(Gr[aux.back()].mr<mr_survive) k++;
            #else
              if(Gr[aux[0]].NumPart>N_part_survive) k++;
              if(Gr[aux.back()].NumPart>N_part_survive) k++;
            #endif

          #else

            if(abs(Rank[aux.back()])>2) k++;
            if(abs(Rank[aux[0]])>2) k++;

          #endif

          contador[k]++;

          fwrite(&k,sizeof(int),1,pfpropiedades);

          k = (int)aux.size();
          fwrite(&k,sizeof(int),1,pfout);
          fwrite(&k,sizeof(int),1,pfpropiedades);

          dux = Gr[aux.back()].Pos[0] - Gr[aux[0]].Pos[0];
          duy = Gr[aux.back()].Pos[1] - Gr[aux[0]].Pos[1];
          duz = Gr[aux.back()].Pos[2] - Gr[aux[0]].Pos[2];

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

       #ifdef CALCULA_MEDIA
      segmentos.push_back(aux);
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

      contador[k]++;
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

  for(i=0;i<3;i++)
    fprintf(stdout,"Segmentos %d %d\n",i,contador[i]);
  fprintf(stdout,"Segmentos %d\n",j);
  free(contador);

  #ifdef CALCULA_MEDIA

  assert(j==(int)segmentos.size());

  BLUE("********** Importante ***********\n");
  sprintf(filename,"Radio      %f Mpc\n",sqrt(rcil2)); RED(filename);
  sprintf(filename,"Largo      %f Mpc\n",(2.0*rlong)); RED(filename);
  sprintf(filename,"Separacion %f Mpc\n",rlong);RED(filename);
  BLUE("**********************************\n");

  i = 0;

  while(!segmentos.empty())
  {
    aux = segmentos.back();    
    segmentos.back().clear();

    npart = 0;
    redmed = 0.0;
    racum = 0.0;
    
    for(k=1;k<(int)aux.size();k++)
    {

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

      dx /= r;
      dy /= r;
      dz /= r;

      if(k==1)
      {

        xc = Gr[aux[0]].Pos[0];
        yc = Gr[aux[0]].Pos[1];
        zc = Gr[aux[0]].Pos[2]; 

        calc_media(xc,yc,zc,dx,dy,dz,&npart,&redmed,rcil2,rlong);

      }else{

        xc = Gr[aux[k-1]].Pos[0];
        yc = Gr[aux[k-1]].Pos[1];
        zc = Gr[aux[k-1]].Pos[2];

        rbin = rlong-racum;
        xc += rbin*dx;
        yc += rbin*dy;
        zc += rbin*dz;   

        #ifdef PERIODIC
        xc = xc >= cp.lbox ? xc-cp.lbox : xc;
        xc = xc <      0.0 ? xc+cp.lbox : xc;
      
        yc = yc >= cp.lbox ? yc-cp.lbox : yc;
        yc = yc <      0.0 ? yc+cp.lbox : yc;
 
        zc = zc >= cp.lbox ? zc-cp.lbox : zc;
        zc = zc <      0.0 ? zc+cp.lbox : zc;
        #endif

        calc_media(xc,yc,zc,dx,dy,dz,&npart,&redmed,rcil2,rlong);

        dux = Gr[aux[k]].Pos[0]-xc;
        duy = Gr[aux[k]].Pos[1]-yc;
        duz = Gr[aux[k]].Pos[2]-zc;

        #ifdef PERIODIC
        dux = dux >= 0.5*cp.lbox ? dux-cp.lbox : dux;
        dux = dux < -0.5*cp.lbox ? dux+cp.lbox : dux;
      
        duy = duy >= 0.5*cp.lbox ? duy-cp.lbox : duy;
        duy = duy < -0.5*cp.lbox ? duy+cp.lbox : duy;
 
        duz = duz >= 0.5*cp.lbox ? duz-cp.lbox : duz;
        duz = duz < -0.5*cp.lbox ? duz+cp.lbox : duz;
        #endif

        r = sqrt(dux*dux+duy*duy+duz*duz);

      }

      lbin = (int)(r/rlong);

      for(bin=0;bin<lbin;bin++)
      {
        xc += rlong*dx;
        yc += rlong*dy;
        zc += rlong*dz;   

        #ifdef PERIODIC
        xc = xc >= cp.lbox ? xc-cp.lbox : xc;
        xc = xc <      0.0 ? xc+cp.lbox : xc;
        
        yc = yc >= cp.lbox ? yc-cp.lbox : yc;
        yc = yc <      0.0 ? yc+cp.lbox : yc;
 
        zc = zc >= cp.lbox ? zc-cp.lbox : zc;
        zc = zc <      0.0 ? zc+cp.lbox : zc;
        #endif

        calc_media(xc,yc,zc,dx,dy,dz,&npart,&redmed,rcil2,rlong);
      }

      dux = xc-Gr[aux[k]].Pos[0];
      duy = yc-Gr[aux[k]].Pos[1];
      duz = zc-Gr[aux[k]].Pos[2];

      #ifdef PERIODIC
      dux = dux >= 0.5*cp.lbox ? dux-cp.lbox : dux;
      dux = dux < -0.5*cp.lbox ? dux+cp.lbox : dux;
      
      duy = duy >= 0.5*cp.lbox ? duy-cp.lbox : duy;
      duy = duy < -0.5*cp.lbox ? duy+cp.lbox : duy;
 
      duz = duz >= 0.5*cp.lbox ? duz-cp.lbox : duz;
      duz = duz < -0.5*cp.lbox ? duz+cp.lbox : duz;
      #endif
      
      racum = sqrt(dux*dux+duy*duy+duz*duz);
      
    }

    calc_media(Gr[aux.back()].Pos[0],Gr[aux.back()].Pos[1],Gr[aux.back()].Pos[2], \
    dx,dy,dz,&npart,&redmed,rcil2,rlong);            

    ///////////////////////////////////////////////////////////////////////////////////

    dx = Gr[aux[1]].Pos[0] - Gr[aux[0]].Pos[0];
    dy = Gr[aux[1]].Pos[1] - Gr[aux[0]].Pos[1];
    dz = Gr[aux[1]].Pos[2] - Gr[aux[0]].Pos[2];

    #ifdef PERIODIC
    if(dx> 0.5*cp.lbox) dx -= cp.lbox;
    if(dx<-0.5*cp.lbox) dx += cp.lbox;
    if(dy> 0.5*cp.lbox) dy -= cp.lbox;
    if(dy<-0.5*cp.lbox) dy += cp.lbox;
    if(dz> 0.5*cp.lbox) dz -= cp.lbox;
    if(dz<-0.5*cp.lbox) dz += cp.lbox;
    #endif

    r = sqrt(dx*dx+dy*dy+dz*dz);

    dx /= r;
    dy /= r;
    dz /= r;

    xc = Gr[aux[0]].Pos[0]+racum*dx;
    yc = Gr[aux[0]].Pos[1]+racum*dy;
    zc = Gr[aux[0]].Pos[2]+racum*dz;

    #ifdef PERIODIC
    xc = xc >= cp.lbox ? xc-cp.lbox : xc;
    xc = xc <      0.0 ? xc+cp.lbox : xc;
    
    yc = yc >= cp.lbox ? yc-cp.lbox : yc;
    yc = yc <      0.0 ? yc+cp.lbox : yc;
 
    zc = zc >= cp.lbox ? zc-cp.lbox : zc;
    zc = zc <      0.0 ? zc+cp.lbox : zc;
    #endif

    //printf("%d %f %d %f\n",i,redmed,npart,redmed/(float)npart);

    calc_media(xc,yc,zc,dx,dy,dz,&npart,&redmed,rcil2,rlong);
    ///////////////////////////////////////////////////////////////////////////////////

    //printf("%d %f %d %f\n",i,redmed,npart,redmed/(float)npart);
    assert(npart!=0);

    segmentos.pop_back();
    i++;
  }

  grid_free();

  #endif
 
  return;

}

#ifdef CALCULA_MEDIA

void calc_media(float xc, float yc, float zc, \
float vx, float vy, float vz, int *numpart, float *red_media, \
float rcil2, float rlong_2)
{
  int i;
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  float Posprima[3];
  float lbox,fac;
  #ifdef PERIODIC
  float lbox2;
  #endif
  float dot,dis;
  int ngrid;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (float)ngrid/lbox;
  #ifdef PERIODIC
  lbox2 = lbox/2.0;
  #endif

  ixc  = (int)(xc*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(yc*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(zc*fac);
  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif
  
      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {
          Posprima[0] = P[i].Pos[0] - xc;
          Posprima[1] = P[i].Pos[1] - yc;
          Posprima[2] = P[i].Pos[2] - zc;

          #ifdef PERIODIC
          if(Posprima[0] >  lbox2) Posprima[0] = Posprima[0] - lbox;
          if(Posprima[1] >  lbox2) Posprima[1] = Posprima[1] - lbox;
          if(Posprima[2] >  lbox2) Posprima[2] = Posprima[2] - lbox;
          if(Posprima[0] < -lbox2) Posprima[0] = Posprima[0] + lbox;
          if(Posprima[1] < -lbox2) Posprima[1] = Posprima[1] + lbox;
          if(Posprima[2] < -lbox2) Posprima[2] = Posprima[2] + lbox;
          #endif

          dot = Posprima[0]*vx+Posprima[1]*vy+Posprima[2]*vz;

          if(point_inside(dot,rlong_2))
          {        

            dis = Posprima[0]*Posprima[0]+Posprima[1]*Posprima[1]+Posprima[2]*Posprima[2]-dot*dot;

            if(dis<rcil2)
            {
              *red_media = *red_media + gal[i].gal.red;
              *numpart   = *numpart + 1;
            }            

          } // cierra dis
           
          i = grid.ll[i];

        } //fin lazo particulas del grid
      } //fin izz
    } //fin iyy
  } //fin ixx

  return;
}

bool point_inside(float dot, float rlong_2)
{

  if(fabs(dot)>rlong_2){

      return false;

  }else{

    return true;

  }

}
#endif
