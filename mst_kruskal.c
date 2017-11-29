#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <functional>
#include <algorithm>

#include "variables.h"
#include "mst_kruskal.h"

int Root(int i, int *Padre)
{

 if(i != Padre[i])
   Padre[i] = Root(Padre[i],Padre);

 return Padre[i];
}

void DLU(int id, int pre, std::vector<std::vector<int> > &adj, int *Padre, int *Rank)
{
  int k,idv;

  Padre[id] = pre;
  Rank[id] = (int)adj[id].size();

  for(k=0;k<Rank[id];k++)
  {
    idv = adj[id][k];
    if(idv != pre)
      DLU(idv,id,adj,Padre,Rank);
  }

  adj[id].clear();

  return;

}
 
void DL(std::vector<std::vector<int> > &vec, int *Padre, int *Rank)
{
  int i;
  
  memset(Padre,-1,ngrup*sizeof(int));
  memset(Rank,0,ngrup*sizeof(int));
  
  for(i=0;i<ngrup;i++)     
    if(vec[i].size()==1)
      DLU(i,Padre[i],vec,Padre,Rank);

  return;
}

bool DFS(int i, std::vector<std::vector<int> > &vec, int cut)
{
  int j,k,id,idv;

  j = 1;
  id = i;
  idv = vec[id][0];

  while(vec[idv].size()==2)
  {

    #ifdef BRANCH_SURVIVE
      #ifdef GAL_LUM
      if(Gr[idv].mr<mr_survive)
      #else
      if(Gr[idv].NumPart>N_part_survive)
      #endif
      {
        j=cut+1;
      }
    #endif    

    k = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];
    id = idv;
    idv = k;
    j++;        

    if(j>cut) break;
  }

  if(j<=cut)
    return true;
  else 
    return false;

}

void Delete_Branch(int i, std::vector<std::vector<int> > &vec)
{
  int k,id,idv;

  id = i;
  idv = vec[id][0];
  vec[id].clear();
 
  while(vec[idv].size()==2)
  {

    k = vec[idv][0]==id ? vec[idv][1] : vec[idv][0];

    vec[idv].clear();

    id = idv;
    idv = k;
  }

  if(vec[idv].size()==1)
    vec[idv].clear();

  return;
}

void Podado(int level, std::vector<std::vector<int> > &vec)
{
  int i,k,cut,N_threads;
  bool itera;
  std::vector<int> aux;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  N_threads = NTHREADS;

  fprintf(stdout,"level podado %d\n",level);
  fflush(stdout);

  cut = 1;

  do{

    itera = false;

    #ifdef GAL_LUM
      #ifdef BRANCH_SURVIVE
        #pragma omp parallel for num_threads(N_threads) schedule(dynamic) default(none) \
        shared(Gr,vec,ngrup,cut,itera,mr_survive) 
      #else
        #pragma omp parallel for num_threads(N_threads) schedule(dynamic) default(none) \
        shared(Gr,vec,ngrup,cut,itera) 
      #endif
    #else
      #ifdef BRANCH_SURVIVE
        #pragma omp parallel for num_threads(N_threads) schedule(dynamic) default(none) \
        shared(Gr,vec,ngrup,cut,itera,N_part_survive) 
      #else
        #pragma pel for num_threads(N_threads) schedule(dynamic) default(none) \
        shared(Gr,vec,ngrup,cut,itera) 
      #endif 
    #endif 
    for(i=0;i<ngrup;i++)
    {
    
      if(vec[i].size()==1)
      {
        
        #ifdef BRANCH_SURVIVE
          #ifdef GAL_LUM
            if(Gr[i].mr<mr_survive) continue;
          #else
            if(Gr[i].NumPart>N_part_survive) continue;
          #endif
        #endif

        if(DFS(i,vec,cut))       
        {
          Delete_Branch(i,vec);
          itera = true;
        }
    
      }

    }
    
    if(itera==false) break;

    #pragma omp parallel for num_threads(N_threads) schedule(dynamic) default(none) \
    private(aux,k) shared(ngrup,vec) 
    for(i=0;i<ngrup;i++)
    {
      if(vec[i].size()>2)
      {
        for(k=0;k<(int)vec[i].size();k++)
        {
          if(!vec[vec[i][k]].empty())
            aux.push_back(vec[i][k]);
        }
  
        swap(aux,vec[i]);
        aux.clear();

      } 
    }

    if(cut==level) break;

    cut++;

  }while(1);
 
  fprintf(stdout,"Termina el podado\n");
  fflush(stdout);

  return;
  
}

void Kruskal(int *Padre, int *Rank, std::vector<std::pair<float,std::pair<int,int> > > &edges, \
std::vector<std::vector<int> > &adjacency_list)
{
  int j,k,id,idv;

  sort(edges.begin(),edges.end(), std::greater<std::pair<float,std::pair<int,int> > >());

  while(!edges.empty()) 
  {
    j = edges.back().second.first;
    k = edges.back().second.second;
    edges.pop_back();

    id = Root(j,Padre);
    idv = Root(k,Padre);

    if(id!=idv)
    {
      if(Rank[id] < Rank[idv])
        Padre[id] = idv;
      else if(Rank[idv] < Rank[id])
        Padre[idv] = id;
      else
      {
        Padre[idv] = id;
        Rank[id]++;
      }

      adjacency_list[j].push_back(k);
      adjacency_list[k].push_back(j);

    } 

  }

}
