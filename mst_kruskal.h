#ifndef KRUSKAL 
#define KRUSKAL 

int  Root(int i, int *Padre);
void DLU(int id, int pre, std::vector<std::vector<int> > &adj, int *Padre, int *Rank);
void DL(std::vector<std::vector<int> > &vec, int *Padre, int *Rank);
void Delete_Branch(int i, std::vector<std::vector<int> > &vec);
bool DFS(int i, std::vector<std::vector<int> > &vec, int cut);
void Podado(int level, std::vector<std::vector<int> > &vec);
void Kruskal(int *Padre, int *Rank, std::vector<std::pair<float,std::pair<int,int> > > &edges, std::vector<std::vector<int> > &adjacency_list);

#endif
