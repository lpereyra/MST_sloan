#ifndef VORONOI 
#define VORONOI 

void locate(double xx[], unsigned long n, float x, unsigned long *j);
bool dist_segment(double fof, int idg, float x, float y, float z);
void Voronoi_Grupos(double fof, std::vector<std::pair<float,std::pair<int,int> > > &edges);

#endif

