#include <assert.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <functional>
#include <algorithm>
#include <vector>

#include "cosmoparam.h"
#include "variables.h"
#include "leesloan.h"
#include "grid.h"
#include "voronoi.h"
#include "voro++.hh"

void locate(double xx[], unsigned long n, float x, unsigned long *j)
{
  unsigned long ju,jm,jl;
  int ascnd;

  jl=0;
  ju=n+1;
  ascnd=(xx[n-1] > xx[0]);
  while(ju-jl > 1)
  {
    jm=(ju+jl) >> 1;
    if((x > xx[jm-1])  == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  *j=jl;
}

bool dist_segment(double fof, int idg, float x, float y, float z)
{
  int i, ibox;
  int ixc, iyc, izc;
  int ixx, iyy, izz;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int xi, yi, zi;
  double fac;
  int itabla,imintabla;
  float rscale,vl,v12;
  float rm12,rmmin12;
  float r[3],dl,dl1,dl2,d12;
  float rint12,rint121,rint122;
  float d1,d2,tita12;
  float red1,delta1,alfa1;
  float tridif,numera,alfa12;
  long unsigned indx;

  r[0] = x + pmin[0];
  r[1] = y + pmin[1];
  r[2] = z + pmin[2];

  d1 = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);                                                     // distancia cosmologica
  locate(dis2red,NTABLADIS,d1,&indx);                                                         
  red1 = (rt[indx]-rt[indx-1])/(dis2red[indx]-dis2red[indx-1])*(d1-dis2red[indx-1])+rt[indx-1]; // redshift
  
  delta1 = atan2(z,sqrt(x*x+y*y));

  if(z>=0.0)
    delta1 += M_PI;
    
  if(x>0)
  {
    if(y>=0)
      alfa1 = atan2(y,x);
    else
      alfa1 = 2*M_PI + atan2(y,x);
  }else if(x<0){
      alfa1 = M_PI + atan2(y,x);
  }else{    
      alfa1 = 0.5*M_PI*(y/fabs(y));
  }

  dl1 = d1*(1.0+red1); // distancia luminosidad

  fac = (double)grid.ngrid/(double)cp.lbox;

  ixc  = (int)((double)x*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)((double)y*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)((double)z*fac);
  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0) ixci = 0;
  if( iyci < 0) iyci = 0;
  if( izci < 0) izci = 0;
  if( ixcf >= (int)grid.ngrid ) ixcf = (int)grid.ngrid - 1;
  if( iycf >= (int)grid.ngrid ) iycf = (int)grid.ngrid - 1;
  if( izcf >= (int)grid.ngrid ) izcf = (int)grid.ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++)
  {
    xi = ixx;
    #ifdef PERIODIC
    if(xi >= ngrid) xi -= grid.ngrid;
    if(xi < 0)      xi += grid.ngrid;
    #endif

    for( iyy = iyci ; iyy <= iycf ; iyy++)
    {

      yi = iyy;
      #ifdef PERIODIC
      if(yi >= ngrid) yi -= grind.ngrid;
      if(yi < 0)      yi += grind.ngrid;
      #endif

      for( izz = izci ; izz <= izcf ; izz++)
      {
        zi = izz;
        #ifdef PERIODIC
        if(zi >= ngrid) zi -= grid.ngrid;
        if(zi < 0)      zi += grid.ngrid;
        #endif

        ibox = (xi * grid.ngrid + yi) * grid.ngrid + zi;

        i = grid.llirst[ibox];

        while(i != -1)
        {

         if(P[i].sub==idg)
         {

          d2   = P[i].Dis;                    // distancia cosmologica
          dl2  = d2*(1.0+gal[i].gal.red);     // distancia luminosidad 

          rm12 = rmaplim-25.0-5.0*log10((dl1+dl2)/2.);
          rmmin12 = rmapmin-25.0-5.0*log10((dl1+dl2)/2.);

          // CALCULA LA INTEGRAL DEL FACTOR DE ESCALA
          itabla=(int)((rm12-mt2min)/dmt2);
          imintabla=(int)((rmmin12-mt1min)/dmt1);

          if(itabla<0 || itabla >= NTABLA)
          {
            fprintf(stdout,"XXXX itabla  = %d\n ",itabla);
            fflush(stdout);
            exit(-33);
          }

          if(imintabla<0 || imintabla >= NTABLA)
          {
            fprintf(stdout,"XXXX imintabla= %d\n ",imintabla);
            fflush(stdout);
            exit(-33);
          }

          rint12  = (rinttabla[imintabla][itabla+1]-rinttabla[imintabla][itabla])/dmt2;
          rint121 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla][itabla];

          rint12  = (rinttabla[imintabla+1][itabla+1]-rinttabla[imintabla+1][itabla])/dmt2;
          rint122 = rint12*(rm12-itabla*dmt2-mt2min)+rinttabla[imintabla+1][itabla];

          rint12 = (rint122-rint121)/dmt1;
          rint12 = rint12*(rmmin12-imintabla*dmt1-mt1min)+rint121;

          rscale = pow(rint12/rintlim,-0.33333);
          dl = fof*rscale;
          vl = v0*rscale;

          alfa12 = gal[i].gal.alfa-alfa1;
          tita12 = sin(gal[i].gal.delta)*sin(delta1)+cos(gal[i].gal.delta)*cos(delta1)*cos(alfa12);

          numera = (cos(gal[i].gal.delta)*sin(alfa12));
          numera *= numera;

          tridif = cos(delta1)*sin(gal[i].gal.delta)-sin(delta1)*cos(gal[i].gal.delta)*cos(alfa12);
          tridif *= tridif;

          numera += tridif;
          numera = sqrt(numera);

          tita12 = atan2(numera,tita12);

          d12 = sin(tita12/2.0)*(d1+d2);
          v12 = fabs(d1-d2);

          if((d12<dl) && (v12<vl))
            return true;

          }

         i = grid.ll[i];

       }

      }

    }

  }

  return false;
}

void Voronoi_Grupos(double fof, std::vector<std::pair<float,std::pair<int,int> > > &edges)
{

  int i, j, idv, id, Ngrid, Tid, itera, N_threads;
  bool xbool, ybool, zbool;
  double x_min, y_min, z_min;	  
  double x_max, y_max, z_max;	  
  double dx,dy,dz,r,d0,r0,frac;
  std::vector<int>  vec;
  voro::voronoicell_neighbor cell;
  float v12,d12;
  float d1,d2,tita12;
  float tridif,numera,alfa12;

  Ngrid = (int)pow((float)ngrup/5.0,1.0/3.0);
  #ifdef PERIODIC
    xbool = ybool = zbool = true;
  #else
    xbool = ybool = zbool = false;
  #endif
  r0 = 1.0e-10; // uso un pequeño gap
  x_min = y_min = z_min = 0.0-r0;	  
  x_max = y_max = z_max = cp.lbox+r0;	  

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  N_threads = NTHREADS;

  #ifdef LEN_FOF_MERCHAN
    d0 = 3.0/(4.0*M_PI*(fof+1)*(0.4*log(10.0)*flfia*rintlim));
  #else
    d0 = 3.0*cp.vol/(4.0*M_PI*(cp.npart*(fof+1.0))); // cbrt(vol/(npart*(fof+1)))
  #endif
  d0 = cbrt(d0);

  r0 = v0;
  if(d0>r0)
  {
    r0=d0;
    fprintf(stdout,"d0 > v0, %f > %f \n",d0,v0);
  }else{
    fprintf(stdout,"v0 > d0, %f > %f \n",v0,d0); 
  }

  /// INICIA LOS PARAMETROS DE LA GRILLA ///
  grid.nobj = cp.npart;
  #ifdef LIM_VOL
    grid.ngrid = (int)(cp.lbox/(r0*pow(intfl(MAGMENOSINF,rmaplim-25.0-5.0*log10(cp.dlummax))/rintlim,-1./3.)));
  #else
    grid.ngrid = (int)(cp.lbox/(r0*pow(intfl(rmapmin-25.0-5.0*log10(cp.dlummax),rmaplim-25.0-5.0*log10(cp.dlummax))/rintlim,-1./3.)));
  #endif
  grid.step = 1;

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  //////////////////////////////////////////
  grid_init();
  grid_build();

  // Creo el container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,Ngrid,Ngrid,Ngrid,xbool,ybool,zbool,8);

  for(i=0;i<ngrup; i++)
    con.put(i,Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2]);	  

  assert(ngrup==con.total_particles());

  voro::c_loop_all clo(con);
  std::vector<std::vector<std::pair<float,std::pair<int,int> > > > lados(N_threads);
  
  if(clo.start()) do if(con.compute_cell(cell,clo))
  {

    id = clo.pid();
    cell.neighbors(vec);

    N_threads = (int)vec.size()>N_threads ? N_threads : (int)vec.size();

    #ifdef GAL_LUM

      d1   = P[Gr[id].id].Dis;                    // distancia cosmologica

    #else

      d1 = sqrt(Gr[id].Pos[0]*Gr[id].Pos[0]+ \
                Gr[id].Pos[1]*Gr[id].Pos[1]+ \
                Gr[id].Pos[2]*Gr[id].Pos[2]); // distancia cosmologica      

    #endif

    #pragma omp parallel for num_threads(N_threads) schedule(static) \
    default(none) private(Tid,j,idv,dx,dy,dz,r,frac,itera, \
    d2,tridif,numera,alfa12,tita12,d12,v12) \
    shared(id,vec,d0,cp,Gr,grid,lados,d1,P,gal)
    for(j=0; j<(int)vec.size(); j++)
    {

      Tid = omp_get_thread_num();
      idv = vec[j];

      if(idv<0) continue;

      if(Gr[id].save!=Gr[idv].save) continue;

      if(Gr[id].id>Gr[idv].id)
      {
        dx = Gr[idv].Pos[0] - Gr[id].Pos[0];
        dy = Gr[idv].Pos[1] - Gr[id].Pos[1];
        dz = Gr[idv].Pos[2] - Gr[id].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;
  
        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif

        r = sqrt(dx*dx+dy*dy+dz*dz);

        itera = (int)(r/d0);
  
        while(itera>0)
        {
          frac = (double)itera*(d0/r);
          if(!dist_segment(d0,Gr[id].save,Gr[id].Pos[0]+frac*dx,Gr[id].Pos[1]+frac*dy,Gr[id].Pos[2]+frac*dz))
            break;
          else  
            itera--; 
        }

        if(itera!=0) continue;

        #ifdef GAL_LUM
    
          d2   = P[Gr[idv].id].Dis;                    // distancia cosmologica

        #else

          d2 = sqrt(Gr[idv].Pos[0]*Gr[idv].Pos[0]+ \
                    Gr[idv].Pos[1]*Gr[idv].Pos[1]+ \
                    Gr[idv].Pos[2]*Gr[idv].Pos[2]); // distancia cosmologica      

        #endif

        #ifdef GAL_LUM

          alfa12 = gal[Gr[idv].id].gal.alfa-gal[Gr[id].id].gal.alfa;

          tita12 = sin(gal[Gr[idv].id].gal.delta)*sin(gal[Gr[id].id].gal.delta)+ \
                   cos(gal[Gr[idv].id].gal.delta)*cos(gal[Gr[id].id].gal.delta)*cos(alfa12);

          numera = (cos(gal[Gr[idv].id].gal.delta)*sin(alfa12));
          numera *= numera;

          tridif = cos(gal[Gr[id].id].gal.delta)*sin(gal[Gr[idv].id].gal.delta)- \
          sin(gal[Gr[id].id].gal.delta)*cos(gal[Gr[idv].id].gal.delta)*cos(alfa12);
          tridif *= tridif;

        #else

          alfa12 = gal[Gr[idv].id].gal.alfa-gal[Gr[id].id].gal.alfa;

          tita12 = sin(gal[Gr[idv].id].gal.delta)*sin(gal[Gr[id].id].gal.delta)+ \
                   cos(gal[Gr[idv].id].gal.delta)*cos(gal[Gr[id].id].gal.delta)*cos(alfa12);

          numera = (cos(gal[Gr[idv].id].gal.delta)*sin(alfa12));
          numera *= numera;

          tridif = cos(gal[Gr[id].id].gal.delta)*sin(gal[Gr[idv].id].gal.delta)- \
          sin(gal[Gr[id].id].gal.delta)*cos(gal[Gr[idv].id].gal.delta)*cos(alfa12);
          tridif *= tridif;

          

        #endif

        numera += tridif;
        numera = sqrt(numera);

        tita12 = atan2(numera,tita12);

        d12 = sin(tita12/2.0)*(d1+d2);
        v12 = fabs(d1-d2);          
        r = d12*d12+v12*v12;

        #ifdef GAL_LUM
          // MAGsun_r = 4.71         // mag_r Hill et al 2010 - https://arxiv.org/pdf/1002.3788.pdf
          r = -(pow(10.0,-0.4*(Gr[id].mr-4.71))*pow(10.0,-0.4*(Gr[idv].mr-4.71)))/r;
        #else
          r = -((float)Gr[id].NumPart*(float)Gr[idv].NumPart)/r;
        #endif

        lados[Tid].push_back(std::make_pair((float)r,std::make_pair(idv,id)));
      }

    }

    vec.clear();

   }while(clo.inc());
  
  con.clear();
  grid_free();
  #ifndef CALCULA_MEDIA
  free(P);
  #endif

  for(i=0;i<N_threads;i++)
  {
    edges.insert(edges.end(),lados[i].begin(),lados[i].end());
    lados[i].clear();
  }

  return;
}


