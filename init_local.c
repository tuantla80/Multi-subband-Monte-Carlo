#include "nrutil.h"
#include "nanowire.h"
#include "petscksp.h"

void initialize_general_constants()
{
  //printf("\n\n At initialize_general_constants()");
  Param  *P = PGetParam() ;
  Energy *E = PGetEnergy() ;

  P->e_si *= e_0 ;
  //printf("\n Dielectric of Si (exact value)P->e_si    = %le",P->e_si);
  
  P->e_ox *= e_0 ;
  //printf("\n Dielectric of Oxide (exact value)P->e_ox = %le",P->e_ox);

  
  E->kT = kB*(P->Temp)/q0 ;
  //printf("\n Thermal energy E->kT             = %f", E->kT);
  
  E->c_si = 0.5*E->g_si ; 
  //printf("\n Conduction band of Si E->c_si    = %f",E->c_si);
  
  E->c_ox = 0.5*E->g_ox ;
  //printf("\n Conduction band of Oxide E->c_ox = %f",E->c_ox);
}

#define TINY 1.0e-10

void set_position_x_direction()
{
  //printf("\n\n At set_position_x_direction()");
  double x1,x2,x3,hx1,hx2,hx3,x,*X,*hx ;
  int k,K,N1,N2,N3,Nx ;
  Dimen D ;
  PoiNum *N ;
  PoiMem *M ;
  Param P ;
  void divide_region() ;
  void SetHamilNx() ;

  D = GetDimen() ;
  N = PGetPoiNum() ;
  M = PGetPoiMem() ;
  P = GetParam() ;

  x1 = D.Lsrc ;
  x2 = x1 + D.Lchannel ;
  x3 = x2 + D.Lsrc ;

  N1 = N->Mx_src ;
  N2 = N->Mx_channel ;
  N3 = N->Mx_src ;
  Nx = N1+N2+N3+1 ;

  X = M->X = dvector(0,Nx-1) ;
  hx = M->hx = dvector(-1,Nx-2) ;

  N->x0 = 0 ;
  N->xa = N->x0+N1 ;
  N->xb = N->xa+N2 ;
  N->x1 = N->xb+N3 ;
  N->x = Nx ;
  //printf("\n N->x0 =%d, N->xa =%d, N->xb =%d, N->x1 =%d, N->x =%d",N->x0,N->xa,N->xb,N->x1,N->x);  

  hx1 = x1/(double)N1 ;
  hx2 = (x2-x1)/(double)N2 ;
  hx3 = (x3-x2)/(double)N3 ;
  //printf("\n hx1 = %f, hx2 = %f, hx3 = %f",hx1,hx2,hx3);

  X[0] = x = 0.0 ;
  K = 1 ;
  divide_region(hx1,x1,X,&K,&x) ;
  divide_region(hx2,x2,X,&K,&x) ;
  divide_region(hx3,x3,X,&K,&x) ;
  
  for ( k=1 ; k<Nx ; k++ ) hx[k-1] = X[k]-X[k-1] ;
  hx[-1] = hx[0] ;

  // Gate

  for ( k=1 ; k<Nx ; k++ ) 
    if ( X[k] > 0.5*(x3-D.Lgate)-TINY ) break ;
  N->xgate0 = k ;
  
  for ( k=1 ; k<Nx ; k++ )
    if ( X[k] > 0.5*(x3+D.Lgate)+TINY ) break ;
  N->xgate1 = k-1 ;

}

void set_position_y_direction()
{
  //printf("\n\n At set_position_y_direction()");
  double y1,y2,y3,hy1,hy2,hy3,y,*Y,*hy ;
  int k,K,N1,N2,N3,Ny ;
  Dimen   D = GetDimen() ;
  Dimen *PD = PGetDimen() ;
  Param   P = GetParam() ;
  PoiNum *N = PGetPoiNum() ;
  PoiMem *M = PGetPoiMem() ;
  void divide_region() ;

  double dev_len_y = 2*D.Wox + D.Wsi ;
  double dev_len_z = D.Tox + D.Tsi + D.Tbox ;

  y1 = D.Wox ;
  y2 = y1 + D.Wsi ;
  y3 = y2 + D.Wox ;

  N1 = N->My_ox ;
  N2 = N->My_si ;
  N3 = N->My_ox ;
  Ny = N1+N2+N3+1 ;

  Y = M->Y = dvector(0,Ny-1) ;
  hy = M->hy = dvector(-1,Ny-2) ;

  N->y0   = 0 ;
  N->ya   = N->y0+N1 ;
  N->yb   = N->ya+N2 ;
  N->y1   = N->yb+N3 ;
  N->y    = Ny ;
  //printf("\n N->y0 =%d, N->ya =%d, N->yb =%d, N->y1 =%d, N->y =%d",N->y0,N->ya,N->yb,N->y1,N->y);

  hy1 = y1/(double)N1 ;
  hy2 = (y2-y1)/(double)N2 ;
  hy3 = (y3-y2)/(double)N3 ;
  //printf("\n hy1 = %f, hy2 = %f, hy3 = %f",hy1,hy2,hy3);

  Y[0] = y = 0.0 ;
  K = 1 ;
  divide_region(hy1,y1,Y,&K,&y) ;
  divide_region(hy2,y2,Y,&K,&y) ;
  divide_region(hy3,y3,Y,&K,&y) ;
  
  for ( k=1 ; k<Ny ; k++ ) hy[k-1] = Y[k]-Y[k-1] ;

  // Determine the location of the gate extended in the y direction

  for ( k=1 ; k<Ny ; k++ ) 
    if ( fabs(Y[k]-D.Wgate)<0.1*hy[k] ) break ;

  N->ygate0 = k ;
  N->ygate1 = N->y0+N->y1-N->ygate0 ;
}

void set_position_z_direction()
{ 
  //printf("\n\n At set_position_z_direction()");   
  double z1,z2,z3,hz1,hz2,hz3,z,R1,R2,*Z,*hz ;
  int k,K,N1,N2,N3,Nz ;
  Dimen D ;
  PoiNum *N ;
  PoiMem *M ;
  void divide_region() ;
  
  D = GetDimen() ;
  N = PGetPoiNum() ;
  M = PGetPoiMem() ;

  z1 = D.Tox ;
  z2 = z1 + D.Tsi  ;
  z3 = z2 + D.Tbox ;

  N1 = N->Mz_ox ;
  N2 = N->Mz_si ;
  N3 = N->Mz_box ;
  Nz = N1+N2+N3+1 ;

  Z = M->Z = dvector(0,Nz-1) ;
  hz = M->hz = dvector(-1,Nz-2) ;

  N->z0   = 0 ;
  N->za   = N->z0+N1 ;
  N->zb   = N->za+N2 ;
  N->z1   = N->zb+N3 ;
  N->z    = Nz ;
  //printf("\n N->z0 =%d, N->za =%d, N->zb =%d, N->z1 =%d, N->z =%d",N->z0,N->za,N->zb,N->z1,N->z);

  hz1 = z1/(double)N1 ;
  hz2 = (z2-z1)/(double)N2 ;
  hz3 = (z3-z2)/(double)N3 ;
  //printf("\n hz1 = %f, hz2 = %f, hz3 = %f",hz1,hz2,hz3);

  Z[0] = z = 0.0 ;
  K = 1 ;
  divide_region(hz1,z1,Z,&K,&z) ;
  divide_region(hz2,z2,Z,&K,&z) ;
  divide_region(hz3,z3,Z,&K,&z) ;
  
  for ( k=1 ; k<Nz ; k++ ) hz[k-1] = Z[k]-Z[k-1] ;

  // Determine the location of the gate extended in the z direction

  for ( k=1 ; k<Nz ; k++ ) 
    if ( Z[k] > D.Tgate+TINY ) break ;
  N->zgate = k-1 ;
}

void set_regions()
{
  int count ;
  Region *R ;
  
  R = PGetR_() ;

  count = 0 ;

  R->_CHANNEL_SILICON = ++count ;
  R->_JUNCTION_SILICON = ++count ;
  R->_OXIDE = ++count ;

  R->_SOURCE_CONTACT = ++count  ;
  R->_DRAIN_CONTACT = ++count ;
  R->_GATE_CONTACT = ++count ;

  R->_SOURCE_WALL = ++count ;
  R->_DRAIN_WALL = ++count ;
  R->_LEFT_WALL = ++count ;
  R->_RIGHT_WALL = ++count ;
  R->_TOP_WALL = ++count ;
  R->_BOTT_WALL = ++count ;

  R->_LEFT_INTERFACE = ++count ;
  R->_RIGHT_INTERFACE = ++count ;
  R->_TOP_INTERFACE = ++count ;
  R->_BOTT_INTERFACE = ++count ;

  R->_CORNERS = ++count ;
}

#define DIV_REG_EPS 1.0e-09

void divide_region(h,xlimit,Pos,K,xv)
     int *K ;
     double h,xlimit,*xv,*Pos ;
{
  int k,k0 ;
  double x ;

  x = *xv ;
  k0 = *K ;

  for ( k=k0 ; ; k++ ) {
    x += h ;
    if ( x < xlimit*(1.0+DIV_REG_EPS) ) Pos[k] = x ;
    else break ;
  }
  
  x -= h ;
  *K = k ;
  *xv = x ;
}

#undef DIV_REG_EPS
