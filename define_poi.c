#include "nrutil.h"
#include "petscksp.h"
#include "nanowire.h"

static PoiNum N ;
static PoiMem M ;
static PoiParam CP ;
static Region R ;

PoiNum *PGetPoiNum() { return(&N) ; }
PoiNum  GetPoiNum() { return(N) ; }
PoiMem *PGetPoiMem() { return(&M) ; }
PoiParam *PGetPoiParam() { return(&CP) ; }
PoiParam  GetPoiParam() { return(CP) ; }
Region *PGetR_() { return(&R) ; }
Region  GetR_() { return(R) ; }
double ***GetPot() { return(M.Phi) ; }
double ***GetPhiKth() { return(M.phi) ; }

double *GetPX() { return(M.X) ; }
double *GetPY() { return(M.Y) ; }
double *GetPZ() { return(M.Z) ; }
double *GetHx() { return(M.hx) ; }
double *GetHy() { return(M.hy) ; }
double *GetHz() { return(M.hz) ; }

int GetPx0() {
  return(N.x0) ;
}

int GetPx1() {
  return(N.x1) ;
}

int GetPy0() {
  return(N.y0) ;
}

int GetPy1() {
  return(N.y1) ;
}

int GetPz0() {
  return(N.z0) ;
}

int GetPz1() {
  return(N.z1) ;
}

int GetPya() {
  return(N.ya) ;
}

int GetPyb() {
  return(N.yb) ;
}

int GetPza() {
  return(N.za) ;
}

int GetPzb() {
  return(N.zb) ;
}

double ***GetN3d() {
  return(M.n3d) ;
}

double ***GetP3d() {
  return(M.p3d) ;
}

static int *n0s,*n1s ;

void SetLocalN(int np, int *I0s, int *I1s)
{
  int i ;
  int x0 = GetPx0() ;

  n0s = ivector(0,np-1) ;
  n1s = ivector(0,np-1) ;

  for ( i=0 ; i<np ; i++ ) {
    n0s[i] = I0s[i] - x0 + 1 ;
    n1s[i] = I1s[i] - x0 + 1 ;
  }
}

void GetLocalN(int node,int *n0, int *n1)
{
  *n0 = n0s[node] ;
  *n1 = n1s[node] ;
}

PCType GetPCType() { 

  if ( !strcmp(CP.pctype,"Jacobi") ) return(PCJACOBI) ;
  else if ( !strcmp(CP.pctype,"Block_Jacobi") ) return(PCBJACOBI) ;
  else if ( !strcmp(CP.pctype,"SOR") ) return(PCSOR) ;
  else if ( !strcmp(CP.pctype,"SOR_Eis") ) return(PCEISENSTAT) ;
  else if ( !strcmp(CP.pctype,"ICC") ) return(PCICC) ;
  else if ( !strcmp(CP.pctype,"ILU") ) return(PCILU) ;
  else if ( !strcmp(CP.pctype,"ASM") ) return(PCASM) ;
  else if ( !strcmp(CP.pctype,"LinearSolv") ) return(PCKSP) ;
  else if ( !strcmp(CP.pctype,"Combinations") ) return(PCCOMPOSITE) ;
  else {
    report_error("PCType Pattern Does Not Match...") ;
    return(PCBJACOBI) ;
  }
}

KSPType GetKSPType() { 
  
  if ( !strcmp(CP.ksptype,"Richardson") ) return(KSPRICHARDSON) ;
  else if ( !strcmp(CP.ksptype,"Chebychev") ) return(KSPCHEBYCHEV) ;
  else if ( !strcmp(CP.ksptype,"CG") ) return(KSPCG) ;
  else if ( !strcmp(CP.ksptype,"BiCG") ) return(KSPBICG) ;
  else if ( !strcmp(CP.ksptype,"GMRES") ) return(KSPGMRES) ;
  else if ( !strcmp(CP.ksptype,"BiCGSTAB") ) return(KSPBCGS) ;
  else if ( !strcmp(CP.ksptype,"CGS") ) return(KSPCGS) ;
  else if ( !strcmp(CP.ksptype,"QMR1") ) return(KSPTFQMR) ;
  else if ( !strcmp(CP.ksptype,"QMR2") ) return(KSPTCQMR) ;
  else if ( !strcmp(CP.ksptype,"CR") ) return(KSPCR) ;
  else if ( !strcmp(CP.ksptype,"LSQR") ) return(KSPLSQR) ;
  else if ( !strcmp(CP.ksptype,"PreOnly") ) return(KSPPREONLY) ;
  else {
    report_error("KSPType Pattern Does Not Match...") ;
    return(KSPGMRES) ; 
  }
}
