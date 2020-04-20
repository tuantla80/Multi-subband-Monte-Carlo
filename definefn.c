#include "nanowire.h"
#include "region.h"

static Param P ;
static Dimen D ;
static Energy E ;
static Density n ;

Param *PGetParam() { 
  return(&P) ; 
}
Param  GetParam() {  
  return(P) ; 
}

Dimen *PGetDimen() { 
  return(&D) ; 
}
Dimen  GetDimen() { 
  return(D) ; 
}

Energy *PGetEnergy() { 
  return(&E) ; 
}
Energy  GetEnergy() { 
  return(E) ; 
}

Density *PGetDensity() { 
  return(&n) ; 
}
Density  GetDensity() { 
  return(n) ; 
}

// EFs

double GetEF() {
  E.F = 0.0 ; // Default
  return(E.F) ;
}

double GetEFnS() {
  return(E.FnS) ; // To be set
}

double GetEFnD() {
  return(E.FnD) ; // To be set
}

double GetESi() { 
  return(P.e_si) ; 
}
double GetEOx() { 
  return(P.e_ox) ; 
}
double GetEc_si() { 
  return(E.c_si) ; 
}
double GetEc_ox() { 
  return(E.c_ox) ; 
}
double GetKT() { 
  return(E.kT) ; 
}
int GetTransSize() { 
  return(P.trans_size) ; 
}

// Phi's

static double Phi_Source ;
static double Phi_Drain ;
static double Phi_Gate ;

void SetPhi_Contact() {

  void assign_built_in_potential() ;
  assign_built_in_potential() ;

  //printf("\n ********** From SetPhi_Contact() *********");

  double Vd = GetDrainVoltage() ;
  //printf("\n Drain Voltage Vd = %f", Vd);
  
  double Vg = GetGateVoltage() ;
  //printf("\n Gate Voltage Vg  = %f", Vg);

  Phi_Source = E.biS ;
  //printf("\n Phi_Source       = %f",Phi_Source);
  
  Phi_Drain  = E.biD + Vd ;
  //printf("\n Phi_Drain        = %f",Phi_Drain);
  
  Phi_Gate   = Vg - P.offset_gate ; // for n-mos;  //Phi_Gate   = -Vg + P.offset_gate + Vd ; // for p-mos with CB_PEMT mode
  //printf("\n Phi_Gate         = %f",Phi_Gate);
  
  //printf("\n ********** Ket thuc SetPhi_Contact() *********");//  getchar();
}

double GetPhi_Source() { 
  return(Phi_Source) ; 
}
double GetPhi_Drain() { 
  return(Phi_Drain) ; 
}
double GetPhi_Gate() {
  return(Phi_Gate) ; 
}

// Ec

double GetEc(int i, double y, double z)
{
  PoiNum N = GetPoiNum() ;
  double *Y = GetPY() ;
  double *Z = GetPZ() ;

  double ya = Y[N.ya]-tiny_06 ;
  double yb = Y[N.yb]+tiny_06 ;
  double za = Z[N.za]-tiny_06 ;
  double zb = Z[N.zb]+tiny_06 ;

  if ( y<ya || y>yb || z<za || z>zb ) return(E.c_ox) ;
  else  return(E.c_si) ;

  return(0.0) ; // dummy
}

// Poisson Linearized

void SetPoissonLinearized(int n) {
  P.poi_lin = n ;
}

int GetPoissonLinearized() {
  return(P.poi_lin) ;
}



