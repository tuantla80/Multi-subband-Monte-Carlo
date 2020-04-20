#include "nanowire.h"

static int PRINT_MidPotFlag ;
static int PRINT_Pot3dFlag ;

void SetPRINT_MidPotFlag(int n) {
  PRINT_MidPotFlag = n ;
}

void SetPRINT_Pot3dFlag(int n ) {
  PRINT_Pot3dFlag = n ;
}

int GetPRINT_MidPotFlag() {
  return(PRINT_MidPotFlag) ;
}

int GetPRINT_Pot3dFlag() {
  return(PRINT_Pot3dFlag) ;
}

void print_output()
{
  // print outputs

  char fn[100] ;
  double vd = GetDrainVoltage() ;
  double vg = GetGateVoltage() ;

  if ( GetPRINT_Pot3dFlag() ) sprintf(fn,"pot3d.r-Vd%.2lf-Vg%.2lf",vd,vg) ;
  else                        sprintf(fn,"pot3d.r") ;
  print_potential_plain(fn) ;
 
  if ( GetPRINT_MidPotFlag() ) {
    sprintf(fn,"mpot.r") ;
    print_midline_potential_plain(fn) ; 
  }
}
