#include "nanowire.h"
#include "region.h"

#define SIGN_CONVFACTOR_DOPINGDENSITY sign*conv_factor*doping_density

static double SourceDoping ;
static double ChannelDoping ;
static double DrainDoping ;

void SetDopingDensityFor(int convflag, char region, char doping_type, double doping_density)
{
  double sign ;
  double conv_factor ;
  
  if ( convflag==1 ) conv_factor = 1.0e06 ;
  else               conv_factor = 1.0 ;
  
  doping_density = fabs(doping_density) ; // necessary for ldos-for-EF routines

  switch ( doping_type ) {
  case 'n' : sign = -1.0 ; break ;
  case 'p' : sign = 1.0  ; break ;
  case 'i' : sign = 0.0  ; break ;
  }

  switch ( region ) {

  case 's' : // Source
    SourceDoping = SIGN_CONVFACTOR_DOPINGDENSITY ;
    break ;

  case 'c' : // Channel
    ChannelDoping = SIGN_CONVFACTOR_DOPINGDENSITY ;
    break ;

  case 'd' : // Drain
    DrainDoping = SIGN_CONVFACTOR_DOPINGDENSITY ;
    break ;
  }
}

double GetDopingDensityFor(char region) 
{
  switch ( region ) {
  case 's' : return(SourceDoping) ;
  case 'c' : return(ChannelDoping) ;
  case 'd' : return(DrainDoping) ;
  }

  return(0.0) ;
}

double GetDoping(int p)
{
  int region ;
  int i,j,k ;
  Region R ;
  PoiNum N ;

  double NL ;
  double *X ;

  N = GetPoiNum() ;
  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;
  X = GetPX() ;

  if ( JUNCTION_SILICON ) {
    if      ( i<N.xa ) return(SourceDoping) ;
    else if ( i>N.xb ) return(DrainDoping) ;
    else 
      report_error("Problem in GetDoping(), junction silicon") ;
  }
  else if ( CHANNEL_SILICON ) 
    return(ChannelDoping) ;
  else
    return(0.0) ;

  return(0.0) ;
}

int signfunc(double x)
{
  if ( x>0.0 ) return(1) ;
  else return(-1) ;
}
