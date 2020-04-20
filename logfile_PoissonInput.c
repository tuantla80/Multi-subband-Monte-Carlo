/* *****************************************************************************
Kiem tra cac dau vao Poisson co dung nhu GS khong

Starting date: June 16, 2010
Latest update: June 16, 2010
****************************************************************************** */
#include <stdio.h>
#include <string.h>
#include "petscksp.h"
#include "nanowire.h"

void logfile_PoissonInput(){
  void logging(FILE *log,const char *fmt,...);
  FILE *logfile = fopen("logfile_PoissonInput.txt","w");

  Param *P ;
  Dimen *D ;
  PoiNum *N ;
  Energy *E ;
  Density *n ;
  PoiParam *CP ;
  P = PGetParam() ;
  D = PGetDimen() ;
  N = PGetPoiNum() ;
  n = PGetDensity() ;
  E = PGetEnergy() ;
  CP = PGetPoiParam() ;
  
  logging(logfile,"\n ************* Physical_Parameters*********************");
  logging(logfile,"\n Temperature                     P->Temp = %f",P->Temp);
  logging(logfile,"\n Band gap of Si                  E->g_si = %f",E->g_si);
  logging(logfile,"\n Band gap of Oxide               E->g_ox = %f",E->g_ox);
  logging(logfile,"\n Dielectric constant of Si       P->e_si = %le",P->e_si);
  logging(logfile,"\n Dielectric constant of Oxide    P->e_ox = %le",P->e_ox);
  logging(logfile,"\n Gate workfunction offset P->offset_gate = %f",P->offset_gate);

  logging(logfile,"\n\n ************* Doping *********************");
  double GetDopingDensityFor(char region);
  logging(logfile,"\n Source Doping GetDopingDensityFor('s')  = %le", GetDopingDensityFor('s')); 
  logging(logfile,"\n Channel Doping GetDopingDensityFor('c') = %le", GetDopingDensityFor('c')); 
  logging(logfile,"\n Drain Doping GetDopingDensityFor('d')   = %le", GetDopingDensityFor('d')); 

  logging(logfile,"\n\n ************* Dimensions(nm)  *********************");
  logging(logfile,"\n Lsource[NoUnit] D->Lsrc     =%f,  Number of points N->Mx_src    =%d",D->Lsrc,N->Mx_src);
  logging(logfile,"\n Lchannel[NoUnit]D->Lchannel =%f, Number of points N->Mx_channel =%d",D->Lchannel,N->Mx_channel);
  logging(logfile,"\n Tox[NoUnit] D->Tox          =%f,   Number of points N->Mz_ox    =%d",D->Tox,N->Mz_ox);
  logging(logfile,"\n Tsi[NoUnit] D->Tsi          =%f,   Number of points N->Mz_si    =%d",D->Tsi,N->Mz_si);
  logging(logfile,"\n Tbox[NoUnit] D->Tbox        =%f,   Number of points N->Mz_box   =%d",D->Tbox,N->Mz_box);
  logging(logfile,"\n Wox[NoUnit] D->Wox          =%f,   Number of points N->My_ox    =%d",D->Wox,N->My_ox);
  logging(logfile,"\n Wsi[NoUnit] D->Wsi          =%f,   Number of points N->My_si    =%d",D->Wsi,N->My_si);
  logging(logfile,"\n Lgate[NoUnit] D->Lgate      =%f",D->Lgate);
  logging(logfile,"\n Tgate[NoUnit] D->Tgate      =%f",D->Tgate);
  logging(logfile,"\n Wgate[NoUnit] D->Wgate      =%f",D->Wgate);

  logging(logfile,"\n\n ***************** Voltages * *********************");
  double GetDrainVoltage();
  double GetGateVoltage();
  logging(logfile,"\n GetDrainVoltage() =%f, GetGateVoltage()=%f",GetDrainVoltage(),GetGateVoltage());

  logging(logfile,"\n\n ************* POISSON::Control_Parameters **** **********");
  logging(logfile,"\n Max Interation CP->MaxIter   =%d",CP->MaxIter);
  logging(logfile,"\n Convergence_Eps CP->ConvEps  =%le",CP->ConvEps);
  logging(logfile,"\n KSPType CP->ksptype          =%s",CP->ksptype);
  logging(logfile,"\n PCType CP->pctype            =%s",CP->pctype);
  logging(logfile,"\n KSP_Rtol CP->ksprtol         =%le",CP->ksprtol);
  logging(logfile,"\n GMRES_Restart CP->gmres_restart =%d",CP->gmres_restart);
  
  fclose(logfile);
  return;

}

#include <stdarg.h>
#include "constants.h"

void logging(FILE *logfile,const char *fmt,...) // fmt la format kieu nhu frintf
{
  
  #ifdef ENABLE_LOG
    va_list ap;
    va_start(ap,fmt);
    vfprintf(logfile, fmt, ap);
    //fprintf(logfile, "########################################\n");
    va_end(ap);
  #endif 
}
