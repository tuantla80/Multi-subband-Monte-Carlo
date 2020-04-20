typedef struct {
  double Temp ;

  double e_si ;
  double e_ox ;
  double offset_gate ;

  int poi_lin ;
  int trans_size ;
} Param ;

typedef struct {
  double Lsrc ;
  double Lchannel ;
  double Lgate ;
  double Tox ;
  double Tsi ;
  double Tbox ;
  double Tgate ;
  double Wox ;
  double Wsi ;
  double Wgate ;
} Dimen ;

typedef struct {
  double kT ;
  double g_si ;
  double g_ox ;
  double c_si ;
  double c_ox ;
  double F ;
  double FnS ;
  double FnD ;
  double biS ;
  double biD ;
} Energy ;

typedef struct {
  double A_zero ;
  double D_zero ;
} Density ;

Param *PGetParam() ;
Param GetParam() ;
Dimen *PGetDimen() ;
Dimen GetDimen() ;
Density *PGetDensity() ;
Density GetDensity() ;
Energy *PGetEnergy() ;
Energy GetEnergy() ;

FILE *GetFout() ;
void SetFout() ;
void   SetDrainVoltage() ;
double GetDrainVoltage() ;
double GetGateVoltage() ;

void SetPhi_Contact() ;
double GetPhi_Source() ;
double GetPhi_Drain() ;
double GetPhi_Gate() ;

void SetEF() ;
void SetEFn() ;
void SetEFnS() ;
void SetEFnD() ;

double GetEF() ;
double GetEFnS() ;
double GetEFnD() ;

double GetESi() ;
double GetEOx() ;
double GetEc_si() ;
double GetEc_ox() ;
double GetEc() ;
double GetKT() ;
double GetN0() ;
double GetP0() ;

int  GetTransSize() ;
void GetLocalN() ;

void SetDrainVoltage() ;
void SetGateVoltage() ;
double GetDrainVoltage() ;
double GetGateVoltage() ;

// Poisson Linearized

void SetPoissonLinearized() ;
int  GetPoissonLinearized() ;

// Doping Density

void SetDopingDensityFor() ;
double GetDopingDensityFor() ;
