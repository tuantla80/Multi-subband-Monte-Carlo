/* **********************************************************************************
This function is to read materials parameter for Silicon from Input file

Starting date: Feb 08, 2010
Update:        Feb 08, 2010
Update:        April 28, 2010 (Them tham so de mapping voi 3D Poisson
Latest update: May 11, 2010: MPI
***************************************************************************** */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "petscksp.h" // for Poisson 3D
#include "nanowire.h" // for Poisson 3D

static double ml, mt; // Effective mass for ellipsoidal
static double constant_acoustic, constant_optical_g_e, constant_optical_g_a;
static double constant_optical_f_e, constant_optical_f_a;
static double Vt, Eg, BandgapOxide, nonparabolicity_factor; //thermal voltage, band gap, 
static double intrinsic_carrier_density,affinity_silicon,gateWF; 
static double emax; // dai nang luong ma ta quan tam doc tu Input file
static double hw0f_phonon,hw0g_phonon;// f-and g-phonon energy can khi tinh Ef o 1D Density of State

void material_param(){
    double T;// temperature
    double eps_Si_constant, eps_Oxide_constant; // La luc Chua nhan voi eps_0
    double eps_sc,eps_oxide; // Da nhan voi eps_0
    double crystal_density,sound_velocity;
    double DefPot_acoustic, DefPot_f,DefPot_g;
    char gateTYPE[30];
    double GateWorkfunctionOffset=0.0;//BandgapSi=Eg luon. mapping at the end of this function
	/* Reading #MATERIAL_PARAMETER */
    FILE *f;
    f=fopen("Input.d","r");
    if(f==NULL) { printf("Can't open file Input.d.\n"); return;}
    //printf("\n Reading MATERIAL_PARAMETER");
    while (!feof(f)) {
        char str[180];
        fscanf(f,"%s",str);
        if(strcmp(str,"#MATERIAL_PARAMETER")==0){
	  fscanf(f,"%*s %*s %le",&T);
          //printf("\n Temperature[K]                 = %f",T);
		  
	  fscanf(f,"%*s %*s %le",&emax);
          //printf("\n Electron energy maximum[eV]    = %f",emax);
          
          fscanf(f,"%*s %*s %le",&Eg);
          //printf("\n Band gap of Silicon[eV]        = %f",Eg);
          
          fscanf(f,"%*s %*s %le",&BandgapOxide);
          //printf("\n Band gap of Oxide[eV]          = %f",BandgapOxide);
         		  
 	  fscanf(f,"%*s %*s %le",&ml);
          //printf("\n Longitudinal effective mass    = %f",ml);

          fscanf(f,"%*s %*s %le",&mt);
          //printf("\n Tranverse effective mass       = %f",mt);

          fscanf(f,"%*s %*s %le ",&nonparabolicity_factor);
          //printf("\n Nonparabolicity factor[1/eV]   = %f",nonparabolicity_factor);
          //nonparabolicity_factor = nonparabolicity_factor/q; // Do cai dau vao la gia tri vi du 0.5[1/eV] doi sang don vi chuan la 0.5/q [1/J]
          // SAI: trong cong thuc cho Nonparablic thi e(1+af*e) thi e [eV] va af[1/eV]          
          fscanf(f,"%*s %*s %le ",&eps_Si_constant);
          //printf("\n Constant dielectric of Silicon = %f",eps_Si_constant);

          fscanf(f,"%*s %*s %le ",&eps_Oxide_constant);
          //printf("\n Constant dielectric of Oxide   = %f",eps_Oxide_constant);

          fscanf(f,"%*s %*s %le ",&intrinsic_carrier_density);
          //printf("\n Intrinsic density[/cm3]        = %le",intrinsic_carrier_density);
          intrinsic_carrier_density=intrinsic_carrier_density*1.0e+6; // cm-3 -> m-3
          
          fscanf(f,"%*s %*s %s",gateTYPE); //(metal_or_nPolysilicon)
          //printf("\n Gate type is                   = %s",gateTYPE);
          
          fscanf(f,"%*s %*s %le ",&affinity_silicon);
          //printf("\n Affinity of Silicon[eV]        = %f",affinity_silicon);

          fscanf(f,"%*s %*s %le ",&gateWF);
          //printf("\n Metal gate workfunction[eV]    = %f",gateWF);
          
          fscanf(f,"%*s %*s %le ",&GateWorkfunctionOffset);
          //printf("\n Gate Workfunction Offset[eV]   = %f",GateWorkfunctionOffset);
          
          fscanf(f,"%*s %*s %le ",&crystal_density);
          //printf("\n Crystal density [kg/m3]      	= %f",crystal_density);
          
          fscanf(f,"%*s %*s %le ",&sound_velocity);
          //printf("\n Sound velocity [m/s]         	= %f",sound_velocity);
          
          fscanf(f,"%*s %*s %le ",&DefPot_acoustic);
          //printf("\n Acoustic_Deformation_Pot[eV] 	= %f",DefPot_acoustic);
          
          fscanf(f,"%*s %*s %le ",&DefPot_f);
          //printf("\n f_phonon_deformation_pot[ev/m] = %le",DefPot_f);
          
          fscanf(f,"%*s %*s %le ",&hw0f_phonon);
          //printf("\n f_phonon_energy[eV]        	= %f",hw0f_phonon);
          
          fscanf(f,"%*s %*s %le ",&DefPot_g);
          //printf("\n g_phonon_deformation_pot[ev/m] = %le",DefPot_g);
          
          fscanf(f,"%*s %*s %le ",&hw0g_phonon);
          //printf("\n g_phonon_energy[eV]         	= %f",hw0g_phonon);

          } // End of if
    }// End of while
     fclose(f);
     /* End of Reading #MATERIAL_PARAMETER */
     
	//Eg = 1.17-4.73*(1.0e-4)*T*T/(T+636);  // [eV] Nhung o day ban dau da lam gi la eV no chi
    // Do mapping voi Poisson nen ta dua vao tu Input   // la 1 con so thoi. vi du: 1.12 Nhung don vi xac dinh la eV
    /* It will be better if we use the fomular for "Temperature dependence of the energy gap"
    Eg= 1.17-4.73*(10e-4)*T*T/(T+636) [eV] where T is temperature in degrees K
    e.g. at T=300K, Eg=1.124519231 [eV]: the same as the value above
    */
        
    Vt=kb*T/q; //thermal voltage 
	//printf("\n Thermal voltage] = %f [eV] ",Vt);
        
    eps_sc = eps_Si_constant*eps_0; // dielectric of Silicon
    eps_oxide = eps_Oxide_constant*eps_0; // dielectric of Oxide
    
    /* Can chu y ve y nghia vat ly. Neu chon muc Fermi reference la 0 thi neu Vs=0 va Vd=1V co nghia la vung window
	   se la mu2-mu1 = Vd-Vs=1eV. Vay thi electron chay tu Source den Drain chi co the nam o cac muc nang luong trong
	   vung window 0 den Vd-Vs ma thoi. Do vay ta chon emax=4eV nhung thuc te electron energy cung chi co 0-1eV ma thoi
	   -> impactionization co the khong ton tai doi voi truong hop nay. Nen co the chon emax=1.1eV la ok neu Vd khong
	   vuot qua 1.1V
	   Neu chon Vs=0 va cho emax = 5 -> Vd co the chon max la 5V
	   */
    
    /* Tham so o day lay tu //Mobility Enhancement of SOI MOSFETs due to Subband
    Modulation in Ultrathin SOI Films - Takagi al., J. Journal of App. Phys.,
    Vol.37, Part1, No.3B, 1998 ****** */
    DefPot_acoustic = DefPot_acoustic*q;    //12 [eV] -> nhan q thanh [J] Deformation potential for acoustic phonon of Si

    DefPot_f=DefPot_f*q;            // [eV/m] tuc la 1.1e9 eV/cm -> nhan q thanh [J/m] Zero-order deformation potential for f process
   	// Doc tu Input file hw0f_phonon = 0.059 ;// [eV] tuc la 59.0meV f-phonon energy
    double w0f_phonon=hw0f_phonon*q/hbar; // wof = hw0f_phonon*q/hbar (nhan cho q de doi ra J, chia cho hbar la theo cong thuc roi
    double nof = 1.0/(exp(hw0f_phonon/Vt)-1.0);  // boi vi Vt=kbtq=kb*tem/q tuc la ta doi ra thanh don vi eV
            // nen hwof_phonon o day cung tinh theo [eV]

    DefPot_g=DefPot_g*q;  // [eV/m] tuc la 8.0e8 eV/cm -> [J/m] Zero-order deformation potential for g process
    // Doc tu input file hw0g_phonon = 0.063; // [eV] g-phonon energy tuc la 63.0 meV
    double w0g_phonon=hw0g_phonon*q/hbar;      //
    double nog = 1.0/(exp(hw0g_phonon/Vt)-1.0);  //For phonon (Bose-Einstein distribution fucntion)

    constant_acoustic=(2.0*pi*DefPot_acoustic*DefPot_acoustic*kb*T)/(hbar*crystal_density*sound_velocity*sound_velocity);
        // La tham so phan dau cua acoustic phonon scattering
        // Duoi day la tham so phan dau cua non-polar optical phonon scattering cho g hoac f process
        // _e va _a la cho emission va absorption
    constant_optical_g_e = (pi*DefPot_g*DefPot_g)*(nog+1.0)/(crystal_density*w0g_phonon); // do Ziv cho g-process la 1
                         // do Ziv cho g-process la 1
    constant_optical_g_a = constant_optical_g_e*nog/(nog+1.0);
    constant_optical_f_e = 2.0*(pi*DefPot_f*DefPot_f)*(nof+1.0)/(crystal_density*w0f_phonon);
                        // do Ziv cho f-process la 2. Xem giai thich trang 103
    constant_optical_f_a = constant_optical_f_e*nof/(nof+1.0);
    //printf("\n constant_acoustic = %le,constant_optical_g_e = %le  ",constant_acoustic,constant_optical_g_e); getchar();

    // Mapping for Poisson 3D 
    Param *P;
    Energy *E;
    
    P = PGetParam();
    E = PGetEnergy();
  
    //printf("\n At material_param(): Mapping for Poisson 3D");
    P->Temp = T; // Temperature
    //printf("\n Temperature                     P->Temp = %f",P->Temp);
    
    E->g_si = Eg;
    //printf("\n Band gap of Si                  E->g_si = %f",E->g_si);
    
    E->g_ox = BandgapOxide;
    //printf("\n Band gap of Oxide               E->g_ox = %f",E->g_ox);
    
    P->e_si = eps_Si_constant;
    //printf("\n Dielectric constant of Si       P->e_si = %f",P->e_si);
    
    P->e_ox = eps_Oxide_constant;
    //printf("\n Dielectric constant of Oxide    P->e_ox = %f",P->e_ox);
    
    P->offset_gate = GateWorkfunctionOffset;
    //printf("\n Gate workfunction offset P->offset_gate = %f",P->offset_gate);

    //Phan printf ra Monitor o node co rank=0
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if(myrank ==0){ 
        printf("\n The MATERIAL PARAMETERs");
	printf("\n Temperature[K]                 = %f",T);
        printf("\n Electron energy maximum[eV]    = %f",emax);
        printf("\n Band gap of Silicon[eV]        = %f",Eg);
        printf("\n Band gap of Oxide[eV]          = %f",BandgapOxide);
        printf("\n Longitudinal effective mass    = %f",ml);
        printf("\n Tranverse effective mass       = %f",mt);
        printf("\n Nonparabolicity factor[1/eV]   = %f",nonparabolicity_factor);
        printf("\n Constant dielectric of Silicon = %f",eps_Si_constant);
        printf("\n Constant dielectric of Oxide   = %f",eps_Oxide_constant);
        printf("\n Intrinsic density[/m3]        = %le",intrinsic_carrier_density);
        printf("\n Gate type is                   = %s",gateTYPE);
        printf("\n Affinity of Silicon[eV]        = %f",affinity_silicon);
        printf("\n Metal gate workfunction[eV]    = %f",gateWF);
        printf("\n Gate Workfunction Offset[eV]   = %f",GateWorkfunctionOffset);
        printf("\n Crystal density [kg/m3]        = %f",crystal_density);
        printf("\n Sound velocity [m/s]           = %f",sound_velocity);
        printf("\n Acoustic_Deformation_Pot[eV]   = %f",DefPot_acoustic/q);// Do sau khi doc vao thi lai nhan cho q nen o day lai chia cho q
        printf("\n f_phonon_deformation_pot[ev/m] = %le",DefPot_f/q);
        printf("\n f_phonon_energy[eV]        	  = %f",hw0f_phonon);
        printf("\n g_phonon_deformation_pot[ev/m] = %le",DefPot_g/q);
        printf("\n g_phonon_energy[eV]         	  = %f",hw0g_phonon);
    }// End of if(myrank ==0)

return;
} // End of void mat_par_initialization()

/* ************************************************************************* 
	Cac ham minh tao ra de TRANH dung bien Global
    Starting date: Feb 09, 2010
	Latest update: Feb 09, 2010
************************************************************************* */
double Get_ml(){
	return ml;
}

double Get_mt(){
	return mt;
}

double Get_constant_acoustic(){
	return constant_acoustic;
}

double Get_constant_optical_g_e(){
	return constant_optical_g_e;
}

double Get_constant_optical_g_a(){
	return constant_optical_g_a;
}

double Get_constant_optical_f_e(){
	return constant_optical_f_e;
}

double Get_constant_optical_f_a(){
	return constant_optical_f_a;
}

double Get_Vt(){
       return Vt;
}

double Get_intrinsic_carrier_density(){
       return intrinsic_carrier_density;
}
 
double Get_affinity_silicon(){
       return affinity_silicon;
}
double Get_gateWF(){
       return gateWF;
}
double Get_Eg(){
       return Eg;
}
double Get_BandgapOxide(){
  return BandgapOxide;
}

double Get_nonparabolicity_factor(){
       return nonparabolicity_factor;
}

double Get_emax(){
       return emax;
}       

double Get_hw0f_phonon(){
       return hw0f_phonon;
       }

double Get_hw0g_phonon(){
       return hw0g_phonon;
       }
       
	
	
	
	
	
	
	
	
	

