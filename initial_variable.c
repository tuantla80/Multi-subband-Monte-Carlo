/* *****************************************************************************
Module nay la de khoi tao cac BIEN va HAM duoc DUNG DI DUNG LAI trong vong lap
for cho (Vg va Vd).
	       
Starting date: Feb 21, 2010
Update:        March 11, 2010
Latest update: May 03, 2010 (Them tham so de mapping voi 3D Poisson)
**************************************************************************** */
#include <stdio.h>
#include "nrutil.h"
#include "constants.h"
#include <time.h>
#include <unistd.h>

static int n_used; // So hat electron ta khoi tao va sau do su dung de chua so hat
                   // electron dang hoat dong trong device
static double **p, *energy;
static int *valley,*subband;// Hat co valley[]=9 (tuong duong voi ip[]=9 o Semi Monte Carlo) thi se bi inactive va sau do la deleted
static double kx, dtau, x_position, electron_energy; 
static int iv, subb;// iv: de chi valley pair index va subb de chi subband index
       // ca bien kx,dtau,x_position,e,iv, subb la de dung o ham emcd() va cac ham duoc goi tu ham do
static long idum; // Cho ham random2()  
static double *****wave, ***eig,***eig_jthSchro; // For 2D Schrodinger Eq. 
static double ***doping,***fai_jthSchro;
static double *pot; // For de chuyen 2D fai(i-hang so,j,k) thanh 1D potential dung cho ham makeh.c trong 2D Schrodinger 
static double **trap_weights; // weighted matrix cho 2D Trapsoidal integration. Day la gia tri co dinh
                               // nen chi can khoi tao 1 lan 
static double ****form_factor; // de tinh cho phonon scattering
static double *****scat_table; // 5D array de luu scattering table. (March 12, 2010 CHUA xet region)
static double *max_gm;// la max cua gamma cho tung section. 
static int iss_out,iss_eli,iss_cre,idd_out,idd_eli,idd_cre;//Tinh so hat VAO/RA de tu do tinh DONG DIEN
static int particle_i_th; // Chi hat THU may 26/03/10 15:13
static double ***electron_density;// la electron density tinh o void electron_density_caculation()
static double *velocity_x_sum, *energy_sum, *current_sum;// Tinh tong cac gia tri cho all time steps
static double *velocity_x_sum_after_transient, *energy_sum_after_transient, *current_sum_after_transient;
static double *pot_x_avg, **pot_yz_avg;// potential average theo phuong x khi y va z CO DINH, va trong mat cross-section yz khi x CO DINH
                                         
void initial_variables(){ // Khoi tao gia tri cho tat ca cac bien o define_functions.c
     // Goi ham
     int Get_nx0(),Get_nx1(),Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
     int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int Get_NSELECT();
     // Cac bien local
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int ny0 = Get_ny0(); int ny1 = Get_ny1(); int nya = Get_nya(); int nyb = Get_nyb();
     int nz0 = Get_nz0(); int nz1 = Get_nz1(); int nza = Get_nza(); int nzb = Get_nzb();
     int NSELECT = Get_NSELECT();
     // Thuc hien
     int i,j,k,s,v,n,m,e_step;
     n_used  = 0;
     
     p       = dmatrix (1,max_electron_number,1,4);
     energy  = dvector (1,max_electron_number);
     valley  = ivector (1,max_electron_number);
     subband = ivector (1,max_electron_number);
     for(i=1; i<=max_electron_number; i++){
            energy[i]  = 0.0;
            valley[i]  = 0;
            subband[i] = 0;
         for(j=1; j<=4; j++){ 
	   p[i][j] = 0.0;
         }
     }
                  
     kx = 0.0; dtau =0.0; x_position = 0.0; electron_energy = 0.0;
     iv =0; subb = 0;
     //idum    = ((long) time(NULL));// * getpid()); //printf("\n idum =%ld", idum); THU KHONG khoi tao xem sao
     
     eig = d3matrix(nx0,nx1,1,3,1,NSELECT);//eig = d3matrix(0,nx_max,1,3,1,NSELECT); ok vi nx0=0 va nx1=nx_max
     eig_jthSchro = d3matrix(nx0,nx1,1,3,1,NSELECT);//d3matrix(0,nx_max,1,3,1,NSELECT);
     int ny = nyb - nya;
     int nz = nzb - nza;
     wave = d5matrix(nx0,nx1,1,3,1,NSELECT,0,ny,0,nz);//cua Loi Silicon. (May 30, 2010) ro rang so diem chia = so khoang +1 
     
     for(s=nx0; s<=nx1; s++){
        for(v=1; v<=3; v++){
            for(i=1; i<=NSELECT; i++){ 
                     eig[s][v][i] = 0.0; 
                     eig_jthSchro[s][v][i]=0.0;
                for(j=0; j<=ny; j++){// Chay tu 0 da the hien so diem chia la + 1 roi con gi 
                   for(k=0; k<=nz; k++){       
                       wave[s][v][i][j][k] = 0.0;
                   }
               }
            }
        }
     }     
    
     doping = d3matrix(nx0,nx1,ny0,ny1,nz0,nz1); 
     for(i=nx0; i<=nx1; i++){
        for(j=ny0; j<=ny1; j++){
            for(k=nz0; k<=nz1; k++){       
               doping[i][j][k] = 0.0;
	    }
        }
     }

     fai_jthSchro = d3matrix(nx0-1, nx1+1, ny0-1, ny1+1, nz0-1, nz1+1);  //NOTE cua GS thi Phi chi so lai them dau la -1, duoi la +1
     for(i=nx0-1; i<=nx1+1; i++){
        for(j=ny0-1; j<=ny1+1; j++){
            for(k=nz0-1; k<=nz1+1; k++){       
               fai_jthSchro[i][j][k] = 0.0;
           }
        }
     }
                
     pot     = dvector(1,(ny+1)*(nz+1));//Do so diem chia o loi silicon se la nyb-nya+1 va nzb-nza+1
     for(i=1; i<=(ny+1)*(nz+1); i++){ pot[i] = 0.0; }
     
     trap_weights = dmatrix (ny0,ny1,nz0,nz1); // Khoi tao bien trap_weights lan dau tien
     for(i=ny0; i<=ny1; i++){ 
         for(j=nz0; j<=nz1; j++){ trap_weights[i][j] = 0.0;
         }
     }
     
     form_factor = d4matrix(nx0,nx1,1,NSELECT,1,NSELECT,1,9);
     for(s=nx0; s<=nx1; s++){
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(v=1; v<=9; v++){      form_factor[s][n][m][v] = 0.0;
                }
             }
         }
     }
          
     scat_table = d5matrix(nx0,nx1,1,NSELECT,1,NSELECT,1,n_lev,1,21);
     for(s=nx0; s<=nx1; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){ 
                   for(v=1; v<=21; v++){         
                        scat_table[s][n][m][e_step][v]=0.0;
                   }
                }
             }
         }
     }
     
     max_gm = dvector(nx0,nx1);

     for(s=nx0; s<=nx1; s++){ 
       max_gm[s] = 0.0;
     }
     
     iss_out = 0;
     iss_eli = 0;
     iss_cre = 0;
     idd_out = 0;
     idd_eli = 0;
     idd_cre = 0;
     
     particle_i_th = 0;
     
     electron_density = d3matrix(nx0,nx1,ny0,ny1,nz0,nz1); //  ny = nyb - nya; nz = nzb - nza;SO KHOANG CHIA trong loi Silicon, nen so diem chia phai tang them 1
     for(i=nx0; i<=nx1; i++){
       for(j=ny0; j<=ny1; j++){// CHay tu da the hien so diem chia la +1 roi con gi 
            for(k=nz0; k<=nz1; k++){       
               electron_density[i][j][k] = 0.0;
             }
         }
     }
  
     velocity_x_sum = dvector(nx0,nx1);
     energy_sum     = dvector(nx0,nx1);
     current_sum    = dvector(nx0,nx1);
     velocity_x_sum_after_transient = dvector(nx0,nx1);
     energy_sum_after_transient     = dvector(nx0,nx1);
     current_sum_after_transient    = dvector(nx0,nx1);
     for(i=nx0; i<=nx1; i++){
         velocity_x_sum[i]=0.0; energy_sum[i]=0.0;
         current_sum[i]=0.0; velocity_x_sum_after_transient[i]=0.0;
         energy_sum_after_transient[i]=0.0; current_sum_after_transient[i]=0.0;
     } 

     pot_x_avg =  dvector(nx0,nx1); 
     for(i=nx0; i<=nx1; i++){ 
       pot_x_avg[i] = 0.0; 
     }

     pot_yz_avg =  dmatrix (ny0,ny1,nz0,nz1);
     for(i=ny0; i<=ny1; i++){ 
         for(j=nz0; j<=nz1; j++){ pot_yz_avg[i][j] = 0.0;
         }
     }
     

  return;
}
// End of void initial_variables(){ // Khoi tao gia tri cho tat ca cac bien o define_functions.c

// For n_used, khoi tao LAN DAU tai void electrons_initialization()
void Set_n_used( int n){
        n_used = n;
 }

int Get_n_used(){
     return n_used;
}// End of for n_used

// The parameters for paticles
double **Get_p(){
       return (p);
       }
       
double *Get_energy(){
       return (energy);
       }
       
int *Get_valley(){
       return (valley);
       }
       
int *Get_subband(){
    return (subband);
}


void Set_kx(double k_momen_x){
     kx = k_momen_x;
}
double Get_kx(){
       return kx;
}

void Set_dtau(double flight_time){
     dtau = flight_time;
}

double Get_dtau(){
       return dtau;
}

void Set_x_position(double position_in_x){
     x_position = position_in_x;
}

double Get_x_position(){
       return x_position;
}

void Set_electron_energy(double ener){
     electron_energy = ener;
}

double Get_electron_energy(){
       return electron_energy;
}

void Set_iv(int valley_pair_index){
     iv = valley_pair_index;
     }
int Get_iv(){
    return iv;
}

void Set_subb(int subband_index){
     subb = subband_index;
}

int Get_subb(){
    return subb;
}


// End of The parameters for paticles

// For random2
long Get_idum(){
     return idum;
     }
     
double *****Get_wave(){
       return (wave);
}     

double ***Get_eig(){
       return (eig);
}

double ***Get_eig_jthSchro(){
       return (eig_jthSchro);
}


double ***Get_doping(){
       return (doping);
}

double ***Get_fai_jthSchro(){
       return (fai_jthSchro);
}

double *Get_pot(){  
       return (pot);
       }
double **Get_trap_weights(){
       return (trap_weights);
}

double ****Get_form_factor(){
       return (form_factor);
}
       
double *****Get_scat_table(){
       return (scat_table);
}

double *Get_max_gm(){
       return (max_gm);
}

void Set_iss_out(int out_p){ // out_p nghia la HAT RA trong tieng Anh
     iss_out = out_p;
}

int Get_iss_out(){
    return iss_out;
}

void Set_iss_eli(int eli_p){ // eli_p nghia la LOAI BO HAT trong tieng Anh
     iss_eli = eli_p;
}

int Get_iss_eli(){
    return iss_eli;
}

void Set_iss_cre(int cre_p){ // cre_p nghia la TAO HAT MOI trong tieng Anh
     iss_cre = cre_p;
}

int Get_iss_cre(){
    return iss_cre;
}

void Set_idd_out(int out_p){ 
     idd_out = out_p;
}

int Get_idd_out(){
    return idd_out;
}

void Set_idd_eli(int eli_p){ 
     idd_eli = eli_p;
}

int Get_idd_eli(){
    return idd_eli;
}

void Set_idd_cre(int cre_p){ // cre_p nghia la TAO HAT MOI trong tieng Anh
     idd_cre = cre_p;
}

int Get_idd_cre(){
    return idd_cre;
}

void Set_particle_i_th( int i){ // i la thu tu hat i_th o emcd() For checking
     particle_i_th = i;
}

int Get_particle_i_th(){
    return particle_i_th;
}

double ***Get_electron_density(){
       return (electron_density);
}     

double *Get_velocity_x_sum(){
       return (velocity_x_sum);
}

double *Get_energy_sum(){
       return (energy_sum);
}

double *Get_current_sum(){
       return (current_sum);
}

double *Get_velocity_x_sum_after_transient(){
       return (velocity_x_sum_after_transient);
}

double *Get_energy_sum_after_transient(){
       return (energy_sum_after_transient);
}
       
double *Get_current_sum_after_transient(){
       return (current_sum_after_transient);
}
       
double *Get_pot_x_avg(){
  return( pot_x_avg);
}

double **Get_pot_yz_avg(){
  return ( pot_yz_avg);
}























            
