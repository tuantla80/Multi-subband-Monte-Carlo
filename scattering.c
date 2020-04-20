/* ******************************************************************************     
   The scattering() function is to
    - select a scattering mechanism by which the particles are scattered
    - calculate the states of particles after scattering 
            (kx, valley(iv), subband, energy, i_region)
   
   Isotropic: is assumed for Acoustic and Non-polar Optical phonon
   
Starting date: March 23, 2010
Latest update: March 25, 2010
********************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "constants.h"

void scattering() {
     // Goi ham 
     void Acoustic_mechanism(int s_th,int valley,int e_step,double ei,double eff_mass,double rr);
     void ZeroOrder_Optical_f_Absorption(int s_th,int val_befo,int val_aft,int e_step,double ei,double eff_mass,double rr,int index);
     void ZeroOrder_Optical_f_Emission(int s_th,int val_befo,int val_aft,int e_step,double ei,double eff_mass,double rr,int index);
     void ZeroOrder_Optical_g_Absorption(int s_th,int val,int e_step,double ei,double eff_mass,double rr,int index); 
     void ZeroOrder_Optical_g_Emission(int s_th,int val,int e_step,double ei,double eff_mass,double rr,int index);
     double *****Get_scat_table();
     double Get_electron_energy(),Get_emax(),Get_x_position(),Get_mesh_size_x();
     double Get_ml(), Get_mt();
     int Get_NSELECT();
     long Get_idum();
     float random2(long *idum); 
     
     // int Get_particle_i_th(); // dua thu tu hat dang simulate va di vao ham scattering
     
     // Cac bien local
     double *****scat_table = Get_scat_table();
     double electron_energy = Get_electron_energy();
     double e_max = Get_emax(); // lay maximum electron energy ta dinh nghia
     double x_position = Get_x_position();
     double mesh_size_x = Get_mesh_size_x();
     double ml = Get_ml();
     double mt = Get_mt();
     int NSELECT = Get_NSELECT();

     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    
     long idum  = Get_idum();
     
// Thuc hien
   int valley_before,valley_after; // before and after scattering
   int index; // chi so cuoi cho mang scat_table[s_th][n][m][ener_step][index]
              // No cung the hien loai scattering nao. Chu y: index=0:Self-scattering 
   double rr; // gia tri random number da duoc chon
   double effective_mass; // phu thuoc vao valley_after la valley nao
   
   // Buoc 1: tu ham dritf() se biet duoc hat o section nao va co energy step la bao nhieu
   int s_th = (int)(round(x_position/mesh_size_x)); // s-th section
       if(s_th < 0)      { s_th = 0;}
       if(s_th > nx_max) {s_th = nx_max;}
   double de = e_max/(double)(n_lev); //energy_interval
   if(electron_energy > e_max)// (May 28, 2010) vi co the electron energy no vuot qua cai nguong ma ta dinh nghia e_max lam thi sao ?
     {
       electron_energy = e_max;
     }
   int ener_step = (int)(electron_energy/de);// Gia tri DONG NANG dau vao thi THUOC energy step THU MAY
       if(ener_step < 1)      { ener_step = 1;}
       if(ener_step > n_lev)  {ener_step = n_lev;}
   // Buoc 2. Chon loai scattering theo index vv=1 den 21. Dua vao bang trang 107 va 122
   double random_number = random2(&idum);
   
   if(random_number > scat_table[s_th][NSELECT][NSELECT][ener_step][21])
        { // self scattering. KHONG thay doi bat cu Parameters nao
          index = 0; // Self-scattering 
          //printf("\n particle %d-th: Self-scattering",Get_particle_i_th());
          return; // ket thuc ham void scattering();. Kiem tra self scattering truoc de tang toc do
        }
        
   else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][1])
        { // Chon Acoustic : Valley 1 den Valley 1
           valley_after = 1; //valley_before = 1; //index=1 TRUNG voi CHI SO valley
           index =1;
           effective_mass = ml*m0;
           rr = random_number; // random_number nay se duoc su dung de chon chinh xac la no o vi tri nao
           Acoustic_mechanism(s_th,valley_after,ener_step,electron_energy,effective_mass,rr);
        }
                      
   else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][2])
        { // Chon Acoustic : Valley 2 den Valley 2
           valley_after = 2; //valley_before = 2;// index=2 TRUNG voi CHI SO valley
           index = 2;
           effective_mass = mt*m0;
           rr = random_number; 
           Acoustic_mechanism(s_th,valley_after,ener_step,electron_energy,effective_mass,rr);
        } 
        
   else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][3])
        { // Chon Acoustic : Valley 3 den Valley 3
           valley_after = 3; //valley_before = 3;//index=3 TRUNG voi CHI SO valley
           index = 3;
           effective_mass = mt*m0;
           rr = random_number; 
           Acoustic_mechanism(s_th,valley_after,ener_step,electron_energy,effective_mass,rr);
        }     
                               
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][4])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 1;
           valley_after  = 2;
           index = 4;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][5])
        { // f-Emission Zero order nonpolar optical
           valley_before = 1;
           valley_after  = 2;
           index = 5;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        } 
    
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][6])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 1;
           valley_after  = 3;
           index = 6;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][7])
        { // f-Emission Zero order nonpolar optical
           valley_before = 1;
           valley_after  = 3;
           index = 7;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        } 
        
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][8])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 2;
           valley_after  = 1;
           index = 8;
           effective_mass = ml*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][9])
        { // f-Emission Zero order nonpolar optical
           valley_before = 2;
           valley_after  = 1;
           index = 9;
           effective_mass = ml*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }  
        
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][10])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 2;
           valley_after  = 3;
           index = 10;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][11])
        { // f-Emission Zero order nonpolar optical
           valley_before = 2;
           valley_after  = 3;
           index = 11;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }  
                                          
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][12])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 3;
           valley_after  = 1;
           index = 12;
           effective_mass = ml*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][13])
        { // f-Emission Zero order nonpolar optical
           valley_before = 3;
           valley_after  = 1;
           index = 13;
           effective_mass = ml*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                           
                         
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][14])
        { // f-Absortion Zero order nonpolar optical
           valley_before = 3;
           valley_after  = 2;
           index = 14;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Absorption(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][15])
        { // f-Emission Zero order nonpolar optical
           valley_before = 3;
           valley_after  = 2;
           index = 15;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_f_Emission(s_th,valley_before,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }  
        
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][16])
        { // g-Absortion Zero order nonpolar optical
           valley_after  = 1; //valley_before = 1; 
           effective_mass = ml*m0;
           index = 16;
           rr = random_number; 
           ZeroOrder_Optical_g_Absorption(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][17])
        { // g-Emission Zero order nonpolar optical
           valley_after  = 1; //valley_before = 1; 
           effective_mass = ml*m0;
           index = 17;
           rr = random_number; 
           ZeroOrder_Optical_g_Emission(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                              
     
     else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][18])
        { // g-Absortion Zero order nonpolar optical
           valley_after  = 2; // valley_before = 2; 
           index = 18;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_g_Absorption(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][19])
        { // g-Emission Zero order nonpolar optical
           valley_after  = 2; //valley_before = 2;
           index = 19;
           effective_mass = mt*m0; 
           rr = random_number; 
           ZeroOrder_Optical_g_Emission(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }
        
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][20])
        { // g-Absortion Zero order nonpolar optical
           valley_after  = 3; // valley_before = 3;
           index = 20;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_g_Absorption(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }                        
                         
    else if(random_number <= scat_table[s_th][NSELECT][NSELECT][ener_step][21])
        { // g-Emission Zero order nonpolar optical
           valley_after  = 3; //valley_before = 3;
           index = 21;
           effective_mass = mt*m0;
           rr = random_number; 
           ZeroOrder_Optical_g_Emission(s_th,valley_after,ener_step,electron_energy,effective_mass,rr,index);
        }  
    
    else{ // Da kiem tra tat ca cac dieu kien roi nen neu o day co nghia la LOI
          printf("\n Some ERRORs in scattering.c");
        }
   //if(index !=0){ // Do self-scattering da xet va printf neu co roi
         //printf("\n particle %d-th with type of scattering mechanism %d",Get_particle_i_th(), index);   
   //  }
 return;
}// End of void scattering() 
//******************************************************************************

/* ******************************************************************************     
   Thuc hien machanism cho Acoustics
   INPUT:   + s-th: section nao
            + valley: o valley nao (before va after la GIONG NHAU cho Acoustic)
            + e_step: o energy step nao?
            + ei: dong nang dau vao la bao nhieu. That ra tu ei cung biet duoc e_step
            + effective mass: Tham so nay duoc DUA vao ham isotropic
            + rr: so random_number nay se chon 1 trong cac kieu scattering cua Acoustic
            + NSELECT: number of subband for each valley
            + eig: eigen value
          
   OUTPUT: Update trang thai sau scattering goi ham isotropic()
  
Starting date: March 23, 2010
Latest update: March 25, 2010
********************************************************************************** */

void Acoustic_mechanism(int s_th,int valley,int e_step,double ei,double effective_mass,double rr)
{
    // Goi cac ham
    void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    double ***Get_eig(),*****Get_scat_table();
    int Get_NSELECT();
    
    // Bien local
    double ***eig = Get_eig();
    double *****scat_table = Get_scat_table();
    int NSELECT = Get_NSELECT();
       
// Thuc hien
    double Ef; // final energy
    int n,m;
    for(n=1; n<=NSELECT; n++){
       for(m=1; m<=NSELECT; m++){
           if(rr <=scat_table[s_th][n][m][e_step][valley]){
              Ef = eig[s_th][valley][n] - eig[s_th][valley][m] + ei;//[eV]
              // valley_final = valley;
              // subband_final = m
              isotropic(valley,m,Ef,effective_mass);
              goto L_end; // Neu chon duoc loai scattering roi thi cho no thoat ra khoi vong lap luon
           }
       }
    }// End of for(n=1; n<=NSELECT; n++)
    
 L_end:
 return;                 
}
// End of void Acoustic_mechanism(
/*********************************************************************************
Thuc hien machanism cho ZeroOrder_Optical_f_Absorption or Emission
   INPUT:   + s-th: section nao
            + valley_before: valley truoc scattering
            + valley_after: valley sau scattering
            + e_step: o energy step nao?
            + ei: dong nang dau vao la bao nhieu. That ra tu ei cung biet duoc e_step
            + effective mass: Tham so nay duoc DUA vao ham isotropic
            + rr: so random_number nay se chon 1 trong cac kieu scattering cua Acoustic
            + NSELECT: number of subband for each valley
            + index: chi so cuoi cung cho mang scat_table
          
   OUTPUT: Update trang thai sau scattering goi ham isotropic()

Starting date: March 23, 2010
Latest update: March 25, 2010
********************************************************************************** */

void ZeroOrder_Optical_f_Absorption(int s_th,int valley_before,int valley_after,int e_step,double ei,double effective_mass,double rr, int index)
{
    // Goi cac ham
    void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    double ***Get_eig(),*****Get_scat_table();
    int Get_NSELECT();
    double Get_hw0f_phonon();
    // Bien local
    double ***eig = Get_eig();
    double *****scat_table = Get_scat_table();
    int NSELECT = Get_NSELECT();
    double hw0f_phonon = Get_hw0f_phonon();
       
// Thuc hien
    double Ef; // final energy
    int n,m;
    for(n=1; n<=NSELECT; n++){
       for(m=1; m<=NSELECT; m++){
           if(rr <=scat_table[s_th][n][m][e_step][index]){
              Ef = eig[s_th][valley_before][n] - eig[s_th][valley_after][m] + ei + hw0f_phonon;//[eV] Absorption
              // valley_final = valley_after;
              // subband_final = m
              isotropic(valley_after,m,Ef,effective_mass);
              goto L_end; // Neu chon duoc loai scattering roi thi cho no thoat ra khoi vong lap luon
           }
       }
    }// End of for(n=1; n<=NSELECT; n++)
    
 L_end:
 return;        
  
} // End of void ZeroOrder_Optical_f_Absorption


void ZeroOrder_Optical_f_Emission(int s_th,int valley_before,int valley_after,int e_step,double ei,double effective_mass,double rr, int index)
{
  // Goi cac ham
    void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    double ***Get_eig(),*****Get_scat_table();
    int Get_NSELECT();
    double Get_hw0f_phonon();
    // Bien local
    double ***eig = Get_eig();
    double *****scat_table = Get_scat_table();
    int NSELECT = Get_NSELECT();
    double hw0f_phonon = Get_hw0f_phonon();
       
// Thuc hien
    double Ef; // final energy
    int n,m;
    for(n=1; n<=NSELECT; n++){
       for(m=1; m<=NSELECT; m++){
           if(rr <=scat_table[s_th][n][m][e_step][index]){
              Ef = eig[s_th][valley_before][n] - eig[s_th][valley_after][m] + ei - hw0f_phonon;//Emission
              // valley_final = valley_after;
              // subband_final = m
              isotropic(valley_after,m,Ef,effective_mass);
              goto L_end; // Neu chon duoc loai scattering roi thi cho no thoat ra khoi vong lap luon
           }
       }
    }// End of for(n=1; n<=NSELECT; n++)
    
 L_end:
 return;        
}
// End of void ZeroOrder_Optical_f_Emission

/*********************************************************************************
Thuc hien machanism cho ZeroOrder_Optical_g_Absorption or Emission
   INPUT:   + s-th: section nao
            + valley: valley truoc va sau scattering GIONG NHAU
            + e_step: o energy step nao?
            + ei: dong nang dau vao la bao nhieu. That ra tu ei cung biet duoc e_step
            + effective mass: Tham so nay duoc DUA vao ham isotropic
            + rr: so random_number nay se chon 1 trong cac kieu scattering cua Acoustic
            + index: chi so cuoi cung cho mang scat_table
            + NSELECT: number of subband for each valley
          
   OUTPUT: Update trang thai sau scattering goi ham isotropic()
       
Starting date: March 23, 2010
Latest update: March 25, 2010
********************************************************************************** */

void ZeroOrder_Optical_g_Absorption(int s_th,int valley,int e_step,double ei,double effective_mass,double rr,int index)
{
   // Goi cac ham
    void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    double ***Get_eig(),*****Get_scat_table();
    int Get_NSELECT();
    double Get_hw0g_phonon();
    // Bien local
    double ***eig = Get_eig();
    double *****scat_table = Get_scat_table();
    int NSELECT = Get_NSELECT();
    double hw0g_phonon = Get_hw0g_phonon();
       
// Thuc hien
    double Ef; // final energy
    int n,m;
    for(n=1; n<=NSELECT; n++){
       for(m=1; m<=NSELECT; m++){
           if(rr <=scat_table[s_th][n][m][e_step][index]){
              Ef = eig[s_th][valley][n] - eig[s_th][valley][m] + ei + hw0g_phonon;//Absorption
              // valley_final = valley;
              // subband_final = m
              isotropic(valley,m,Ef,effective_mass);
              goto L_end; // Neu chon duoc loai scattering roi thi cho no thoat ra khoi vong lap luon
           }
       }
    }// End of for(n=1; n<=NSELECT; n++)
    
 L_end:
 return;        
} // End of void ZeroOrder_Optical_f_Absorption


void ZeroOrder_Optical_g_Emission(int s_th,int valley,int e_step,double ei,double effective_mass,double rr,int index)
{
    // Goi cac ham
    void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    double ***Get_eig(),*****Get_scat_table();
    int Get_NSELECT();
    double Get_hw0g_phonon();
    // Bien local
    double ***eig = Get_eig();
    double *****scat_table = Get_scat_table();
    int NSELECT = Get_NSELECT();
    double hw0g_phonon = Get_hw0g_phonon();
       
// Thuc hien
    double Ef; // final energy
    int n,m;
    for(n=1; n<=NSELECT; n++){
       for(m=1; m<=NSELECT; m++){
           if(rr <=scat_table[s_th][n][m][e_step][index]){
              Ef = eig[s_th][valley][n] - eig[s_th][valley][m] + ei - hw0g_phonon;//Emission
              // valley_final = valley;
              // subband_final = m
              isotropic(valley,m,Ef,effective_mass);
              goto L_end; // Neu chon duoc loai scattering roi thi cho no thoat ra khoi vong lap luon
           }
       }
    }// End of for(n=1; n<=NSELECT; n++)
    
 L_end:
 return;      

}// End of void ZeroOrder_Optical_g_Emission(
//**********************************************************************************
/*********************************************************************************
Thuc hien isotropic
   INPUT:  + valley_final: valley after scattering
           + subband_final: subband after scattering
           + Ef: final energy after scattering
           + effective mass: la ml HAY mt TUY THUOC vao valley_final
           
   OUTPUT: Update trang thai sau scattering if Ef>0
           + o valley nao: iv
           + o subband nao: m
           + gia tri energy la bao nhieu: electron_energy (xet khi Ef>0)
           + gia tri kx la bao nhieu: kx

Starting date: March 25, 2010
Latest update: March 26, 2010
********************************************************************************** */

void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass){
   // Goi ham
   void Set_iv(int valley_pair_index),Set_subb(int subband_index);
   void Set_kx(double k_momen_x),Set_electron_energy(double ener);
   double Get_nonparabolicity_factor();
   long Get_idum();
   float random2(long *idum); 
   
   // Bien local
   double af = Get_nonparabolicity_factor();
   long idum  = Get_idum();
   
// Thuc hien
   double kx,k_update,random_number;
   if(energy_final <=0.0){
          return; // khong thuc hien viec gi ca
       }
   else{ // energy_final >0.0
         // Update carrier wavevector
         k_update=sqrt(2.0*effective_mass*q)*sqrt(energy_final*(1+af*energy_final))/hbar;
         random_number = random2(&idum);
         
         if(random_number <=0.5)
              { 
                kx = k_update; // chon forward process
              }
         else { // random_number > 0.5
                kx = -k_update; // chon backward process
              }
              
         // Update TAT CA ket qua cua HAT SAU SCATTERING
         Set_iv(valley_final);             // valley
         Set_subb(subband_final);          // subband
         Set_kx(kx);                       // momentum
         Set_electron_energy(energy_final);// energy 
         //printf("\n val_final=%d, subb_final=%d,kx=%le, energy=%le",valley_final,subband_final,kx,energy_final);
   } // End of else{ // energy_final >0.0
return;     
} // End of void isotropic(double Ef)
