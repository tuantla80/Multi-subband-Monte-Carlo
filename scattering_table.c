/* ****************************************************************
 TO CREATE THE SCATTERING TABLE and then NORMALIZE
   Scattering mechanism:
     - Acoustic phonons
     - Zero-order Non-polar optical phonons 
     - First-order Non-polar optical phonons  (not yet)
     - Coulomb scattering (not yet)
     - Interface roughness scattering (not yet)
CHU Y 1:       
  flag_mech = 1 ==> Isotropic scattering process (acoustic or g-process) vi o g-process
                    no xay ra o SAME axis nen cung xep vao loai nay
  flag_mech = 2 ==> Isotropic cho f-process (xay ra o DIFFRENT axis)
  flag_mech = 3 ==> Anistropic (vi du Coulomb scattering)
  
CHU Y 2: Ro rang scattering se KHAC nhau o region KHAC nhau (vi du Coulomb scattering)
         Vay o MANG scat_table can 1 index cho region. (March 12, 2010: CHUA them vao)
 
Starting date: March 11, 2010
Update:        March 12, 2010
Update:        March 23, 2010. SUA lai dang g-process
Latest update: April 15, 2010. Gop ham  Density_of_State_1D() vao ham scattering_table()
               De khong can su dung mang DOS1D => giam memory can dung khi chay chuong trinh
****************************************************************** */
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "constants.h"

void scattering_table(){
     // Goi ham
     double *****Get_scat_table();
     double ****Get_form_factor();
     double ***Get_eig();
     double Get_ml(),Get_mt(),Get_emax();
     int Get_NSELECT(); // chay cho section va subband
     double Get_nonparabolicity_factor(),Get_hw0f_phonon(),Get_hw0g_phonon();
         // Goi cac ham flag cho scattering 
     int Get_flag_acoustic(),Get_flag_zero_order_optical();
     int Get_flag_first_order_optical(),Get_flag_Coulomb(),Get_flag_Surface_roughness();
         // Goi cac ham hang so
     double Get_constant_acoustic();
     double Get_constant_optical_g_e(),Get_constant_optical_g_a();
     double Get_constant_optical_f_e(),Get_constant_optical_f_a();
     
     // Bien local
     double *****scat_table = Get_scat_table();
     double ****form_factor = Get_form_factor();
     double ***eig = Get_eig();
     double ml = Get_ml(); // chi la gia tri ml, CHUA nhan voi m0
     double mt = Get_mt();
     double emax = Get_emax();

     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
     
     int NSELECT = Get_NSELECT();
     double af = Get_nonparabolicity_factor(); //printf("\n Nonparabolicity factor[1/eV]   = %f",af); getchar();
     double hw0f_phonon = Get_hw0f_phonon();
     double hw0g_phonon = Get_hw0g_phonon();
          // Cac bien cho flag: ADD cac loai Scattering NAO ? 
     int flag_acoustic = Get_flag_acoustic();
     int flag_zero_order_optical = Get_flag_zero_order_optical();
     int flag_first_order_optical = Get_flag_first_order_optical();
     int flag_Coulomb = Get_flag_Coulomb();
     int flag_Surface_roughness = Get_flag_Surface_roughness();
         // Cac bieb cho hang so
     double constant_acoustic    = Get_constant_acoustic();
     double constant_optical_g_e = Get_constant_optical_g_e();
     double constant_optical_g_a = Get_constant_optical_g_a();
     double constant_optical_f_e = Get_constant_optical_f_e();
     double constant_optical_f_a = Get_constant_optical_f_a(); 

// Cu moi lan tinh scattering table thi ta reset no cho chac chan
     int s,v,n,m,e_step;
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

// Thuc hien
     double de=emax/(double)(n_lev); // energy interval
     double *E; E=dvector(1,n_lev); // [eV] Dinh nghia dai nang luong minh quan tam E di tu 1*de den "emax"
                                   // Neu di tu 0 thi se co loi. The hien phan kinetic cua electron truoc scattering
     
     
     for(e_step=1; e_step<=n_lev; e_step++){ 
       E[e_step] = e_step*de;
     }
     double Ef = 0.0; // Ef trong cong thuc tinh DOS 
     double theta = 0.0; // la de tinh ham step function trong cong thuc tinh DOS; 
     double DOS_temp = 0.0;//DOS_temp la bien trung gian tinh DOS tranh dung mang nhu o version truoc
     double const_DOS_ml, const_DOS_mt; // Thanh phan dau tien cua DOS1D. Xem cong thuc 1 trang 93 
     const_DOS_ml = sqrt(ml*m0)/(2.0*sqrt(2.0)*pi*hbar*sqrt(q));//printf("\n const_DOS_ml = %le", const_DOS_ml);
     const_DOS_mt = sqrt(mt*m0)/(2.0*sqrt(2.0)*pi*hbar*sqrt(q));//printf("\n const_DOS_mt = %le", const_DOS_mt); getchar();
    
// I. Acoustic phonons scattering (Intra-Valley (Inter- and Intra-subband))
if(flag_acoustic==1){ // Include: Acoustic phonon <--Tinh gia tri scattering nay

    //m* = ml: Cap valley pair 1. Xem page 84
    //v =1; // Cap valley pair 1, const_DOS_ml = sqrt(ml*m0)/(2.0*sqrt(2.0)*pi*hbar*sqrt(q));
    for(s=0; s<=nx_max; s++){ // Cho moi section slice
        for(n=1; n<=NSELECT; n++){ // Scattering tu subband n-th   
            for(m=1; m<=NSELECT; m++){ // den subband m-th 
               for(e_step=1; e_step<=n_lev; e_step++){// cac buoc nang luong
		  Ef = eig[s][1][n] - eig[s][1][m] + E[e_step];// [eV] printf("\n Ef = %le", Ef); // For acoustic phonon
		 
                   if(Ef > 0.0) { //ham step function -> Ef<0 thi theta=0, Ef>0 thi theta=1
                     theta = 1.0;
                     DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                     scat_table[s][n][m][e_step][1] = constant_acoustic*form_factor[s][n][m][1]*DOS_temp;
		     //printf("\n scat_table[s][n][m][e_step][1]  = %le",scat_table[s][n][m][e_step][1]); getchar();
                     }
	               else {// Ef <0.0
                         //theta =0.0; DOS_temp =0.0
                         scat_table[s][n][m][e_step][1] = 0.0;
                      }  
               }// End of for(e_step=1; e_step<=n_lev; e_step++)
             } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(int s=0; s<=nx_max; s++){ // Cho moi section slice 
   
    //m* = mt: Cap valley pair 2 va 3. Xem page 84, const_DOS_mt = sqrt(mt*m0)/(2.0*sqrt(2.0)*pi*hbar*sqrt(q));
    for(s=0; s<=nx_max; s++){ // Cho moi section slice
        for(v=2; v<=3; v++){// v=2 va v=3: 2 cap valley pairs co cung mt
            for(n=1; n<=NSELECT; n++){ // Scattering tu subband n-th   
                for(m=1; m<=NSELECT; m++){ // den subband m-th 
                   for(e_step=1; e_step<=n_lev; e_step++){// cac buoc nang luong
		       Ef = eig[s][v][n] - eig[s][v][m] + E[e_step];// [eV]  // For acoustic phonon
		       
                       if(Ef > 0.0){ 
                          theta = 1.0;
                          DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                          scat_table[s][n][m][e_step][v] = constant_acoustic*form_factor[s][n][m][v]*DOS_temp;
                        }
	                   else{// theta =0.0; DOS_temp = 0.0
                             scat_table[s][n][m][e_step][v] = 0.0;
                           }
                 }// End of for(e_step=1; e_step<=n_lev; e_step++)
             } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
       }// for(v=2; v<=3; v++)
     } // End of for(s=0; s<=nx_max; s++){ // Cho moi section slice 
} // End of if(flag_acoustic==1)  
// End of Part I. Acoustic phonons scattering rate                     

/* **************************************************************************** 
    Part II. Non-polar Optical phonon (Inter-Valley (Inter- and Intra-subband))
    f-process: tu Valley NAY den Valley KIA  KHAC truc AXIS
    g-process: tu Valley NAY den Valley KIA  CUNG truc Axis
*******************************************************************************/
if(flag_zero_order_optical ==1){
  // II.1. f-process
  // II.1.1. Tu valley pair 1 (v=1) --> Valley pair 2 (v=2): m* la mt
     // Chi so cho form factor la 4, chi so cho scat_table la tu 4-5
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][1][n] - eig[s][2][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0) { 
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][4] = constant_optical_f_a*
                                                       form_factor[s][n][m][4]*
                                                       DOS_temp;
                      //printf("\n scat_table[s][n][m][e_step][4]  = %le",scat_table[s][n][m][e_step][4]); getchar();
                      }
	               else { //theta =0.0; DOS_temp = 0.0
                          scat_table[s][n][m][e_step][4] = 0.0;
                        } 
                 
                   //f-process EMISSION
                   Ef = eig[s][1][n] - eig[s][2][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){ 
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][5] = constant_optical_f_e*
                                                       form_factor[s][n][m][4]*
                                                       DOS_temp;
                      }
	               else{ //theta =0.0; DOS_temp = 0.0
                         scat_table[s][n][m][e_step][5] = 0.0;
                        } 
                  
              }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(int s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.1.1. Tu valley 1 (v=1) --> Valley 2 (v=2): m* la mt

// II.1.2. Tu valley 1 (v=1) --> Valley 3 (v=3): m* la mt
    // Chi so cho form factor la 5, chi so cho scat_table la tu 6-7
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][1][n] - eig[s][3][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][6] = constant_optical_f_a*
                                                       form_factor[s][n][m][5]*
                                                       DOS_temp;
                      }
	               else{//theta =0.0;DOS_temp=0.0
	                    scat_table[s][n][m][e_step][6]=0.0;        
                       }
              
                   //f-process EMISSION
                   Ef = eig[s][1][n] - eig[s][3][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][7] = constant_optical_f_e*
                                                       form_factor[s][n][m][5]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][7]=0.0;
                        }
               }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(int s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.1.2. Tu valley 1 (v=1) --> Valley 3 (v=3): m* la mt

// II.1.3. Tu valley 2 (v=2) --> Valley 1 (v=1): m* la ml
   // Chi so cho form factor la 6, chi so cho scat_table la tu 8-9
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][2][n] - eig[s][1][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][8] = constant_optical_f_a*
                                                       form_factor[s][n][m][6]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
	                    scat_table[s][n][m][e_step][8]=0.0;
                       }
   
                   //f-process EMISSION
                   Ef = eig[s][2][n] - eig[s][1][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][9] = constant_optical_f_e*
                                                       form_factor[s][n][m][6]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
	                    scat_table[s][n][m][e_step][9]=0.0;
                      } 
               }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.1.3. Tu valley 2 (v=2) --> Valley 1 (v=1): m* la ml

// II.1.4. Tu valley 2 (v=2) --> Valley 3 (v=3): m* la mt
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][2][n] - eig[s][3][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][10] = constant_optical_f_a*
                                                        form_factor[s][n][m][7]*
                                                        DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
	                    scat_table[s][n][m][e_step][10]=0.0;
                       }
 
                   //f-process EMISSION
                   Ef = eig[s][2][n] - eig[s][3][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][11] = constant_optical_f_e*
                                                        form_factor[s][n][m][7]*
                                                        DOS_temp;
                       }
	               else{ //theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][11]=0.0;
                       } 
               }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(int s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.1.4. Tu valley 2 (v=2) --> Valley 3 (v=3): m* la mt

// II.1.5. Tu valley 3 (v=3) --> Valley 1 (v=1): m* la ml
   // Chi so cho form factor la 8, chi so cho scat_table la tu 12-13
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][3][n] - eig[s][1][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][12] = constant_optical_f_a*
                                                        form_factor[s][n][m][8]*
                                                        DOS_temp;
                       }
	               else{ //theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][12]=0.0;
                       } 
          
                   //f-process EMISSION
                   Ef = eig[s][3][n] - eig[s][1][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][13] = constant_optical_f_e*
                                                        form_factor[s][n][m][8]*
                                                        DOS_temp;
                       }
	               else{// theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][13]=0.0;    
                       } 
               }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.1.5. Tu valley 3 (v=3) --> Valley 1 (v=1): m* la ml

// II.1.6. Tu valley 3 (v=3) --> Valley 2 (v=2): m* la mt
    // Chi so cho form factor la 9, chi so cho scat_table la tu 14-15
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                   // f-process ABSORPTION         
                   Ef = eig[s][3][n] - eig[s][2][m] + E[e_step]+ hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][14] = constant_optical_f_a*
                                                        form_factor[s][n][m][9]*
                                                        DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][14]=0.0;
                       } 
         
                   //f-process EMISSION
                   Ef = eig[s][3][n] - eig[s][2][m] + E[e_step] - hw0f_phonon;
                   if(Ef > 0.0){
                      theta = 1.0;
                      DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                      scat_table[s][n][m][e_step][15] = constant_optical_f_e*
                                                        form_factor[s][n][m][9]*
                                                        DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][15]=0.0;
                       } 
              }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(s=0; s<=nx_max; s++){ // Cho moi section slice   
 // End of II.1.6. Tu valley 3 (v=3) --> Valley 2 (v=2): m* la mt
 // End of II.1
 
// II.2. g-process
// II.2.1. Tu valley 1 (v=1) --> Valley 1' (v=1): m* la ml
    // Chi so cho form factor la 1, chi so cho scat_table la tu 16-17
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                 //g-process ABSORPTION
                   Ef = eig[s][1][n] - eig[s][1][m] + E[e_step] + hw0g_phonon;
                   if(Ef > 0.0){
                     theta = 1.0;
                     DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                     scat_table[s][n][m][e_step][16] = constant_optical_g_a*
                                                       form_factor[s][n][m][1]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp = 0.0
                        scat_table[s][n][m][e_step][16]=0.0;
                        } 
              
                   //g-process EMISSION
                   Ef = eig[s][1][n] - eig[s][1][m] + E[e_step] - hw0g_phonon;
                   if(Ef > 0.0){
                       theta = 1.0;
                       DOS_temp = const_DOS_ml*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                       scat_table[s][n][m][e_step][17] = constant_optical_g_e*
                                                         form_factor[s][n][m][1]*
                                                         DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
                        scat_table[s][n][m][e_step][17]=0.0;
                       } 
              }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for(s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.2.1. Tu valley 1 (v=1) --> Valley 1' (v=1): m* la ml

// II.2.2. Tu valley 2 (v=2) --> Valley 2' (v=2): m* la mt
    // Chi so cho form factor la 2, chi so cho scat_table la tu 18-19
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                 //g-process ABSORPTION
                   Ef = eig[s][2][n] - eig[s][2][m] + E[e_step] + hw0g_phonon;
                   if(Ef > 0.0){
                     theta = 1.0;
                     DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                     scat_table[s][n][m][e_step][18] = constant_optical_g_a*
                                                       form_factor[s][n][m][2]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
                        scat_table[s][n][m][e_step][18]=0.0;
                        } 
              
                   //g-process EMISSION
                   Ef = eig[s][2][n] - eig[s][2][m] + E[e_step] - hw0g_phonon;
                   if(Ef > 0.0){
                       theta = 1.0;
                       DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                       scat_table[s][n][m][e_step][19] = constant_optical_g_e*
                                                         form_factor[s][n][m][2]*
                                                         DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
                        scat_table[s][n][m][e_step][19]=0.0;
                        } 
              }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for( s=0; s<=nx_max; s++){ // Cho moi section slice   
// End of II.2.2. Tu valley 2 (v=2) --> Valley 2' (v=2): m* la mt

// II.2.3. Tu valley 3 (v=3) --> Valley 3' (v=3): m* la mt
    for(s=0; s<=nx_max; s++){ 
        for(n=1; n<=NSELECT; n++){ 
            for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){
                             
                 //g-process ABSORPTION
                   Ef = eig[s][3][n] - eig[s][3][m] + E[e_step] + hw0g_phonon;
                   if(Ef > 0.0){
                     theta = 1.0;
                     DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                     scat_table[s][n][m][e_step][20] = constant_optical_g_a*
                                                       form_factor[s][n][m][3]*
                                                       DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
                        scat_table[s][n][m][e_step][20]=0.0;
                        } 
              
                   //g-process EMISSION
                   Ef = eig[s][3][n] - eig[s][3][m] + E[e_step] - hw0g_phonon;
                   if(Ef > 0.0){
                       theta = 1.0;
                       DOS_temp = const_DOS_mt*theta*((1+2.0*af*Ef)/sqrt(Ef*(1+af*Ef)));
                       scat_table[s][n][m][e_step][21] = constant_optical_g_e*
                                                         form_factor[s][n][m][3]*
                                                         DOS_temp;
                       }
	               else{//theta =0.0; DOS_temp=0.0
                        scat_table[s][n][m][e_step][21]=0.0;
                        } 
              }  // End of for(e_step=1; e_step<=n_lev; e_step++)
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     } // End of for( s=0; s<=nx_max; s++){ // Cho moi section slice   
   // End of II.2.1. Tu valley 3 (v=3) --> Valley 3' (v=3): m* la mt
 // End of Part II.2
// End of Part II
} // End of if(flag_zero_order_optical ==1)
//*****************************************************************************
  // Free local variables
  free_dvector(E,1,n_lev);   
  return;
}// End of void scattering_table()
/* ********************************************************************************* */

/* ****************************************************************
 De SAVE Scattering rate vao trong 1 file. Chi Save tai 1 section thoi
 INPUT:  + s_th: s-th section. Chi save tai section s_th thoi - so luong tu 0 den nx_max
         + valley pair - so luong tu 1 den 3
         + NSELECT: subband index n (1 den NSELECT) den m (1 den NSELECT)
         + n_lev: So enerfy steps: tu constants.h
         + scat_table[s][n][m][e_step][scatName]
                  
 OUTPUT: + SAVE scat_table[s][n][m][e_step][scatName] tai 1 section TONG 3 valley pair xac dinh 
           Vao file Scattering_rate.dat
          gom Acoustic va Nonpolar Optical Phonon
        
Starting date: March 13, 2010
Update:        March 13, 2010
Latest update: March 23, 2010 De kiem tra scattering rate tai cac valley
****************************************************************** */
void save_scattering_table(int s_th, char *fn){
     // Goi ham
     double *****Get_scat_table(), Get_emax();
     int Get_NSELECT();
     // Local bien 
     double *****scat_table = Get_scat_table();
     double emax = Get_emax();// printf("\n emax=%f, n_lev=%d",emax, n_lev);
     
     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    
     int NSELECT = Get_NSELECT();
     
     // Thuc hien
     if (s_th < 0){ s_th = 0;}// Kiem tra dau vao s_th de khong vuot qua chi so cua cross-section
     if (s_th > nx_max) { s_th = nx_max; }
     
     double *Scat_Acous_All_Valley = dvector(1,n_lev); // bien trung gian tinh Acoustic la tong cac o cac subband va valley
     double *Scat_Acous_Valley_1 = dvector(1,n_lev); //tinh Acoustic la tong cac o cac subband o Valley pair 1
     double *Scat_Acous_Valley_2 = dvector(1,n_lev); //tinh Acoustic la tong cac o cac subband o Valley pair 2
     double *Scat_Acous_Valley_3 = dvector(1,n_lev); //tinh Acoustic la tong cac o cac subband o Valley pair 3
     
     double *Scat_Optical_All_Valley = dvector(1, n_lev);
     //double *Scat_Optical_Valley_1 = dvector(1, n_lev);
     //double *Scat_Optical_Valley_2 = dvector(1, n_lev);
     //double *Scat_Optical_Valley_3 = dvector(1, n_lev);
     
     int i; 
     for(i=1; i<=n_lev; i++){// KHONG HIEU TAI sao neu khong khoi tao thi no lai lay gia tri cua DOS
         Scat_Acous_All_Valley[i]=0.0;   Scat_Acous_Valley_1[i]=0.0;
         Scat_Acous_Valley_2[i]=0.0;     Scat_Acous_Valley_3[i]=0.0;
         Scat_Optical_All_Valley[i]=0.0; //Scat_Optical_Valley_1[i]=0.0;
         //Scat_Optical_Valley_2[i]=0.0; //  Scat_Optical_Valley_3[i]=0.0;
        } 
     
     double de=emax/(double)(n_lev); // energy interval
     FILE *f; f = fopen(fn,"w"); // "w" neu tep ton tai no se bi xoa
     //FILE *f; f = fopen("Scat_table.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file at save_scattering_table.c"); return ; }
     fprintf(f,"\n #ener_step  AcouVall_1 AcouVall_2 AcouVall_3 AcouAllValley OpticalAllValley  s_th \n");
     //fprintf(f,"\n ener_step  AcouVall_1 AcouVall_2 AcouVall_3 AcouAllValley OptVall_1 OptVall_2 OptVall_3 OptAllVall  s_th \n");
     int n,m,e_step,v;//
     for(e_step=1; e_step<=n_lev; e_step++){
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                 // Tinh cho Acoustic           
                 Scat_Acous_Valley_1[e_step]=Scat_Acous_Valley_1[e_step]+scat_table[s_th][n][m][e_step][1]; 
                 Scat_Acous_Valley_2[e_step]=Scat_Acous_Valley_2[e_step]+scat_table[s_th][n][m][e_step][2];
                 Scat_Acous_Valley_3[e_step]=Scat_Acous_Valley_3[e_step]+scat_table[s_th][n][m][e_step][3];
                 for(v=1; v<=3; v++){// Acoustic                   
                    Scat_Acous_All_Valley[e_step]= Scat_Acous_All_Valley[e_step]+ scat_table[s_th][n][m][e_step][v];
                 }  // Nghia la Tong o 3 cap valley
                
                // Tinh cho Non-polar Optical Phonon
                for(v=4; v<=21; v++){ // Optical. KINH NGHIEM: it nhat gia tri 21 nay can 1 cai TEN
                    Scat_Optical_All_Valley[e_step] = Scat_Optical_All_Valley[e_step] + scat_table[s_th][n][m][e_step][v];            
                 } 
            } // End of for(m=1; m<=NSELECT; m++){
         } // End of for(n=1; n<=NSELECT; n++)
     }// End of for(e_step=1; e_step<=n_lev; e_step++)
    
    for(e_step=1; e_step<=n_lev; e_step++){ // Ghi vao file
        fprintf(f,"%le %le %le %le %le %le  %d  \n",
                        e_step*de,
                        Scat_Acous_Valley_1[e_step],
                        Scat_Acous_Valley_2[e_step],
                        Scat_Acous_Valley_3[e_step],
                        Scat_Acous_All_Valley[e_step],
                        Scat_Optical_All_Valley[e_step],
                        s_th); 
        }
        
   fclose(f); 
  
    // Free local matrix
    free_dvector(Scat_Acous_Valley_1,1,n_lev);
    free_dvector(Scat_Acous_Valley_2,1,n_lev);
    free_dvector(Scat_Acous_Valley_3,1,n_lev);
    free_dvector(Scat_Acous_All_Valley,1,n_lev);
    free_dvector(Scat_Optical_All_Valley,1,n_lev);
    
    return; 
} // End of void save_scattering_table(int s_th){

/* *********************************************************************************
    To normalize the scattering table for EACH SECTION s-th
       - choose the maximum of scattering rate
       - normalize of scattering rate
       
Starting date: March 12, 2010
Latest update: March 16, 2010
*********************************************************************************** */
void normalize_table(){
     // Goi ham 
     double *****Get_scat_table();
     double *Get_max_gm();
     int  Get_NSELECT();
          
     // Cac bien local
     double *****scat_table = Get_scat_table();
     double *max_gm = Get_max_gm();
     
     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    
     int NSELECT = Get_NSELECT();

// Reset max_gm truoc moi lan chay
     int s;
     for(s=nx0; s<=nx1; s++){ 
       max_gm[s] = 0.0;
     }
         
     // Thuc hien
     double **gm = dmatrix(0,nx_max,1,n_lev);// gm la gamma day. Do ta se gan GIAN TIEP gm nen CAN KHOI TAO
     int i,j; 
     for(i=0; i<=nx_max; i++){
         for(j=1; j<=n_lev; j++){ 
	   gm[i][j] = 0.0; 
	 }
     }
                  
     double temp = 0.0; // la bien de tinh tong
     double *****temp_scat = d5matrix(0,nx_max,1,NSELECT,1,NSELECT,1,n_lev,1,21);// KINH NGHIEM it nhat thi gia tri 21 nay CAN CO 1 CAI TEN
     // La BIEN TRUNG GIAN cho scat_table NO PHAI MATCH HOAN TOAN voi bien scat_table
     int i1,i2,i3,i4,i5;
     for(i1=0; i1<=nx_max; i1++){
	     for(i2=1; i2<=NSELECT; i2++){
             for(i3=1; i3<=NSELECT; i3++){
                 for(i4=1; i4<=n_lev; i4++){
                     for(i5=1; i5<=21; i5++){
                         temp_scat[i1][i2][i3][i4][i5]=0.0;
                  }
                }
              }
           }
       } // End of for(i1=0; i1<=nx_max; i1++)
          
     int v,n,m,e_step;
     // NOTE: Luc da chon IT NHAT 1 loai SCATTERING thi mechanism cua SCATTERING DO luon lon hon 1.
     // Vi du: Acoustic cung co 3 index cho 3 kieu intra-Valley
     
     // Buoc 1: tinh tong scattering rate tai s-th section va energy_step
      for(s=0; s<=nx_max; s++){ 
         for(e_step=1; e_step<=n_lev; e_step++){ 
             temp = 0.0; 
             for(v=1; v<=21; v++){//v=1 den 3: Acoustic; v=4 den 21: Zero-order Optical 
                 // NOTE: phai chay cho v TRUOC roi moi den n va m thi moi tao ra bang trang 107 duoc 
                 for(n=1; n<=NSELECT; n++){ 
                     for(m=1; m<=NSELECT; m++){ 
                         gm[s][e_step] = gm[s][e_step] + scat_table[s][n][m][e_step][v];
                         temp = temp + scat_table[s][n][m][e_step][v];
                         temp_scat[s][n][m][e_step][v]= temp;   
                     }
                  }
              }
         }
      } // End of for(s=0; s<=nx_max; s++)
      
      // Buoc 2: Tinh gia tri gm max tai moi section
      for(s=0; s<=nx_max; s++){ 
          max_gm[s] = gm[s][n_lev]; // gan cho gia tri tai diem energy cuoi cung
          for(e_step=1; e_step<=n_lev-1; e_step++){ 
              if(max_gm[s]< gm[s][e_step]){
                  max_gm[s] = gm[s][e_step];
                }
           } 
          //printf("\n Maximum scatering rate at %d section is %le",s,max_gm[s]); //Da check: Ket qua DUNG
          
      }// End of for(s=0; s<=nx_max; s++) 
          
      // Buoc 3: Normalized scattering table with max_gm[s]
     for(s=0; s<=nx_max; s++){ 
        for(e_step=1; e_step<=n_lev; e_step++){     
            for(n=1; n<=NSELECT; n++){ 
                for(m=1; m<=NSELECT; m++){ 
                   for(v=1; v<=21; v++){
                      scat_table[s][n][m][e_step][v]=temp_scat[s][n][m][e_step][v];//Matchinh scat_table lay tu temp_scat
                      scat_table[s][n][m][e_step][v] = scat_table[s][n][m][e_step][v]/max_gm[s];
		      //printf("\n Normalized scatering rate at %d section is %f",s,scat_table[s][n][m][e_step][v]); getchar(); 
                      }
                  }
              }
          }
      } // End of for(s=0; s<=nx_max; s++)
 // Free cac bien local 
    free_dmatrix(gm,0,nx_max,1,n_lev);
    free_d5matrix(temp_scat,0,nx_max,1,NSELECT,1,NSELECT,1,n_lev,1,21);
    return;
}

/* ****************************************************************
 De SAVE normalize table  vao trong 1 file. Chi Save tai tat ca section thoi
 INPUT:  + s_th: s-th section. Chi save tai section s_th thoi - so luong tu 0 den nx_max
         + valley pair - so luong tu 1 den 3
         + NSELECT: subband index n (1 den NSELECT) den m (1 den NSELECT)
         + n_lev: So enerfy steps: tu constants.h
         + scat_table[s][n][m][e_step][scatName] : cai da Normalized NHE
                  
 OUTPUT: + SAVE scat_table[s][n][m][e_step][scatName] tai 1 section TONG 3 valley pair xac dinh 

 Chi show neu gia tri cua no    1.0e-6 > gia tri >= 1.
 Y dinh cua ta la biet duoc tai moi section, o subband nao den subband nao va energy step la bao nhieu thi gia tri normalized la 1 
       
Starting date: June 03, 2010
Latest update: June 03, 2010 De kiem tra 
****************************************************************** */
void save_normalized_table(){
     // Goi ham
     double *****Get_scat_table();// Bang scat da DUOC normalized
     double *****scat_table = Get_scat_table();

     double Get_emax();
     double emax = Get_emax();
     
     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; 
    
     int Get_NSELECT();
     int NSELECT = Get_NSELECT();
        
     FILE *f; f = fopen("normalized_table.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file normalized_table.dat"); return ; }
     fprintf(f,"\n # Section_sth Subband_n Subband_m  e_step normalized  \n");
     
     int s, n,m,e_step,v;
     for(s=0; s<=nx_max; s++){ 
       for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(e_step=1; e_step<=n_lev; e_step++){ 
                   for(v=1; v<=21; v++){  
		     // NOTE if 
		     if((scat_table[s][n][m][e_step][v]<1.0e-6)||(scat_table[s][n][m][e_step][v]>=1))//Chung to xem no co o trong KHOANG (0,1)
		       { 
			 fprintf(f,"%d %d %d %d %le \n ",s, n,m, e_step, scat_table[s][n][m][e_step][v]);
		       }
                   }
                }
             }
         }
     }
   
   fclose(f); 
  return; 
} // End of void save_normalized_table(int s_th)
