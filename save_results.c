/* ****************************************************************
 De SAVE ket qua <---KIEM TRA duoc KET QUA DUNG hay SAI
 
Starting date: March 11, 2010
Latest update: March 12, 2010
****************************************************************** */
#include <stdio.h>
#include "nrutil.h"
#include "constants.h"

void save_results(){
     int Get_nx0(),Get_nx1(),Get_NSELECT();
     // Cac ham save 
         //void save_eig_wave( NEN de o ham Solved_2D_Schro_for_MSMC() vi no lien quan den so valley can dung, vv
     void save_doping_potential_initialization();
     void save_form_factor_calculation(int s_th); // tai 1 section cu the
     void save_scattering_table(int s_th); // // tai 1 section cu the
         // Cac ham flag de SAVE hay KHONG
     int Get_save_doping_potential_init();
     int Get_save_form_factor();
     int Get_save_scattering_rate();
     
     // Cac bien local
     int nx0 = Get_nx0();
     int nx1 = Get_nx1();
     int NSELECT =  Get_NSELECT();
     int save_doping_potential_init = Get_save_doping_potential_init();
     int save_form_factor = Get_save_form_factor(); 
     int save_scattering_rate = Get_save_scattering_rate();
     
     // Thuc hien
     int s_th = (int)((nx1-nx0)/2); // chon section o GIUA 
     // 1. Save Doping and Potential Initialization  
        if(save_doping_potential_init==1){
	    save_initial_potential_doping();
            printf("\n Doping and Potential Initialization are Saved");
        }
     
     // 2.  Save Form factor  
        if(save_form_factor==1){ 
             save_form_factor_calculation(s_th);                    
             printf("\n Form factor is Saved");
         }
     
     // 3. SAVE Scattering rate    
        if(save_scattering_rate==1){
	     save_scattering_table(s_th);                        
             printf("\n Scattering rate is Saved");
        }
    
   return;
} // End of void save_results(){

     
