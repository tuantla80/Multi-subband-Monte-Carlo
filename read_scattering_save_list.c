/* ************************************************************************
  To read parameters from the SCATTERING_LIST and SAVE_LIST in input.dat file

Starting date: March 11, 2010
Latest update: March 11, 2010 
*************************************************************************** */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "constants.h"
#include "nrutil.h"

int static num_scat_used =0; // Tinh so luong scattering mechanism su dung 
// For #SCATTERING_LIST
int static flag_ballistic_transport,flag_acoustic,flag_zero_order_optical;
int static flag_first_order_optical,flag_Coulomb, flag_Surface_roughness;

// For #SAVE_LIST
int static save_eig_wavefuntion;
int static save_doping_potential_init;
int static save_electron_initlization;
int static save_form_factor;  
int static save_scattering_rate; // tuc la scattering table 
int static save_electron_population; // Save N_s,v,i
int static save_electron_density;    //[1/m3]
int static save_VelEnerCurr_all_time_steps;// ca velocity energy va current
int static save_VelEnerCurr_after_transient;//  ca velocity energy va current
int static save_init_free_flight;
//******************************************************************************

void read_scattering_save_list(){
 
    FILE *f;
//  Read #SCATTERING_LIST
    f=fopen("Input.d","r");
    if(f==NULL){printf("Error: can't open file Input.d.\n"); return;}
    else { // printf("\n\n Reading #SCATTERING_LIST"); 
         }
      
    while(!feof(f)) {
         char str[180]; 
         fscanf(f,"%s",str);
         if(strcmp(str,"#SCATTERING_LIST")==0){
             fscanf(f,"%*s %*s %d ",&flag_ballistic_transport); 
             fscanf(f,"%*s %*s %d ",&flag_acoustic); 
             fscanf(f,"%*s %*s %d ",&flag_zero_order_optical);
             fscanf(f,"%*s %*s %d ",&flag_first_order_optical);
             fscanf(f,"%*s %*s %d ",&flag_Coulomb);
             fscanf(f,"%*s %*s %d ",&flag_Surface_roughness);
          } // End of if
    }// End of while
    fclose(f);
 
    //If flag_ballistic_transport=1 Ballistic; =0: Diffusive; //DEFAULT: Ballistic=0=>Diffusive Model
	if(flag_ballistic_transport==1) { //=1 Ballistic (NO Scattering)
	   // O ham emcd.c cung phan ra hai truong hop nhu the nay
           //printf("\n We use BALLISTIC TRANSPORT MODEL");
           flag_acoustic = 0; // set lai tat ca cac scattering flags bang 0 het
           flag_zero_order_optical = 0;
           flag_first_order_optical = 0;
           flag_Coulomb = 0;
           flag_Surface_roughness = 0;
         }
    else {// ballistic_transport !=1 -> Diffusive: Chon loai scattering
           //printf("\n We use DIFFUSIVE TRANSPORT MODEL");// Gia tri cac scattering_flags tu Input.d file
           if(flag_acoustic==1){
                  num_scat_used = num_scat_used +1;// so luong scattering mechanims tang len 1              
                  //printf("\n Acoustic phonon scattering is Included");
                }
           else { //printf("\n Acoustic phonon scattering is NOT included");
                }
           
           if(flag_zero_order_optical==1){
                  num_scat_used = num_scat_used +1;                        
                  //printf("\n Zero-order Non-polar Optical phonon is Included");
                }
           else { //printf("\n Zero-order Non-polar Optical phonon is NOT included");
                }
           
           if(flag_first_order_optical==1){
                  num_scat_used = num_scat_used +1;                         
                  //printf("\n First-order Non-polar Optical phonon is Included");
                }
           else { //printf("\n First-order Non-polar Optical phonon is NOT included");
                }
           
           if(flag_Coulomb==1){
                  num_scat_used = num_scat_used +1;             
                  //printf("\n Coulomb scattering is Included");
                }
           else { //printf("\n Coulomb scattering is NOT included");
                }
           
           if(flag_Surface_roughness==1){
                  num_scat_used = num_scat_used +1;                       
                  //printf("\n Surface roughness scattering is Included");
                }
           else { //printf("\n Surface roughness scattering is NOT included");
                }
           
           //printf("\n Number of Scattering Mechanims in used = %d",num_scat_used);
           if(num_scat_used > n_scatt_max){
              printf("\n Number of scattering mechanisms exceeds maximum number");
              printf("\n You can change n_scatt_max in constants.h file \n");
              nrerror("scattering mechanisms exceeds maximum number");
           }
      //else { printf("\n Number of scattering mechanisms using is = %d",i_count);}
   // }// End of if (doping_density[i_g]!=0.0)
     
         }// End of  else {// ballistic_transport !=1 -> Diffusive: Chon loai scattering
// END of  Read SCATTERING_LIST   

// Read #SAVE_LIST
    f=fopen("Input.d","r");
    if(f==NULL){printf("Error: can't open file Input.d.\n"); return;}
    else { //printf("\n Reading #SAVE_LIST");
         }
      
    while(!feof(f)) {
         char str[180]; 
         fscanf(f,"%s",str);
         if(strcmp(str,"#SAVE_LIST")==0){
             fscanf(f,"%*s %*s %d ",&save_eig_wavefuntion); 
             fscanf(f,"%*s %*s %d ",&save_doping_potential_init);
             fscanf(f,"%*s %*s %d ",&save_electron_initlization);
             fscanf(f,"%*s %*s %d ",&save_form_factor);
             fscanf(f,"%*s %*s %d ",&save_scattering_rate);
             fscanf(f,"%*s %*s %d ",&save_electron_population);
             fscanf(f,"%*s %*s %d ",&save_electron_density);
             fscanf(f,"%*s %*s %d ",&save_VelEnerCurr_all_time_steps);
             fscanf(f,"%*s %*s %d ",&save_VelEnerCurr_after_transient);
	     fscanf(f,"%*s %*s %d ",&save_init_free_flight);
         } // End of if
    }// End of while
      
    fclose(f);
    //printf("\n save_eig_wavefuntion             = %d",save_eig_wavefuntion);
    //printf("\n save_doping_potential_init       = %d",save_doping_potential_init);
    //printf("\n save_electron_initlization       = %d",save_electron_initlization);
    //printf("\n save_form_factor                 = %d",save_form_factor);
    //printf("\n save_scattering_rate             = %d",save_scattering_rate);
    //printf("\n save_electron_population         = %d",save_electron_population);
    //printf("\n save_electron_density            = %d",save_electron_density);
    //printf("\n save_VelEnerCurr_all_time_steps  = %d",save_VelEnerCurr_all_time_steps);
    //printf("\n save_VelEnerCurr_after_transient = %d",save_VelEnerCurr_after_transient);
    //printf("\n save_init_free_flight = %d",save_init_free_flight);

    //Phan printf ra Monitor o node co rank=0
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if(myrank ==0){ 
      printf("\n\n Reading #SCATTERING_LIST"); 
      if(flag_ballistic_transport==1) { 
	printf("\n We use BALLISTIC TRANSPORT MODEL");
      }
      else {// ballistic_transport !=1 -> Diffusive: Chon loai scattering
	printf("\n We use DIFFUSIVE TRANSPORT MODEL");// Gia tri cac scattering_flags tu Input.d file
	if(flag_acoustic==1){
	  printf("\n Acoustic phonon scattering is Included");
	}
	else { printf("\n Acoustic phonon scattering is NOT included");}
           
	if(flag_zero_order_optical==1){
	  printf("\n Zero-order Non-polar Optical phonon is Included");
	}
	else { printf("\n Zero-order Non-polar Optical phonon is NOT included");}
           
	if(flag_first_order_optical==1){
	  printf("\n First-order Non-polar Optical phonon is Included");
	}
	else { printf("\n First-order Non-polar Optical phonon is NOT included");}
           
	if(flag_Coulomb==1){
	  printf("\n Coulomb scattering is Included");
	}
	else { printf("\n Coulomb scattering is NOT included");}
           
	if(flag_Surface_roughness==1){
	  printf("\n Surface roughness scattering is Included");
	}
	else { printf("\n Surface roughness scattering is NOT included");}
           
	printf("\n Number of Scattering Mechanims in used = %d",num_scat_used);
	if(num_scat_used > n_scatt_max){
	  printf("\n Number of scattering mechanisms exceeds maximum number");
	  printf("\n You can change n_scatt_max in constants.h file \n");
	  nrerror("scattering mechanisms exceeds maximum number");
	}
      }
     }//End of if(myrank ==0)
    return;
// End of Read #SAVE_LIST
}// End of void read_scattering_save_list() 
// *****************************************************************************
int Get_num_scat_used(){
    return num_scat_used;
}

int Get_flag_ballistic_transport(){
    return flag_ballistic_transport;
}

int Get_flag_acoustic(){
    return flag_acoustic;
}

int Get_flag_zero_order_optical(){
    return flag_zero_order_optical;
}

int Get_flag_first_order_optical(){
    return flag_first_order_optical;
}

int Get_flag_Coulomb(){
    return flag_Coulomb;
}

int Get_flag_Surface_roughness(){
    return flag_Surface_roughness;
}

int Get_save_eig_wavefuntion(){
    return save_eig_wavefuntion;
}

int Get_save_doping_potential_init(){
    return save_doping_potential_init;
}

int Get_save_electron_initlization(){
    return save_electron_initlization;
}

int Get_save_form_factor(){
    return save_form_factor;
}

int Get_save_scattering_rate(){
    return save_scattering_rate;
}

int Get_save_electron_population(){
    return save_electron_population;
}

int Get_save_electron_density(){
    return save_electron_density;
}

int Get_save_VelEnerCurr_all_time_steps(){
    return save_VelEnerCurr_all_time_steps;
}

int Get_save_VelEnerCurr_after_transient(){
    return save_VelEnerCurr_after_transient;
}

int Get_save_init_free_flight(){
  return save_init_free_flight;
}



