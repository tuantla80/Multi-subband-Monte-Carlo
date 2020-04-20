/* ************************************************************************
    To read parameters in the SIMULATION_LIST part from input.dat file

Starting date: Feb 23, 2010
Latest update: Feb 23, 2010 
*************************************************************************** */
#include <stdio.h>
#include <string.h>
#include "mpi.h"

double static dt,tot_time,transient_time,length_plotted,depth_plotted,width_plotted;
int static NSELECT;// Number of subbands for each valley
int static NTimeSteps; // Number of time steps to solve 2D Schrodinger Eq. De giam thoi gian cho 2D Schrodinger 
double static artificial_factor_cell_volume;

void read_simulation_list(){
          
      /* ******** Read simulation list from Input.d file ********* */
      FILE *f; //char s[180];
      f=fopen("Input.d","r");
      if(f==NULL){printf("Error: can't open file Input.d \n"); return;}
      //printf("\n\n Reading #SIMULATION_LIST");
      
   while (!feof(f)) {
      char str[180]; 
      fscanf(f,"%s",str);
      if(strcmp(str,"#SIMULATION_LIST")==0){
         fscanf(f,"%*s %*s %le",&dt); 
         //printf("\n dt...MC_time_step[fs]= %3.2f",dt);
         dt = dt*1.0e-15; // chuyen [fs] -> [s]
      
         fscanf(f,"%*s %*s %le ",&tot_time);
         //printf("\n Total time[ps]= %3.8f",tot_time); 
         tot_time = tot_time*1.0e-12; // chuyen [ps] -> [s]
      
         fscanf(f,"%*s %*s %le ",&transient_time);
         //printf("\n Transient time after which results are calculated[ps]= %3.8f",transient_time); 
         transient_time = transient_time*1.0e-12; // chuyen [ps] -> [s]
         
         fscanf(f,"%*s %*s %le ",&length_plotted);
         //printf("\n Length where the 3D plots are plotted[nm]= %3.2f",length_plotted);
         length_plotted = length_plotted*1.0e-9; // chuyen [nm] -> [m]
      
         fscanf(f,"%*s %*s %le ",&depth_plotted);
         //printf("\n Depth where the 3D plots are plotted[nm]= %3.2f",depth_plotted);
         depth_plotted = depth_plotted*1.0e-9;  // chuyen [nm] -> [m]
         
         fscanf(f,"%*s %*s %le ",&width_plotted);
         //printf("\n Width where the 3D plots are plotted[nm]= %3.2f",width_plotted);
         width_plotted = width_plotted*1.0e-9;  // chuyen [nm] -> [m] 
         
         fscanf(f,"%*s %*s %d ",&NSELECT);
         //printf("\n Number of Subbands wanted for 2D Schrodinger = %d",NSELECT);
         
         fscanf(f,"%*s %*s %d ",&NTimeSteps);
         //printf("\n Number of time steps to solve 2D Schrodinger = %d",NTimeSteps);
         
         fscanf(f,"%*s %*s %le ",&artificial_factor_cell_volume);
         //printf("\n Artificial factor of cell volume = %le",artificial_factor_cell_volume);
                
     } // End of if
  }// End of while
     fclose(f);

   //Phan printf ra Monitor o node co rank=0
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if(myrank ==0){ 
       printf("\n\n Reading #SIMULATION_LIST");
       printf("\n dt...MC_time_step[fs]= %3.2f",dt/1.0e-15);
       printf("\n Total time[ps]= %3.8f",tot_time/1.0e-12); 
       printf("\n Transient time after which results are calculated[ps]= %3.8f",transient_time/1.0e-12); 
       printf("\n Length where the 3D plots are plotted[nm]= %3.2f",length_plotted/1.0e-9);
       printf("\n Depth where the 3D plots are plotted[nm]= %3.2f",depth_plotted/1.0e-9);
       printf("\n Width where the 3D plots are plotted[nm]= %3.2f",width_plotted/1.0e-9);
       printf("\n Number of Subbands wanted for 2D Schrodinger = %d",NSELECT);
       printf("\n Number of time steps to solve 2D Schrodinger = %d",NTimeSteps);
       printf("\n Artificial factor of cell volume = %le",artificial_factor_cell_volume);
    }// End of if(myrank ==0)

return;
} 
// ********* End of void read_simulation_list() *************
  
double Get_dt(){
       return dt;
       }
       
double Get_tot_time(){
       return tot_time;
       }

double Get_transient_time(){
       return transient_time;
       }
       
double Get_length_plotted(){
       return length_plotted;
       }
       
double Get_depth_plotted(){
       return depth_plotted;
       }
       
double Get_width_plotted(){
       return width_plotted;
       }
       
int Get_NSELECT(){
    return NSELECT;
}

int Get_NTimeSteps(){
    return NTimeSteps;
}

double Get_artificial_factor_cell_volume(){
       return artificial_factor_cell_volume;
}
