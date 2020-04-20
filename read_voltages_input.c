/* ***********************************************************************
 To read voltages from an input file
 
 Latest update January 19, 2009 
************************************************************************* */
#include <stdio.h>
#include <string.h>
#include "mpi.h"

static double V_source_input,Vd_start,Vd_end,Vd_step,Vg_start,Vg_end,Vg_step;
static double DrainVoltage, GateVoltage; //static double Vd,Vg;Dung ten bien kieu GS

void read_voltages_input(){
          
      /* ******** Read voltages from an input file ********* */
      FILE *f;
      char s[180];
      f=fopen("Input.d","r"); //rewind(f); // Dua co tro ve dau tep
      if(f==NULL){printf("Error:can't open file Input.d \n");return;}
      //printf("\n\n Reading VOLTAGES_INPUT ");
      
   while (!feof(f)) {
      char str[180]; 
      fscanf(f,"%s",str); //fgets(str, 180, f); 
      if(strcmp(str,"#VOLTAGES_INPUT")==0){
                                           
      fscanf(f,"%*s %*s %lf",&V_source_input); 
      //printf("\n V_source_input[V]= %5.2lf",V_source_input);
      
      fscanf(f,"%*s %*s %lf %*s %lf %*s %lf ",&Vd_start,&Vd_end,&Vd_step); 
      //printf("\n Vd_start[V]= %3.2lf Vd_end[V]= %3.2lf Vd_step= %3.2lf",Vd_start,Vd_end,Vd_step);
      
      fscanf(f,"%*s %*s %lf %*s %lf %*s %lf ",&Vg_start,&Vg_end,&Vg_step); 
      //printf("\n Vg_start[V]= %3.2lf Vg_end[V]= %3.2lf Vg_step= %3.2lf",Vg_start,Vg_end,Vg_step);
      } // End of if(strcmp(str,"#VOLTAGES_INPUT")==0){
   } // End of while (!feof(f)) 
     fclose(f);

   // Phan in ra o monitor khi rank=0
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if(myrank ==0){ 
       printf("\n\n Reading VOLTAGES_INPUT ");
       printf("\n V_source_input[V]= %5.2lf",V_source_input);
       printf("\n Vd_start[V]= %3.2lf Vd_end[V]= %3.2lf Vd_step= %3.2lf",Vd_start,Vd_end,Vd_step);
       printf("\n Vg_start[V]= %3.2lf Vg_end[V]= %3.2lf Vg_step= %3.2lf",Vg_start,Vg_end,Vg_step); 
    }// End of if(myrank ==0)
} 
/* ******** END OF Read voltages from file ********* */

double Get_V_source_input(){
       return V_source_input;
}
       
double Get_Vd_start(){
       return Vd_start;
}
       
double Get_Vd_end(){
       return Vd_end;
}
       
double Get_Vd_step(){
       return Vd_step;
}
       
double Get_Vg_start(){
       return Vg_start;
}

double Get_Vg_end(){
       return Vg_end;
}
       
double Get_Vg_step(){
       return Vg_step;
}

// For Vd
void SetDrainVoltage(double v) {
  DrainVoltage = v ;
}

double GetDrainVoltage() {
  return(DrainVoltage);
}

// For Vg
void SetGateVoltage(double v) {
  GateVoltage = v;
}

double GetGateVoltage() {
  return(GateVoltage);
} 
















