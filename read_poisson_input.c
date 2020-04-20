#include <stdio.h>
#include <string.h>
#include "petscksp.h"
#include "nanowire.h"

void read_poisson_input(){
  char cvar1, cvar2, cvar3 ;
    Param *P ;
    PoiParam *CP ;
    P = PGetParam() ;
    CP = PGetPoiParam() ;

    FILE *f;
    f = fopen("Input.d","r");
    if(f == NULL){printf("Can't open file Input.d\n"); return;}
    //printf("\n\n Reading #POISSON_INPUT");
    
while (!feof(f)){
    char str[180]; 
    fscanf(f,"%s",str);
    if (strcmp(str,"#POISSON_INPUT")==0){
	   fscanf(f,"%*s %*s %d",&(CP->MaxIter)) ;
       //printf("\n Max Interation   =%d",CP->MaxIter);
  
       fscanf(f,"%*s %*s %lf",&(CP->ConvEps)) ;
       //printf("\n Convergence_Eps  =%le",CP->ConvEps);
  
       fscanf(f,"%*s %*s %s",CP->ksptype) ;
       //printf("\n KSPType          =%s",CP->ksptype);
  
       fscanf(f,"%*s %*s %s",CP->pctype) ;
       //printf("\n PCType           =%s",CP->pctype);
  
       fscanf(f,"%*s %*s %le",&(CP->ksprtol)) ; 
       //printf("\n KSP_Rtol         =%le",CP->ksprtol);
  
       fscanf(f,"%*s %*s %d",&(CP->gmres_restart)) ;
       //printf("\n GMRES_Restart    =%d",CP->gmres_restart);

       fscanf(f,"%*s %*s %c", &cvar1);
       SetPoissonLinearized(ConvertBooleanCharToInt(cvar1)) ;
       //printf("\n Pois_Linea(n/y)  =%c",cvar1);

       fscanf(f,"%*s %*s %d",&(P->trans_size)) ;
       //printf("\n Trans_Data_Siz   =%d",P->trans_size);
 
       fscanf(f,"%*s %*s %c", &cvar2);
       SetPRINT_MidPotFlag(ConvertBooleanCharToInt(cvar2)) ;
       //printf("\n Print_Midline_Potential(y/n) =%c",cvar2);
    
       fscanf(f,"%*s %*s %c", &cvar3);
       SetPRINT_Pot3dFlag(ConvertBooleanCharToInt(cvar3)) ;
       //printf("\n Print_3D_Potential_at_Each_Bias_Point(y/n) =%c \n",cvar3);

   }// End of if(strcmp(s,"#POISSON_INPUT")==0)
 } // End of while
 fclose(f) ;
 
 //Phan printf ra Monitor o node co rank=0
 int myrank;
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

 if(myrank ==0){
    printf("\n\n Reading #POISSON_INPUT");
    printf("\n Max Interation   =%d",CP->MaxIter);
    printf("\n Convergence_Eps  =%le",CP->ConvEps);
    printf("\n KSPType          =%s",CP->ksptype);
    printf("\n PCType           =%s",CP->pctype);
    printf("\n KSP_Rtol         =%le",CP->ksprtol);
    printf("\n GMRES_Restart    =%d",CP->gmres_restart);
    printf("\n Pois_Linea(n/y)  =%c",cvar1);
    printf("\n Trans_Data_Siz   =%d",P->trans_size);
    printf("\n Print_Midline_Potential(y/n) =%c",cvar2);
    printf("\n Print_3D_Potential_at_Each_Bias_Point(y/n) =%c \n",cvar3);
 }// End of if(myrank ==0)

   return;
}// End of void read_poisson_input()
