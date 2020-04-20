/* ***********************************************************************
To delete extra particles based on their valley[]
 Neu valley[tai hat i-th] = 9 <--Can delete HAT nay

Starting date: March 28, 2010
Latest update: March 29, 2010 
************************************************************************* */
#include <stdbool.h> // De khai bao kieu bool
#include <stdio.h>

void delete_particles(){ 
     // Goi ham
     double **Get_p();
     double *Get_energy();
     int *Get_valley();
     int *Get_subband();
     int Get_n_used(); // Tinh SO HAT dang LUU valley=1,2,3 HAY 9
     void Set_n_used( int n);
     int count_used_particles();
     
     // Ca bien local
     double **p     = Get_p();
     double *energy = Get_energy();
     int *valley = Get_valley();
     int *subband = Get_subband();
     int n_used = Get_n_used(); // TAT ca cac HAT valley=1,2,3 HAY 9
     
// Thuc hien    
     int n_fix = 0;  
     bool flag_conv;  
     //int ne = count_used_particles();
     //printf("\n Number of electrons from count_used_particles() = %d ",ne); // phai bang so hat sau delete
     //printf("\n Number of electrons Before delete = %d ",n_used);  
     int i,j;   
     for(i = 1; i<=n_used; i++){
        if(valley[i] == 9){ // Particle i-th has iv=9 need to delete
           flag_conv = false; // =0 means FALSE
           n_fix = n_used + 1; // Cong them 1 hat vao n_used
            
           while(!flag_conv){
                if(valley[n_fix] != 9){ // Neu hat cong them vao do co iv khac 9
                                        //thi ta tien hanh chuyen tat ca cac parameters cua hat do
                                        //cho hat thu i va giam gia tri di 1 
                    valley[i] = valley[n_fix];  // move valley index
                    energy[i] = energy[n_fix];  // move energy
                    subband[i] = subband[n_fix];// move subband
                    for(j = 1; j <= 4; j++)
                       {                     // j=1: move mometum; j=2: move position
                         p[i][j]=p[n_fix][j];// j=3: move free flight time hien tai;
                       }                     // j=4: move i-region> Hien tai (March 29) CHUA DUNG 
          
                    valley[n_fix] = 9; // inactive particle n_fix -> will be delete it
                    n_fix = n_fix - 1; //
                    flag_conv = true; // =1 means TRUE
                 } // End of if(valley[n_fix] != 9)
               
                // Below, means ip[n_fix]==9
                n_fix = n_fix - 1; // reduce n_fix and check ip index again 
                if(n_fix < i) { flag_conv = true; } //  =1 means TRUE              
           } // End of while(!flag_conv){
        }// End of if(valley[i] == 9)
     } //End of for (int i=1; i<=n_used; i++){
  
    // Gio day no da NOI LEN ca roi thi ta moi GIAM   
    while (valley[n_used]==9){
           n_used = n_used - 1;
          }
    Set_n_used(n_used);
    //printf("\n Number of particles used (after deleting) = %d, Kieu 2 =%d",n_used,Get_n_used());
  return;
}

