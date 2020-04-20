/* ***********************************************************************
               To count the number of used electrons
               
Starting date: Feb 22, 2010
Latest update: Feb 22, 2010 
************************************************************************* */
#include <stdio.h>
#include "constants.h"

int count_used_particles(){ // ne here is the number of electrons using in the simulator
    // Goi ham
    int *Get_valley(); // Ta su dung bien valley[] thay cho bien ip[] o Semiclassical MC
    
    // Bien local
    int *valley = Get_valley();
    
    //Thuc hien
    int i, ne = 0;
    for(i=1; i<=max_electron_number; i++){
        if(valley[i] != 9) { // KHAC 9 co nghia la hat dang duoc SU DUNG
              ne = ne + 1; 
           }
       }
    //printf("\n Number of electrons using = %d ",ne);
   
    return ne;
 }

