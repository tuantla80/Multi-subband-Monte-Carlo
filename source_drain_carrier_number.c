/* *****************************************************************************
+ To calculate number of carriers in the source and drain regions based on doping 
  (cua cac cell cuoi 2 dau mut ay).
+ Do la ohmic contact tai cac dau mut cua S va D nen ta co the xet no la charge neutrality
  cai luon luon in thermal equilibrium tham chi la khi current is flowing -> Keep number of particles
  la constant tai cac dau mut do.
    
Starting date: Feb 17, 2010
Update:        Feb 18, 2010
Update:        March 29, 2010. TAI SAO lai thay doi cach luu tru tu ** thanh dang hang so
       + Vi thuc te la VI TRI cua particle chi phu thuoc vao x-direction. Nen cho du
         co UPDATE vi tri cua hat theo phuong y va z thi NO CUNG KHONG THAY DOI trong suot
         qua trinh simulate
Update:        May 03, 2010 (Mo rong ra vung Oxide de mapping voi 3D Poisson)
               Nhung so hat chi tinh o vung contact LOI Silicon
Latest update: May 05, 2010 su dung phuong phap artificial_factor_of_cell_volume
****************************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "nrutil.h"
static int nsource_side_carriers,ndrain_side_carriers; 

void source_drain_carrier_number(){
    // Goi cac ham 
    double Get_artificial_factor_cell_volume();
    double Get_mesh_size_x(), Get_mesh_size_y(), Get_mesh_size_z(); 
    int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
    double *Get_doping_density();
    // Cac bien local
    double artificial_factor_cell_volume = Get_artificial_factor_cell_volume();
    double mesh_size_x = Get_mesh_size_x();
    double mesh_size_y = Get_mesh_size_y();
    double mesh_size_z = Get_mesh_size_z();
    int nya = Get_nya(); int nyb = Get_nyb();
    int nza = Get_nza(); int nzb = Get_nzb();
    double *doping = Get_doping_density();//1,2 or 3 ->S,D or Channel
    // Khoi tao 2 bien ngoai
    nsource_side_carriers = 0;
    ndrain_side_carriers  = 0; 
    
    int n_number = 0;// So hat tai moi adjacent cell thu i day
    double cell_volume = 0.5*mesh_size_x*mesh_size_y*mesh_size_z;
        /* Fig 4.10, "Numerical Simulation of Submicron Semiconductor Devices" by Tomizawa */
      // Tinh cho 1/2 cell ngay canh Source va Drain
    artificial_factor_cell_volume = artificial_factor_cell_volume/mesh_size_z;
    cell_volume = cell_volume*artificial_factor_cell_volume;
    
    double denn = 0.0,factor = 0.0;
    int j,k;// i
    
    // Ta xet SIDE CONTACT  
    //i = nx0; // At the Source edge contact phan LOI Silicon
    for(j = nya; j<=nyb; j++){
       for(k=nza; k<=nzb; k++){
         denn = doping[1]; //real value of doping
         if((j==nya)||(j==nyb)||(k==nza)||(k==nzb)) 
           { 
            denn = 0.5*denn; //Left Topmost or Left Bottomost: 1/4 cell only
           }
         factor = denn*cell_volume;// printf("\n So hat thuc te = %le",factor); //Number of particles
         n_number = (int)(factor+0.5);// NOTE no luon rat nho so voi 1 NEU KHONG dung artificial_factor_cell_volumn
         //if(n_number<1){n_number = 1;}// printf("\n n_number = %d",n_number);
         nsource_side_carriers += n_number;// Used for check_source_drain_contacts() function
       } // End of for(k=nza; k<=nzb; k++){         
     } // End of for(j = nya; j<=nyb; j++){
         
    //i = nx1; // At the Drain edge phan LOI Silicon
    for(j = nya; j<=nyb; j++){
       for(k=nza; k<=nzb; k++){
           denn = doping[2];
           if((j==nya)||(j==nyb)||(k==nza)||(k==nzb)) 
              { 
               denn = 0.5*denn; //Right Topmost or Right Bottomost: 1/4 cell only
              }
          factor = denn*cell_volume; //printf("\n So hat thuc te = %le",factor);  //Number of particles
          n_number = (int)(factor+0.5);// NOTE no luon rat nho so voi 1
          //if(n_number<1){n_number = 1;}// printf("\n n_number = %d",n_number);
          ndrain_side_carriers += n_number;// Used for check_source_drain_contacts() function 
       } // End of for(k=nza; k<=nzb; k++){                   
     } // End of for(j = nya; j<=nyb; j++){
          
    //printf("\n nsource_side_carriers = %d",nsource_side_carriers);
    //printf("\n ndrain_side_carriers  = %d",ndrain_side_carriers); 
    if((nsource_side_carriers==0)||(ndrain_side_carriers==0)){
        printf("\n Report error at source_drain_carrier_number()");
        system("pause");
     }  
return;
}
/****** End of void source_drain_carrier_number() ****/

int Get_nsource_side_carriers(){
    return (nsource_side_carriers);
}

int Get_ndrain_side_carriers(){
    return (ndrain_side_carriers);
}

