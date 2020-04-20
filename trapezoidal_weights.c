/* ****************************************************************
 Gia dinh la chia tren truc y deu la dy va chia tren truc z deu la dz
 INPUT:  + ny, nz chinh la ny_max va nz_max: intervals tren truc y va z
         + dy, dz la step tren truc y va z.  
 OUTPUT: + **trap_weights la mang 2 chieu chua gia tri trapezoidal_weights

Xem cong thuc 1, trang 91, 92 quyen vo cua ta        
Starting date: March 05, 2010
Update:        March 05, 2010
Latest update: May 07, 2010 de mapping voi 3D Poisson
****************************************************************** */
#include <stdio.h>

void trapezoidal_weights(){
     // Goi ham 
     double **Get_trap_weights();
     double Get_mesh_size_y(),Get_mesh_size_z();
     int Get_nya(),Get_nyb(); // int Get_ny_max(), Get_nz_max();
     int Get_nza(),Get_nzb();
     // Bieb local
     double **trap_weights = Get_trap_weights();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
     int nya = Get_nya(); // Hien nay ny_max=nyb-nya
     int nyb = Get_nyb(); //int ny_max = Get_ny_max();
     int nza = Get_nza(); // nz_max=nzb-nza
     int nzb = Get_nzb(); //int nz_max  = Get_nz_max();
     
     // Thuc hien
     int ny_max=nyb-nya;
     int nz_max=nzb-nza;
     int i,j;
     double const_area = (mesh_size_y/1.0e-9)*(mesh_size_z/1.0e-9)/4.0; // Xem hang so phan dau cua cong thuc 2 trang 92
     //do thuc te cong thuc (2) trang 92 la khong co don vi. Nen ta can chia cho 1.0e-9 <--Khong CO don vi
     //printf("\n const_area = %f",const_area);
     double weight;
     for(i=0; i<=ny_max; i++){// Do chay tu 0 nen da xet DUNG so DIEM CHIA (May 30, 2010)
         // printf("\n");  // for checking only   
         for(j=0; j<=nz_max; j++){
            // Xet 4 goc va 4 canh
            if (i==0){ // Canh thu nhat I
                if((j==0)||(j==nz_max)) { weight = 1.0;}
                else                    { weight = 2.0;}
              }
              
            else if(i==ny_max){// Canh thu hai II
                 if((j==0)||(j==nz_max)) { weight = 1.0;}
                 else                    { weight = 2.0;}
              } 
            else if((j==0)&&((i!=0)||(i!=ny_max)))     { weight = 2.0;} // canh thu ba III   
            else if((j==nz_max)&&((i!=0)||(i!=ny_max))){ weight = 2.0;} // canh thu ba III
            else { // Interior
                 weight = 4.0; 
                 }
            trap_weights[i][j] = weight; // Luu vao mang trap_weights, la mang weight cong thuc 1 trang 92
            trap_weights[i][j] = const_area*trap_weights[i][j]; // dua hang so o cong thuc 2 vao trong
            //printf("\n trap_weights[%d][%d] = %le", i,j,trap_weights[i][j]);
           // printf(" %le",trap_weights[i][j]);
       }// End of for(j=0; j<=nz; j++)
     } // End of for(i=0; i<=ny; i++){           
}// End of void trapezoidal_weights 


