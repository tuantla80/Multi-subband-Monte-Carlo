/* ****************************************************************
 Ham de thuc hien viec tinh Form Factor cho tat ca cross-section, 3 cap valley pairs
 va tu subband n den subband m
 INPUT:  + nx: s-th section - so luong tu 0 den nx_max
         + subband index n (1 den NSELECT) den m (1 den NSELECT)
         + Wave_s,v,n_(y,z) va Wave_s,v,m_(y,z)
         + So diem chia theo phuong y (tu 0 den ny_max) va z ( tu 0 den nz_max)
 
 OUTPUT: + FormFactor_s,n,m,v  : 4D array. NOTE chi so v la 1 kieu ket hop cua Intra va Inter_Valley   

RAT CHU Y: Cai ta tinh duoc la CHUA normalize theo DON VI. De tinh theo DON VI
           ta CAN NHAN them 1.0e+18. Xem trang 92
            
Starting date: March 05, 2010
Update:        March 10, 2010 (Big change)
Latest update: March 15, 2010: Y nghia cua don vi dang 1/m2 hoac 1/cm2
****************************************************************** */
#include <stdio.h>

void form_factor_calculation(){
     // Goi ham
     double *****Get_wave(); // Lay wave function tu 2D Schrodinger
     double **Get_trap_weights();// Lay phan cong thuc 1 va hang so cua cong thuc 2 trang 92
     double ****Get_form_factor();// chua form factor
    
     // Local bien
     double *****wave   = Get_wave();
     double **trap_weights = Get_trap_weights();
     double ****form_factor = Get_form_factor();

     int Get_nx0(),Get_nx1(), Get_nya(),Get_nyb(),Get_nza(),Get_nzb(); 
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 

     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
     int ny_max = nyb - nya;// De thay doi it nhat co the
     int nz_max = nzb - nza; 

     int Get_NSELECT(); 
     int NSELECT = Get_NSELECT();

// Thuc hien
     int i,j; // chay cho so diem tren truc y, z
     int s,v,n,m;// section, valley, tu subband n den subband m
     double sum = 0.0;
   
     // Reset form_factor truoc moi lan vao
     for(s=nx0; s<=nx1; s++){
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(v=1; v<=9; v++){      form_factor[s][n][m][v] = 0.0;
                }
             }
         }
     }
/* ***************************************************************************** 
Part I. + Acoustic phonon (Intra-Valley (Inter- and Intra-subband))
        + Or g-process zero-order optical: do MAC DU no la Inter-Valley nhung
          xay ra o SAME axis nen subband energy va wave cung GIONG voi truong hop
          acoustic.            
 ******************************************************************************* */
     for(s=0; s<=nx_max; s++){ // chay cho tung section
        for(v=1; v<=3; v++){// chay cho 3 cap valley. Xem trang 96
            for(n=1; n<=NSELECT; n++){ // chay cho subband n
                for(m=1; m<=NSELECT; m++){ // chay cho subband m
		  sum =0.0;
                   for(i=0; i<=ny_max; i++){  // chay cho diem tren truc y
                       for(j=0; j<=nz_max; j++){// chay cho diem tren truc z 
                          sum = sum + 
                          trap_weights[i][j]*wave[s][v][n][i][j]*wave[s][v][n][i][j]*wave[s][v][m][i][j]*wave[s][v][m][i][j];
                       } // End of for(j=0; j<=nz_max; j++)
                   }// End of for(i=0; i<=ny_max; i++)
                  
               form_factor[s][n][m][v] = sum*1.0e+18; // Nhan them 1.0e+18 Xem trang 92
               sum = 0.0; // de quay ve vong lap ke tiep
             }// End of for(m=1; m<=NSELECT;m++)
           }// End of for(n=1; n<=NSELECT; n++)
       }// End of for(v=1; v<=3; v++)
     }// End of for(s=0; s<=nx_max; s++)

/* **************************************************************************** 
    Part II. f-process Non-polar Optical phonon 
            (Inter-Valley (Inter- and Intra-subband)) (trang 99)
 *******************************************************************************/
// II.1. Tu valley pair 1 (v=1) --> Valley pair 2
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ // chay cho subband n
             for(m=1; m<=NSELECT; m++){ // chay cho subband m
                
                 for(i=0; i<=ny_max; i++){  // chay cho diem tren truc y
                     for(j=0; j<=nz_max; j++){// chay cho diem tren truc z 
                        sum = sum + 
                        trap_weights[i][j]*wave[s][1][n][i][j]*wave[s][1][n][i][j]*wave[s][2][m][i][j]*wave[s][2][m][i][j];
                       } // End of for(j=0; j<=nz_max; j++)
                   }// End of for(i=0; i<=ny_max; i++)
                  
               form_factor[s][n][m][4] = sum*1.0e+18;// Nhan 1.0e+18 xem trang 92. De normalized lai theo m2
               sum = 0.0; // de quay ve vong lap ke tiep
               
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx_max; s++)
// End of II.1. Tu valley 1 (v=1) --> Valley 2 (v=2)

// II.2. Tu valley 1 (v=1) --> Valley 3 (v=3)
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(i=0; i<=ny_max; i++){  
                     for(j=0; j<=nz_max; j++){
                        sum = sum + 
                        trap_weights[i][j]*wave[s][1][n][i][j]*wave[s][1][n][i][j]*wave[s][3][m][i][j]*wave[s][3][m][i][j];
                       } // End of for(j=0; j<=nz_max; j++)
                   }// End of for(i=0; i<=ny_max; i++)
                  
               form_factor[s][n][m][5] = sum*1.0e+18;
               sum = 0.0; // de quay ve vong lap ke tiep
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx_max; s++)
// End of II.2. Tu valley 1 (v=1) --> Valley 3 (v=3)

// II.3. Tu valley 2 (v=2) --> Valley 1 (v=1)
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(i=0; i<=ny_max; i++){  
                     for(j=0; j<=nz_max; j++){
                        sum = sum + 
                        trap_weights[i][j]*wave[s][2][n][i][j]*wave[s][2][n][i][j]*wave[s][1][m][i][j]*wave[s][1][m][i][j];
                       } // End of for(j=0; j<=nz; j++)
                 }// End of for(i=0; i<=ny; i++)
                  
               form_factor[s][n][m][6] = sum*1.0e+18;
               sum = 0.0; // de quay ve vong lap ke tiep
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx; s++)
// End of II.3. Tu valley 2 (v=2) --> Valley 1 (v=1)

// II.4. Tu valley 2 (v=2) --> Valley 3 (v=3)
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(i=0; i<=ny_max; i++){  
                     for(j=0; j<=nz_max; j++){
                        sum = sum + 
                        trap_weights[i][j]*wave[s][2][n][i][j]*wave[s][2][n][i][j]*wave[s][3][m][i][j]*wave[s][3][m][i][j];
                       } // End of for(j=0; j<=nz; j++)
                 }// End of for(i=0; i<=ny; i++)
                  
               form_factor[s][n][m][7] = sum*1.0e+18;
               sum = 0.0; // de quay ve vong lap ke tiep
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx; s++)
// End of II.4. Tu valley 2 (v=2) --> Valley 3 (v=3)

// II.5. Tu valley 3 (v=3) --> Valley 1 (v=1)
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(i=0; i<=ny_max; i++){  
                     for(j=0; j<=nz_max; j++){
                        sum = sum + 
                        trap_weights[i][j]*wave[s][3][n][i][j]*wave[s][3][n][i][j]*wave[s][1][m][i][j]*wave[s][1][m][i][j];
                       } // End of for(j=0; j<=nz; j++)
                 }// End of for(i=0; i<=ny; i++)
                  
               form_factor[s][n][m][8] = sum*1.0e+18;
               sum = 0.0; // de quay ve vong lap ke tiep
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx; s++)
// End of II.5. Tu valley 3 (v=3) --> Valley 1 (v=1)

// II.6. Tu valley 3 (v=3) --> Valley 2 (v=2)
     sum = 0.0;
     for(s=0; s<=nx_max; s++){ 
         for(n=1; n<=NSELECT; n++){ 
             for(m=1; m<=NSELECT; m++){ 
                for(i=0; i<=ny_max; i++){  
                     for(j=0; j<=nz_max; j++){
                        sum = sum + 
                        trap_weights[i][j]*wave[s][3][n][i][j]*wave[s][3][n][i][j]*wave[s][2][m][i][j]*wave[s][2][m][i][j];
                       } // End of for(j=0; j<=nz; j++)
                 }// End of for(i=0; i<=ny; i++)
                  
               form_factor[s][n][m][9] = sum*1.0e+18;
               sum = 0.0; // de quay ve vong lap ke tiep
           }// End of for(m=1; m<=NSELECT;m++)
         }// End of for(n=1; n<=NSELECT; n++)
      }// End of for(s=0; s<=nx; s++)
return;
}//End of void form_factor_calculation
//************************************************************************

/* ****************************************************************
 De in form factor
 INPUT:  + nx: s-th section - so luong tu 0 den nx_max
         + subband index n (1 den NSELECT) den m (1 den NSELECT)
         + FormFactor_s,v,n->m
         
 OUTPUT: + In ra FormFactor_s,n,m,v  tren man hinh
        
Starting date: March 05, 2010
Latest update: March 10, 2010
****************************************************************** */
void printf_form_factor_calculation(int nx, int NSELECT){
     // Goi ham
     double ****Get_form_factor();// chua form factor
     // Local bien 
     double ****form_factor = Get_form_factor();
     // Thuc hien
     int s,v,n,m;// section, valley, tu subband n den subband m
     for(s=0; s<=nx; s++){ // chay cho tung section
        for(n=1; n<=NSELECT; n++){ // chay cho subband n
            for(m=1; m<=NSELECT; m++){ // chay cho subband m
                for(v=1; v<=9; v++){// 1 den 3: Acoustic; 4 den 9: Optical Phonon                  
                    printf("\n form_factor[%d][%d][%d][%d]=%le",s,n,m,v,form_factor[s][n][m][v]);
               }
            }
         }
     }
} // End of void printf_form_factor_calculation
//************************************************************************

/* ****************************************************************
 De SAVE form factor vao trong 1 file
 INPUT:  + nx: s-th section - so luong tu 0 den nx_max
         + subband index n (1 den NSELECT) den m (1 den NSELECT)
         + FormFactor_s,n,m,v
         
 OUTPUT: + SAVE FormFactor_s,n,m,v  Vao File Form_factor.dat
        
Starting date: March 06, 2010
Latest update: March 10, 2010
****************************************************************** */
void save_form_factor_calculation(int s_th, char *fn){
    
     double ****Get_form_factor();// chua form factor
     double ****form_factor = Get_form_factor();

     int Get_NSELECT();
     int NSELECT = Get_NSELECT();
  // Thuc hien
     // To store Form factor

     FILE *f;
     f = fopen(fn,"w"); // "w" neu tep ton tai no se bi xoa
     //f = fopen("Form_factor.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file at save_form_factor_calculation.c"); return ; }
     fprintf(f,"\n #s_th n m v form_factor[1/m2] \n");
     int v,n,m;// section, valley, tu subband n den subband m
         for(n=1; n<=NSELECT; n++){ // chay cho subband n
             for(m=1; m<=NSELECT; m++){ // chay cho subband m 
                 for(v=1; v<=9; v++){// Kieu ket hop Intra- va Inter-Valley                
                    fprintf(f,"%d    %d %d %d %le \n",s_th,n,m,v,form_factor[s_th][n][m][v]);
               }
            }
         }
     fclose(f);
return; 
} // End of void save_form_factor_calculation
