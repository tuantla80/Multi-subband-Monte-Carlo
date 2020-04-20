/* *****************************************************************************
   First order correction for eigen energy using time-independent perturbation
   Reference: Applied Quantum Mechanics, Levi A.F.J., 2003. Chapter 10.
   INPUT: StepsUsedCorrection: buoc j_th time steps ta su dung eigen energy correction
   OUTPUT:
          eigen energy da duoc correction qua khoang step thu j-th
   
Starting date: April 19, 2010
Update:        April 20, 2010.
Latest update: May 17, 2010
******************************************************************************** */
#include <stdio.h>
#include "nrutil.h"
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"


void eigen_energy_correction(int StepsUsedCorrection){
   
    void save_eig_correction(int StepsUsedCorrection);
   
    double *****Get_wave();
    double *****wave = Get_wave(); // de lay wave

    double ***Get_eig(),***Get_eig_jthSchro();
    double ***eig = Get_eig(); //eigen energy se thay doi sau nhung khoang thoi gian dt
    double ***eig_jthSchro = Get_eig_jthSchro();// giu gia tri eig tai jth Schrodinger solution
    
    double ***Phi = GetPot();

    double ***Get_fai_jthSchro(),*Get_pot();//chuyen fai(s-th,j,k) thanh pot 1D (potential 1D)
    double ***fai_jthSchro = Get_fai_jthSchro();// giu gia tri cua potential tai jth Schrodinger solution
    double *pot   = Get_pot();
  
    int Get_nx0(),Get_nx1(),Get_ny0(),Get_nya(),Get_nyb(),Get_ny1(),Get_nz0(),Get_nza(),Get_nzb(),Get_nz1(); 
    int nx0 = Get_nx0(); int nx1 = Get_nx1();
    int ny0 = Get_ny0(); int nya = Get_nya(); int nyb = Get_nyb(); int ny1 = Get_ny1(); 
    int nz0 = Get_nz0(); int nza = Get_nza(); int nzb = Get_nzb(); int nz1 = Get_nz1(); 
    int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    int ny_max = nyb - nya;// De thay doi it nhat co the
    int nz_max = nzb - nza;// DO sau do ta chay tu 0 den ny_max va nz_max nen do DIEM CHIA ls ok (May 30, 2010)
    
    int Get_NSELECT(),Get_NTimeSteps();
    int NSELECT = Get_NSELECT();
    int NTimeSteps = Get_NTimeSteps();
  
// Thuc hien
   int s,v,i,j,k;
   //double temp = 0.0;
   //Buoc 1.
   // Lan dau tien goi ham eigen_energy_correction se phai giu gia tri eigen energy 
   // va potential vi no chinh la gia tri tu j-th Schrodinger solution  
   if(StepsUsedCorrection %  NTimeSteps ==1){
       // Buoc 1.1. Gan eig tai jth Schrodinger solution cho eig_jthSchro                   
       for(s=0; s<=nx_max; s++){
          for(v=1; v<=3; v++){
              for(i=1; i<=NSELECT; i++){ 
                  eig_jthSchro[s][v][i] = eig[s][v][i];
              }
          }
       }// End of for(s=0; s<=nx_max; s++)
       
      // Buoc 1.2. Gan potential tai jth Schrodinger solution cho eig_jthSchro                    
      for(s=0; s<=nx_max; s++){
	for(j=ny0; j<=ny1; j++){
	  for(k=nz0; k<=nz1; k++){
             fai_jthSchro[s][j][k] = Phi[s][j][k];
           }
         }
      }// End of for(i=0; i<=nx_max; i++)
      
      // Buoc 1.3. Khong can gan wave boi vi minh khong update wave neu khong giai Schrodinger
   } // End of if(StepsUsedCorrection %  NTimeSteps ==1)                
    
    
    // Buoc 2. Tinh eigen energy duoc correction
    // Buoc 2.1. Tinh delta_fai
    double ***delta_fai=d3matrix(0,nx_max,ny0,ny1,nz0,nz1);
    for(s=0; s<=nx_max; s++){
	for(j=ny0; j<=ny1; j++){
	  for(k=nz0; k<=nz1; k++){
              // CHECK LAI XEM LA DAU - hay +. Check theo bai bao cua ho thi del_fai=Vcu(eV)-Vmoi(eV)
               delta_fai[s][j][k]= fai_jthSchro[s][j][k]-Phi[s][j][k];
           }
         }
      }// End of for(i=0; i<=nx_max; i++)
    // Buoc 2.2 Tinh eigen energy
    double delta_ener; // CHU Y: index cua wave la o trong vung ma ta giai 2D Schrodinger
    for(s=0; s<=nx_max; s++){// cung chinh la i o tren
          for(v=1; v<=3; v++){
              for(i=1; i<=NSELECT; i++){
                       
                 delta_ener = 0.0;
		 for(j=0; j<=ny_max; j++){// chay tu 0 da the hien so diem chia can +1 (May 30, 2010)
                    for(k=0; k<=nz_max; k++){ 
                        delta_ener +=  wave[s][v][i][j][k]*delta_fai[s][j][k]*wave[s][v][i][j][k];
                      }
                  }// end of for(j=0; j<=ny_max; j++)
                                    
                  eig[s][v][i] = eig_jthSchro[s][v][i] + delta_ener;
                 
              }// End of for(i=1; i<=NSELECT; i++)
          }//End of for(v=1; v<=3; v++)
       }// End of for(s=0; s<=nx_max; s++)
 // Save ca xem sao
    save_eig_correction(StepsUsedCorrection);
 // Free local array
 free_d3matrix(delta_fai,0,nx_max,ny0,ny1,nz0,nz1);    
}// End of void eigen_energy_correction()
//******************************************************************************
/* *****************************************************************************
Starting date: April 21, 2010
Latest update: April 21, 2010.
******************************************************************************** */
void save_eig_correction(int StepsUsedCorrection){
     // Goi ham
      double ***Get_eig();
      int Get_NSELECT();
      // bien local
      double ***eig = Get_eig();
      int NSELECT = Get_NSELECT();

      int Get_nx0(),Get_nx1();
      int nx0 = Get_nx0(); int nx1 = Get_nx1();
      int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
      
// Thuc hien
      int s=(int)(nx_max/2+0.5);// lay section o giua  
      int i,v;
      FILE *f; f = fopen("eig_correction.dat","a"); // "w" neu tep ton tai no se bi xoa
      if(f == NULL) { printf("\n Cannot open file eig_correction.dat"); return ; }
      fprintf(f,"\n # s_th v_th i_subband time eigen_correction[eV] \n");
      for(v=1; v<=3; v++){
         for(i=1; i<=NSELECT; i++){ // do m o ham dsyevx tu diasym la bang NSELECT
	     fprintf(f,"%d      %d       %d     %d        %le \n",s,v,i,StepsUsedCorrection,eig[s][v][i]);
          }
     }                                      
      fclose(f);
      return;
}// End of void void save_eig_correction()
