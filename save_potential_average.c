/* ******************************************************************************
  To write potential average after transient time 
  NOTE: La tinh trung binh tai 1 gia tri (depth va width) hoac lenght doc tu Input file
  
   INPUT: time: thoi gian ma simulation dang chay (tu dt den total_time) 
          iter_reference: if ==0: tinh trung binh cac gia tri (luc dt da chay den diem cuoi cung cua total_time)
                          else: tinh tong cac gia tri.  
        
   OUTPUT: To write average related parameters
   NOTE: Tinh toan va chay cho 1 node thoi   
Starting date: May 24, 2010
Latest Update: June 15, 2010          
        
****************************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "nrutil.h"
#include "constants.h"

void save_potential_average(FILE *f1, FILE *f2, double time, int iter_reference) {//f1 ve cho x, f2 ve cho yz

  double ***Phi = GetPot();
 
  double *Get_pot_x_avg();
  double *pot_x_avg = Get_pot_x_avg();

  double **Get_pot_yz_avg();
  double **pot_yz_avg = Get_pot_yz_avg();

  int Get_nx0(),Get_nx1(), Get_ny0(),Get_ny1(), Get_nz0(),Get_nz1();
  int nx0 = Get_nx0(); int nx1 = Get_nx1(); 
  int ny0 = Get_ny0(); int ny1 = Get_ny1();
  int nz0 = Get_nz0(); int nz1 = Get_nz1();
  int nx_max = nx1-nx0;
  int ny_max = ny1-ny0;
  int nz_max = nz1-nz0;

  double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
  double mesh_size_x = Get_mesh_size_x();
  double mesh_size_y = Get_mesh_size_y();
  double mesh_size_z = Get_mesh_size_z();

  double Get_Eg(); // Band gap
  double Eg = Get_Eg();
  double delta_Ec = Eg/2.0;//  double delta_Ec = 0.0;// NOTE ta ve dang -Phi(i,j,k) thoi

  double Get_length_plotted(), Get_depth_plotted(), Get_width_plotted();
  double length = Get_length_plotted();// doc vao o read_simulation_list.c da chuyen sang dang [m]
  double depth = Get_depth_plotted();
  double width = Get_width_plotted();
  
  double Get_transient_time();
  double transient_time = Get_transient_time(); // Chon bang BAO NHIEU can KINH NGHIEM

  int i,j,k;
  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

if(myrank ==0){// NOTE chi tinh toan tren 1 node thoi

  /* ***************************************************************************************
              time > transient_time: Tinh tong cac gia tri 
  **************************************************************************************** */
  if(time > transient_time){// Neu thoi gian time lon hon gia tri transient time duoc doc tu input file thi ta se tinh tong 
   
    // Tinh trung binh theo phuong x voi gia tri co dinh y va z doc tu input file 
    j = (int)(depth/mesh_size_y + 0.5); 
    k = (int)(width/mesh_size_z + 0.5);
                             
    for (i=0; i<=nx_max; i++){
      pot_x_avg[i] += (delta_Ec - Phi[i][j][k]);
    } 
           
    // Tinh trung binh theo yz khi x la co dinh
    i = (int)(length/mesh_size_x + 0.5);
    for(j=0; j<=ny_max; j++){
      for(k=0; k<=nz_max; k++){
	pot_yz_avg[j][k] += (delta_Ec - Phi[i][j][k]);
      }
    }
      
  } // End of if(time > transient_time) nghia la het phan cong cac gia tri
  /* ***************************************************************************************
           End of  time > transient_time: Tinh tong cac gia tri 
  **************************************************************************************** */

  /* ***************************************************************************************
      iter_reference ==0: Calculating time averaging excluding initial transient
  **************************************************************************************** */
  // Den phan tinh trung binh cong cac gia tri theo phuong x HOAC mat yz
  double x_pos = 0.0, y_pos = 0.0, z_pos = 0.0;

  if(iter_reference == 0){  //iter_reference==0: tinh trung binh cac gia tri (luc time da chay den diem cuoi cung cua total_time) 
   double Get_dt(),Get_tot_time();
   double dt = Get_dt(); // Observation time step
   double tot_time = Get_tot_time(); // Tong thoi gian chay Simulation
   int n_time_steps_after_transient = (int)((tot_time-transient_time)/dt);//so khoang step cho dt tinh tu sau transient_time

   //FILE *f;  //f = fopen("average_potential_along_x.dat","w");// Do la "w" nen neu chay nhieu cap (Vd,Vg) thi no la CAP CUOI CUNG               
   //if(f == NULL){printf("\n Cannot open file average_potential_along_x.dat");return; }
   //fprintf(f,"\n #x[nm] Pot_x[eV] \n");

    x_pos = 0.0;
    if(n_time_steps_after_transient > 0){
      for(i=0;i<=nx_max;i++){
	x_pos = (double)(i)*mesh_size_x/1.0e-9;// [m] -> [nm]
	pot_x_avg[i] = pot_x_avg[i]/((double)(n_time_steps_after_transient));// Vay ghi vao file la da ghi gia tri trung binh roi nhe !         
	fprintf(f1, "%le %le \n", x_pos,pot_x_avg[i]);
      }
    } // End of if(n_time_steps_after_transient > 0)
  
   //fclose(f);
   
   //FILE *f_avg_yz;   //f_avg_yz = fopen("average_potential_yz_plane.dat","w");                  
   //if(f_avg_yz == NULL){printf("\n Cannot open file average_potential_yz_plane.dat");return; }
   //fprintf(f_avg_yz,"\n #y[nm] z[nm] Pot_avg_yz[eV] \n");  
   
    y_pos = 0.0, z_pos = 0.0;
    if(n_time_steps_after_transient > 0){
      for(j=0;j<=ny_max;j++){
        y_pos = (double)(j)*mesh_size_y/1.0e-9;// [m] -> [nm]
	for(k=0; k<=nz_max; k++){
	  z_pos = (double)(k)*mesh_size_z/1.0e-9;
	    
	  pot_yz_avg[j][k] =  pot_yz_avg[j][k]/((double)(n_time_steps_after_transient));
	  fprintf(f2," %le %le %le \n",y_pos,z_pos, pot_yz_avg[j][k]);// potential[eV]

	} //End of for(k=0; k<=nz_max; k++)
      } // End of for(j=0;j<=ny_max;j++){
    }// End of  if(n_time_steps_after_transient > 0)
 
    //fclose(f_avg_yz); 

    //NOTE: reset; Thuc te thi phan nay chi tinh 1 lan duy nhat khi time da chay den total time Nhung do no chay cho cac gia tri Vd va Vg khac nhau nen can reset
    for(i=nx0; i<=nx1; i++){ 
       pot_x_avg[i] = 0.0; 
    }

    for(i=ny0; i<=ny1; i++){ 
        for(j=nz0; j<=nz1; j++){ pot_yz_avg[i][j] = 0.0;
	}
     }

   } // End of if(iter_reference ==0){
    /* ***************************************************************************************
       End of iter_reference ==0: Calculating time averaging excluding initial transient
    **************************************************************************************** */ 
}// End of if(myrank ==0){

return;
}// End of void save_potential_average(double time, int iter_reference) 
