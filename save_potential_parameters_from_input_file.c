/* ****************************************************************     
 To save parameters:  Electrostatic potential Sau khi chay cuoi cung cua time 
     flag = 1: xyz 
     flag = 0: along x_direction (yz) or y_direction (xz) or z_direction (xy)
     flag = 2: co dinh y va z VE THEO PHUONG x only

Starting date: May 24, 2010
Latest Update: May 25, 2010     
****************************************************************** */
#include <stdio.h>
#include <math.h>
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "nrutil.h"    
#include "constants.h"

void save_potential_parameters_from_input_file(int flag){

    double ***Phi = GetPot();

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
    double delta_Ec = Eg/2.0;

    double x_pos = 0.0, y_pos = 0.0, z_pos = 0.0;// Can khoi tao se tranh nham bien toan cuc
    int i,j,k;
/* ****************** flag == 1 ve cac tham so tren xyz  ************ */
 if(flag==1){ // Ve cac tham so tren theo ca phuong x,y, va z
    FILE *f;
    f=fopen("potential_xyz.dat","w");//"w" nen neu chay cho nhieu (Vd,Vg) point thi no la cho CAP CUOI CUNG
    fprintf(f,"\n #x[nm]  y[nm] z[nm]  pot_xyz[eV] \n");

    x_pos = 0.0, y_pos = 0.0, z_pos=0.0;
    for(i=0;i<=nx_max;i++){
       x_pos = (double)(i)*mesh_size_x/1.0e-9;// [m] -> [nm]
       for(j=0;j<=ny_max;j++){
          y_pos = (double)(j)*mesh_size_y/1.0e-9;// [m] -> [nm]
	  for(k=0; k<=nz_max; k++){
	    z_pos = (double)(k)*mesh_size_z/1.0e-9;
            fprintf(f,"%le %le %le %le \n",x_pos,y_pos,z_pos, (delta_Ec - Phi[i][j][k]));// potential[eV]
                            
         } //End of for(k=0; k<=nz_max; k++)
       } // End of for(j=0;j<=ny_max;j++){
    }// End of for(i=0;i<=nx_max;i++){
  
  fclose(f);                            
}// End of  if(flag ==1){ // Ve cac tham so tren theo ca phuong x,y,z
/* **************** End of flag ==1 ve cac tham so tren xy plane ************ */


/* ***************** flag==0 ve cac tham so theo mat phang xy, yz and xz************ */
if(flag == 0){ 
  double Get_length_plotted(), Get_depth_plotted(), Get_width_plotted();
  double length =  Get_length_plotted();// doc vao o read_simulation_list.c da chuyen sang dang [m]
  double depth =  Get_depth_plotted();
  double width = Get_width_plotted();
   // **************************Cat theo phuong x ****************************   
    FILE *f1;
    f1 = fopen("potential_yz_plane.dat","w");// "w" ghi de
    if(f1==NULL){ printf("\n Cannot open file potential_yz_plane.dat");return; }
    fprintf(f1,"\n #y[nm] z[nm] pot_yz[eV] \n");
    
    i = (int)(length/mesh_size_x + 0.5);
    y_pos = 0.0, z_pos = 0.0;
    for(j=0;j<=ny_max;j++){
        y_pos = (double)(j)*mesh_size_y/1.0e-9;// [m] -> [nm]
	for(k=0; k<=nz_max; k++){
	    z_pos = (double)(k)*mesh_size_z/1.0e-9;
            fprintf(f1," %le %le %le \n",y_pos,z_pos,(delta_Ec - Phi[i][j][k]));// potential[eV]
        } //End of for(k=0; k<=nz_max; k++)
     } // End of for(j=0;j<=ny_max;j++){

    fclose(f1);
   // ********************** End of Cat theo phuong x  ****************************       
 
  
   /* ************************** Cat theo phuong y **************************** */ 
  
    FILE *f2;
    f2 = fopen("potential_xz_plane.dat","w");// Kieu mo "w" ghi de
    if(f2==NULL){ printf("\n Cannot open file potential_xz_plane.dat"); return; }
    fprintf(f2,"\n #x[nm] z[nm] pot_xz[eV] \n");
   
    j = (int)(depth/mesh_size_y + 0.5); 
    x_pos = 0.0, z_pos = 0.0;
    for(i=0;i<=nx_max;i++){
       x_pos = (double)(i)*mesh_size_x/1.0e-9;// [m] -> [nm]
       for(k=0; k<=nz_max; k++){
	  z_pos = (double)(k)*mesh_size_z/1.0e-9;
          fprintf(f2," %le %le %le \n",x_pos,z_pos,(delta_Ec - Phi[i][j][k]));// potential[eV]
      } //End of for(k=0; k<=nz_max; k++)
    }// End of for(i=0;i<=nx_max;i++){
  
    fclose(f2);
   // ********************** End of Cat theo phuong y  **************************** 

    // **************************Cat theo phuong z ****************************   
    FILE *f3;
    f3 = fopen("potential_xy_plane.dat","a");
    if(f3==NULL){ printf("\n Cannot open file potential_xy_plane.dat");return; }
    fprintf(f3,"\n #x[nm] y[nm] pot_xy[eV] \n");
    
    k = (int)(width/mesh_size_z + 0.5);
    x_pos = 0.0, y_pos = 0.0;

    for(i=0;i<=nx_max;i++){
       x_pos = (double)(i)*mesh_size_x/1.0e-9;// [m] -> [nm]
       for(j=0;j<=ny_max;j++){
          y_pos = (double)(j)*mesh_size_y/1.0e-9;// [m] -> [nm]
	  fprintf(f3," %le %le %le \n",x_pos,y_pos, (delta_Ec - Phi[i][j][k]));// potential[eV]
       } // End of for(j=0;j<=ny_max;j++){
    }// End of for(i=0;i<=nx_max;i++){
   fclose(f3);
   // ********************** End of Cat theo phuong z  ****************************  

} // End of if(flag ==0)
/* *************** End of  flag==0 ve cac tham so theo mat phang xy, yz and xz  *********** */

/* ***************** flag==2 VE THEO x thoi con y va z duoc co dinh va doc tu Input file************ */
if(flag == 2){ 
  double Get_depth_plotted(), Get_width_plotted();
  double depth =  Get_depth_plotted();// doc vao o read_simulation_list.c da chuyen sang dang [m]
  double width = Get_width_plotted();
   // **************************Ve theo phuong x only ****************************   
    FILE *f4;
    f4 = fopen("potential_x_direction.dat","a");
    if(f4==NULL){ printf("\n Cannot open file potential_x_direction.dat");return; }
    fprintf(f4,"\n #x[nm]  pot_x[eV] \n");

    j = (int)(depth/mesh_size_y + 0.5); 
    k = (int)(width/mesh_size_z + 0.5);
    x_pos = 0.0;

    for(i=0;i<=nx_max;i++){
       x_pos = (double)(i)*mesh_size_x/1.0e-9;// [m] -> [nm]
       fprintf(f4,"%le %le \n",x_pos, (delta_Ec - Phi[i][j][k]));// potential[eV]
    }// End of for(i=0;i<=nx_max;i++){
    
    fclose(f4);
 }
   // ********************** End of Ve theo phuong x only   ****************************      
 
return;
}


