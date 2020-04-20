/* ***********************************************************************
 To initialize free flight time for each electron dua vao tham so x-direction cua hat
 thi se biet no o SECTION nao -> Su dung max_gm[s] tai SECTION do
  Matrix  p[][] (parameters for particle)
          p[n][1] = kx (is stored kx)
          p[n][2] = x
          p[n][3] = is to store ts (ts or tc is free flight time)
          p[n][4] = i_region (which region we are considering?)
 
 NOTE: Phai dat ham nay SAU ham normalize_table()    
  
Starting date: March 17, 2010
Update:        March 17, 2010
Latest update: May 19, 2010
************************************************************************** */
#include <math.h>
#include <stdio.h>
#include "mpi.h"
#include "constants.h"
void init_free_flight(){
     // Goi cac ham
     double **Get_p(); // THAM SO cua HAT
     double *Get_max_gm(); // Gamma_max cho tung SECTION
     int count_used_particles(); // Tinh SO HAT dang SU DUNG
     double Get_mesh_size_x();
     long Get_idum();
     float random2(long *idum);
     
     // Cac bien local
     long idum  = Get_idum();
     double **p     = Get_p();
     double *max_gm = Get_max_gm(); // printf("\n Maximum scatering rate at 0 section is %le",max_gm[0]);
     int ne = count_used_particles();
     double mesh_size_x = Get_mesh_size_x();
     // Thuc hien
     //printf("\n Number of electrons using before init_free_flight()  = %d ",ne); getchar();
     int i,s; // s la chi so cho SECTION
     double x_position;
     
     for(i=1;i<=ne;i++){ // Calculate for each particle 
         x_position = p[i][2]; // lay POSITION cua HAT
         
         // Tim duoc vi tri SECTION cua HAT
         s = (int)(round(x_position/mesh_size_x)); //printf("\n Hat thu %d co vi tri SECTION %d",i,s); getchar();
         
         // Tinh free-flight time cho HAT o SECTION s-th
         p[i][3] = -log(random2(&idum))/max_gm[s]; // [s]
         //printf("\n Free-flight time la %le[s] cho HAT tai SECTION %d",p[i][3],s);  getchar();
         
         }

     // Phan save vao file
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
     int Get_save_init_free_flight();
     int save_init_free_flight_or_not =  Get_save_init_free_flight();
     FILE *f;
     if(myrank ==0){ 
       if(save_init_free_flight_or_not==1){ // Co save
	 f = fopen("init_free_flight.dat","w"); // "w" neu tep ton tai no se bi xoa
	 if(f==NULL){printf("\n Cannot open file init_free_flight.dat"); return ; }
	 fprintf(f,"\n #particle_i_th InSection free_flight_time \n");
	 for(i=1;i<=ne;i++){
	   x_position = p[i][2]; // lay POSITION cua HAT
	   s = (int)(round(x_position/mesh_size_x)); 
	   
	   fprintf(f,"      %d        %d     %le\n",i,s,p[i][3]);
	 }
	 fclose(f);
       }
     }
    
     return;
} // End of void init_free_flight()
