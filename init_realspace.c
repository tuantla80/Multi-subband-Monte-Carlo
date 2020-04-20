/* ***********************************************************************
 To initialise position for each particle (electron) based on node location
         1D transport direction only nen chi can tham so i la du
  
Starting date: Feb 22, 2010
Update:        Feb 22, 2010
Update:        March 28, 2010 (Giai thich TAI SAO khong update cho y va z)
       + Vi hat chi chay tren x direction ma thoi
       + Hat co scattering theo phuong x, NHUNG no chi CO THE lam thay doi momentum,
         energy, valley hay subband ma thoi
       + Cach tinh charge density cung KHAC voi Semiclassical MC
Latest update: May 06, 2010 (De phu hop voi 3D Poisson)
************************************************************************* */
#include <math.h>
#include <stdio.h>
    
void init_realspace(int ne,int i){
     // Goi cac ham 
     double **Get_p();
     int Get_nx0(),Get_nx1();
     long Get_idum();
     float random2(long *idum);
     double Get_mesh_size_x();
     
     // Cac bien local
     double **p = Get_p();
     int nx0 = Get_nx0();
     int nx1 = Get_nx1(); 
     long idum  = Get_idum();
     double mesh_size_x = Get_mesh_size_x();
     
     // Thuc hien
     double device_length = ((double)(nx1-nx0))*mesh_size_x;  
     double x_position=0.0,factor = 0.0;// Define particle coordinates based on grid location
     if((i!=nx0) && (i!= nx1)){
        x_position = ((double)(i)+random2(&idum)-0.5) *mesh_size_x;
        goto L11; 
       }
     else if(i==nx0){
        x_position = 0.5*random2(&idum)*mesh_size_x;
        goto L11;
       }
     else { // if(i==nx1){
        x_position = device_length-0.5*random2(&idum)*mesh_size_x;
       }
       
     L11: //Label: lam cac lenh tiep theo
        
     if(x_position <0.0) {x_position=-x_position;}
       
     if(x_position > device_length){
        factor = fabs(x_position - device_length);
        x_position = device_length - factor;
     }
  
    //    Map particle atributes
    p[ne][2] = x_position; // [m]
    
    //printf("\n Init real space x position=%le",x_position);
    //if(p[ne][2]>1.0e-8){printf("\n Ket qua: x=%le ",p[ne][2]); getch();}           
    return;
}


