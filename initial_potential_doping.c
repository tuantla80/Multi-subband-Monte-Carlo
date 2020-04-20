/* ****************************************************************************
   De initial potential
   - Yeu cau: dang potential chinh xac
   - De HOI TU nhanh

Gop 2 ham doping_potential_initialization() va apply volgate() lai.
NOTE: Cach khoi tao potential da advance RAT NHIEU
Starting date: May 30, 2010
Latest update: May 31, 2010
***************************************************************************** */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "mpi.h"
#include "nrutil.h"
#include "constants.h"

void initial_potential_doping(){
   double ***Get_doping();
   double ***doping = Get_doping();// la donor thi dau + ; con acceptor thi dau tru 
          
   // Toa do
   int Get_nx0(),Get_nxa(),Get_nxb(),Get_nx1();
   int Get_ny0(),Get_nya(),Get_nyb(),Get_ny1();
   int Get_nz0(),Get_nza(),Get_nzb(),Get_nz1();
   int Get_nxgate0(),Get_nxgate1();
   int nx0 = Get_nx0(); int nxa = Get_nxa(); int nxb = Get_nxb(); int nx1 = Get_nx1(); 
   int ny0 = Get_ny0(); int nya = Get_nya(); int nyb = Get_nyb(); int ny1 = Get_ny1();
   int nz0 = Get_nz0(); int nza = Get_nza(); int nzb = Get_nzb(); int nz1 = Get_nz1();
   int nxgate0 = Get_nxgate0();// gate
   int nxgate1 = Get_nxgate1();// gate

   double Get_mesh_size_x(),Get_mesh_size_y(), Get_mesh_size_z();
   double mesh_size_x = Get_mesh_size_x();
   double mesh_size_y = Get_mesh_size_y();
   double mesh_size_z = Get_mesh_size_z();

   /* ************ For Doping part ********************************** */
   int find_region(int i,int j,int k);
   double *Get_doping_density();
   double *dop_term = Get_doping_density();//1,2 or 3 ->S,D or Channel
   //printf("\n Source doping </cm3> = %5.1le",dop_term[1]); printf("\n Drain doping </cm3>  = %5.1le",dop_term[2]);  printf("\n Channel doping </cm3>= %5.2le",dop_term[3]);
   
   int i,j,k,n=0;
   
   for(i=nx0; i<=nx1; i++){
      for(j=ny0; j<=ny1; j++){
         for(k=nz0; k<=nz1; k++){
            n = find_region(i,j,k); // n=1,2,3 or 0 ->o region nao S,D,Channel or outside (Oxide)
            if((n==1)||(n==2)||(n==3)){// S,D,Channel; //printf("\n n = %d",n);
              doping[i][j][k] = dop_term[n]; // [1/m3]
	    } 
            else { // n=0 Oxide region
	      doping[i][j][k] =1.0; //1.0e+10; // ;/1.0;// 1.0e+5;// Cho 1 gia tri rat nho. Ho cho la 1.0e+5 cho cm3
	    }
	 } // End of  for(k=nz0; k<=nz1; k++){  
      }// End of for(j=ny0; j<=ny1; j++){
   }// End of  for(i=nx0; i<=nx1; i++){
   /* ************************End of For Doping part****************** */
  return;
}// End of void initial_potential()

//*******************************************************************************
// Ham tinh Fermi
#include <math.h>
double Fermihalf (double x){ 
 // Fermi Dirac integral of order 1/2. Cong thuc (10) in Notes on Fermi-Dirac Integrals - Lundstrom
 double a = x*x*x*x + 33.6*x*(1.0-0.68*exp(-0.17*(x+1.0)*(x+1.0)))+50;
 
  double y = 1.0/( 0.75*sqrt(3.14159)*pow(a,-3.0/8.0) + exp(-x));
 return y;
}
//************************************************************************************

void save_initial_potential_doping(){
    
     int Get_nx0(),Get_nx1(),Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int ny0 = Get_ny0(); int ny1 = Get_ny1();
     int nz0 = Get_nz0(); int nz1 = Get_nz1();

     double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     double ***Get_doping();
     double ***doping = Get_doping();
     
     double ***Phi = GetPot();
     
     int i,j,k;
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank ==0){ 
       FILE *f; f = fopen("potential_doping_initial_yz_at_midpoint_x.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f==NULL){ printf("\n Cannot open file potential_doping_initial_yz_at_midpoint_x.dat"); return ; }
       fprintf(f,"\n # y[nm] z[nm]  doping[1/m3] potential[eV] \n");
     
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       for(j=ny0; j<=ny1; j++){
	 for(k=nz0; k<=nz1; k++){
	    fprintf(f,"%le  %le %le %le \n",j*mesh_size_y/1.0e-9,k*mesh_size_z/1.0e-9,doping[i][j][k],-Phi[i][j][k]);
	 } // NOTE: ta bieu dien o dang -Phi nhe                                                        
       }
       fclose(f);
     
       // Thu va tuong duong voi mpot.r0 cua GS
       FILE *f1; f1 = fopen("potential_initial_x_at_midpoint_yz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f1==NULL){ printf("\n Cannot open file potential_initial_x_at_midpoint_yz.dat"); return ; }
       fprintf(f1,"\n # x[nm] potential_x[eV] \n");
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc z
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f1);

       FILE *f2; f2 = fopen("potential_doping_initial_xy_at_midpoint_z.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f2==NULL){ printf("\n Cannot open file potential_doping_initial_xy_at_midpoint_z.dat "); return ; }
       fprintf(f2,"\n # x[nm] y[nm]  doping[1/m3] potential[eV] \n");
     
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           for(j=ny0; j<=ny1; j++){
	    fprintf(f2,"%le  %le %le %le \n",i*mesh_size_x/1.0e-9,j*mesh_size_y/1.0e-9,doping[i][j][k],-Phi[i][j][k]);
	 } // NOTE: ta bieu dien o dang -Phi nhe                                                        
       }
       fclose(f2);
       
        
       FILE *f3; f3 = fopen("potential_initial_y_at_midpoint_xz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f3==NULL){ printf("\n Cannot open file potential_initial_y_at_midpoint_xz.dat"); return ; }
       fprintf(f3,"\n # y[nm] potential_y[eV] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(j=ny0; j<=ny1; j++){
           fprintf(f3,"%le %le \n",j*mesh_size_y/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f3);

       FILE *f4; f4 = fopen("potential_initial_z_at_midpoint_xy.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f4==NULL){ printf("\n Cannot open file potential_initial_z_at_midpoint_xy.dat"); return ; }
       fprintf(f4,"\n # z[nm] potential_z[eV] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc y
       
       for(k=nz0; k<=nz1; k++){
          fprintf(f4,"%le %le \n",k*mesh_size_z/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f4);

     }// End of if(myrank ==0)
     return;
} // End of void initial_potential_doping()
//************************************************************************************************************************

void save_initial_charge_density_for_Poisson(){
     
     int Get_nx0(),Get_nx1(),Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int ny0 = Get_ny0(); int ny1 = Get_ny1();
     int nz0 = Get_nz0(); int nz1 = Get_nz1();

     double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     double  ***Get_electron_density(); 
     double ***electron_density = Get_electron_density();

     int i,j,k;
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank ==0){ 
       FILE *f; f = fopen("charge_initial_for_Poisson_yz_at_midpoint_x.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f==NULL){ printf("\n Cannot open file charge_initial_for_Poisson_yz_at_midpoint_x.dat "); return ; }
       fprintf(f,"\n # y[nm] z[nm]  electron[1/m3] hole[1/m3] \n");
     
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       for(j=ny0; j<=ny1; j++){
	 for(k=nz0; k<=nz1; k++){
	    fprintf(f,"%le  %le %le \n",j*mesh_size_y/1.0e-9,k*mesh_size_z/1.0e-9,electron_density[i][j][k]);
	 }                                                      
       }
       fclose(f);
     
       
       FILE *f1; f1 = fopen("charge_initial_for_Poisson_x_at_midpoint_yz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f1==NULL){ printf("\n Cannot open file charge_initial_for_Poisson_x_at_midpoint_yz.dat "); return ; }
       fprintf(f1,"\n # x[nm] electron[1/m3] \n");
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc z
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f1);

       FILE *f2; f2 = fopen("charge_initial_for_Poisson_xy_at_midpoint_z.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f2==NULL){ printf("\n Cannot open file charge_initial_for_Poisson_xy_at_midpoint_z.dat "); return ; }
       fprintf(f2,"\n # x[nm] y[nm] electron[1/m3] \n");
     
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           for(j=ny0; j<=ny1; j++){
	    fprintf(f2,"%le  %le %le \n",i*mesh_size_x/1.0e-9,j*mesh_size_y/1.0e-9,electron_density[i][j][k]);
	 }                                                        
       }
       fclose(f2);
       
        
       FILE *f3; f3 = fopen("charge_initial_for_Poisson_y_at_midpoint_xz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f3==NULL){ printf("\n Cannot open file charge_initial_for_Poisson_y_at_midpoint_xz.dat"); return ; }
       fprintf(f3,"\n # y[nm] electron[1/m3] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(j=ny0; j<=ny1; j++){
           fprintf(f3,"%le %le \n",j*mesh_size_y/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f3);

       FILE *f4; f4 = fopen("charge_initial_for_Poisson_z_at_midpoint_xy.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f4==NULL){ printf("\n Cannot open file charge_initial_for_Poisson_z_at_midpoint_xy.dat"); return ; }
       fprintf(f4,"\n # z[nm] electron[1/m3] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc y
       
       for(k=nz0; k<=nz1; k++){
          fprintf(f4,"%le %le \n",k*mesh_size_z/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f4);

     }// End of if(myrank ==0)
     return;
} // End of void save_initial_charge_density_for_Poisson()

//**********************************************************************************************
//************************************************************************************************************************

#include "mpi.h"
void save_charge_density_for_Poisson(){
     
     int Get_nx0(),Get_nx1(),Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int ny0 = Get_ny0(); int ny1 = Get_ny1();
     int nz0 = Get_nz0(); int nz1 = Get_nz1();

     double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     double  ***Get_electron_density(); 
     double ***electron_density = Get_electron_density();

     int i,j,k;
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank ==0){ 
       FILE *f; f = fopen("charge_density_for_Poisson_yz_at_midpoint_x.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f==NULL){ printf("\n Cannot open file charge_density_for_Poisson_yz_at_midpoint_x.dat "); return ; }
       fprintf(f,"\n # y[nm] z[nm]  electron[1/m3] hole[1/m3] \n");
     
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       for(j=ny0; j<=ny1; j++){
	 for(k=nz0; k<=nz1; k++){
	    fprintf(f,"%le  %le %le \n",j*mesh_size_y/1.0e-9,k*mesh_size_z/1.0e-9,electron_density[i][j][k]);
	 }                                                      
       }
       fclose(f);
     
       
       FILE *f1; f1 = fopen("charge_density_for_Poisson_x_at_midpoint_yz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f1==NULL){ printf("\n Cannot open file charge_density_for_Poisson_x_at_midpoint_yz.dat "); return ; }
       fprintf(f1,"\n # x[nm] electron[1/m3] \n");
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc z
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f1);

       FILE *f2; f2 = fopen("charge_density_for_Poisson_xy_at_midpoint_z.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f2==NULL){ printf("\n Cannot open file charge_density_for_Poisson_xy_at_midpoint_z.dat "); return ; }
       fprintf(f2,"\n # x[nm] y[nm] electron[1/m3] \n");
     
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           for(j=ny0; j<=ny1; j++){
	    fprintf(f2,"%le  %le %le \n",i*mesh_size_x/1.0e-9,j*mesh_size_y/1.0e-9,electron_density[i][j][k]);
	 }                                                        
       }
       fclose(f2);
       
        
       FILE *f3; f3 = fopen("charge_density_for_Poisson_y_at_midpoint_xz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f3==NULL){ printf("\n Cannot open file charge_density_for_Poisson_y_at_midpoint_xz.dat"); return ; }
       fprintf(f3,"\n # y[nm] electron[1/m3] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(j=ny0; j<=ny1; j++){
           fprintf(f3,"%le %le \n",j*mesh_size_y/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f3);

       FILE *f4; f4 = fopen("charge_density_for_Poisson_z_at_midpoint_xy.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f4==NULL){ printf("\n Cannot open file charge_density_for_Poisson_z_at_midpoint_xy.dat"); return ; }
       fprintf(f4,"\n # z[nm] electron[1/m3] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc y
       
       for(k=nz0; k<=nz1; k++){
          fprintf(f4,"%le %le \n",k*mesh_size_z/1.0e-9,electron_density[i][j][k]);
       }
       fclose(f4);

     }// End of if(myrank ==0)
     return;
} // End of void initial_potential_doping()

//**********************************************************************************************


/* ****************************************************************************
   De ve potential tai bat cu time step nao ma ta DAT HAM NAY o do

Starting date: June 03, 2010
Latest update: June 03, 2010
***************************************************************************** */

void save_potential(){
    
     int Get_nx0(),Get_nx1(),Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int ny0 = Get_ny0(); int ny1 = Get_ny1();
     int nz0 = Get_nz0(); int nz1 = Get_nz1();

     double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
    
     double ***Phi = GetPot();

     int i,j,k;
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank ==0){ 
       FILE *f; f = fopen("potential_yz_at_midpoint_x.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f==NULL){ printf("\n Cannot open file potential_yz_at_midpoint_x.dat"); return ; }
       fprintf(f,"\n # y[nm] z[nm] potential[eV] \n");
     
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       for(j=ny0; j<=ny1; j++){
	 for(k=nz0; k<=nz1; k++){
	    fprintf(f,"%le  %le %le \n",j*mesh_size_y/1.0e-9,k*mesh_size_z/1.0e-9,-Phi[i][j][k]);
	 } // NOTE: ta bieu dien o dang -Phi nhe                                                        
       }
       fclose(f);
     
       // Thu va tuong duong voi mpot.r0 cua GS
       FILE *f1; f1 = fopen("potential_x_at_midpoint_yz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f1==NULL){ printf("\n Cannot open file potential_x_at_midpoint_yz.dat"); return ; }
       fprintf(f1,"\n # x[nm] potential_x[eV] \n");
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc z
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f1);

       FILE *f2; f2 = fopen("potential_xy_at_midpoint_z.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f2==NULL){ printf("\n Cannot open file potential_xy_at_midpoint_z.dat "); return ; }
       fprintf(f2,"\n # x[nm] y[nm] potential[eV] \n");
     
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(i=nx0; i<=nx1; i++){
           for(j=ny0; j<=ny1; j++){
	    fprintf(f2,"%le  %le %le \n",i*mesh_size_x/1.0e-9,j*mesh_size_y/1.0e-9,-Phi[i][j][k]);
	 } // NOTE: ta bieu dien o dang -Phi nhe                                                        
       }
       fclose(f2);
       
        
       FILE *f3; f3 = fopen("potential_y_at_midpoint_xz.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f3==NULL){ printf("\n Cannot open file potential_y_at_midpoint_xz.dat"); return ; }
       fprintf(f3,"\n # y[nm] potential_y[eV] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       k=(int)((nz1-nz0)/2.0 + 0.5);// lay o giua truc z
       for(j=ny0; j<=ny1; j++){
           fprintf(f3,"%le %le \n",j*mesh_size_y/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f3);

       FILE *f4; f4 = fopen("potential_z_at_midpoint_xy.dat","w"); // "w" neu tep ton tai no se bi xoa
       if(f4==NULL){ printf("\n Cannot open file potential_z_at_midpoint_xy.dat"); return ; }
       fprintf(f4,"\n # z[nm] potential_z[eV] \n");
       i=(int)((nx1-nx0)/2.0 + 0.5);// lay o giua truc x
       j=(int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc y
       
       for(k=nz0; k<=nz1; k++){
          fprintf(f4,"%le %le \n",k*mesh_size_z/1.0e-9,-Phi[i][j][k]);// NOTE: ta bieu dien o dang -Phi nhe 
       }
       fclose(f4);

     }// End of if(myrank ==0)
     return;
} // End of void initial_potential_doping()
//************************************************************************************************************************

    
