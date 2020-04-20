/* ***********************************************************************
 Cach tinh: 
   + Buoc 1. Tinh electron_population[s][v][i]: N_s,v,i 
          - Phuong phap 1: Nearest grid points cho huong thep phuong x
          - Phuong phap 2: NEC (Nearest-Element-Center)
   + Buoc 2. Tinh Ns la tong cua N_s,v,i theo v va i
   + Buoc 3. Tinh Psi_s,v,i(y,z) square
   + Buoc 4 va 5. Tinh electron density [1/m3]
 electron density se la dau vao cho Poisson  

Starting date: March 30, 2010
Update:        March 31, 2010 
Update:        May 20, 2010 MOT LOI rat CO BAN. Can reset gia tri electron_density cho lan chay ke tiep. Vi no la tinh tong va lap di lap lai nhieu lan theo thoi gian
Latest update: June 2, 2010. Da mo rong electron density ra ca vung Oxide nhung KHONG TINH o BIEN
************************************************************************* */
#include <stdio.h>
#include "mpi.h"
#include "nrutil.h"

void electron_density_caculation(){
     // Goi ham
     double ***Get_electron_density();
     int Get_save_electron_population();
     int Get_save_electron_density();
     void save_electron_population(double ***electron_population, int nx_max, int v_max, int NSELECT);
     
     int Get_NSELECT();
     int Get_n_used();// Vi ham nay NGAY SAU ham delete_particle() nen n_used o day CHI BAO GOM iv=1,2,3
     double **Get_p(),Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     int *Get_valley(),*Get_subband();
     double *****Get_wave(); // Lay wave function tu 2D Schrodinger
         
     // Cac bien local
     double ***electron_density = Get_electron_density();
     int electron_population_save_or_not = Get_save_electron_population();// Co save vao file khong ?
     int electron_density_save_or_not    = Get_save_electron_density();

     int Get_nx0(),Get_nx1(), Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();// da thay doi thanh electron_density = d3matrix(nx0,nx1,0,ny,0,nz);ny=nyb-nya, nz=nzb-nza
     int nya = Get_nya(); int nyb = Get_nyb();// Nhung do wave function chi lay o LOI Silicon nen ta chi tinh trong doan do thoi
     int nza = Get_nza(); int nzb = Get_nzb(); 

     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
     int ny_max = nyb - nya;// De thay doi it nhat co the
     int nz_max = nzb - nza; 
    
     int NSELECT = Get_NSELECT();
     int n_used = Get_n_used(); //printf("\n n_used at electron_density_caculation() = %d",n_used);
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
     double **p   = Get_p();
     int *valley  = Get_valley();
     int *subband = Get_subband();
     double *****wave   = Get_wave();
     
// Thuc hien
    int s,v,i,n,j,k;
    double x_position;
    double ***electron_population;// N_s,v,i: electron population cho MOI section s-th, valley v-th va subband i-th
    double *Ns; // Bien trung gian la tong cua N_s,v,i theo v va i
    double *****wave_square; // bien trung gian de tinh wave * wave
    electron_population = d3matrix(0,nx_max,1,3,1,NSELECT);// Can khoi tao
    Ns = dvector(0,nx_max); // chay cho tat ca cac section
    wave_square = d5matrix(0,nx_max,1,3,1,NSELECT,0,ny_max,0,nz_max);
    
    for(s=0; s<=nx_max; s++){ Ns[s] = 0.0;
        for(v=1; v<=3; v++){     
            for(i=1; i<=NSELECT; i++){ electron_population[s][v][i]=0.0;
               for(j=0; j<=ny_max; j++){  // chay cho diem tren truc y
                  for(k=0; k<=nz_max; k++){// chay cho diem tren truc z  
                     wave_square[s][v][i][j][k] = 0.0;
                  }  
               }
            }
        }
     } // End of for(s=0; s<=nx_max; s++){
   
   // Buoc 1: Tinh electron population  
   /*  //Cach 1. Nearest grid points          
   for(n = 1; n<=n_used; n++){
       if(valley[n] != 9){
          x_position = p[n][2];
          s = (int)(x_position/mesh_size_x + 0.5);// vi tri section cua hat n-th. PHAI CONG THEM 0.5 nhe 
          if(s < 0)       { s = 0;}
          if(s > nx_max)  { s = nx_max;} // De khong vuot qua chi so mang
          v = valley[n];  // valley cua hat thu n-th
          i = subband[n]; // subband cua hat n-th
          electron_population[s][v][i] += 1; // tai cac vi tri s,v,i thi tang len 1
       }// End of if(valley[n] != 9)
    }// End of for(n = 1; n<=n_used; n++)
    */ // Ket thuc cach 1
    
   ///* // Cach 2: NEC (Nearest-Element-Center)
   for(n = 1; n<=n_used; n++){
      if(valley[n] != 9){// chi xet hat co iv=1,2,3 ma thoi. Thuc ra khong can lenh nay. Nhung lam cho chac
         x_position = p[n][2];
         s = (int)(x_position/mesh_size_x + 0.5); // vi tri section cua hat n-th
         if(s < 0)    { s=0;}
         if(s>=nx_max){ s=nx_max-1;}// do phuong phap NEC yeu cau
         v = valley[n];  // valley cua hat thu n-th
         i = subband[n]; // subband cua hat n-th
         electron_population[s][v][i] += 0.5; // 
         electron_population[s+1][v][i] += 0.5;
      }// End of if(valley[n] != 9)
    }// End of for(n = 1; n<=n_used; n++)
      // Trick de tranh tao ra electron density qua bien doi. Co can dung TRICK khong? CHUA BIET 01/04/10 14:28
      //electron_population[0][v][i] = 2.0*electron_population[0][v][i];//((s==0)||(s==nx_max)) thi nhan cho 2
      //electron_population[nx_max][v][i] = 2.0*electron_population[nx_max][v][i];
    //*/ // Ket thuc cach 2  
      
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
     if(myrank ==0){ // Save tai mot node thoi
        //SAVE electron population or NOT
        if (electron_population_save_or_not ==1){ // =1: Co SAVE
             int v_max=3;// la 3 valley pairs vao the nay cho chuyen nghiep
             save_electron_population(electron_population,nx_max,v_max,NSELECT);
        }
     }

    // Buoc 2. Tinh Ns la tong cua N_s,v,i the v va i
    for(s=0; s<=nx_max; s++){ 
        for(v=1; v<=3; v++){     
            for(i=1; i<=NSELECT; i++){ 
               Ns[s] += electron_population[s][v][i];
             }
         }
     }// End of for(s=0; s<=nx_max; s++)

    //for(s=0; s<=nx_max; s++){
    //printf("\n Number of partice at s-th %d =%f",s,Ns[s]);     
        //if(Ns[s]== 0.0){
	  //printf("\n Ns[%d]=%f",s,Ns[s]);        
          // Ns[s] = 1.0e-1;} //(May 28, 2010) (SAI: Co the) trick de tranh chia cho 0 luc tinh electron_density. Xay ra luc KIEM TRA nhu the nay       
    //}
        
    // Buoc 3: Tinh Psi_s,v,i(y,z) square
    for(s=0; s<=nx_max; s++){
        for(v=1; v<=3; v++){     
            for(i=1; i<=NSELECT; i++){
	      for(j=0; j<=ny_max; j++){ // Chay tu 0 da the hien so diem chia la + 1 roi (May 30, 2010) 
                  for(k=0; k<=nz_max; k++){
                     wave_square[s][v][i][j][k] = wave[s][v][i][j][k]*wave[s][v][i][j][k];
                  }  
               }
            }
        }
     } // End of for(s=0; s<=nx_max; s++){ 
         
    // Buoc 4. Tinh electron_density [1/m3]
    // Truoc khi tinh electron_density can reset no vi neu khong lai tinh ca lan truoc do -  no la tinh TONG va duoc lap di lap lai theo thoi gian
    // Buoc 4.1. Reset. NHUNG CHI O PHAN LOI SILICON THOI (June 02, 2010)
    for(s=0; s<=nx_max; s++){
       for(j=nya; j<=nyb; j++){  
           for(k=nza; k<=nzb; k++){
	     electron_density[s][j][k]=0.0;// Phan Loi silicon cho ve 0
	   }
       }
    }
    // Buoc 4.2. Tinh toan gia tri CAN DICH HE TOA DO
     for(s=0; s<=nx_max; s++){
       if(Ns[s] != 0.0){// Vi neu khong co hat nao o trong section thi co nghia la electron_population=0 do do electron_density =0 (May 28, 2010)
        for(v=1; v<=3; v++){     
            for(i=1; i<=NSELECT; i++){
               for(j=0; j<=ny_max; j++){  
                  for(k=0; k<=nz_max; k++){
		    
                      electron_density[s][j+nya][k+nza] += (electron_population[s][v][i]/Ns[s])*wave_square[s][v][i][j][k];
                    // printf("\n electron_density = %le",electron_density[s][j][k]);
		    }
                  }
               }
            }
        }// End of if(Ns[s] =! 0.0){ 
      }// End of for(s=0; s<=nx_max; s++) 
   
     
    // Buoc 5: Tinh electron density theo THE TICH [1/m3]. Chi tinh tai loi Silicon
    //double volume = mesh_size_x*mesh_size_y*mesh_size_z;//the tich. SAI hat bay gio theo phuong x la classical con theo y va z la quantum
     double volume = mesh_size_x*mesh_size_y*ny_max*mesh_size_z*nz_max;
    for(s=0; s<=nx_max; s++){
        for(j=0; j<=ny_max; j++){  
            for(k=0; k<=nz_max; k++){   
                electron_density[s][j+nya][k+nza] = electron_density[s][j+nya][k+nza]/volume;// [1/m3]
            }
        }
    }// End of for(s=0; s<=nx_max; s++)
     
    // Free local variables
    free_d3matrix(electron_population,0,nx_max,1,3,1,NSELECT);
    free_dvector(Ns,0,nx_max);
    free_d5matrix(wave_square,0,nx_max,1,3,1,NSELECT,0,ny_max,0,nz_max);
 return;
} // End of electron_density_caculation

/* *****************************************************************************
 De SAVE N_s,v,i: electron population cho MOI section s-th, valley v-th va subband i-th
 INPUT:  + electron_population[s][v][i] 
         + s:  so section di tu 0 den nx_max
         + v=1,2,or 3: 3 cap valley pairs 
         + i=1 den NSELECT so subbands cho moi valley
         
 OUTPUT: + SAVE electron_population[s][v][i]
        
Starting date: March 30, 2010
Latest update: March 30, 2010
***************************************************************************** */
void save_electron_population(double ***electron_population,int nx_max,int v_max,int NSELECT){
     
     FILE *f; f = fopen("Electron_population.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file Electron_population.dat"); return ; }
     fprintf(f,"\n #section valley subband Number_of_Electron \n");
     int s,v,i;
     for(s=0; s<=nx_max; s++){ 
        for(v=1; v<=v_max; v++){     
            for(i=1; i<=NSELECT; i++){ 
               fprintf(f,"   %d      %d       %d       %f\n",s,v,i,electron_population[s][v][i]);     
            }
         }
     }
  fclose(f);
  return; 
} // End of void save_Density_of_State_1D(int s_th)

/* *****************************************************************************
 De SAVE electron_density[s][y][z] co nghia la electron_density[x][y][z]
 INPUT:  + electron_density[s-th][j][k] 
         + s-th: dinh save o section nao? Cu the 1 section ma thoi tu 0 den nx_max
         + j: diem chay theo phuong y di tu 0 den ny_max (ny_max = nyb - nya)
         + k: diem chay theo phuong z di tu 0 den nz_max (nz_max = nzb - nza)
       
 OUTPUT: + SAVE electron_density[s-th][y][z] tai s-th section
        
Starting date: March 31, 2010
Latest update: June 15, 2010
***************************************************************************** */
void save_electron_density(char *fn)
  {
     double ***Get_electron_density();
     double ***electron_density = Get_electron_density();

     double Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     int Get_nx0(),Get_nx1(), Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nya = Get_nya(); int nyb = Get_nyb();
     int nza = Get_nza(); int nzb = Get_nzb(); 

     int nx_max = nx1-nx0; 
     int ny_max = nyb - nya;// De thay doi it nhat co the
     int nz_max = nzb - nza; 

     FILE *f;f = fopen(fn,"w"); // "w" neu tep ton tai no se bi xoa
     //f = fopen("Electron_Density.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file at save_electron_density.c"); return ; }
     fprintf(f,"\n   #y           z          electron_density[1/m3] section_th \n");

     int myrank, i,j, s_th;
     double y,z; // x
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank == 0){ // Save tai mot node thoi
       s_th =(int)(nx_max/2 + 0.5); // Luon save cho section o GIUA // NOTE: nen o nhung time step dau tien co the no se bang 0
	                             // ngay time step dau tien thi voi s_th=3 Ket qua electron_density rat dep
       for(i=0; i<=ny_max; i++){ 
	 y = (double)(i)*mesh_size_y/1.0e-9;// doi tu [m] sang [nm]       
	 for(j=0; j<=nz_max; j++){
            z = (double)(j)*mesh_size_z/1.0e-9;// doi tu [m] sang [nm] 
            fprintf(f," %le %le  %le %d \n",y,z,electron_density[s_th][i+nya][j+nza],s_th);      
	 }
       }
     }
     fclose(f);   
  return; 
} // End of void save_Density_of_State_1D(int s_th)

/* *****************************************************************************
 De SAVE electron_density[s][y][z] co nghia la electron_density[x][y][z]
 INPUT:  + electron_density[s-th][j][k] 
         + s-th: save for all sections
         + j: diem chay theo phuong y di tu 0 den ny_max (ny_max = nyb - nya)
         + k: diem chay theo phuong z di tu 0 den nz_max (nz_max = nzb - nza)
       
 OUTPUT: + SAVE electron_density[s][y][z] tai all sections
 Note: For checing only   
Starting date: June 04, 2010
Latest update: June 04, 2010
***************************************************************************** */
void save_electron_density_all_sections(char *fn)
  {
     double ***Get_electron_density();
     double ***electron_density = Get_electron_density();
     int Get_nx0(),Get_nx1(), Get_nya(),Get_nyb(), Get_nza(),Get_nzb();// int Get_nya(), Get_nza();// Can dich he toa do
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nya = Get_nya(); int nyb = Get_nyb();
     int nza = Get_nza(); int nzb = Get_nzb(); 
     
     int ny_max = nyb - nya;// De thay doi it nhat co the
     int nz_max = nzb - nza; 
     
     double Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
   
// Thuc hien
     FILE *f; f = fopen(fn,"w"); // "w" neu tep ton tai no se bi xoa
     //f = fopen("Electron_Density_All_Sections.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL) { printf("\n Cannot open file at save_electron_density_all_sections.c"); return ; }
     fprintf(f,"\n   #y           z          electron_density[1/m3] section_th \n");
     
     int s,i,j,myrank;
     double y,z; // x
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

     if(myrank ==0){ // Save tai mot node thoi
       for(s=nx0; s<=nx1; s++){
	 for(i=0; i<=ny_max; i++){ 
	   y = (double)(i)*mesh_size_y/1.0e-9;// doi tu [m] sang [nm]       
	   for(j=0; j<=nz_max; j++){
	     z = (double)(j)*mesh_size_z/1.0e-9;// doi tu [m] sang [nm] 
	     fprintf(f," %le %le  %le %d \n",y,z,electron_density[s][i+nya][j+nza],s);      
	   }
	 }
       }
     }
  fclose(f);
  return; 
} // End of void save_Density_of_State_1D(int s_th)
