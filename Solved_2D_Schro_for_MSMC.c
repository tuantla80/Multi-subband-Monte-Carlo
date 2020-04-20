/* ****************************************************************
 Ham de thuc hien viec giai 2D Schrodinger cho tat ca cross-section
 va 3 cap valley pairs.
 INPUT:  + s: s-th section
         + v: v-th valley pair
         + m1,m2: la hai tham so cho cac valley pairs. Vi du (m1=0.19, m2=0.916)
         + potential Phi(i,j,k) tu viec guest potential va sau do tu 3D Poisson
 
 OUTPUT: + E_s,v,i (eigen energy for each cross-section, valley pair and subband
         + Wave_s,v,i(y,z): wave tuong ung       
        
Starting date: Feb 23, 2010
Update:        Feb 24, 2010
Update:        May 09, 2010 (Them tham so de mapping voi 3D Poisson)
Update:        May 10, 2010 (Parallel)
Latest update: May 30: SUA SAI o SO DIEM CHIA. THUC TE SO DIEM CHIA o cross-section CAN TANG THEM 1
****************************************************************** */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "nrutil.h"
#include "constants.h"
#include "mpi.h"

void Solved_2D_Schro_Parallel_for_MSMC(int *argc, char ***argv){
     // Goi cac ham
     void Schrodinger2D(int s, int v, double m1, double m2);
   
     // Cac bien local
     int Get_nx0(),Get_nx1(); // int Get_nx_max();
     int nx0 = Get_nx0(); 
     int nx1 = Get_nx1();
     int nx_max = nx1-nx0; //printf("\n nx0=%d, nx1=%d,nx_max=%d",nx0,nx1,nx_max);// de thay doi it nhat
     int NXMAX = nx_max+1;// Do ta chay section tu 0 den nx_max nen tong so lan can chay la NXMAX
     double Get_ml(), Get_mt();
     double ml = Get_ml();
     double mt = Get_mt(); //printf("\n ml=%f, mt=%f", ml,mt);
     int Get_NSELECT();
     int NSELECT = Get_NSELECT();
     int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 
     int ny_max = nyb - nya;
     int nz_max = nzb - nza; 
     int ny = ny_max + 1;// (May 30, 2010) int ny = ny_max; Do so diem chia tren loi Silicon can +1 
     int nz = nz_max + 1; //int nz = nz_max;
          
     // Thuc hien
     int s,v; // s-th section, v-th valley pair
     double m1, m2;
     int size,rank;// So luong processor va rank tuong ung
     int ierr; // thong bao loi
 
     // Buoc 1. Lay rank va size
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);// Lay number of processor;  // Initialize MPI o ham main() roi
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);// lay rank cua this process
     
     if ( rank==0 ){ // master node
         printf ( "\n" );
 	     printf ( "The number of processes is %d\n", size );
	     if (size==1){ 
	       //printf("\n KHONG CHAY SONG SONG o Master node: STOP");
	       //exit(1);//Se dung lenh nay khi KHONG su dung master node
	     }
      }// End of if ( rank==0 )

     //NOTE.chay tai rank=0,1,2,3,4.
      printf("\n Chay cho cac valley");

      // Buoc 2. Chay SONG SONG cho cac node. Moi node ung voi rank cua minh se giai 2D Schrodinger o cac section khac nhau
      // Dau tien phai cho eig va wave ve 0. Neu chay 1 node thi khong can nhung ta lai chay nhieu node va moi node thi ghi 1 phan du lieu vao. Nen can
      // xoa du lieu cu
      double ***Get_eig(), *****Get_wave();
      double ***eig      = Get_eig();
      double *****wave   = Get_wave();
      int i,j,k;
      for(s=nx0; s<=nx1; s++){
        for(v=1; v<=3; v++){
            for(i=1; i<=NSELECT; i++){ eig[s][v][i] = 0.0; 
	      for(j=0; j<=ny_max; j++){// Chay tu 0 den ny_max = nyb-nya da the hien so diem chia can cong them 1 
                   for(k=0; k<=nz_max; k++){       
                       wave[s][v][i][j][k] = 0.0;
                   }
               }
            }
        }
     }
      
      // Valley pair 1.
      m1 = mt;  m2 = mt; v  = 1;
      for(s=NXMAX*(rank)/(size); s<=(NXMAX*(rank+1))/(size)-1; s++){// Cong thuc khi chay ca Master node
            //printf("\n ***********************");
            Schrodinger2D(s,v,m1,m2);
	    //printf("\n The section %d with valley %d are solved",s,v);
        }
      // Valley pair 2
      m1 = ml;  m2 = mt; v  = 2;
      for(s=NXMAX*(rank)/(size); s<=(NXMAX*(rank+1))/(size)-1; s++){// Cong thuc khi chay ca Master node
            //printf("\n ***********************");
            Schrodinger2D(s,v,m1,m2);
	    //printf("\n The section %d with valley %d are solved",s,v);
        } 
     // Valley pair 3. // chi can giai khi no khong doi xung. Tinh Sau
      m1 = mt;  m2 = ml; v  = 3;   
      for(s=NXMAX*(rank)/(size); s<=(NXMAX*(rank+1))/(size)-1; s++){// Cong thuc khi chay ca Master node
            //printf("\n ***********************");
            Schrodinger2D(s,v,m1,m2);
	    //printf("\n The section %d with valley %d are solved",s,v);
        } 

      // Buoc 3. Sau khi chay song song thi cac MANG eig va wave o cac node chi luu cac gia tri ma no SOLVE thoi
      // LAP DAY cac gia tri cua eig va wave o cac node. Ta KHOI TAO mang eig va wave giong nhau ve KICH THUOC cho tat ca cac node.
      void convert_3Dmatrix_to_array(int M, int N,int P, double ***mtx, double *a);// Su dung cho eigen energy
      void convert_array_to_3Dmatrix(int M, int N,int P, double *a, double ***mtx);
      void convert_5Dmatrix_to_array(int M, int N,int P,int Q,int R, double *****mtx, double *a);// su dung cho wave
      void convert_array_to_5Dmatrix(int M, int N,int P,int Q,int R, double *a,double *****mtx);
      int Get_save_eig_wavefuntion(); // Flag chi ra co SAVE hay KHONG?
   
      // Buoc 3.1. Thuc hien cho eigen
      int count1 = NXMAX*3*NSELECT;//number of elements cua eigen energy to be sent. PHU THUOC section, valley, subband
      double *eig_array = dvector(0,count1);//do ***eig la mang 3 chieu ta can chuyen no sang mang 1 chieu
      double *eig_sum =  dvector(0,count1);//tong eigen-energy tu TAT CA cac node: LAP DAY gia tri eigen TUONG UNG cho chay NOI TIEP
      for(i=0; i<=count1; i++){ eig_array[i]=0.0; eig_sum[i]=0.0; }

      convert_3Dmatrix_to_array(nx_max,3,NSELECT,eig,eig_array);// chuyen 3D eigen thanh 1D eigen
      MPI_Allreduce(eig_array,eig_sum,count1+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);// Su dung ham MPI_Allreduce de CONG va PHAN BO cho tat ca cac node
      convert_array_to_3Dmatrix(nx_max,3,NSELECT, eig_sum,eig);// Chuyen lai 1D eigen thanh 3D eigen

      // Buoc 3.2. Thuc hien cho wave (May 13, 2010)
     
      int count2 = NXMAX*3*NSELECT*ny*nz;//number of elements cua wave to be sent. PHU THUOC section, valley, subband, y va z
                                         // vi index 0 den ny_max tu la co ny=ny_max+1 diem               
      double *wave_array = dvector(0,count2);//do *****wave la mang 5 chieu ta can chuyen no sang mang 1 chieu
      double *wave_sum =  dvector(0,count2);//tong wave tu TAT CA cac node: LAP DAY gia tri wave TUONG UNG cho chay NOI TIEP
      for(i=0; i<=count2; i++){ wave_array[i]=0.0; wave_sum[i]=0.0; }
      
      convert_5Dmatrix_to_array(nx_max,3,NSELECT,ny_max,nz_max,wave,wave_array);//chuyen 5D wave thanh 1D wave

      MPI_Allreduce(wave_array,wave_sum,count2+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);// chay 3 section 5x5 thi thoi gian khong dang ke

      convert_array_to_5Dmatrix(nx_max,3,NSELECT,ny_max,nz_max,wave_sum,wave);

    
     //  */ 
      /*	Se dung cai nay
     if (rank !=0){// Ta khong chay o master node. Ta chi chay tai rank=1,2,3,4..thoi
        // Valley pair 1
        m1 = mt;  m2 = mt; v  = 1;
	    for(s=NXMAX*(rank-1)/(size-1); s<=(NXMAX*rank)/(size-1)-1; s++)// Cong thuc khi KHONG DUNG master node
          {
            //printf("\n ***********************");//printf("\n The section %d with valley %d are solving",s,v);
	        Schrodinger2D(s,v,m1,m2);
          }
        
        // Valley pair 2
        m1 = ml;  m2 = mt; v  = 2;
        for(s=NXMAX*(rank-1)/(size-1); s<=(NXMAX*rank)/(size-1)-1; s++)
          {
            //printf("\n ***********************"); //printf("\n The section %d with valley %d are solving",s,v);      
            Schrodinger2D(s,v,m1,m2);
          }
       
        // Valley pair 3. // chi can giai khi no khong doi xung. Tinh Sau
        m1 = mt;  m2 = ml; v  = 3;
        for(s=NXMAX*(rank-1)/(size-1); s<=(NXMAX*rank)/(size-1)-1; s++)
          { 
           //printf("\n ***********************");  //printf("\n The section %d with valley %d are solving",s,v);     
           Schrodinger2D(s,v,m1,m2);
          }
      } // End of if (rank !=0){
  */         
     // Free local variables
      free_dvector(eig_array,0,count1);
      free_dvector(eig_sum,0,count1);
      free_dvector(wave_array,0,count2);
      free_dvector(wave_sum,0,count2);
 return;
     
} // End of Solved_2D_Schro_Parallel_for_MSMC(int *argc, char ***argv)

void convert_3Dmatrix_to_array(int M, int N,int P, double ***mtx, double *a)
{//eig = d3matrix(0,nx1,1,3,1,NSELECT);
  int i,j,p;
  int k = 0 ;
 for(p=1; p<=P; p++){
   for ( j=1 ; j<=N ; j++ ){
     for ( i=0 ; i<=M ; i++ ){
      a[k++] = mtx[i][j][p];
     }
   }
 }
}

void convert_array_to_3Dmatrix(int M, int N,int P, double *a, double ***mtx)
{
  int i,j,p;
  int k = 0 ;
for(p=1; p<=P; p++){
  for ( j=1 ; j<=N ; j++ ){
    for ( i=0 ; i<=M ; i++ ){
      mtx[i][j][p] = a[k++];
    }
  }
 }
}

void convert_5Dmatrix_to_array(int M, int N,int P,int Q,int R, double *****mtx, double *a)
{//wave = d5matrix(nx0,nx1,1,3,1,NSELECT,0,ny,0,nz);//cua Loi Silicon; nx0=0;// xem ham initial_variable()
  int i,j,p,qq,r;
  int k = 0 ;
  for(r=0; r<=R; r++){//cho z
    for(qq=0; qq<=Q; qq++){// cho y
	for(p=1; p<=P; p++){// subband
	  for ( j=1 ; j<=N ; j++ ){//valley
	    for ( i=0 ; i<=M ; i++ ){//section
	      a[k++] = mtx[i][j][p][qq][r];
	    }
	  }
	}
    }
  }
}

void convert_array_to_5Dmatrix(int M, int N,int P,int Q,int R,double *a, double *****mtx)
{
  int i,j,p,qq,r;
  int k = 0 ;
  for(r=0; r<=R; r++){//cho z
    for(qq=0; qq<=Q; qq++){// cho y
	for(p=1; p<=P; p++){// subband
	  for ( j=1 ; j<=N ; j++ ){//valley
	    for ( i=0 ; i<=M ; i++ ){//section
	       mtx[i][j][p][qq][r]= a[k++];
	    }
	  }
	}
    }
  }
}
/* ****************************************************************
 Ham de thuc hien viec giai 2D Schrodinger cho tat ca cross-section
 va 3 cap valley pairs.
 INPUT:  + s: s-th section
         + v: v-th valley pair
         + m1,m2: la hai tham so cho cac valley pairs. Vi du (m1=0.19, m2=0.916)
         + potential Phi(i,j,k) tu viec guest potential va sau do tu 3D Poisson
 
 OUTPUT: + E_s,v,i (eigen energy for each cross-section, valley pair and subband
         + Wave_s,v,i(y,z): wave tuong ung       
        
Starting date: Feb 23, 2010
Latest update: Feb 24, 2010
Latest update: May 09, 2010 (Them tham so de mapping voi 3D Poisson)
****************************************************************** */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "constants.h"

void Solved_2D_Schro_for_MSMC(){
     // Goi cac ham
     void Schrodinger2D(int s, int v, double m1, double m2);
     int Get_nx0(),Get_nx1(); // int Get_nx_max();
     double Get_ml(), Get_mt();
     
     // Cac bien local
     int nx0 = Get_nx0(); 
     int nx1 = Get_nx1();
     int nx_max = nx1-nx0; //printf("\n nx0=%d, nx1=%d,nx_max=%d",nx0,nx1,nx_max);// de thay doi it nhat
     double ml = Get_ml();
     double mt = Get_mt(); //printf("\n ml=%f, mt=%f", ml,mt);
          
     // Thuc hien
     int s,v; // s-th section, v-th valley pair
     double m1, m2;
     ///*
     int Get_NSELECT();
     int NSELECT = Get_NSELECT();
     int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 
     int ny_max = nyb - nya;
     int nz_max = nzb - nza; 
     double ***Get_eig(), *****Get_wave();
     double ***eig      = Get_eig();
     double *****wave   = Get_wave();
     int i,j,k;
     // Reset eigen va wave moi lan goi 2D Schrodinger
     for(s=nx0; s<=nx1; s++){
        for(v=1; v<=3; v++){
            for(i=1; i<=NSELECT; i++){ eig[s][v][i] = 0.0; 
	      for(j=0; j<=ny_max; j++){// Chay tu 0 den ny_max = nyb-nya da the hien so diem chia can cong them 1 
                   for(k=0; k<=nz_max; k++){       
                       wave[s][v][i][j][k] = 0.0;
                   }
               }
            }
        }
     }    
     
     // Valley pair 1
     m1 = mt;  m2 = mt; v  = 1;
     for(s=0; s<=nx_max; s++){ 
         //printf("\n ***********************");
         Schrodinger2D(s,v,m1,m2);
	 //printf("\n The section %d with valley %d is solved",s,v);
      }
          
     // Valley pair 2
     m1 = ml;  m2 = mt; v = 2;
     for(s=0; s<=nx_max; s++){
         //printf("\n ***********************");
         Schrodinger2D(s,v,m1,m2);
	 //printf("\n The section %d with valley %d is solved",s,v);
       }
     
     // Valley pair 3 chi can giai khi no khong doi xung. Tinh Sau
     m1 = mt;  m2 = ml; v = 3;
     for(s=0; s<=nx_max; s++){ 
         //printf("\n ***********************");
         Schrodinger2D(s,v,m1,m2);
	 //printf("\n The section %d with valley %d is solved",s,v);
       }
    
      return;
     
} // End of void Solved_2D_Schro_for_MSMC()
//****************************************************************************

/* The program is to find eigenfunctions and eigenvalues of 2-dimensional Schroedinger Eq
   in k-space.  Su dung ham DSYEVX Lappack

  INPUT:  + s: s-th section
          + v: v-th valley pairs
          + m1,m2: la hai tham so cho cac valley pairs. Vi du (m1=0.19, m2=0.916) 
            m1 ung voi truc y. m2 ung voi truc z
          + Cac tham so khac

  OUTPUT: + eigenvalues: eig(1,NSELECT)
          + eigenvetors: wave(1,NSELECT,1,nz,1,ny); da CHUAN HOA
 
  STARTING DATE:  January 08, 2010
  LASTEST UPDATE: Feb 23, 2010
*************************************************************************************** */
void Schrodinger2D(int s, int v, double m1, double m2){
     // Goi cac ham
     void get_potential(int s, int ny, int nz);
     void makeh(int ny,int nz,double ly,double lz,double *potential,double **ham,double m1, double m2);
     double *Get_pot();// Sau khi goi get_potential(int s) thi no se chuyen Phi(i_HANGSO,j,k) thanh 1D pot
     void diasym(int s, int v,int ny, int nz,double ly,double lz,int NSELECT, double **ham, double ***eig, double *****wave);
     void normalize_wave(int s,int v,int ny,int nz,int NSELECT, double *****wave);
     double ***Get_eig(), *****Get_wave();
     int Get_nya(),Get_nyb(), Get_nza(),Get_nzb();
     double Get_mesh_size_y(), Get_mesh_size_z();
     int Get_NSELECT();
     int Get_save_eig_wavefuntion(); // Flag chi ra co SAVE hay KHONG?
     
     // Cac bien local
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 
     int ny_max = nyb - nya;// De thay doi it nhat co the
     int nz_max = nzb - nza; //printf("\n nya=%d, nyb=%d,ny_max=%d, nza=%d,nzb=%d,nz_max=%d",nya,nyb,ny_max,nza,nzb,nz_max);
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
     int NSELECT = Get_NSELECT();
     
     double ***eig      = Get_eig();
     double *****wave   = Get_wave();
     int save_eig_wavefuntion = Get_save_eig_wavefuntion();
     
     // Thuc hien
     int ny = ny_max + 1; // Do so KHOANG CHIA la ny_max va nz_max NEN SO DIEM CHIA la ny_max+1 va nz_max+1 (May 30, 2010)
     int nz = nz_max + 1; // int ny = ny_max; int nz = nz_max;
     int i,j;
     double **ham; // theo ta xay dung "ham" trong ham Schrodinger la hop ly
     ham=dmatrix(1,ny*nz,1,ny*nz); // Hamiltonian di tu chi so 1 cho phu hop voi ham makeh.
     for(i=1; i<=ny*nz; i++){// Noi chung neu khoi tao array thi nen khoi tao gia tri dau cho no. KHONG se co RAT NHIEU LOI
       for(j=1; j<=ny*nz; j++){ // (June, 03, 2010)
	 ham[i][j] =0.0;
       }
     }

     //printf("\n Get potential"); 
     get_potential(s,ny,nz);// Get potential. Sau ham nay duoc array chua potential
     double *pot = Get_pot();// Vi sau ham get_potential() thi Get_pot moi thay doi
     
     //printf("\n Make Hamiltonian");
     double ly = mesh_size_y*ny_max/1.0e-9;// do ly la DOAN nen nhan voi ny_max la HOP LY  printf("\n ly =%f",ly);//NOTE: [nm] Kkong phai la m
     double lz = mesh_size_z*nz_max/1.0e-9;
     makeh(ny,nz,ly,lz,pot,ham,m1,m2);// To contruct the Hamiltonian

     //printf("\n Calling diasym");
     diasym(s,v,ny,nz,ly,lz,NSELECT,ham,eig,wave);// Call Lapack routine // da check thu roi, eig va wave o dau ra la ok

     //printf("\n Normalize wave");
     normalize_wave(s,v,ny_max,nz_max,NSELECT,wave);// Normalize wave

     /* Free variables */
     free_dmatrix(ham,1,ny*nz,1,ny*nz) ;

}// End of void Schrodinger2D(double m1, double m2)
//***********************************************************************************
/* *****************************************************************************     
  De lay potential Phi(i,j,k)
  Nhung do tai giai cho tung Cross-section nen cai ma ta lay la o dang Phi(0,j,k)
  roi Phi(1,j,k), ...., Phi(nx_max, j,k)
  
  INPUT: + s: s-th section
         + ny va nz: la so diem chia theo phuong y va z
         + Phi(i,j,k)
         
  OUTPUT: potential[] chuyen tu matrix sang dang array. potential[] la bien global
 
 Starting date: April 19,2010 
 Update:        April 19,2010 
 Latest update: May 09, 2010: Ta tinh theo dung cong thuc tinh potential
 Ta phai chu y toa do cua Potential ma ta dinh nghia.
 Ham nay la thay doi mang tinh CO BAN 
********************************************************************************* */
//#include <stdio.h>
//#include "nrutil.h"

void get_potential(int s, int ny, int nz) { // s-th section
   
     double ***Phi = GetPot();

     double  *Get_pot();//chuyen Phi(s-th,j,k) thanh pot 1D (potential 1D)
     double *pot   = Get_pot();
     
     double Get_Eg(); // Band gap
     double Eg = Get_Eg(); //printf("\n Eg = %f",Eg);// E->c_si = Eg/2.0
     //double dEc = 0; // Jun 02, 2010

     int Get_nya(),Get_nza();
     int nya = Get_nya(); // Ro rang toa do cua potential bi dich 1 doan la nya va nza
     int nza = Get_nza(); 
     
     int j_k=0,j,k;
     for(j=1; j<=ny; j++){// ny = nyb-nya+1; tuong tu cho nz
         for(k=1; k<=nz; k++){      
             j_k=j_k+1;
             pot[j_k]= Eg/2.0 - Phi[s][j+nya-1][k+nza-1];//potential o loi Silicon  //pot[j_k]= dEc - Phi[s][j+nya-1][k+nza-1];
            
	 }
     }
    
    return;
} // End of void get_potential(int s)
//********************************************O******************************************************

/* *****************************************************************************
  Tinh ma tran Hamiltonian
  
  INPUT: + ny va nz: la so diem chia theo phuong y va z 
	     + ly va lz: la do dai cua rectangular theo phuong y va z
	     + *potential la mang potential.CO THEM beta; Can sua lai cho hop ly voi potential cua ta
	     + m1,m2: la hai tham so cho cac valley pairs. Vi du (m1=0.19, m2=0.916)
	       m1 ung voi truc y. m2 ung voi truc z

  OUTPUT: + **ham: Hamiltonia matrix (1,ny*nz,1,ny*nz)

 Note: Su dung open boundary tuc la dang hard wall.

 Starting date: December 21, 2009
 Latest update: January 15, 2010
********************************************************************************/
//#include <stdio.h>
//#include <math.h>
//#include "constants.h"
//#include "nrutil.h"

void makeh(int ny,int nz,double ly,double lz,double *potential,double **ham,double m1, double m2)
  { 
    double energy0(int ky,int kz,double ly,double lz,double m1, double m2);
    void multiply_matrix(int INDEX,double alpha, double **matrixA, double **matrixB, double **matrixC);

    int nynz = ny*nz;
    double **U  = dmatrix(1,nynz,1,nynz); // Chua Uk
    double **VU = dmatrix(1,nynz,1,nynz); // Chua V x Uk
    double **Vk = dmatrix(1,nynz,1,nynz); // la gia tri cua hangso*U_t*VU
    int i,j;
    for(i=1; i<=ny*nz; i++){// Noi chung neu khoi tao array thi nen khoi tao gia tri dau cho no. KHONG se co RAT NHIEU LOI
       for(j=1; j<=ny*nz; j++){ // (June, 03, 2010)
	 U[i][j]=0.0; VU[i][j]=0.0; Vk[i][j]=0.0;
       }
     }

    int m,n,pp,qq; // do q da duoc dinh nghia o constant, p da dinh nghia la bien TOAN CUC roi
    double ym,zn;
    // Dau tien tinh thanh phan diagonal
    for(i=1; i<=nynz; i++){
      n=1+(i-1)%nz; // chinh la kz cua ta
      m=1+(i-1)/nz;// chinh la ky cua ta
      ham[i][i]= energy0(m,n,ly,lz,m1,m2);// thanh phan diagonal
    }
    // Tinh thanh phan Uk va VrUk duoc tao doi potential
    for(i=1; i<=nynz; i++){
      n=1+(i-1)%nz; 
      m=1+(i-1)/nz;
      for(j=1; j<=nynz; j++){
	   qq=1+(j-1)%nz;
	   pp=1+(j-1)/nz;
	   ym=(double)(m)/(double)(ny+1);
	   zn=(double)(n)/(double)(nz+1);
	   U[i][j]=sin(pp*pi*ym)*sin(qq*pi*zn); //chinh la Uk //printf(" U[%d][%d]= %le",i,j,U[i][j]);
	   VU[i][j]=beta*potential[i]*U[i][j]; //chinh la V*Uk;/printf("\n potential[%d] = %f",i,potential[i]);//check potential
       }
     }
 
    //Thuc hien viec nhan
    double alpha = 4.0/((double)((ny+1)*(nz+1)));
    multiply_matrix(nynz,alpha,U,VU,Vk);// dau ra Vk chinh la thanh phan Hamiltonian do potential tao ra 
  
    // Cong 2 thanh phan tao nen ham(i,j) lai voi nhau
    for(i=1; i<=nynz; i++){
       for(j=1; j<=nynz; j++){
	 ham[i][j] = ham[i][j]+Vk[i][j];
       }
    }
    
    // Free_local vector
    free_dmatrix(U,1,nynz,1,nynz); free_dmatrix(VU,1,nynz,1,nynz); free_dmatrix(Vk,1,nynz,1,nynz);

} // End of void makeh(int ny,int nz,double ly,double lz,double *potential,double **ham,double m1, double m2)
//*********************************************************************************************
/* *****************************************************************************
  De tinh gia tri betaEkykz

  INPUT: + ky va kz la chia theo k-space, ky di tu 1 den ny, kz di tu 1 den nz
         + ly va lz: la do dai cua rectangular theo phuong y va z
	     + m1,m2: la hai tham so cho cac valley pairs. Vi du (m1=0.19, m2=0.916)
	       m1 ung voi truc y. m2 ung voi truc z
 
   OUTPUT:+ Gia tri energy E0 ma ta tim duoc. Don vi beta*E

 Note: Su dung open boundary tuc la dang hard wall.
  
 Starting date: December 28, 2009
 Latest update: December 28, 2009
*********************************************************************************/
//#include <stdio.h>
//#include "constants.h"

double energy0(int ky,int kz,double ly,double lz,double m1, double m2){
	
    double pi2  = (pi*pi)/2.0;
    double kyly = (double)(ky)/ly;
    double kzlz = (double)(kz)/lz;
    double energy;
    energy = pi2*((1.0/m1)*kyly*kyly +(1.0/m2)* kzlz*kzlz);
    //printf("\n Energy0 = %le",energy); 
    return energy;
	
} // End of energy0(int ky,int kz,double ly,double lz,double m1, double m2)

/* *****************************************************************************
  Nhan 2 ma tran su dung dgemm routine tu Lapack 
  
  INPUT: + Ma tran 2 chieu matrixA va matrixB 
         + INDEX la chi so mang day la ma tran vuong
	     + alpha: la so vo huong

  OUTPUT: + Ma tra 2 chieu matrixC =  alpha*matrixA_t*matrixB (t la transpose)+beta*matrixC nhung beta=0.0

 Starting date: January 11, 2010
 Latest update: January 15, 2010
********************************************************************************/

void multiply_matrix(int INDEX,double alpha, double **matrixA, double **matrixB, double **matrixC){

   void convert_mtx_to_array(int M, int N, double **mtx, double *a);
   void convert_array_to_mtx(int M, int N, double *a, double **mt);
   // **matrixA=dmatrix(1,INDEX,1,INDEX),**matrixB=dmatrix(1,INDEX,1,INDEX), **matrixC=dmatrix(1,INDEX,1,INDEX);   
   double *A=dvector(1,INDEX*INDEX),*B=dvector(1,INDEX*INDEX),*C=dvector(1,INDEX*INDEX);
   int i,j,n;
   for(i=1; i<=INDEX*INDEX; i++){// Noi chung neu khoi tao array thi nen khoi tao gia tri dau cho no. KHONG se co RAT NHIEU LOI
     A[i]=0.0; B[i]=0.0; C[i]=0.0;// (June, 03, 2010)
   }

   // Chuyen ma tran A thanh array A de phu hop voi ham dgemm, tuong tu cho ma tran B
   convert_mtx_to_array(INDEX,INDEX,matrixA,A);   
   convert_mtx_to_array(INDEX,INDEX,matrixB,B);   
   // Ket qua duoc mang C 
 
   // Matrix-Matrix Multiply       
   n = INDEX; 
   double beta_beta=0.0; // do beta la hang so dinh nghia o constant.h
   //printf("\n Multiply matrix");
   dgemm("t","n",&n,&n,&n,&alpha,A,&n,B,&n,&beta_beta,C,&n);
   // "t' dau nghia la A la transpose, "n" sau nghia la B la ok
   // Chuyen ket qua tu array sang ma tran
   convert_array_to_mtx(INDEX,INDEX,C,matrixC);

   //Free-local vector
   free_dvector(A,1,INDEX*INDEX); free_dvector(B,1,INDEX*INDEX); free_dvector(C,1,INDEX*INDEX);
   return;
}// End of void multiply_matrix(double **A, double **B, double **C)
//***********************************************************************************
void convert_mtx_to_array(int M, int N, double **mtx, double *a)
{
  int i,j ;
  int k = 0 ;

  for ( j=1 ; j<=N ; j++ )
    for ( i=1 ; i<=M ; i++ )
      a[k++] = mtx[i][j] ;
}

void convert_array_to_mtx(int M, int N, double *a, double **mtx)
{
  int i,j ;
  int k = 0 ;

  for ( j=1 ; j<=N ; j++ )
    for ( i=1 ; i<=M ; i++ )
      mtx[i][j] = a[k++];
}
//****************************************************************************************
/****************************************************************************************
   - Su dung ham DSYEVX Lappack
  INPUT: + s: s-th section
         + v: v-th valley pair
         + ny va nz: la so diem chia theo phuong y va z
         + ly va lz: la do dai cua rectangular theo phuong y va z
         + Nselect: so subband ta can lay
         + Hamitonian matrix ham[n][n]

  OUTPUT: + eigenvalues: eig(1,NSELECT)
          + eigenvetors: wave(1,NSELECT,1,ny,1,nz);

  STARTING DATE:  January 04, 2010
  LASTEST UPDATE: Feb 24, 2010

*************************************************************************************** */
//#include <stdio.h>
//#include <string.h>
//#include "nrutil.h"

void diasym(int s, int v,int ny, int nz,double ly,double lz,int NSELECT, double **ham, double ***eig, double *****wave){
     
     double psirealspace(double y,double z,int ny,int nz,double ly,double lz,double **vec, int index);
     void convert_mtx_to_array(int M, int N, double **mtx, double *a);
     void print_matrix( char* desc, int m, int n, double* a, int lda );
     
     int i,j;
     //  SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,IFAIL, INFO )
     int N = ny*nz; // la bien trung gian de tranh bi thay doi khi vao ham dsyevx
     int LDA = N, LDZ = N; //  NSELECT = 10 (for example);// global variable khai bao tu ham main
     /* Locals */
     int n = N, il, iu, m, lda = LDA, ldz = LDZ, info, lwork;
     double abstol, vl, vu;
     double wkopt;
     double* work;
     /* Local arrays */
     int *iwork = ivector(0,5*N-1); for(i=0; i<=5*N-1; i++){ iwork[i] = 0;}      // iwork dimension should be at least 5*n 
     int *ifail = ivector(0,N-1); for(i=0; i<=N-1; i++){ ifail[i] = 0; }
     double *w = dvector(0,N-1); for(i=0; i<=N-1; i++){ w[i]=0.0; }  //chinh la eigvalues ta tim duoc
     double *z = dvector(0,LDZ*NSELECT-1); for(i=0; i<= LDZ*NSELECT-1; i++){ z[i]=0.0;}   // chua eigenvectors ma ta can tim
     double *a = dvector(0,LDA*N-1); for(i=0; i<=LDA*N-1; i++){ a[i]=0.0;}// chinh la array de chua ma tran Hamiltonian

     convert_mtx_to_array(LDA,N,ham,a);   // do LDA=N=ny*nz
     
     /* Executable statements */
     //printf( "\n Running DSYEVX " );
     abstol = -1.0; /* negative abstol means using the default value */
     il = 1; /* Set il, iu to compute NSELECT smallest eigenvalues */
     iu = NSELECT;
     // Query and allocate the optimal workspace 
     lwork = -1;
     dsyevx("Vectors","Indices","Upper", &n, a, &lda, &vl, &vu, &il, &iu,&abstol, &m, w, z, &ldz, &wkopt, &lwork, iwork, ifail, &info);
     lwork = (int)wkopt; // Chu y goi dsyevx o tren la &wkopt day nhe !
     work = dvector(0,lwork-1);
     /* Solve eigenproblem */
     dsyevx("Vectors","Indices","Upper", &n, a, &lda, &vl, &vu, &il, &iu,&abstol, &m, w, z, &ldz,  work,  &lwork, iwork, ifail, &info );
     /* Check for convergence */
     if( info > 0 ) {
         printf( "The algorithm dsyevx failed to compute eigenvalues.\n" );
         exit( 1 );
       }
     /* Print the number of eigenvalues found */
     //printf( "\n Number of subbands for each valley:  %2i", m );// m=iu-il+1 = NSELECT
     // Print eigenvalues //  print_matrix( "Selected eigenvalues from DSYEVX", 1, m, w, 1 );
     // Tinh eigen_value tuong ung voi so NSELECT ma ta chon
     int k=0; // w chinh la mang chua eigen value
     for( i = 0; i <1; i++ ) {// chay theo 0. 
	    for( j = 0; j < m; j++ ) { // m chinh la so NSELECT
	       k=k+1;
	       eig[s][v][k]= w[i+j*1];//eigvalue o day la tu DSYEVX ma ra. Ta se di tu chi so 1 cho lowest eigenvalue
	       eig[s][v][k]=eig[s][v][k]/beta; // [eV]  theo cong thuc tinh toan Schrodinger cua ta
	   }
      }

      // Tinh eigenvectors. //print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
      double **vec=dmatrix(1,N,1,NSELECT);  // n cung chinh bang N 
      for( i = 0; i < N; i++ ) {// di tu 1 den ny*nz
	  for( j = 0; j < m; j++ ) {//m = NSELECT so luong tran thai state can luu. 
	    vec[i+1][j+1]= z[i+j*ldz] ; // vec lay cot cua mang z; do chi so mang vec la di tu 1
	  }
	}
      // To store state psi to wave array de su dung tinh scattering
      int index;
      double y_point, z_point; // Gia tri de xac dinh so diem chia thi khong quan trong nhung no can matching voi index o wave
      int ny_max = ny -1; // (May 30, 2010)
      int nz_max = nz -1;// do ny=nyb-nya+1; nen ny_max moi la index chinh xac o wave, no di tu 0 den ny_max
        for(i=0;i<=ny_max;i++){ // so diem chia theo phuong y. vi du chia thanh 100 diem
	       y_point=ly*(double)(i)/(double)(ny_max); // diem chia theo y di tu 0 den ly
	       for(j=0;j<=nz_max;j++){ // so diem chia theo phuong z. vi du chia thanh 100 diem
	          z_point=lz*(double)(j)/(double)(nz_max); // diem chia theo z di tu 0 den lz
	          for(index=1; index<=NSELECT; index++){
	             wave[s][v][index][i][j]= psirealspace( y_point,z_point,ny,nz,ly,lz,vec,index);
	   }// End of for(index=1; index<=NSELECT; index++)
	 } // End of for(i=0;i<=ny_max;i++){ 
     }// End of for(j=0;j<=nz_max;j++){ 
    
     // Free local variables
     free_dvector(work,0,lwork-1);
     free_ivector(iwork,0,5*N-1), 
     free_ivector(ifail,0,N-1);
     free_dvector(w,0,N-1);//chinh la eigvalues ta tim duoc
     free_dvector(z,0,LDZ*NSELECT-1);// chua eigenvectors ma ta can tim
     free_dvector(a,0,LDA*N-1);      // chinh la array de chua ma tran Hamiltonian
     free_dmatrix(vec,1,N,1,NSELECT);

     //printf("\n Ket thuc ham diasym");
 } // End of void diasym(double **ham, double *eig, double ***wave)

/**************************************************************************************** */
/* *****************************************************************************
  To calculate the wavefunction in real space
   psi = Sum(Ck*|k>)
   
  INPUT: + y,z la diem vao de tinh gia tri psi
  		   Hai bien nay chi can luc ve psi thoi. Thuc ra lay ny va nz de ve cung 
  		   ok. Nhung co the ta can nhieu diem chia hon
         + ny va nz: la so diem chia theo phuong y va z
         + ly va lz: la do dai cua rectanggular theo phuong y va z
	 + vec[][] la he so Ck day. Chu y he so thek la tong hop cua kx va ky
	 + index: la chi so mang thu 2 cua vec de chi viec ta dinh lay psi ung voi eigen thu may
	          nen index se di tu 1 den NSELECT
		 
  OUTPUT:+ psirealspace: gia tri psi trong real space

 Note: Su dung open boundary tuc la dang hard wall.

 Starting date: December 29, 2009
 Latest update: December 31, 2009
**************************************************************************************************/
//#include <stdio.h>
//#include <math.h>
//#include "constants.h"

double psirealspace(double y,double z,int ny,int nz,double ly,double lz,double **vec, int index){
      int i,j,ky,kz;
      double psi,phiy,phiz;
      psi=0.0;
      for(i=1; i<=ny*nz; i++){// la gia tri k chay tu 1 den N=ny*nz
	  //printf( "\n %6.8f",vec[i] );
	  kz = 1+(i-1)%nz;    // chay theo phuong ky
	  ky = 1+(i-1)/nz;    // roi chay theo phuong kz
	  phiy = sin((double)(ky)*pi*y/ly); // la ham sin theo ky
	  phiz = sin((double)(kz)*pi*z/lz);   // la ham sin theo kz
	  psi = psi + vec[i][index]*phiy*phiz; // Xem cong thuc tren
      } // End of for
   
      psi = 2.0*psi/(sqrt(ly*lz)); // Tai sao thieu so 2 nhi ?
      //psi = psi/(sqrt(ly*lz));// Thie so 2 thi co nghia chua normalized. Sau do se normalized. Nhung nhan 2 o day van hay hon
      //printf( "\n psi = %6.8f",psi );

   return psi;
} // End of double psirealspace(double y,double z,int ny,int nz,double ly,double lz,vec[]){


// Chon kieu nay thi n chinh la NSELECT
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
  int i, j;
        printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
	  for( j = 0; j < n; j++ ) {
	    printf( " %6.8f", a[i+j*lda] );
	    }
          printf( "\n" );
        }
}
//*******************************************************************************
/****************************************************************************************
   - To normalize wave function su dung trapezoidal
  INPUT: + s: s-th section
         + v: v-th valley pairs
         + ny,nz la so diem chia tren 2 truc ly va lz
         + Nselect: so subband ta can lay
         + trapezoidal_weights[][] Truoc (April 08, 2010) su dung cach don gian hon
         + wave(1,NSELECT,0,ny,0,nz) CHUA normalize
         
  OUTPUT: +  wave(1,NSELECT,0,ny,0,nz) DA DUOC normalize
 
  STARTING DATE:  April 08, 2010
  LASTEST UPDATE: April 08, 2010
*************************************************************************************** */
void normalize_wave(int s,int v,int ny,int nz,int NSELECT, double *****wave){
  // Goi ham
  double **Get_trap_weights();// Lay phan cong thuc 1 va hang so cua cong thuc 2 trang 92
  // Local bien
  double **trap_weights = Get_trap_weights();
  
// Thuc hien
  double mag_psi;
  int index,i,j;
  for(index=1; index<=NSELECT; index++){// Chay cho tung subband
    mag_psi = 0.0;
    for(i=0; i<=ny; i++){
      for(j=0; j<=nz; j++){ 
     	  mag_psi += trap_weights[i][j]*wave[s][v][index][i][j]* wave[s][v][index][i][j];
      }
    }
    mag_psi = 1.0/sqrt(mag_psi);
    //printf("\n Trapzoidal Tham so A can tim de normalization = %le", mag_psi);
     for(i=0; i<=ny; i++){
      for(j=0; j<=nz; j++){ 
   	     wave[s][v][index][i][j]=  mag_psi*wave[s][v][index][i][j];// Normalize wave
      }
    }
   
  }// End of for( index=1; index<=NSELECT; index++){// Chay cho tung subband

}// End of void normalize_wave(){
//*************************************************************************************

/****************************************************************************************
  INPUT: + s: s-th section
         + v: v-th valley pair
         + ny va nz: la so diem chia theo phuong y va z
         + ly va lz: la do dai cua rectangular theo phuong y va z
         + Nselect: so subband ta can lay
         + eigenvalues: eig(1,NSELECT)
         + eigenvetors: wave(1,NSELECT,1,ny,1,nz);

  OUTPUT: + eig.dat chua eigen energy
          + wave.dat chua psi da duoc CHUAN HOA

  STARTING DATE:  January 08, 2010
  LASTEST UPDATE: Feb 24, 2010

*************************************************************************************** */
//#include <stdio.h>
//#include "nrutil.h"

void save_eig_wave(int s,int v, char *fneigen, char *fnwave){

     double ***Get_eig(), *****Get_wave();
     double ***eig      = Get_eig();
     double *****wave   = Get_wave();
   
     int Get_NSELECT();
     int NSELECT = Get_NSELECT();
      
     int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 
     int ny = nyb - nya; // int ny = ny_max = nyb - nya;
     int nz = nzb - nza; // int nz = nz_max = nzb - nza; 

     double Get_mesh_size_y(), Get_mesh_size_z();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();
     double ly = mesh_size_y*ny/1.0e-9;//Do ly la DOAN  printf("\n ly =%f",ly);//NOTE: [nm] Kkong phai la m
     double lz = mesh_size_z*nz/1.0e-9;
    
// Thuc hien
      int index,i,j;
      // To store eigen value
      FILE *f; f = fopen(fneigen,"a"); // "w" neu tep ton tai no se bi xoa
      if(f == NULL) { printf("\n Cannot open file at save_eig_wave.c"); return ; }
      fprintf(f,"\n # s_th v_th i_subband  eigen_value[eV] \n");
      for( i=1; i<=NSELECT; i++){ // do m o ham dsyevx tu diasym la bang NSELECT
	  fprintf(f,"%d      %d       %d         %le \n",s,v,i, eig[s][v][i]);
	}                                      
      fclose(f);
	
      // To store wave (psi)
      FILE *f1; 
      f1 = fopen(fnwave,"a"); if(f1== NULL) { printf("\n Cannot open file save_eig_wave.c"); return;}
      fprintf(f1,"\n # s-th v_th    y_point    z_point    psi_square \n");
      double y_point, z_point;
 // Chi in ra file wave khi NSELECT = 10
 if(NSELECT ==10){
      for(index=1; index<=NSELECT; index++){// chay thu tu tung wave
          for(i=0;i<=ny;i++){ // so diem chia theo phuong y. 
	         y_point=ly*(double)(i)/(double)(ny); // diem chia theo y di tu 0 den ly
	         for(j=0;j<=nz;j++){ // so diem chia theo phuong z. 
	            z_point=lz*(double)(j)/(double)(nz); // diem chia theo z di tu 0 den lz
	            fprintf(f1,"%d     %d     %le      %le      %le       %le %le %le %le %le %le %le %le %le \n",
		              s,v,z_point,y_point,
			      wave[s][v][1][i][j]*wave[s][v][1][i][j],
		              wave[s][v][2][i][j]*wave[s][v][2][i][j],
		              wave[s][v][3][i][j]*wave[s][v][3][i][j],
		              wave[s][v][4][i][j]*wave[s][v][4][i][j],
		              wave[s][v][5][i][j]*wave[s][v][5][i][j],
		              wave[s][v][6][i][j]*wave[s][v][6][i][j],
		              wave[s][v][7][i][j]*wave[s][v][7][i][j],
		              wave[s][v][8][i][j]*wave[s][v][8][i][j],
		              wave[s][v][9][i][j]*wave[s][v][9][i][j],
		              wave[s][v][10][i][j]*wave[s][v][10][i][j]);	// Binh phuong wave	
	     
	        }// End of for(j=0;j<=nz;j++){ 
         } // End of for(i=0;i<=ny;i++){ 
      }//  End of for(index=1; index<=NSELECT; index++)
}// End of  if(NSELECT ==10)
     fclose(f1);
}// End of void save_eig_wave( double *eig, double ***wave){
  
/****************************************************************************************
  INPUT: + nx_max: so diem theo phuong x
         + ny va nz: la so diem chia theo phuong y va z
         + ly va lz: la do dai cua rectangular theo phuong y va z
         + Nselect: so subband ta can lay
         + eigenvalues: eig(1,NSELECT)
 
  OUTPUT: + subband.dat chua eigen energy

  STARTING DATE:  May 26, 2010
  LASTEST UPDATE: May 27, 2010

*************************************************************************************** */
void save_subband_energy(char *fnsubband){// con cho theory thi can gi doi ten

  double ***Get_eig();
  double ***eig      = Get_eig();
    
  int Get_NSELECT();
  int NSELECT = Get_NSELECT();

  int Get_nx0(),Get_nx1();
  int nx0 = Get_nx0(); int nx1 = Get_nx1();
  int nx_max = nx1-nx0;
  
  int s; // E_s,v,i
      
  if(NSELECT >= 4){// Chi save khi so subband lon hon 4 tro len
     FILE *f; f = fopen(fnsubband,"w"); // "w" neu tep ton tai no se bi xoa
     if(f == NULL) { printf("\n Cannot open file subband.dat"); return ; }
     fprintf(f,"\n # s_th Es_1_1 E_s_1_2 E_s_1_3 E_s_1_4 Es_2_1 E_s_2_2 E_s_2_3 E_s_2_4 Es_3_1 E_s_3_2 E_s_3_3 E_s_3_4 \n");
     
     double Get_mesh_size_x(); double mesh_size_x = Get_mesh_size_x();
     for( s=0; s<=nx_max; s++){ 
	fprintf(f,"%le %le %le %le %le %le %le %le %le %le %le %le %le\n", (double)(s)*mesh_size_x/1.0e-9, 
		     eig[s][1][1], eig[s][1][2], eig[s][1][3], eig[s][1][4],
		     eig[s][2][1], eig[s][2][2], eig[s][2][3], eig[s][2][4],
		     eig[s][3][1], eig[s][3][2], eig[s][3][3], eig[s][3][4]);
     }// End of for( s=0; s<=nx_max; s++)                                      
     fclose(f);
  }// End of if(NSELECT >= 4){
  
// Tinh energy spacing theo ly thuyet
  int Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
  int nya = Get_nya(); int nyb = Get_nyb(); 
  int nza = Get_nza(); int nzb = Get_nzb();
  
  double Get_mesh_size_y(), Get_mesh_size_z();
  double mesh_size_y = Get_mesh_size_y(); 
  double mesh_size_z = Get_mesh_size_z();

  double Tsi = (double)(nyb - nya)*mesh_size_z; // chu y la theo phuong  nhe
  double Wsi = (double)(nzb - nza)*mesh_size_y;

  double Get_ml(), Get_mt();
  double ml = Get_ml();
  double mt = Get_mt();
  
  int y_interger = 5 , z_interger = 5; // integer number for y and z direction. Chon tuy y  // ta chon bang 5 la du roi
  double my, mz, **E1, **E2, **E3;
  E1=dmatrix(1,y_interger,1,z_interger), E2=dmatrix(1,y_interger,1,z_interger), E3=dmatrix(1,y_interger,1,z_interger);
  int i,j;
  for (i=1; i<=y_interger; i++){
    for(j=1; j<=z_interger; j++){
      E1[i][j]=0.0; E2[i][j]=0.0; E3[i][j]=0.0;
    }
  }
  
  //Valley pair 1: my=mt, mz=mt
  my = mt, mz=mt;
  for (i=1; i<=y_interger; i++){
    for(j=1; j<=z_interger; j++){
      E1[i][j] = ((double)(i*i)*pi*pi*hbar*hbar)/(2.0*my*m0*Tsi*Tsi) + ((double)(j*j)*pi*pi*hbar*hbar)/(2.0*mz*m0*Wsi*Wsi);
    }
  }

  // Valley pair 2:
  my = ml, mz = mt;
  for (i=1; i<=y_interger; i++){
    for(j=1; j<=z_interger; j++){
      E2[i][j] = ((double)(i*i)*pi*pi*hbar*hbar)/(2.0*my*m0*Tsi*Tsi) + ((double)(j*j)*pi*pi*hbar*hbar)/(2.0*mz*m0*Wsi*Wsi);
    }
  }
  // Valley pair 3
  my = mt, mz = ml;  
  for(j=1; j<=z_interger; j++){// Muon chung to giong valley 2 khi doi xung thi can nhu the nay
    for (i=1; i<=y_interger; i++){
    //for(j=1; j<=z_interger; j++){
      E3[j][i] = ((double)(i*i)*pi*pi*hbar*hbar)/(2.0*my*m0*Tsi*Tsi) + ((double)(j*j)*pi*pi*hbar*hbar)/(2.0*mz*m0*Wsi*Wsi);
    }
  }
  
  FILE *f1; f1 = fopen("subband_spacing_theory.dat","w"); // "w" neu tep ton tai no se bi xoa
  if(f1 == NULL) { printf("\n Cannot open file subband__spacing_theory.dat"); return ; }
  fprintf(f1,"\n # Valley1_spacing[ev] Valley2_spacing[eV] Valley3_spacing[eV] NOT shown in the ascending order  \n");
  fprintf(f1,"%f %f %f \n",(E1[1][2]-E1[1][1])/q,(E2[1][2]-E2[1][1])/q,(E3[1][2]-E3[1][1])/q); // first spacing if my=mz
  fprintf(f1,"%f %f %f \n",(E1[2][2]-E1[1][2])/q,(E2[2][2]-E2[1][2])/q,(E3[2][2]-E3[1][2])/q); // second spacing if my=mz
  fprintf(f1,"%f %f %f \n",(E1[1][3]-E1[2][2])/q,(E2[1][3]-E2[2][2])/q,(E3[1][3]-E3[2][2])/q); // thirsd spacing if my=mz
  fprintf(f1,"%f %f %f \n",(E1[2][3]-E1[1][3])/q,(E2[2][3]-E2[1][3])/q,(E3[2][3]-E3[1][3])/q); // forth spacing if my=mz
  fprintf(f1,"%f %f %f \n",(E1[1][4]-E1[2][3])/q,(E2[1][4]-E2[2][3])/q,(E3[1][4]-E3[2][3])/q); // 5th spacing if my=mz
  fclose(f1);

  FILE *f2; f2 = fopen("subband_theory.dat","w"); // "w" neu tep ton tai no se bi xoa
  if(f2 == NULL) { printf("\n Cannot open file subband__theory.dat"); return ; }
  fprintf(f2,"\n # i j  Valley1[ev] Valley2[eV] Valley3[eV]  \n");
  for (i=1; i<=y_interger; i++){
    for(j=1; j<=z_interger; j++){
      fprintf(f2," %d %d %f %f %f \n",i,j, E1[i][j]/q, E2[i][j]/q, E3[i][j]/q); 
    }
  }
  fclose(f2);
  
  // Free local
  free_dmatrix(E1,1,y_interger,1,z_interger), free_dmatrix(E2,1,y_interger,1,z_interger), free_dmatrix(E3,1,y_interger,1,z_interger);
  return;

}// End of void save_subband_energy(char *fnsubband)
//**********************************************************************************************

/* *************************************************************************************************  
  Save tai midpoint cua cac phuong x or ( y va z)
  Chi save o phan loi Silicon

  Starting date: May 25, 2010         
  Latest update: May 27, 2010:
************************************************************************************************* */
void save_potential_used_for_2D_Schrodinger(char *fn_x, char *fn_yz){//For checking only (May 25, 2010)
     
     double ***Phi = GetPot();

     double Get_Eg(); // Band gap
     double Eg = Get_Eg(); 
     //double dEc = 0.0;

     int Get_nx0(),Get_nx1(),Get_nya(),Get_nyb(), Get_nza(),Get_nzb();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0;
    
     int nya = Get_nya(); int nyb = Get_nyb(); 
     int nza = Get_nza(); int nzb = Get_nzb(); 
     int ny = nyb - nya +1;// (May 30, 2010) // int ny = nyb - nya;
     int nz = nzb - nza +1; // int nz = nzb - nza; 

     double Get_mesh_size_x(),Get_mesh_size_y(),Get_mesh_size_z();
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     int j,k;
     double x_pos = 0.0, y_pos = 0.0, z_pos = 0.0;
     int midpoint_x = (int)(nx_max/2.0 + 0.5);
     int s = midpoint_x;
     FILE *ff;
     ff = fopen(fn_yz,"w");
     //ff = fopen("potential_yz_at_midpoint_x_in2DSchrodinger.dat","w");
     if(ff==NULL){ printf("\n Cannot open file at save_potential_used_for_2D_Schrodinger.c ");return; }
     fprintf(ff,"\n # y[nm] z[nm] pot_yz[eV] \n");
     for(j=1; j<=ny; j++){
         y_pos = (double)(j-1)*mesh_size_y/1.0e-9;// [m] -> [nm]; vi tri se di tu 0 den nyb-nya
         for(k=1; k<=nz; k++){  
	    z_pos = (double)(k-1)*mesh_size_z/1.0e-9;
            fprintf(ff," %le %le %le \n",y_pos,z_pos,( Eg/2.0 - Phi[s][j+nya-1][k+nza-1]));// potential[eV] 
	    //fprintf(ff," %le %le %le \n",y_pos,z_pos,( dEc - Phi[s][j+nya-1][k+nza-1]));// potential[eV] 
         }
     }// End of for(j=1; j<=ny; j++){
     fclose(ff);
     
   // Ve theo phuong x khi dat y va z rai midpoint (May 26, 2010)
   int Get_ny0(),Get_ny1(),Get_nz0(),Get_nz1();
   int ny0 = Get_ny0(); int ny1 = Get_ny1();
   int nz0 = Get_nz0(); int nz1 = Get_nz1();
   int midpoint_y = (int)((ny1-ny0)/2.0 + 0.5);// lay o giua truc y do tinh doi xung cua oxide
   int midpoint_z = (int)((nz1-nz0)/2.0 + 0.5);
   FILE *f1; 
   f1 = fopen(fn_x,"w");
   //f1 = fopen("potential_x_at_midpoint_yz_in2DSchrodinger.dat","w"); // "w" neu tep ton tai no se bi xoa
   if(f1==NULL){ printf("\n Cannot open file at save_potential_used_for_2D_Schrodinger.c"); return ; }
   fprintf(f1,"\n # x[nm] potential_x[eV] \n");
   int i;
   for(i=nx0; i<=nx1; i++){
     fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,Eg/2.0 - Phi[i][midpoint_y][midpoint_z]);
     // fprintf(f1,"%le %le \n",i*mesh_size_x/1.0e-9,dEc - Phi[i][midpoint_y][midpoint_z]);
   }
   fclose(f1);
   return;
}
