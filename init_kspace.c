/* ***********************************************************************
 To INITIALIZE CARRIER ENERGY AND WAVEVECTOR ACCORDING TO THE
                                 MAXWELL-BOLTZMANN STATISTICS
   - Initial electron energy using thermal energy in 1D do o day ta dung 1D transport. Day la DONG NANG thoi
   - Initial momentum of electrons
   - Initial valley and subband where the particle is residing for the first time
    
   Matrix  p[][] (parameters for particle)
           p[n][1] = kx (is stored kx)
           p[n][2] = x
           p[n][3] = is to store ts (ts or tc is free flight time)
           p[n][4] = i_region (which region we are considering?)
             
   Array  valley[]: hat thu n-th nam o vallye pair nao 
          valley[n] =1 : valley pair 1
          valley[n] =2 : valley pair 2 
          valley[n] =3 : valley pair 3 
          
   Matran subband[n]: de chua subband index cua hat thu n-th 
   
   Array energy[n]: electron energy of the n-th particle
   
NOTE: ta hoan toan co the de cac gia tri valley, subband, energy cua n-th particle vao mang 
      p[][] luon. Nhung tach ra nhu ta de cho no truc quan hon
   
Starting date: Feb 21, 2010
Update:        Feb 22, 2010
Update:        March 29, 2010. Hien tai chua dung i_region va neu dung thi CAU HOI
       se la co CAN THIET update ca vi tri theo phuong y va z khong?
Latest update: May 06, 2010 (De phu hop voi 3D Poisson)
************************************************************************** */
#include <math.h>
#include <stdio.h>
#include "constants.h"
   
void init_kspace(int ne,int i){// March 29, 2010
     //void init_kspace(int ne,int i,int j,int k){// ne-th electron to inilitialize at cell(i,j,k)
     // Goi cac ham
     double **Get_p(),*Get_energy();
     int *Get_valley(),*Get_subband();
     double Get_Vt(),Get_ml(),Get_mt(),Get_nonparabolicity_factor();
     long Get_idum();
     float random2(long *idum);
     int Get_nx0(),Get_nx1();
          
     // Cac bien local
     double **p     = Get_p();
     double *energy = Get_energy();
     int *valley    = Get_valley();
     int *subband   = Get_subband();
     double Vt = Get_Vt();
     double ml = Get_ml();
     double mt = Get_mt();
     double af = Get_nonparabolicity_factor();
     long idum = Get_idum();
     int nx0 = Get_nx0();
     int nx1 = Get_nx1(); 
          
     // Thuc hien 
     double k_momen = 0.0,kx = 0.0, electron_energy = 0.0; // k_momen la bien cho momentum TRANH nham voi k la bien chay cho truc z
     int iv = 0, subband_index = 0; // iv: valley index; i_region = 0
      
     // Initial particle energy (using Boltzman statistics)
     double rr = 0.0;
     do {
          rr=random2(&idum);
         }
     while ((rr<=1.0e-6)||(rr>1.0));
    
     // Day chi la DONG NANG thoi
     electron_energy = -(0.5*Vt)*log(rr); //[eV]// printf("\n electron energy =%f",electron_energy); // electron energy in 1D
     
     //Khoi tao valley index - Tu do khoi tao momentum va subband tuong ung trong valley do 
	 rr = 3.0*random2(&idum); //iv valley index luc initial thi chon theo random number	
     if(rr <=1.0){ 
           iv=1; //valley index =1 <- valley pair 1, m*=ml
           // Initial wavevector
           k_momen = sqrt(2.0*m0*ml)*sqrt(q)/hbar*sqrt(electron_energy*(1+af*electron_energy)); // Don vi chuan
                     // Nhan manh: sqrt(q) la de chuyen energy sang don vi J. Nhung de ra ngoai vi cong thuc
                     // sqrt(electron_energy(1+af*electron_energy))thi electron_energy phai la eV va af[1/eV]
	             // NHAN MANH: electron_energy o day la DONG NANG
           kx      = k_momen; // 1D transport
           // Khoi tao subband trong valley do
           subband_index = 1; // assumed o subband DAU TIEN
      } 
     else if(rr <=2.0){
           iv=2; //valley index =2 <- valley pair 2, m*=mt
           // Initial wavevector
           k_momen = sqrt(2.0*m0*mt)*sqrt(q)/hbar*sqrt(electron_energy*(1+af*electron_energy));
           kx      = k_momen;
           // Khoi tao subband trong valley do
           subband_index = 1; // assumed o subband DAU TIEN
        } 
     else  { // if(rr <=3.0) ; 
           iv=3; //valley index =3 <- valley pair 3, m*=mt
           // Initial wavevector
           k_momen = sqrt(2.0*m0*mt)*sqrt(q)/hbar*sqrt(electron_energy*(1+af*electron_energy));
           kx      = k_momen;
           // Khoi tao subband trong valley do
           subband_index = 1; // assumed o subband DAU TIEN
        } 
      
     // Khoi tao region dua vao toa do (i,j,k) cua hat
     //i_region = find_region(i,j,k); // tim duoc hat dang o region nao: 1:S, 2:D, 3:Channel 

     // Check boundaries
     if((i==nx0)&&(kx<0.0))     { kx = -kx; }     //reflection
     if((i==nx1)&&(kx>0.0))     { kx = -kx; }

     // Vi chua co bang scattering nen chua the Initial free-flight
      
     // Map particle atributes
     p[ne][1]    = kx;  // Don vi chuan            
     // p[ne][2] = x; // in init_realspace() function
     // p[ne][3] = ts or tc; //in init_free_flight() phai sau khi tinh bang scattering
     //p[ne][4]    = i_region; //29/03/10 15:44. Hien tai chua dung
     valley[ne]  = iv;
     subband[ne] = subband_index;
     energy[ne]  = electron_energy;// [eV] electron energy of particle neth
	 //printf("\n Ket qua: kx=%le region=%d valley=%d subband=%d energy=%le",kx, i_region,iv,subband_index,electron_energy);
     //printf("\n Ket qua: kx=%le region=%d valley=%d subband=%d energy=%le",
                        //p[ne][1],(int)(p[ne][4]),valley[ne],subband[ne],energy[ne] );
     //printf("\n Tu init k-space kx=%le",kx);
     //getchar();
     
     return;
}

