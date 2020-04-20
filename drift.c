/* *****************************************************************************
        PERFORM THE K-SPACE AND REAL-SPACE MOTION OF THE CARRIERS
This function is to calculate the momentums and positions of the particle
which drifts during time (tau)
   + F tinh theo gradient cua eigen energy
   + The change of the momentum under electric field F in the x-direction is
                dkx = [(-e*F)/hbar]*tau  
   + The new position is obtained from
         x(at t+tau)=x(at t)+(hbar/2m*)*[kx(at t)+kx(at t+tau)]*tau
                    = x(at t)+(hbar/m*)*[kx(at t)+0.5*dkx]*tau
          
   + NOTE: For spherical non-parabolic conduction band
           E(1+alpha*E)=(hbar*hbar)*(k*k)/(2m*) =gamma(k)
     so the group velocity of electron (at x-direction) will be
           vg=(hbar/m*)*kx/sqrt[1+4*alpha*gamma(k)] that is why in the below code we have
           the parameter sq = sqrt[1+4*alpha*gamma(k)]  
           
Starting date: March 19, 2010
Latest update: March 20, 2010 
***************************************************************************** */
#include <math.h>
#include <stdio.h>
#include "constants.h"

void drift(double tau){
  // Goi cac ham
    void check_boundary();
    int Get_iv(),Get_subb();
    double Get_kx(),Get_ml(), Get_mt(), Get_nonparabolicity_factor();
    double ***Get_eig(), Get_mesh_size_x(),Get_x_position();
    void Set_electron_energy(double ener),Set_x_position(double position_in_x),Set_kx(double k_momen_x);
    
    //Cac bien local
    int v = Get_iv();    // v-th valley
    int i  = Get_subb(); // i-th subband
    double kx = Get_kx(); //printf("\n Before kx=%le",kx);
    double ml = Get_ml();
    double mt = Get_mt();
    double af = Get_nonparabolicity_factor();// printf("\n af = %le", af);//xem da o dang 0.5/eV CHUA?
    double ***eig = Get_eig();
    
// Thuc hien    
    // Tim hat dang o SECTION nao ?
    double mesh_size_x = Get_mesh_size_x();
    double x_position = Get_x_position(); //printf("\n Before x position=%le",x_position);
    int s = (int)(round(x_position/mesh_size_x)); // s-th section 
    
    // Ban dau KHONG cho HAT VUOT QUA vi tri cho phep
    int Get_nx0(),Get_nx1();
    int nx0 = Get_nx0(); int nx1 = Get_nx1();
    int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    if(s < 0)      { s = 0;}
    if(s > nx_max) {s = nx_max;}
        
    // Tinh electric field theo phuong x
    double Fx;
    if(s >= nx_max) { s = nx_max-1; } // Do eig ta chi co den gia tri la nx_max ma thoi
    Fx = (eig[s+1][v][i]-eig[s][v][i])/mesh_size_x;//KHONG nhan q de doi ra don vi chuan. Xem chu y 1 trang 113
                                                    // Cong thuc (2) page 112
    // Tinh dkx
    double dkx;
    //dkx = -(q/hbar)*Fx*tau;// [don vi chuan] // cong thuc (1) trang 112
    dkx = (q/hbar)*Fx*tau;// [don vi chuan] // cong thuc (1) trang 112
    
    // Update energy
    double E_parab; // energy tinh theo parabolic
    switch (v)
          {
             case 1: // v=1: valley pair 1 thi m* la ml
                  E_parab = (hbar*hbar)*(kx+dkx)*(kx+dkx)/(2.0*ml*m0); // [J]
                  break;
             case 2: // v=2: valley pair 2 thi m* la mt
                  E_parab = (hbar*hbar)*(kx+dkx)*(kx+dkx)/(2.0*mt*m0); // [J]
                  break;
             case 3: // v=3: valley pair 3 thi m* la mt
                  E_parab = (hbar*hbar)*(kx+dkx)*(kx+dkx)/(2.0*mt*m0); // [J]
                  break;
             default:
	       {
                  printf("\n Case iv=9 was considered in emcd(). Wrong! CHECK dritf.c");
		  exit(1);
	       }
          } // End of switch (v)

    // Chuyen E_parab tu [J] sang [eV] 
    E_parab = E_parab/q; // [eV]  
    double E_Nonparab; // energy tinh theo Non-parabolic: TA DUNG
    E_Nonparab = (sqrt(1+4.0*af*E_parab) -1.0) /(2.0*af); // [eV]
    
    Set_electron_energy(E_Nonparab);// [eV] // Update gia tri cho bien energy
    
    // Update position
    switch (v)
          {
             case 1: // v=1: valley pair 1 thi m* la ml
                  x_position=x_position+hbar*tau*(kx+0.5*dkx)/((ml*m0)*sqrt(1+4.0*af*E_parab));//[m]
                  break;
             case 2: // v=2: valley pair 2 thi m* la mt
                  x_position=x_position+hbar*tau*(kx+0.5*dkx)/((mt*m0)*sqrt(1+4.0*af*E_parab));//[m]
                  break;
             case 3: // v=3: valley pair 3 thi m* la mt
                  x_position=x_position+hbar*tau*(kx+0.5*dkx)/((mt*m0)*sqrt(1+4.0*af*E_parab));//[m]
                  break;
             default:
                  printf("\n Case iv=9 was considered in emcd(). Wrong! CHECK dritf.c");
          } // End of switch (v)
  
    Set_x_position(x_position);// [m] // Update gia tri cho bien x_position
    
    // Update momentum
    kx = kx +dkx; // [1/m]
    Set_kx(kx); // Update gia tri cho bien kx   
    
    //printf("\n tau=%le    dkx = %le",tau,dkx);
    //printf("\n section=%d Electric field Fx = %le",s, Fx); 
    //printf("\n Energy=%le, x position=%le, kx=%le",E_Nonparab,x_position,kx);
            
    // Goi ham check_boundary() DE KIEM TRA dieu kien boundary     
    check_boundary(); // Sau ham nay duoc gia tri iv MOI: Neu iv KHAC 9 thuc hien o day
                      //                                  Neu iv=9 xem o ham emcd()
    /* Chu thay y nghia cua doan nay    
    int ix;
    int iv_update = Get_iv();
    if (iv_update !=9){ // HAT DUOC  DUNG do co iv KHAC 9
       ix = (int)(x_position/mesh_size_x);          
       if(ix > nx_max) { ix = nx_max; }
       if(ix < 0){ ix = 0; }
    }// End of if (iv_update !=9)
    */

  return;
} // End of void drift(double tau)
/*******************************************************************************
/* *****************************************************************************
    + To check boundaries: 
       => Source and drain regions
 
Starting date: March 21, 2010
Latest update: March 21, 2010 
*********************************************************************************/

void check_boundary(){
  
     int Get_nx0(),Get_nx1();
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; 

     double Get_x_position(),Get_mesh_size_x();
     double x_position = Get_x_position(); // Lay VI TRI cua HAT
     double device_length=(double)(nx_max)*Get_mesh_size_x();// printf("\n device_length=%le",device_length);//Lay DO DAI cua Device    
     
// Thuc hien
     int iv_update,iss_out_update,idd_out_update;
     int Get_iss_out(),Get_idd_out();
     void Set_iv(int valley_pair_index),Set_iss_out(int out_p),Set_idd_out(int out_p);
     // HAT co the OUT neu VUOT QUA GIOI HAN o S va D. luon la Side Contact cho MSMC
     if(x_position <= 0.0){// Absorption this particle o vung Source contact
        iv_update = 9; // ngu y rang HAT nay se bi LOAI
        Set_iv(iv_update);// Update gia tri ve valley index cho bien iv
                
        // 3 lenh sau day chi TUONG DUONG voi 1 lenh iss_out = iss_out + 1; neu su dung GLOBAL variable
        iss_out_update = Get_iss_out();//lay so luong hat iss_out DA CO
        iss_out_update = iss_out_update + 1;// Them 1 hat nua bi LOAI
        Set_iss_out(iss_out_update); // Update gia tri ve SO HAT bi LOAI cho bien iss_out
         
        //printf("\n Vuot qua Souce Contact iss_out_update=%d",iss_out_update);
     	// return; // KHONG CAN thuc hien cac lenh sau nua NEU da thuc hien cac lenh tren
      } // End of if(x_position <= 0.0){
     
     if(x_position >= device_length){// Absorption this particle o vung Drain contact
        iv_update = 9;
        Set_iv(iv_update);
        idd_out_update = Get_idd_out();
        idd_out_update = idd_out_update + 1;
        Set_idd_out(idd_out_update);
        //printf("\n Vuot qua Drain Contact idd_out_update=%d",idd_out_update);
            
     }// End of if(x_position >= device_length)

return;
}
