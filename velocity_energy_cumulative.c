/* *****************************************************************************
To calculate cumulative velocities, energies and current density (su dung
   average velocity) of the particles
      
 INPUT:+ time:
            - Neu khong co dieu kien gi: tinh tu first dt cho den total time
            - Neu co DIEU KIEN time > transient_time: chi tinh sau transient time
            NOTE: chi bang cach tinh nhu the nay thi moi BIET chon transient time
            bao nhieu la vua
       + n_time_steps_average=(tot_time-transient_time)/dt -> so lan tinh trung binh
       + iteration_reference: 
              =1: tinh tong cac gia tri theo tung khoang thoi gian dt  
              =0: luc dt da chay den diem cuoi cung cua total_time ->trung binh cac gia tri
                         
 OUTPUT:
       + Write averages in a file (1 file cho all time, 1 file cho after transient time)
           x-position, velocity_x_average, energy_average, current_sum

Starting date: April 1, 2010
Latest update: June 15, 2010
****************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "nrutil.h"
#include "constants.h"
   
void velocity_energy_cumulative(FILE *f1, FILE *f2, double time,int iteration_reference){//f1 cho All, f2 cho Transient
   // Goi cac ham
   int Get_n_used(),*Get_valley();// lay so hat dang co trong device co iv=1,2,3 va co the 9
   double **Get_p(),*Get_energy();// tham so hat
   double Get_nonparabolicity_factor(),Get_ml(), Get_mt(); 
   double Get_dt(),Get_tot_time(),Get_transient_time(); // Xac dinh tu Input.d
   double Get_mesh_size_x();
   double *Get_velocity_x_sum();
   double *Get_energy_sum();
   double *Get_current_sum();
   double *Get_velocity_x_sum_after_transient();
   double *Get_energy_sum_after_transient();
   double *Get_current_sum_after_transient();
   int Get_save_VelEnerCurr_all_time_steps(); // ca velocity, energy va current
   int Get_save_VelEnerCurr_after_transient();// ca velocity, energy va current
   void save_VeloEnerCurrent_all_time_steps(FILE *f1, double *velo,double *ener,double *curr,int steps,int nx_max,double mesh_size_x);
   void save_VeloEnerCurrent_after_transient(FILE *f2,double *velo,double *ener,double *curr,int n_ti,int nx_max,double mesh_size_x);
   // Cac bien local
   int n_used = Get_n_used(); //printf("\n n_used at velocity_energy_cumulative() = %d",n_used);
   int *valley = Get_valley();
   double **p   = Get_p();
   double *energy= Get_energy();
   double af = Get_nonparabolicity_factor();//[1/eV]
   double ml = Get_ml(); // chi la gia tri ml, CHUA nhan voi m0
   double mt = Get_mt();
   double dt = Get_dt(); // Observation time step
   double tot_time = Get_tot_time(); // Tong thoi gian chay Simulation
   double transient_time = Get_transient_time(); // Chon bang BAO NHIEU can KINH NGHIEM
   
   int Get_nx0(),Get_nx1(); 
   int nx0 = Get_nx0(); int nx1 = Get_nx1();
   int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
  
   double mesh_size_x = Get_mesh_size_x();
   double *velocity_x_sum = Get_velocity_x_sum();
   double *energy_sum = Get_energy_sum();
   double *current_sum = Get_current_sum();
   double *velocity_x_sum_after_transient = Get_velocity_x_sum_after_transient();
   double *energy_sum_after_transient = Get_energy_sum_after_transient();
   double *current_sum_after_transient = Get_current_sum_after_transient();
   int save_velo_ener_current_all_time_steps = Get_save_VelEnerCurr_all_time_steps();
   int save_velo_ener_current_after_transient = Get_save_VelEnerCurr_after_transient();
   
// Thuc hien
   int total_time_steps  =  (int)(tot_time/dt); // tong so khoang thoi gian dt
   int n_time_steps_after_transient = (int)((tot_time-transient_time)/dt);//so khoang step cho dt tinh tu sau transient_time
   double x_position = 0.0;
   int n_valley = 3;// number of valley pairs = 3
   double **velx_sum,**sum_ener,*el_number; // Bien temperary luc tinh TAT CA time steps
   double **velx_sum_after_transient,**sum_ener_after_transient,*el_number_after_transient; // neu chon chi tinh sau transient time
   velx_sum = dmatrix(0,nx_max,1,n_valley);// Van toc_x o vi tri i va valley index j
   sum_ener = dmatrix(0,nx_max,1,n_valley);// Nang luong cua electron o vi tri i va valley index j
   el_number = dvector(0,nx_max); // So luong electron o vi tri i (chinh la cho moi section)
   velx_sum_after_transient  = dmatrix(0,nx_max,1,n_valley);// Chu "sum" nghia la sume cho cac time step
   sum_ener_after_transient  = dmatrix(0,nx_max,1,n_valley);// do vay sau do phai chia cho tong time step
   el_number_after_transient = dvector(0,nx_max);           // ma no sum
   int i,j;
   for(i=0; i<=nx_max; i++){ el_number[i]=0.0; el_number_after_transient[i]=0.0;
       for(j=1; j<=n_valley; j++){
           velx_sum[i][j]=0.0; velx_sum_after_transient[i][j]=0.0;
           sum_ener[i][j]=0.0; sum_ener_after_transient[i][j]=0.0;
       }
    }// End of for(i=0; i<=nx_max; i++)
    
// Buoc 1. CALCULATE THE AVERAGE VALUE IF iteration_reference==0 ********* 
if(iteration_reference==0){ //Ban dau gan iteration_reference=1 sau do neu j_iter=iter_total
                            // thi ta cho iter_reference=0; de tinh trung binh cho tat ca time steps ma ta cong
                               
      // Buoc 1.1. Co SAVE cho all time steps KHONG?
      if(save_velo_ener_current_all_time_steps == 1){//=1 co SAVE
         // Goi ham de save                          
	save_VeloEnerCurrent_all_time_steps(f1,velocity_x_sum,energy_sum,current_sum,total_time_steps,nx_max,mesh_size_x);
         //Sau khi ghi vao file xong thi minh refresh no lai de dung cho lan sau tuc la cap (Vg va Vd khac)
         for(i=0; i<=nx_max; i++){
             velocity_x_sum[i] = 0.0;
             energy_sum[i]     = 0.0;
             current_sum[i]    = 0.0;
          } // End of for
     }// End of if( save_velo_ener_curr_cumu_all_time_steps_or_not ==1)
     // ***** End of Buoc 1.1.
     
      // Buoc 1.2. Co SAVE after transient time KHONG?  
      // Can 2 DIEU KIEN: Co save after transient va time > transient time
      if((save_velo_ener_current_after_transient==1)&&(time > transient_time)){ // =1 co SAVE  
        // Goi ham de SAVE
	save_VeloEnerCurrent_after_transient(f2, velocity_x_sum_after_transient,energy_sum_after_transient,
                               current_sum_after_transient,n_time_steps_after_transient,nx_max,mesh_size_x);
        //Sau khi ghi vao file xong thi minh refresh no lai de dung cho lan sau tuc la cap (Vg va Vd khac)
        for(i=0; i<=nx_max; i++){
             velocity_x_sum_after_transient[i] = 0.0;
             energy_sum_after_transient[i]     = 0.0;
             current_sum_after_transient[i]    = 0.0;
        } // End of for(i=0; i<=nx_max; i++)
     }// End of if((save_velo_ener_current_after_transient==1)&&(time > transient_time))
     // ***** End of Buoc 1.2.
}// End of  if(iter_reference==0)
//**** End of Buoc 1. CALCULATE THE AVERAGE VALUE IF iteration_reference==0       
 
//Buoc 2. CALCULATE TONG velocity, energy, current NEU iter_reference !=0
         // thuc hien tinh cac gia tri tren theo tung time step dt, tu dt, 2dt,3dt,...
  int n=0,valley_index=0;
  double ee=0.0, velx=0.0, kx=0.0;
  
  // Buoc 2.1. ****************************************************************
  if(save_velo_ener_current_all_time_steps == 1){//=1 co SAVE. Yeu cau SAVE thi moi CAN TINH 
     for(n=1; n<=n_used; n++){ // Chay cho tung hat tu 1 den n_used
         if(valley[n] != 9){
            i = (int)(p[n][2]/mesh_size_x+0.5);
            if(i <0)      { i=0; }
            if(i >nx_max) { i=nx_max;}
                 
            valley_index = valley[n]; // co the la 1, 2 hoac 3
            ee = energy[n]; // [ev] lay energy cua hat
            kx = p[n][1]; // [1/m]lay momentum
          
            switch (valley_index)
              {
               case 1: // valley_index=1: valley pair 1 thi m* la ml
                     velx = (hbar*kx)/(ml*m0*(1.0+2.0*af*ee));//[m/s]
                     break;
               case 2: // valley_index=2: valley pair 2 thi m* la mt
                     velx = (hbar*kx)/(mt*m0*(1.0+2.0*af*ee));//[m/s]
                     break;
               case 3: // valley_index=3: valley pair 3 thi m* la mt
                     velx = (hbar*kx)/(mt*m0*(1.0+2.0*af*ee));//[m/s]
                     break;
               default:
                     printf("\n Something Wrong! CHECK velocity_energy_cumulative.c");
             } // End of switch (valley_index)
   
             velx_sum[i][valley_index] += velx; // cong cho tat ca cac hat o cung section
             sum_ener[i][valley_index] += ee;                       
             el_number[i] += 1;    
        }//End of if(valley[n] != 9)     
      } // End of for(int n=1;n<=n_used;n++)
       
      double temp_velo=0.0, temp_ener=0.0;
      for(i=0; i<=nx_max; i++){ // Chay theo phuong x chi tung section
          if(el_number[i]!=0){
             temp_velo = velx_sum[i][1]+velx_sum[i][2]+velx_sum[i][3];//tong velocity tren 3 valleys
             velocity_x_sum[i] += temp_velo/el_number[i]; //[m/s]
             current_sum[i]    += temp_velo*q/mesh_size_x; //[A/m] Tinh theo cach 2
             temp_ener = (sum_ener[i][1]+sum_ener[i][2]+sum_ener[i][3])/el_number[i];
             energy_sum[i] += temp_ener;
          }// End of if(el_number[i]!=0)  
      }// End of for(i=0; i<=nx_max; i++)
  }// End of if(save_velo_ener_current_all_time_steps == 1)
  // End of Buoc 2.1. **********************************************************
   free_dmatrix(velx_sum,0,nx_max,1,n_valley);
   free_dmatrix(sum_ener,0,nx_max,1,n_valley);
   free_dvector(el_number,0,nx_max); // Cac bien KHONG CAN DUNG thi nen Free luon tranh nham lan
    
  // Buoc 2.2. ****************************************************************
  // Can 2 DIEU KIEN: Co save after transient va time > transient time
  if((save_velo_ener_current_after_transient==1)&&(time > transient_time)){// =1 co SAVE
     for(n=1; n<=n_used; n++){ 
        if(valley[n] != 9){
           i = (int)(p[n][2]/mesh_size_x+0.5);
           if(i <0)      { i=0; }
           if(i >nx_max) { i=nx_max;}
           valley_index = valley[n]; 
           ee = energy[n]; 
           kx = p[n][1]; 
           switch (valley_index)
             {
               case 1: // valley_index=1:  ml
                     velx = (hbar*kx)/(ml*m0*(1.0+2.0*af*ee));
                     break;
               case 2: // valley_index=2: mt
                     velx = (hbar*kx)/(mt*m0*(1.0+2.0*af*ee));
                     break;
               case 3: // valley_index=3: mt
                     velx = (hbar*kx)/(mt*m0*(1.0+2.0*af*ee));//[m/s]
                     break;
               default:
                     printf("\n Something Wrong! CHECK velocity_energy_cumulative.c");
            } // End of switch (valley_index)
   
            velx_sum_after_transient[i][valley_index] += velx; 
            sum_ener_after_transient[i][valley_index] += ee;                       
            el_number_after_transient[i] += 1;    
       }//End of if(valley[n] != 9)     
     } // End of for(int n=1;n<=n_used;n++)
       
    double temp_velo=0.0, temp_ener=0.0;
    for(i=0; i<=nx_max; i++){ // Chay theo phuong x cho tung section
        if(el_number_after_transient[i]!=0){
           temp_velo=velx_sum_after_transient[i][1]+velx_sum_after_transient[i][2]+velx_sum_after_transient[i][3];
           velocity_x_sum_after_transient[i] += temp_velo/el_number_after_transient[i]; //[m/s]
           current_sum_after_transient[i]    += temp_velo*q/mesh_size_x; //[A/m] Tinh theo cach 2
           temp_ener = (sum_ener_after_transient[i][1]+sum_ener_after_transient[i][2]+sum_ener_after_transient[i][3])/el_number_after_transient[i];
           energy_sum_after_transient[i] += temp_ener;
         }// End of if(el_number[i]!=0)  
     }// End of for(i=0; i<=nx_max; i++)
  }//End of if((save_velo_ener_current_after_transient==1)&&(time > transient_time) 
  // End of // Buoc 2.2. *******************************************************
//End of buoc 2
   free_dmatrix(velx_sum_after_transient,0,nx_max,1,n_valley);
   free_dmatrix(sum_ener_after_transient,0,nx_max,1,n_valley);
   free_dvector(el_number_after_transient,0,nx_max);  
  // Chu y neu cac bien cuc bo ma ta khoi tao trong 1 ham bang phuong phap trong quyen sach phuong phap so thi
  //thi cuoi ham ta can free no, neu khong no se gay ra loi tran bo nho khi ta dung ham nay lap di lap lai trong vong while 
  //chang han
     
  return;
}
//*****************************************************************************
/* *****************************************************************************
 De save_velo_ener_current_all_time_steps vao trong 1 file. 
 INPUT:  + velocity_x_sum[i]: i chay tu 0 den nx_max
         + energy_sum[i]: i chay tu 0 den nx_max
         + current_sum[i]: i chay tu 0 den nx_max
         + total_time_step
         + nx_max
         + mesh_size_x
                 
 OUTPUT: + velocity_x_sum[i], energy_sum[i] va current_sum[i] da CHIA CHO total_time_steps
         de lay duoc gia tri trung binh
         Save vao file ten la Velo_ener_current_averages_all_time_steps.dat
        
Starting date: April 02, 2010
Latest update: April 02, 2010
****************************************************************************** */
void save_VeloEnerCurrent_all_time_steps(FILE *f1, double *velocity_x_sum,double *energy_sum,double *current_sum,int total_time_steps,int nx_max,double mesh_size_x)
{
    double x_position = 0.0;
    //FILE *f;
    //f=fopen("Velo_ener_current_averages_all_time_steps.dat","a"); // Chu y kieu mo la "a" co the GHI THEM
    //if(f==NULL){ printf("\n Cannot open file Velo_ener_current_averages_all_time_steps.dat");return;} 
    //fprintf(f,"\n   #x[nm]    velocity_x[m/s]  energy[eV]   current\n");
   
    int i , myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank ==0){
      for(i=0; i<=nx_max; i++){
	x_position = (double)(i)*mesh_size_x/1.0e-9;// doi tu [m] sang [nm] 
	fprintf(f1,"%f   %le    %le  %le \n",
                     x_position, 
                     velocity_x_sum[i]/((double)(total_time_steps)), // [m/s]
                     energy_sum[i]/((double)(total_time_steps)), // [eV]   
                     current_sum[i]/((double)(total_time_steps))
                 );     
      } // End of for(i=0;i<=nx_max;i++)
    }
    //fclose(f);
}// End of void save_VeloEnerCurrent_all_time_steps(

/* *****************************************************************************
 De save_velo_ener_current_after_transient vao trong 1 file. 
 INPUT:  + velocity_x_sum_after_transient[i]: i chay tu 0 den nx_max
         + energy_sum_after_transient[i]: i chay tu 0 den nx_max
         + current_sum_after_transient[i]: i chay tu 0 den nx_max
         + n_time_steps_after_transient
         + nx_max
         + mesh_size_x
                 
 OUTPUT: + velocity_x_sum_after_transient[i], energy_sum_after_transient[i] va
           current_sum_after_transient[i] da CHIA CHO n_time_steps_after_transient
           de lay duoc gia tri trung binh
           Save vao file ten la Velo_ener_current_averages_after_transient.dat
        
Starting date: April 02, 2010
Latest update: April 02, 2010
****************************************************************************** */
void save_VeloEnerCurrent_after_transient(FILE *f2, double *velocity_x_sum_after_transient,double *energy_sum_after_transient,
      double *current_sum_after_transient,int n_time_steps_after_transient,int nx_max,double mesh_size_x)
{
  //FILE *f;
  //f=fopen("Velo_ener_current_averages_after_transient.dat","a"); // Chu y kieu mo la "a"
  //if(f==NULL){ printf("\n Cannot open file Velo_ener_current_averages_after_transient.dat");return;} 
  //fprintf(f,"\n #x[nm] velocity_x[m/s] energy[eV] current\n");
    double x_position = 0.0;
    int i, myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank ==0){
      for(i=0; i<=nx_max; i++){
        x_position = (double)(i)*mesh_size_x/1.0e-9;// doi tu [m] sang [nm] 
        fprintf(f2,"%f %le %le %le \n",
                  x_position, 
                  velocity_x_sum_after_transient[i]/((double)(n_time_steps_after_transient)), // [m/s]
                  energy_sum_after_transient[i]/((double)(n_time_steps_after_transient)), // [eV]   
                  current_sum_after_transient[i]/((double)(n_time_steps_after_transient))
             );     
      } // End of for(i=0;i<=nx_max;i++)
    }
      //fclose(f);
}// End of void save_VeloEnerCurrent_after_transient



