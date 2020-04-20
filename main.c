/* *****************************************************************************
        To solve MSMC

Starting date: Feb 8, 2010
Update:        April 6, 2010
Latest update: May 14, 2010
****************************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "functions.h"
#include "constants.h"
#include "nanowire.h"
#include "petscksp.h"

main(int argc, char **argv)
{
//Buoc 1. Tinh thoi gian bat dau
    time_t time_start, time_finish;
    time_start = time(NULL);// Lay gio he thong
    //printf("The start time is %s\n",asctime(localtime(&time_start)));// xem o buoc 3
    
// Buoc 2. Initializa MPI 
    MPI_Init(&argc,&argv);// Sau ham MPI_Init thi cac processor BAT DAU Lam viec SONG SONG
    
// Buoc 3. Lay rank va size
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    double time_tic, time_toc; // thoi gian bat dau va ket thuc 1 chu trinh nao do
    FILE *logfile = fopen("logfile_MSMC.txt","w");
    if(myrank ==0){ 
      time_tic = MPI_Wtime(); //  void logging(FILE *log,const char *fmt,...);
      
      printf("The start time is %s\n",asctime(localtime(&time_start)));
      }

// Buoc 4.Phan chay cho TAT CA cac processor. Ta da check rat nhieu lan roi. Bang cach in ket qua khi rank khac 0
    material_param();
    device_structure();
    read_voltages_input();
    read_simulation_list();
    read_scattering_save_list();
    initial_variables();// Can khoi tao sau ham device_structure() va read_simulation_list();
	                    // va truoc ham doping_potential_initialization(); 

    Initialize(&argc,&argv); // for 3D Poisson . logfile_PoissonInput(); // For checking

    initial_potential_doping();// Can dat ham nay truoc ham source_drain_carrier_number() de co doping ma thoi. o day do chua set Vs,Vd,Vg nen no bang 0
    source_drain_carrier_number();

    if(myrank ==0){ 
      logging(logfile,"\n Number of particles for Source and Drain contact"); 
      logging(logfile,"\n nsource_side_carriers = %d",Get_nsource_side_carriers());
      logging(logfile,"\n ndrain_side_carriers  = %d",Get_ndrain_side_carriers());  
     
      printf("\n nsource_side_carriers = %d",Get_nsource_side_carriers());
      printf("\n ndrain_side_carriers  = %d",Get_ndrain_side_carriers());       
     }
 
    trapezoidal_weights();// tranh chay nhieu lan cho ham nay //save_results();
      
// Buoc 5. Open Current-Voltage files // Mo kieu nay co nghia o 1 node thoi
    FILE *f_IdVg, *f_IdVd;

    if(myrank ==0){ // Mo 2 file chua I-V tai 1 node thoi

      f_IdVg=fopen("IdVg.dat","a");// Mo ra de ghi de, vi co the se chay nhieu vong lap cho Drain va Gate voltage
      if(f_IdVg==NULL){printf("\n Cannot open file IdVg.dat");return 0;}
      fprintf(f_IdVg,"\n #V_drain_input_fixed  V_gate_input  Id \n");
      
      f_IdVd=fopen("IdVd.dat","a"); 
      if(f_IdVd==NULL){printf("\n Cannot open file IdVd.dat");return 0;}
      fprintf(f_IdVd,"\n #V_gate_input_fixed  V_drain_input  Id \n");
    }

    if(myrank ==0){ 
      time_toc = MPI_Wtime(); //  void logging(FILE *log,const char *fmt,...);
      logging(logfile,"\n Time to initialize all functions before Applied Voltage Bias loop =%f [s]",time_toc - time_tic); 
      logging(logfile,"\n ******************************************************************************* ");
    }
    
 // Buoc 6.
// ************************* APPLIED VOLTAGE BIAS LOOP *************************
// Solve 1D Multi Subband MC, 2D Schrodinger and 3D Poisson self-consistently for each point Vd and Vg
    // Lay cac gia tri voltage tu Input file
    double Vg_start = Get_Vg_start();
    double Vg_end   = Get_Vg_end();
    double Vg_step  = Get_Vg_step();
    double Vd_start = Get_Vd_start();
    double Vd_end   = Get_Vd_end();
    double Vd_step  = Get_Vd_step(); 

    if(myrank ==0){ 
      logging(logfile,"\n    Vgate start = %f", Vg_start);
      logging(logfile,"\n    Vgate end   = %f", Vg_end );
      logging(logfile,"\n    Vgate step  = %f", Vg_step);
      logging(logfile,"\n    Vdrain start= %f", Vd_start);
      logging(logfile,"\n    Vdrain end  = %f", Vd_end);
      logging(logfile,"\n    Vdrain step = %f", Vd_step); 
      logging(logfile,"\n*********************************************************\n\n");
    }

    double Id=0.0,cur_av=0.0,V_gate_input=0.0,V_drain_input=0.0;

for(V_gate_input=Vg_start; V_gate_input<=Vg_end; V_gate_input=V_gate_input+Vg_step){
    for(V_drain_input=Vd_start; V_drain_input<=Vd_end; V_drain_input=V_drain_input+Vd_step){

// Buoc 6.1. 
      /* ************************************************************************************* 
	 Solve Poisson equation LAN DAU TIEN for electron initialization
	 - Do initial state of the system cang gan steady state thi convergence cang nhanh
	 -> Initial guess su dung 1 phuong phap co ve ADVANCE hon
	 -> To achieve a smoother initial potential and electron density, which leads to faster
	 convergence 
      ************************************************************************************* */
     
       SetDrainVoltage(0.0);// NOTE cho Vd=0.0 o day
       SetGateVoltage(V_gate_input);//printf("\n Vd=%f, Vgate_input=%f",GetDrainVoltage(),GetGateVoltage());
                         
       void SetPhi_Contact(); // o definefn.c. No dung ham void assign_built_in_potential(). Ta DA SUA ham nay mot ti roi
       SetPhi_Contact();//Tat ca cac node.  Phai dat ham nay TRUOC ham guess_potential(), void set_initial_potential() Neu khong No chua co gia tri Vs va Vd nen =0 tat ca.

       void guess_potential();// Da them tham so vao ham guess_potential() de chay tot
       guess_potential();     //Tat ca cac node. void set_initial_potential();// Duoc Phi(i,j,k) cua GS
                             //set_initial_potential();// Ta khong doc tu file pot.r nen NO SE la guess_potential(); Cai nay la tu GS neu no chuyen truc tiep vao Poisson luon
             
       void get_initial_charge_density_for_Poisson();// Tinh cho 1 node va distribute-> Tat ca cac node
       get_initial_charge_density_for_Poisson();// tinh duoc electron density lan dau tien no cung chinh la n3d lan dau tien
            //Co save cai charge dua vao Poisson lan dau tien khong ?
            //save_initial_charge_density_for_Poisson();
     	   
       if ( myrank==0 ){
	 printf("\n Solving Poisson the first time to inilialize electrons");
       }

       if (GetPoissonLinearized()) { SolvePoisson = solve_3d_poisson_linearized;}
       else                        { SolvePoisson = solve_3d_poisson;}

       (*SolvePoisson)();// Chon kieu gi o tren thi solve Poisson theo kieu do 

       void mpi_gather_potential_at_node0();// Sau khi giai 3D Poisson xong thi can gather chu potential lai.RAT NOTE day. Vi potential ban dau dang rai rac
       mpi_gather_potential_at_node0();     // Neu chon Print_Midline_Potential(y/n) = y hoac print_output(); thi no da gather lai cho minh roi.

       void mpi_distribute_potential_from_node0();//RAT NOTE Sau ham mpi_gather_potential_at_node0() ta duoc 1 potential hoan hao va ta gui no den cac node khac
       mpi_distribute_potential_from_node0(myrank, mysize);// node = myrank, np = mysize
      
       electrons_initialization();// tu potential o Poisson_0 ta initialization cho electrons. Da kiem tra rang cac node da tao ra so electrons giong nhau

       if(myrank ==0){ 
	 logging(logfile,"\n For Equilibrium Vg =%f, Vd = %f",GetGateVoltage(),GetDrainVoltage());
	 logging(logfile,"\n Number of initialized particles using Equilibrium above = %d",count_used_particles()); 
	 logging(logfile,"\n Maximum number of electrons allowed (in constants.h)    = %d", max_electron_number);

	 printf("\n Number of particles in use = %d",count_used_particles());
       }
      
            // Save parameters cua electron lan dau tien khong     
            //char fn[100];
            //sprintf(fn,"ElectronInitParameters-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //save_electron_parameters(fn);
            
            //char fnSubband[100], fnEnergy[100];
            //sprintf(fnSubband,"ElectronInitDistributionSubbands-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //sprintf(fnEnergy,"ElectronInitDistributionEnergy-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //save_electron_distribution(fnSubband,fnEnergy);
	
       /*End of solving equilibrium Poisson equation for electron initialization */
       /* ***************************************************************************** */

 // Buoc 6.2.
       /* **************************************************************************************       
                             For Non-Equilibrium case
                       Apply all actual biases at the contacts
	*********************************************************************************** */
       SetDrainVoltage(V_drain_input); // Vs=0.0; Vd=Vd_input
       SetGateVoltage(V_gate_input);// Vg=Vgate_input. Da set o tren roi

       if(myrank ==0){ 
	 logging(logfile,"\n");
	 logging(logfile,"\n Solving  MSMC, Schrodinger and Poisson for the bias");
	 logging(logfile,"\n    Source Voltage [V] = %3.3f", Get_V_source_input());
	 logging(logfile,"\n    Drain Voltage  [V] = %3.3f", GetDrainVoltage()); 
	 logging(logfile,"\n    Gate Voltage   [V] = %3.3f", GetGateVoltage());

	 printf("\n");// Dang in cho node dau tien
	 printf("\n ****************************************************");
	 printf("\n   Solving  MSMC, Schrodinger and Poisson for the bias");
	 printf("\n        Source Voltage [V] = %3.3f", Get_V_source_input()); 
	 printf("\n        Drain Voltage  [V] = %3.3f", GetDrainVoltage());  
	 printf("\n        Gate Voltage   [V] = %3.3f", GetGateVoltage());
	 printf("\n ****************************************************\n");
       }
                   
      SetPhi_Contact();// Phai dat ham nay TRUOC ham guess_potential(), void set_initial_potential() Neu khong No chua co gia tri Vs va Vd nen =0 tat ca.

      guess_potential();
       
      get_initial_charge_density_for_Poisson();// tinh duoc electron density lan dau tien no cung chinh la n3d lan dau tien
            //Co save cai charge dua vao Poisson lan dau tien khong ?
            //save_initial_charge_density_for_Poisson(); logfile_PoissonInput(); // For checking
       
      if (GetPoissonLinearized()) { SolvePoisson = solve_3d_poisson_linearized;}
      else                        { SolvePoisson = solve_3d_poisson;}

      (*SolvePoisson)();// Chon kieu gi o tren thi solve Poisson theo kieu do 

      mpi_gather_potential_at_node0();// DU co print_output() thi van nen de no o dang tuong minh nhu the nay
      mpi_distribute_potential_from_node0(myrank, mysize);// node = myrank, np = mysize

      print_output();
     
      if ( myrank==0 ){
	 printf("\n Solving 2D Schrodinger the first time");
      }
      //Solved_2D_Schro_Parallel_for_MSMC(&argc, &argv);// Su dung parallel
      Solved_2D_Schro_for_MSMC();// Giai Schrodinger LAN DAU se goi ham  get_potential(int s, int ny, int nz);
            /*
            //Save potential cai su dung de giai 2D Schrodinger
            char fn_x_FirstTime[100], fn_yz_FirstTime[100];
            sprintf(fn_x_FirstTime,"potential_x_at_midpoint_yz_in2DSchrodingerFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            sprintf(fn_yz_FirstTime,"potential_yz_at_midpoint_x_in2DSchrodingerFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            save_potential_used_for_2D_Schrodinger(fn_x_FirstTime,fn_yz_FirstTime);

	    //Save subband_energy cho lan dau tien. Neu o get_potential() cho Schrodinger lay pot=-Phi va show ket qua cung chi la -Phi thi TRUNG KHIT
            char fnsubbandFirstTime[100];
	    sprintf(fnsubbandFirstTime,"subbandFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	    save_subband_energy(fnsubbandFirstTime);
            */          
            //SAVE wave va eigen neu muon cho lan dau tien
            //int save_eig_wavefuntion = Get_save_eig_wavefuntion();
            //char fneigenFirstTime[100], fnwaveFirstTime[100];
            //sprintf(fneigenFirstTime,"eigenFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //sprintf(fnwaveFirstTime,"waveFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //int NumSections = Get_nx1();
            //int s,v;
            //if ( myrank==0 ){ 
            // if(save_eig_wavefuntion==1){//Save eigen function and wave function ?
	    //   for(v=1; v<=3; v++){// 3 valleys
                  ////for(s=0; s<=NumSections;s++){// Tat ca cac sections
                 // s=3;
                 // save_eig_wave(s,v,fneigenFirstTime,fnwaveFirstTime);//neu de file dang "a" thi khi ve wave lai bi chong chap can chu y
                // }
	       // }
            // }                 
            //}// End of if ( myrank==0 )
      
       form_factor_calculation();// Tinh Form factor LAN DAU. // Tinh o TAT CA cac node
            // Co save form factor lan dau khong ?
            //char fnformfactorFirstTime[100];
            //sprintf(fnformfactorFirstTime,"FormfactorFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	    //int NumSections = Get_nx1();
            //int FormfactoSavedAtsth = (int)(NumSections/2 +0.5);// lay o giua
            //save_form_factor_calculation(FormfactoSavedAtsth, fnformfactorFirstTime);
            
       
       scattering_table(); //LAN DAU // save_results(); // truoc khi normalized scattering // Tinh o TAT CA cac node
            // Save scattering lan dau tien khong ?
            //char fnscatFirstTime[100];
            //sprintf(fnscatFirstTime,"ScatTableFirstTime-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	    //int NumSections = Get_nx1();
            //int ScatSavedAtsth = (int)(NumSections/2 +0.5);// lay o giua
            //save_scattering_table(ScatSavedAtsth,fnscatFirstTime);
           
       normalize_table(); //LAN DAU // //save_results(); getchar();// La save sau khi normalize table. KHONG CO Y NGHIA GI cho scattering table ca
                          // Tinh o TAT CA cac node
                           //void save_normalized_table();
                          //save_normalized_table();
             
       init_free_flight(); //LAN DAU // Tinh o TAT CA cac node

      
       multi_subbands_MC(logfile, &argc,&argv);// Vao MSMC o day

            // Save electron final parameter sau khi chay xong multi_subbands_MC() khong ?
            //void save_electron_parameters(char *fn);// Dat o dau thi save o day
            //char fnFinal[100];
            //sprintf(fnFinal,"ElectronFinalParameters-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //save_electron_parameters(fnFinal);
           
            //void save_electron_distribution(char *fnSubband,char *fnEnergy);
            // char fnSubbandFinal[100], fnEnergyFinal[100];
            //sprintf(fnSubbandFinal,"ElectronFinalDistributionSubbands-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //sprintf(fnEnergyFinal,"ElectronFinalDistributionEnergy-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
            //save_electron_distribution(fnSubbandFinal,fnEnergyFinal);
           
       
       current_calculation(&cur_av); // Tinh current cho cap (Vd va Vg) nay
             //if((V_drain_input >= 0)&&(cur_av < 0.0)) { cur_av = 0.0;}// (May 24 2010 KHONG Dung lenh nay voi de kiem tra xem the nao )

       if(myrank ==0){ // ghi file I-V vao 1 node
	 fprintf(f_IdVd,"%f   %f   %le \n",GetGateVoltage(),GetDrainVoltage(),cur_av);
       }
       
       if(myrank ==0){ 
	 logging(logfile,"\n    The results for the bias:");
	 logging(logfile,"\n    Vd (V) = %f", GetDrainVoltage());
	 logging(logfile,"\n    Vg (V) = %f", GetGateVoltage());
	 logging(logfile,"\n    Id (A) = %le\n",cur_av);

	 printf("\n The results for the bias:");
	 printf("\n Vd (V) = %f", GetDrainVoltage());
	 printf("\n Vg (V) = %f", GetGateVoltage());
	 printf("\n Id (A) = %le\n",cur_av);// De dua ra micro A thi NHAN voi 1.0e+6
       }

       // For checking
       //save_potential(0);// 0:xy, yz, xz;//save_potential(1);// 1:xyz save the nay la save cho tung cap Vd va Vg day
       //save_potential(2);// ve theo x only
      if(myrank ==0){ 
	time_toc = MPI_Wtime(); //  void logging(FILE *log,const char *fmt,...);
	//logging(logfile,"\n Vg =%f, Vd = %f",GetGateVoltage(),GetDrainVoltage());
	//logging(logfile,"\n    Time to run  =%f [s]",time_toc - time_tic); 
	double time_to_run = time_toc - time_tic;
	logging(logfile,"\n    Time to simulate %d hours %d minutes %d seconds ",
                             (int)(time_to_run/3600), (int)(((int)time_to_run%3600)/60), ((int)time_to_run%3600)%60);
	logging(logfile,"\n\n *****************************************************");
      }
       

    }//End of Single bias for(V_drain_input=Vd_start;V_drain_input< n

    if(myrank ==0){ // ghi file I-V vao 1 node
      fprintf(f_IdVg,"%f  %f  %le \n",GetDrainVoltage(),GetGateVoltage(),cur_av);
    }

}// END for ALL BIASES (Gate voltage bias)

  
    fclose(f_IdVd);
    fclose(f_IdVg);
      
  // Tinh thoi gian ket thuc chuong trinh
    time_finish=time(NULL);
    //printf("\n The finish time is %s",asctime(localtime(&time_finish)));
    double time_to_simulate = difftime(time_finish,time_start);
    if(myrank ==0){ 
        printf("\n Time to simulate %d hours %d minutes %d seconds \n", 
           (int)(time_to_simulate/3600), (int)(((int)time_to_simulate%3600)/60), ((int)time_to_simulate%3600)%60);
	printf("\n The end !");
       logging(logfile,"\n     The end !");
    }
    fclose(logfile);
    return 0;
}
