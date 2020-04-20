/* *****************************************************************************
De giai MSMC
Starting date: April 3, 2010
Update:        April 5, 2010
Latest update: May 16, 2010

****************************************************************************** */
#include <stdio.h>
#include "functions.h"
#include "nanowire.h"
#include "petscksp.h"
#include "region.h"

void multi_subbands_MC(FILE *logfile, int *argc, char ***argv){
    // Goi cac ham
    void emcd();
    void check_source_drain_contacts();
    void delete_particles();
    void electron_density_caculation();
    void eigen_energy_correction(int StepsUsedCorrection);
    void Solved_2D_Schro_for_MSMC();
    void Solved_2D_Schro_Parallel_for_MSMC(int *argc, char ***argv);
    void form_factor_calculation();
    void scattering_table();
    void save_results();
    void normalize_table();
    int Get_iss_out(),Get_iss_eli(),Get_iss_cre(),Get_idd_out(),Get_idd_eli(),Get_idd_cre();
    double Get_dt(),Get_tot_time(),Get_transient_time(); // Xac dinh tu Input.d
    int Get_NTimeSteps();
    void mpi_gather_potential_at_node0(); 
    void mpi_distribute_potential_from_node0();
    void save_potential();
    
    // Cac bien local
    double dt = Get_dt(); // Observation time step
    double tot_time = Get_tot_time(); // Tong thoi gian chay Simulation
    double transient_time = Get_transient_time(); // Chon bang BAO NHIEU can KINH NGHIEM
    int NTimeSteps = Get_NTimeSteps();
    
    int Get_nx1();
    int NumSections = Get_nx1();

// Thuc hien
    int StepsUsedCorrection=0; // tai buoc j_th using time independent perturbation for eigen energy correction

    // Buoc 1.1. Lay rank, size va  Mo 2 file de tinh current
       int myrank, mysize;
       MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
       MPI_Comm_size(MPI_COMM_WORLD,&mysize);

       FILE *ff, *ff_transient; 

       if(myrank ==0){// Chi mo cac files TAI 1 NODE thoi
	 ff=fopen("charge_source_drain.dat","w"); // Chu y dieu kien chon la "w" day !
	 if(ff==NULL){printf("Error: can't open file charge_source_drain.dat.\n");return;}
	 // Thu tu time | source_factor | drain_factor cai nay la tinh tu thoi diem dau tien tro di 
	 // Ta khong nen ghi comment vao file nay vi file nay se duoc doc o ham current_calculation() 
	 ff_transient=fopen("charge_source_drain_transient.dat","w"); // Chu y dieu kien chon la "w" day !
	 if(ff_transient==NULL){printf("Error: can't open file charge_source_drain_transient.dat\n");return;}
	 // Thu tu time | source_factor | drain_factor cai nay la tinh tu thoi diem transient tro di 
       }

    // Buoc 1.2. Mo cac files de dung trong vong lap while
       // Cho ham void velocity_energy_cumulative
       // Cho ham void save_VeloEnerCurrent_all_time_steps(FILE *fname....
       FILE *fAllTimeSteps, *fAfterTransient, *fPotetialAveragex,*fPotentialAvewrageyz;

       if(myrank ==0){
	 fAllTimeSteps=fopen("Velo_ener_current_averages_all_time_steps.dat","w"); // Chu y kieu mo la "a" co the GHI THEM
	 if(fAllTimeSteps==NULL){ printf("\n Cannot open file Velo_ener_current_averages_all_time_steps.dat");return;} 
	 fprintf(fAllTimeSteps,"\n   #x[nm]    velocity_x[m/s]  energy[eV]   current\n");

	 // Cho ham void save_VeloEnerCurrent_after_transient(FILE *fname.....
	 //FILE *fAfterTransient;
	 fAfterTransient=fopen("Velo_ener_current_averages_after_transient.dat","w"); // Chu y kieu mo la "a"
	 if(fAfterTransient==NULL){ printf("\n Cannot open file Velo_ener_current_averages_after_transient.dat");return;} 
	 fprintf(fAfterTransient,"\n #x[nm] velocity_x[m/s] energy[eV] current\n");

	 // Cho ham void save_potential_average(FILE *f1, FILE *f2,double time, int iter_reference)
	 //FILE *fPotetialAveragex;
	 fPotetialAveragex = fopen("average_potential_along_x.dat","w");// Do la "w" nen neu chay nhieu cap (Vd,Vg) thi no la CAP CUOI CUNG               
	 if(fPotetialAveragex == NULL){printf("\n Cannot open file average_potential_along_x.dat");return; }
	 fprintf(fPotetialAveragex,"\n #x[nm] Pot_x[eV] \n");

	 //FILE *fPotentialAvewrageyz;
	 fPotentialAvewrageyz = fopen("average_potential_yz_plane.dat","w");                  
	 if(fPotentialAvewrageyz == NULL){printf("\n Cannot open file average_potential_yz_plane.dat");return; }
	 fprintf(fPotentialAvewrageyz,"\n #y[nm] z[nm] Pot_avg_yz[eV] \n");
       }
      
    // Buoc 2. Khoi tao cac bien dung trong vong lap while
        int issOut=0, issEli=0, issCre=0, iddOut=0, iddEli=0, iddCre=0;
        double drain_factor=0.0, source_factor=0.0;
	double time = 0.0; // thoi gian chay cua hat

	int iter_total=(int)(tot_time/dt+0.5); 

	if(myrank ==0){
	  logging(logfile,"\n    Total number of iterations = %d",iter_total); 
	  printf("\n Total number of iterations = %d",iter_total);
	}

	// Ta khong the save tat ca duoc. Chi co the save tai 1 vai iteration hoac o last interation ma thoi.
	int TimeStepsToSave = iter_total; // Save tai last iteration ma thoi
    
     // Buoc 3. Start MSMC phase
	int iter_reference=1;
	int j_iter=0;

   while(j_iter < iter_total){ // Chay tung lan cho tat ca cac iterations. Neu ma co dau = thi so luong j_iter la total number of interation +1
        
         j_iter = j_iter + 1;// Lenh nay quan trong lam day
          
         time=dt*(double)(j_iter); // thoi gian o interation thu j_iter la bao nhieu
          
         emcd();// EMC: The function drift() and scattering() are used trong 1 khoang dt for all particles. Can SONG SONG
                  
         check_source_drain_contacts();// To check the charge neutrality of the source and drain ohmic contacts
                                       //call initial_kspace and initial_realspace if needed. KHONG can SONG SONG
         delete_particles(); // Neu can.  KHONG can SONG SONG
                         
         electron_density_caculation();// la tu MSMC. KHONG can SONG SONG
	      // Save electron density de kiem tra khong ?
	      if(myrank ==0){// Chi  TAI 1 NODE thoi
		if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
		  char fnElectronDensity[100];
		  sprintf(fnElectronDensity,"Electron_Density-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
		  save_electron_density(fnElectronDensity);

		  char fnAllSections[100];
		  sprintf(fnAllSections,"Electron_Density_All_Sections-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
		  save_electron_density_all_sections(fnAllSections);// For Checking only
	      	}
	      }

	 //PHAN GIAI 3D POISSON
	 get_charge_density_for_Poisson();// chuyen charge density tu MSMC sang Poisson. No chuyen den TAT CA cac node
	      // Co save charge dua vao Poisson tai cac Interim khong
	      //save_charge_density_for_Poisson();

	 //logfile_PoissonInput(); // For checking
	 
	 if (GetPoissonLinearized()) { SolvePoisson = solve_3d_poisson_linearized;}
	 else                        { SolvePoisson = solve_3d_poisson;}

	 (*SolvePoisson)();// Chon kieu gi o tren thi solve Poisson theo kieu do 

	 //void mpi_gather_potential_at_node0();// Sau khi giai 3D Poisson xong thi can gather chu potential lai.RAT NOTE day
         mpi_gather_potential_at_node0();     // Neu chon Print_Midline_Potential(y/n) = y hoac print_output(); thi no da gather lai cho minh roi.
	 mpi_distribute_potential_from_node0(myrank, mysize); 
           
        // KET THUC giai 3D Poisson
	
         if(myrank ==0){// Chi  TAI 1 NODE thoi
	    if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
	      print_output();// Chi de cuoi cung thoi. Neu de o day se gay loi sau khoang 502 vong lap
	                     //Nhung do GS con gather potential roi moi lam nen. Chua hieu lam ???
              save_potential();//printf("\n Dung tai  Solved_2D_Schro_for_MSMC() lan thu 2");
	    }
	  }
	

         iter_reference=1; // Ban dau gan iter_reference =1 sau do neu j_iter=iter_total
                           //(tuc bien chay j-iter di den diem cuoi cung)thi ta cho iter_reference=0; 
         if(j_iter == iter_total) { 
	   iter_reference=0;
	 }
         
         // Tinh cac gia tri trung binh
	 save_potential_average(fPotetialAveragex, fPotentialAvewrageyz, time,iter_reference);//Trong source code minh tinh toan va chay cho 1 node

	 if(myrank ==0){// Chi  TAI 1 NODE thoi

	   velocity_energy_cumulative(fAllTimeSteps,fAfterTransient, time,iter_reference);
		     
	   // Lay cac gia tri cua HAT VAO (CREATION) va RA (OUT or ELIMINATE)
	   issOut = Get_iss_out();
	   issEli = Get_iss_eli();
	   issCre = Get_iss_cre();
	   iddOut = Get_idd_out();
	   iddEli = Get_idd_eli();
	   iddCre = Get_idd_cre();
	   source_factor =(double)(issOut+issEli-issCre); // tong so hat vao ra
	   drain_factor = (double)(iddOut+iddEli-iddCre); // tong so hat vao ra
       	   fprintf(ff,"%le %le %le \n",time,source_factor,drain_factor);//Ghi vao file charge_source_drain.dat

	   if(time >transient_time){// Ghi vao file charge_source_drain_transient.dat
	     fprintf(ff_transient,"%le %le %le \n",time,source_factor,drain_factor);
	   }
         }// End of  if(myrank ==0)


         // NEED to solve 2D Scchrodinger or NOT
         if(j_iter % NTimeSteps !=0){// chi can dung phuong phap correction
                StepsUsedCorrection = j_iter;   
                eigen_energy_correction(StepsUsedCorrection);// First order correction for energy using time-independent pertubation
	 }
         else { // j_iter % NTimeSteps ==0 can giai 2D Schrodinger.

	   if(myrank ==0){ 
	     printf("\n Solved 2D Schrodinger at interation =%d",j_iter);
	   }// End of if(myrank ==0)
		
	   //Solved_2D_Schro_Parallel_for_MSMC(argc,argv);
           Solved_2D_Schro_for_MSMC(); // LAN TIEP
	          //if(myrank ==0){
		     //Save potential cai su dung de giai 2D Schrodinger
		     //if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
		          //char fn_x_Interim[100], fn_yz_Interim[100];
		          //sprintf(fn_x_Interim,"potential_x_at_midpoint_yz_in2DSchrodingerInterim-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
		          //sprintf(fn_yz_Interim,"potential_yz_at_midpoint_x_in2DSchrodingerInterim-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
		         //save_potential_used_for_2D_Schrodinger(fn_x_Interim,fn_yz_Interim);
		     

		     //Save subband_energy 
		     //char fnsubbandInterim[100];
		     //sprintf(fnsubbandInterim,"subbandInterim-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
		     //save_subband_energy(fnsubbandInterim); 
		     // }// End of if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
	           // }

                form_factor_calculation();// LAN TIEP
		     // Co save form factor cac lan TRUNG GIAN de Kiem tra khong ?
		     if(myrank ==0){// Chi  TAI 1 NODE thoi
		       if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
			 char fnformfactor[100];
			 sprintf(fnformfactor,"FormfactorLastInteration-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
			 int FormfactoSavedAtsth = (int)(NumSections/2 +0.5);// lay o giua
			 save_form_factor_calculation(FormfactoSavedAtsth, fnformfactor);
		       }
		     }


                scattering_table(); //LAN TIEP
		     // Save scattering cho nhung lan trung gian de kiem tra khong ?
		     if(myrank ==0){// Chi  TAI 1 NODE thoi
		       if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
			 char fnscat[100];
			 sprintf(fnscat,"ScatTableLastIteration-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
			 int ScatSavedAtsth = (int)(NumSections/2 +0.5);// lay o giua
			 save_scattering_table(ScatSavedAtsth,fnscat);
		       }
		     }
		     
	        normalize_table(); //LAN TIEP de quay lay emcd()

	 }//End of else { // j_iter % NTimeSteps ==0 can giai 2D Schrodinger
        
       if(myrank ==0){// Chi  TAI 1 NODE thoi
	 if(j_iter % TimeStepsToSave ==0){// Save. TimeStepsToSave = iter_total nen save tai last
	   void save_electron_parameters(char *fn);// Dat o dau thi save o day
	   char fn[100];
	   sprintf(fn,"ElectronParametersLastIteration-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	   save_electron_parameters(fn);
      
	   void save_electron_distribution(char *fnSubband,char *fnEnergy);
	   char fnSubbandInter[100], fnEnergyInter[100];
	   sprintf(fnSubbandInter,"ElectronDistributionSubbandsLastIteration-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	   sprintf(fnEnergyInter,"ElectronDistributionEnergyLastIteration-Vd%0.2f-Vg%0.2f.dat",GetDrainVoltage(),GetGateVoltage());
	   save_electron_distribution(fnSubbandInter,fnEnergyInter);
	 }
       }

  }  // End of the while(j_iter < iter_total) Time Loop

  fclose(ff); // Close file ff= ("charge_source_drain.dat","w"); 
  fclose(ff_transient);  // Close file ff_transient= ("charge_source_drain_transient.dat","w"); 
  fclose(fAllTimeSteps);
  fclose(fAfterTransient);
  fclose(fPotetialAveragex);
  fclose(fPotentialAvewrageyz);
  
}// End of multi_subbands_MC()
