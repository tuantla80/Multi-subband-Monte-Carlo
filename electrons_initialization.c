/* ***********************************************************************
 - To calculate the number of initial electrons
 - Co 2 cach
      + Cach 1: Don gian. Nhung HAN CHE la TAO RA SO HAT RAT LON
      + Cach 2: neu guest potential ma ok thi se dung 
       
 Initialises parameters for each particle (electron)
   - Initial electron energy 
   - Initial position and momentum of electrons
   - Initial valley and subband where the particle is residing for the first time
  
Starting date: Feb 19, 2010
Update:        Feb 19, 2010
Update:        May 05, 2010 (Giai 3D Poisson at equilibrium roi sau do moi khoi tao electron)
Latest update: June 02, 2010: Cach moi
************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "nrutil.h"
#include "constants.h"

// ******************* Cach 2: Dung khi guest potential la NGON *************** 
void electrons_initialization(){
     // Goi cac ham
     void init_kspace(int ne,int i);
     void init_realspace(int ne,int i);
     void Set_n_used( int n);
     int Get_n_used();
     int find_region(int i,int j,int k);
     
     double Get_artificial_factor_cell_volume();
     double artificial_factor_cell_volume = Get_artificial_factor_cell_volume();

     int Get_save_electron_initlization();
     int save_electron_initlization_or_not = Get_save_electron_initlization();

     double Get_mesh_size_x(), Get_mesh_size_y(), Get_mesh_size_z(); 
     double mesh_size_x = Get_mesh_size_x();
     double mesh_size_y = Get_mesh_size_y();
     double mesh_size_z = Get_mesh_size_z();

     int Get_nx0(),Get_nxa(),Get_nxb(),Get_nx1(), Get_nya(),Get_nyb(),Get_nza(),Get_nzb();
     int nx0 = Get_nx0(); int nxa = Get_nxa(); int nxb = Get_nxb(); int nx1 = Get_nx1(); 
     int nya = Get_nya(); int nyb = Get_nyb();
     int nza = Get_nza(); int nzb = Get_nzb();

     double ***Phi = GetPot();

     double  Get_intrinsic_carrier_density(),Get_Vt();
     double intrinsic_carrier_density = Get_intrinsic_carrier_density();
     double Vt = Get_Vt(); 

     int *Get_valley();
     int *valley   = Get_valley();
     
// Thuc hien
     int i,j,k;          
     double cell_volume = mesh_size_x*mesh_size_y*mesh_size_z; //volume of a cell
     artificial_factor_cell_volume = artificial_factor_cell_volume/mesh_size_z;// Do ngay xua ta lay device width la 1.0e-6 m
     cell_volume = cell_volume*artificial_factor_cell_volume; // Neu khong thi cell_volume qua nho dan den so hat trong cell luon bang 0
      
     int np_ijk = 0; // number of electrons in a cell(i,j,k)
     double denn = 0.0;
     int ne = 0; // Number of electrons which we initialize
     
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
         
     FILE *f;
     if(myrank ==0){ 
       if(save_electron_initlization_or_not==1){ // Co save
	 f = fopen("electron_init.dat","w"); // "w" neu tep ton tai no se bi xoa
	 if(f==NULL){printf("\n Cannot open file electron_init.dat"); return ; }
	 fprintf(f,"\n #i j k  Electron_number \n");
       }
     }

     double Fermihalf (double x);
     double Nc_Si = 3.22*1.0e+25;//[m3] // Tai room temperature Nc cua Silicon la 3.22e+19/cm3.Semiconductor Device Fundamental - F.Pierret (page 51)
                                         //Can GIONG VOI void get_initial_charge_density_for_Poisson()
     double Ni = Nc_Si;

     double GetDrainVoltage();
     double Vd = GetDrainVoltage(); //printf("\n Vd = %f", Vd); getchar();

     Energy *E;
     E= PGetEnergy();

     int n =0;
     for(i=nx0; i<=nx1; i++){
       for(j=nya; j<=nyb; j++){ // chi doan co the khoi tao thoi. Tuc la LOI Silicon
           for(k=nza; k<=nzb; k++){ 
	      n = find_region(i,j,k);
	      if(n==1){ //S
	       // Boltzmann statistics
	       denn = Ni*exp((-Phi[i][j][k]+E->biS)/Vt);// vi tri toa de cua Phi va n3d la TUONG DUONG NHAU mac du Phi co trai rong index re 1 ti
	       //Fermi-Dirac statistics
	       //denn = Ni*Fermihalf((-Phi[i][j][k]+E->biS)/Vt);	       
	     }

	     else if (n==2){//D
	       //Boltzmann statistics
	       denn = Ni*exp((-Phi[i][j][k]+Vd+ E->biD)/Vt);
	       //Fermi-Dirac statistics
	       //n3d[i][j][k] = Ni*Fermihalf((-Phi[i][j][k]+Vd+ E->biD)/Vt);
	     }

	     else if (n==3){//Ch
	       //Boltzmann statistics
	       denn = Ni*exp((-Phi[i][j][k])/Vt);
	       //Fermi-Dirac statistics
	       //n3d[i][j][k] = Ni*Fermihalf((-Phi[i][j][k])/Vt);
	     }
	    
               if(i==nx0||i==nx1) { denn = denn/2.0; } // at the boundary the denn is only 1/2
               if(j==nya||j==nyb) { denn = denn/2.0; }
               if(k==nza||k==nzb) { denn = denn/2.0; }
           
	           np_ijk = (int)(denn*cell_volume+0.5); // number of electrons in a cell(i,j,k)
	                 // tai vi int(0.7)=0 nen ta CAN cong them 0.5 thi chinh xac hon
		   if(myrank ==0){ 
		     if ((np_ijk >=1)&&(save_electron_initlization_or_not==1)){// thi save
		       fprintf(f,"%d %d %d  %d \n",i,j,k, np_ijk);
		     }
		   }
               if (np_ijk==0) { goto L20; } // comeback to the loops
               //else: np_ij is different 0 -> initialize the parameters for each particle (electron) 
               while( np_ijk > 0){
         	     ne = ne + 1;
                     if(ne > max_electron_number) {
                        printf("\n Number of initial electrons = %d",ne);
			printf("\n Problem at electrons_initialization.c");
			printf("\n You can change max_electron_number in constants.h Xem denn o Tren \n");
			nrerror("Number of particles exceed max_electron_number");//Dung ham thong bao loi trong nr.h
	              }
                
                    init_kspace(ne,i); //initial momentum and energy           
                    init_realspace(ne,i);  // Initial position
             
                    np_ijk = np_ijk - 1;
	       }// End of while( np_ijk > 0)
       L20:
       continue;
         } // End of for(k=nza; k<=nzb; k++)         
      }// End of for(j=nya; j<=nyb; j++){
   }//End of for(i=nx0; i<=nx1; i++) 
   
   Set_n_used(ne); // Tuong duong voi cau lenh int n_used = ne;
   int n_used = Get_n_used();
   if(myrank ==0){ 
     printf("\n Maximum number of electrons allowed = %d",max_electron_number);
     printf("\n Number of electrons initialized     = %d",n_used);
     if (n_used < 5000){// initial co so hat nho hon 5,000 hat. Can them so hat
     	printf("\n In electrons_initialization()");
	printf("\n Need change something");
	nrerror("Number of particles initialized quite small");
     }
   }

   for(i=n_used+1; i<=max_electron_number; i++){
            valley[i] = 9; // inactive these particles, nhung hat ma vuot ngoai gia tri n_used thi phai inactive
   } 
     
   if(myrank ==0){ // close file da mo
     if(save_electron_initlization_or_not==1){ // Co save
       fclose(f);
     }
   }


  return;
}
// Cach 1 xem phan cuoi cua ham
//***********************************************************************************************************
/* ***********************************************************************
De save cac tham so cua electron de kiem tra

Starting date: June 06, 2010
Latest update: June 06, 2010)

NOTE: for checking only. Dat o dau thi kiem tra o day
************************************************************************* */
#include "mpi.h"

void save_electron_parameters(char *fn){

  int Get_n_used();// lay so hat dang su dung
  int n_used = Get_n_used();

  double  **Get_p(),*Get_energy();
  int *Get_valley(),*Get_subband();// lay kx,x,ts, valley, subband va energy cua hat thu i
  double *energy = Get_energy();
  double **p     = Get_p();
  int *valley    = Get_valley();
  int *subband   = Get_subband();
  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  int i;

  FILE *f;
  if(myrank ==0){ 
    f = fopen(fn,"w"); // "w" neu tep ton tai no se bi xoa
    //f = fopen("electron_parameters.dat","w"); // "w" neu tep ton tai no se bi xoa
    if(f==NULL){printf("\n Cannot open file save_electron_parameters.c"); return ; }
    fprintf(f,"\n# ne_thPaticle x_position[nm]       kx[1/m]       FreeFlightTime[fs] Valley Subband  Energy[eV] \n");
    for (i=1; i<=n_used; i++){
      fprintf(f,"        %d        %3.5f          %le    %3.5f              %d      %d       %3.5f \n",i, p[i][2]/1.0e-9, p[i][1], p[i][3]/1.0e-15, valley[i], subband[i], energy[i]);
    }

  }// End of  if(myrank ==0){ 
  return;
}// End of void save_electron_parameters(){
//***************************************************************************************
/* ***********************************************************************
De save cac tham so cua electron de kiem tra
Gom 2 phan

Starting date: June 07, 2010
Latest update: June 07, 2010

NOTE: for checking only. Dat o dau thi kiem tra o day
************************************************************************* */
#include "mpi.h"

void save_electron_distribution(char *fnSubband,char *fnEnergy){

  int Get_n_used();// lay so hat dang su dung
  int n_used = Get_n_used();

  double  **Get_p(),*Get_energy();
  int *Get_valley(),*Get_subband();// lay kx,x,ts, valley, subband va energy cua hat thu i
  double *energy = Get_energy();
  double **p     = Get_p();
  int *valley    = Get_valley();
  int *subband   = Get_subband();

  int Get_NSELECT();
  int NSELECT = Get_NSELECT();

//************************************************************************************
  // Phan 1. Electron distribution theo Valley va Subband
  // Cac bien local
  int i, j;
  int *ElectronDistributionValley1, *ElectronDistributionValley2,*ElectronDistributionValley3;// Tinh so electron o Valley va subband tuong ung
  ElectronDistributionValley1 = ivector(1, NSELECT);
  ElectronDistributionValley2 = ivector(1, NSELECT);
  ElectronDistributionValley3 = ivector(1, NSELECT);
  for(i=1; i<=NSELECT; i++){ ElectronDistributionValley1[i]=0; ElectronDistributionValley2[i]=0;ElectronDistributionValley3[i]=0; }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  FILE *f, *f1;
  // Tinh electron dictribution in each valley and subbands
  if(myrank ==0){ 
    for (i=1; i<=n_used; i++){// Chay cho tat ca cac hat

      if(valley[i]==1){// Valley 1
	for(j=1; j<=NSELECT; j++){
	  if(subband[i]==j){// Chu y la subband tai hat thu i nhe
	    ElectronDistributionValley1[j] += 1; // tang them 1 hat
	  }
	}
      }// end of if(valley[i]==1)

      else if (valley[i]==2){// Valley 2
	for(j=1; j<=NSELECT; j++){
	  if(subband[i]==j){// Chu y la subband tai hat thu i nhe
	    ElectronDistributionValley2[j] += 1; // tang them 1 hat
	  }
	}
      }//end of else if (valley[i]==2)

      else if (valley[i]==3){// Valley 3
	for(j=1; j<=NSELECT; j++){
	  if(subband[i]==j){// Chu y la subband tai hat thu i nhe
	    ElectronDistributionValley3[j] += 1; // tang them 1 hat
	  }
	}
      }//end of else if (valley[i]==3)

    }// End of for (i=1; i<=n_used; i++){
    
    f = fopen(fnSubband,"w"); // "w" neu tep ton tai no se bi xoa
    //f = fopen("electron_distributionSubbands.dat","w"); // "w" neu tep ton tai no se bi xoa
    if(f==NULL){printf("\n Cannot open file at save_electron_distribution for Subband "); return ; }
    fprintf(f,"\n# Subband Valley1 Valley2 Valley3 \n");
    for(j=1; j<=NSELECT; j++){// chay cho NSELECT subband
      fprintf(f," %d        %d %d %d \n",j, ElectronDistributionValley1[j],ElectronDistributionValley2[j],ElectronDistributionValley3[j]);
    }
    fclose(f);
  }// End of  if(myrank ==0){ 
//************************************************************************************
  // KET THUC Phan 1. Electron distribution theo Valley va Subband



//************************************************************************************
  // Phan 2. Electron distribution theo Valley va Energy
  // Tinh electron distribution in energy range
  double Get_emax();
  double emax = Get_emax();
  double de=emax/(double)(n_lev); // energy interval
  double *E; E=dvector(0,n_lev); // [eV] Dinh nghia dai nang luong E di tu 1*de den "emax" 
                                // Neu di tu 0 thi se co loi neu tinh scattering. The hien phan kinetic cua electron
  for(i=0; i<=n_lev; i++){ E[i] = i*de; }
     
  int *ElecDistributionEnergyValley1, *ElecDistributionEnergyValley2,*ElecDistributionEnergyValley3, *ElecDistributionEnergyAllValley;
  ElecDistributionEnergyValley1 = ivector(0,n_lev); 
  ElecDistributionEnergyValley2 = ivector(0,n_lev);
  ElecDistributionEnergyValley3 = ivector(0,n_lev);
  ElecDistributionEnergyAllValley = ivector(0,n_lev);

  for(i=0; i<=n_lev; i++){
    ElecDistributionEnergyValley1[i]=0; ElecDistributionEnergyValley2[i]=0;
    ElecDistributionEnergyValley3[i]=0; ElecDistributionEnergyAllValley[i]=0;
  }
  
   if(myrank ==0){ 

     for (i=1; i<=n_used; i++){// Chay cho tat ca cac hat

       if(valley[i]==1){// Valley 1
	 for(j=0; j<=n_lev; j++){
	   if ( (energy[i]>= E[j]) && (energy[i]<=E[j+1])){ // Chu y chi so cua energy la i cua E la j
	     ElecDistributionEnergyValley1[j] += 1;// Tang them 1 hat
	   }
	 }
       } // end of if(valley[i]==1)

       else if(valley[i]==2){// Valley 2
	 for(j=0; j<=n_lev; j++){
	   if ( (energy[i]>= E[j]) && (energy[i]<=E[j+1])){ // Chu y chi so cua energy la i cua E la j
	     ElecDistributionEnergyValley2[j] += 1;// Tang them 1 hat
	   }
	 }
       } // end of if(valley[i]==2)

       else if(valley[i]==3){// Valley 3
	 for(j=0; j<=n_lev; j++){
	   if ( (energy[i]>= E[j]) && (energy[i]<=E[j+1])){ // Chu y chi so cua energy la i cua E la j
	     ElecDistributionEnergyValley3[j] += 1;// Tang them 1 hat
	   }
	 }
       } // end of if(valley[i]==3)
 
     }// end of  for (i=1; i<=n_used; i++)
     
     // For all valleys
     for (i=1; i<=n_used; i++){// Chay cho tat ca cac hat
      for(j=0; j<=n_lev; j++){
	   if ( (energy[i]>= E[j]) && (energy[i]<=E[j+1])){ // Chu y chi so cua energy la i cua E la j
	     ElecDistributionEnergyAllValley[j] += 1;// Tang them 1 hat
	   }
      }
     }
          
    f1 = fopen(fnEnergy,"w"); // "w" neu tep ton tai no se bi xoa
    //f1 = fopen("electron_distributionEnergyrange.dat","w"); // "w" neu tep ton tai no se bi xoa
    if(f1==NULL){printf("\n Cannot open file Cannot open file at save_electron_distribution for Energy.dat"); return ; }
    fprintf(f1,"\n# EnergyStep Valley1 Valley2 Valley3 \n");
    for(j=1; j<=n_lev; j++){// chay cho tung energy step
      fprintf(f1," %f  %d %d %d %d \n",j*de,ElecDistributionEnergyValley1[j],ElecDistributionEnergyValley2[j],
	                                    ElecDistributionEnergyValley3[j],ElecDistributionEnergyAllValley[j]);
    }
    fclose(f1);

   }// End of  if(myrank ==0){ 

  // Free local variables
   free_ivector(ElectronDistributionValley1, 1, NSELECT);
   free_ivector(ElectronDistributionValley2, 1, NSELECT);
   free_ivector(ElectronDistributionValley3, 1, NSELECT);
   free_dvector(E,0,n_lev);
   free_ivector(ElecDistributionEnergyValley1, 0,n_lev); 
   free_ivector(ElecDistributionEnergyValley2, 0,n_lev);
   free_ivector(ElecDistributionEnergyValley3, 0,n_lev);
   free_ivector(ElecDistributionEnergyAllValley, 0,n_lev);
  
  return;
}//


/*
//   Cach 1: Don gian hien tai dang dung Feb 19, 2010 
void electrons_initialization(){
     // Goi cac ham
     void init_kspace(int ne,int i);
     void init_realspace(int ne,int i); 
     void Set_n_used( int n);
     int Get_n_used();
     int Get_nx_max(), Get_ny_max(), Get_nz_max(),Get_n_s(),Get_n_d();
     int *Get_valley();
     // Cac bien local
     int nx_max    = Get_nx_max();
     int ny_max    = Get_ny_max();
     int nz_max    = Get_nz_max();
     int n_s       = Get_n_s();
     int n_d       = Get_n_d();
     int *valley   = Get_valley();
     // Thuc hien
     int i,j,k;          
     int ne = 0; // Number of electrons which we initialize
     // Tai Source
     for(i=0; i<=n_s; i++){
       for(j=0; j<=ny_max; j++){ 
           for(k=0; k<=nz_max; k++){
               ne = ne+1;
               if(ne > max_electron_number){
                   printf("\n Number of initial electrons = %d",ne);
                   printf("You can change max_electron_number in constants.h file \n");
                   nrerror("Number of particles exceed max_electron_number");
                  }
                //Khoi tao position, valley, subband, momentum, energy cho hat do
                init_kspace(ne,i);      //initial momentum and energy           
                init_realspace(ne,i);  // Initial position
           }
       }
     } 
         
     // Tai Drain
     for(i=n_d+1; i<=nx_max; i++){
       for(j=0; j<=ny_max; j++){ 
           for(k=0; k<=nz_max; k++){
               ne = ne+1;
               if(ne > max_electron_number){
                   printf("\n Number of initial electrons = %d",ne);
                   printf("You can change max_electron_number in constants.h file \n");
                   nrerror("Number of particles exceed max_electron_number");
                  }
                //Khoi tao position, valley, subband, momentum, energy cho hat do
                init_kspace(ne,i);      //initial momentum and energy           
                init_realspace(ne,i);  // Initial position
           }
       }
     } 
     
     Set_n_used(ne); // Tuong duong voi cau lenh int n_used = ne;
     int n_used = Get_n_used();
     printf("\n Number of electrons initialized     = %d",n_used);
     printf("\n Maximum number of electrons allowed = %d",max_electron_number);
     
     for(i=n_used+1; i<=max_electron_number; i++){
            valley[i] = 9; // inactive these particles, nhung hat ma vuot ngoai gia tri n_used thi phai inactive
        }
    
 return;
} // End of void electrons_initialization()
*/

        

