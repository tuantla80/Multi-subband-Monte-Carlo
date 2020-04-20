/* *****************************************************************************
Ensemble Monte Carlo (EMC) algorithm for a device.
 - based on the succesive and simultaneous motions of many particles during a small time increment dt
 - The function drift() and scattering() are used in EMC during dt for all particles.
 -    t-------------t+dt ; ts: free-flight time is determined by a random number.
        + (1) if (ts>t+dt) means t---------(t+dt)-----ts -> the particle drifts during dt 
                                                         -> call drift(tau) where tau=(t+dt)-t
        + (2) if (ts<=t+dt) means t--------(ts)------t+dt-> the particle first drifts during (ts-t) 
        and then scattered at ts -> call drift(tau=ts-t) and then scattering().
        => Then generates a new scattering time ts (determined by another random number) 
        and we check again ts is larger than t+dt or not -> comes to (1) or (2) -> continue to do
        so to the end of the simulation   
        
Starting date: June 20, 2010
Latest update: June 20, 2010 
***************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "mpi.h"

void emcd_parallel(){ //Ensemble Monte Carlo for the Device 
      
    void drift(double tau); //void Set_particle_i_th( int i); //Set (Lay) thu tu hat dang simulate For checking
    void scattering();
    
    // Cac bien local
    double **Get_p();
    double **p     = Get_p();

    double *Get_energy(); 
    double *energy = Get_energy();

    int *Get_valley(),*Get_subband();// THAM SO cua HAT
    int *valley = Get_valley();
    int *subband = Get_subband();

    double Get_mesh_size_x();
    double mesh_size_x = Get_mesh_size_x();

    int Get_nx0(),Get_nx1();
    int nx0 = Get_nx0(); int nx1 = Get_nx1(); 
    int nx_max = nx1 - nx0;

    float random2(long *idum); 
    long Get_idum();
    long idum  = Get_idum();

    double Get_dt(), *Get_max_gm(); // max_gm: Gamma_max cho tung SECTION
    double dt = Get_dt();// Time step: observation time Ta TU DINH NGHIA doc tu Input.d
    double *max_gm = Get_max_gm();
   
    
// Thuc hien
    int i, iv_check, ix_position;
    double dte=0.0, dt2=0.0, dte2=0.0, rr=0.0, dt3=0.0, dtp=0.0, tau=0.0;
    double x_update=0.0; // position update
    int s_update=0;// section update

    int Get_flag_ballistic_transport();
    int flag_ballistic_transport;
    flag_ballistic_transport = Get_flag_ballistic_transport();//=1 la BALLISTIC
    
    // Reset the electron number
    void Set_iss_out(int out_p),Set_iss_eli(int eli_p), Set_iss_cre(int cre_p);
    void Set_idd_out(int out_p),Set_idd_eli(int eli_p), Set_idd_cre(int cre_p);
    Set_iss_out(0);// Nghia la iss_out = 0;
    Set_iss_eli(0);
    Set_iss_cre(0);
    Set_idd_out(0);
    Set_idd_eli(0);
    Set_idd_cre(0);//printf("\n idd_out=%d, idd_eli=%d, idd_cre=%d",Get_idd_out(),Get_idd_eli(),Get_idd_cre()); getchar();
    
    int Get_n_used(); // so hat dang su dung
    int n_used = Get_n_used();    

    void Set_kx(double k_momen_x),Set_dtau(double flight_time),Set_x_position(double position_in_x);
    void Set_electron_energy(double ener),Set_iv(int valley_pair_index),Set_subb(int subband_index);
    double Get_kx(),Get_dtau(),Get_x_position(),Get_electron_energy();
    int Get_iv(),Get_subb();

    double *ArrayParam = dvector(1,6*n_used);//co 6 tham so can la kx,x,tau,iv,subb,energy. Can khoi tao mang nay sau n_used = Get_n_used()
    for(i=1; i<=6*n_used; i++){ ArrayParam[i] = 0.0;}// Dung cho MPI

    int np; //number of particles per node. Do da khoi tao MPI tu dau chuong trinh nen bien nay cung duoc KHOI TAO TAI TAT CA cac nodes

    int myrank=0, mysize=1;// Khoi tao thoi
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if((n_used % mysize) !=0){//Neu SO HAT KHONG chia het cho SO Processors thi (mysize-1) processors se co so hat bang nhau va BANG np
                              // con last processor se co SO HAT = n_used - np*(mysize-1);
      np = (int)(n_used/(mysize-1));
      if(myrank==(mysize-1)){//last processor
	np = n_used - np*(mysize-1);
      }
    }
    else{
      np = (int)(n_used/(mysize));// Neu SO HAT chia het cho SO Processor thi np o cac node la giong nhau
    }
    printf("\n np = %d", np);//Checking only

    MPI_Barrier(MPI_COMM_WORLD);//NOTE Ta can Blocks until all processes in the communicator have reached this routine
    int count =0; // Bien dem de tinh sp phan tu cho ArrayParam

    for(i=1+myrank*np; i<=(myrank+1)*np; i++){// Chay tu 1 den n_used nhung o cac processors khac nhau. Chu y ung voi moi rank KHAC NHAU thi np CO THE khac nhau
                                              // for(i=1; i<=n_used; i++){
        if(i > n_used){// Thuc ra thi ko xay ra nhung van kiem tra
	  goto Label_End;
	}

      // Inverse mapping of particle atributes
        Set_kx(p[i][1]);         // kx = p[i][1]
        Set_x_position(p[i][2]); // x_position = p[i][2]
        Set_dtau(p[i][3]);       // dtau = p[i][3]
        // i_region = p[n][4]: Hien tai khong dung tham so nay
        Set_iv(valley[i]);       // iv = valley[i]
        Set_electron_energy(energy[i]);// electron_energy = energy[i]
        Set_subb(subband[i]);    // subb= subband[i]
        
        // Checking, postion va  iv Neu ==9 thi HAT nay bi LOAI khong can thuc hien nua
	ix_position = (int)(Get_x_position()/mesh_size_x + 0.5);
	if(ix_position <0 ){
	  ix_position = 0;
	}
	if(ix_position > nx_max){
	  ix_position = nx_max;
	}

        iv_check = Get_iv();
        if(iv_check==9) {
	  goto L403;// Bat dau 1 vong moi cua chu trinh ma khong thuc hien tiep cac lenh dang sau
	}
        
        // BALLISTIC TRANSPORT or DIFFUSIVE TRANSPORT ?
        if( flag_ballistic_transport==1)//ballistic
           {
             drift(dt);// ballistic nen trong khoang thoi gian dt no chi drift thoi
             iv_check = Get_iv();// Do trong qua trinh Drift no CO THE THAY DOI iv
             if(iv_check==9){
	       goto L403;
	     }
             goto L403; 
           }
        // Neu khong thi no se la 
        // DIFFUSIVE TRANSPORT case
        tau = Get_dtau(); //printf("\n Drift DAU TIEN of particle %d with drift time Tau= %le while dt=%le",i,tau,dt); getchar();
        dte = tau;
        if(dte >= dt)
            { 
              dt2=dt; 
            }
        else{ // dte < dt
              dt2=dte;
            }
        // -> dt2 is a drift time
        drift(dt2); // Call drift() function during dt2
        
        iv_check = Get_iv();// Do trong qua trinh Drift no CO THE THAY DOI iv
        if(iv_check==9) { goto L403;}
      
        while(dte < dt){// do = da xet o tren 
             // Free-flight and scatter part
             dte2 = dte; // Da xet drift o tren roi

             // Lay thu tu hat dang simulate ma se co scattering
             //Set_particle_i_th(i); // Hat thu i dang simulate se co scattering

             scattering(); // Goi ham scattering. KHONG SONG SONG
             
             do {
                  rr=random2(&idum);
                }
             while ((rr<=1.0e-6)||(rr>=1.0));
             
             // Do no da drift 1 doan nen co the vi tri cua HAT da sang SECTION KHAC 
             x_update = Get_x_position();
             s_update = (int)(round(x_update/mesh_size_x));
             dt3 = -log(rr)/max_gm[s_update]; // [s]

             dtp = dt - dte2;	// remaining time to scatter in dt-interval
         
             if(dt3<=dtp) {
	       dt2 = dt3;
	     }
             else {
	       dt2 = dtp;
	     }
      
             drift(dt2); //printf("\n Drift SAU KHI goi scattering() of particle %d with drift time Tau= %le",i,dt2);
             
             iv_check = Get_iv();// Do trong qua trinh Drift no CO THE THAY DOI iv
             if(iv_check==9) {
	       goto L403;
	     }
             	   
             // Update times
             dte2 = dte2 + dt3;
             dte = dte2;
         } // End of the "while" loop
           
        // Meaning dte >= dt  // after "while"  loops
       dte = dte - dt;
       tau = dte;
       Set_dtau(tau); // tau la gia tri cho khoang thoi gian drift moi
       
   L403:
        /* Note: khong the dung  lenh continue vi continue trong C/C++ va trong
          Fortran co mot vai diem khac nhau.
          - Trong C/C++ neu gap continue trong vong for thi no se quay tro lai vong lap for
            voi chi so ke tiep ma khong lam cac lenh sau do
          - Trong Fortran thi lenh continue la "chi dan rang" lam cac lenh dang sau continue
          Do do neu o day dung lenh continue thi no se khong lam cac lenh trong doan
          Map particle atributes -> cac gia tri iss_out, idd_eli, iss_cre, idd_out
          idd_eli, idd_cre se bang 0 het -> sai */ 
  
  // Map particle atributes

       count = count +1;// Do Ta simulate cho tung hat, va tung hat lai co thu tu tu kx,x,tau,electron_energy,iv,subb
       ArrayParam[count] = Get_kx(); // p[i][1] = Get_kx();

       count = count +1;
       ArrayParam[count] = Get_x_position();// p[i][2] = Get_x_position();

       count = count +1;
       ArrayParam[count] = Get_dtau();// p[i][3] = Get_dtau();
       
       count = count +1;
       ArrayParam[count] = Get_electron_energy();// energy[i] = Get_electron_energy(); 

       count = count +1;
       ArrayParam[count] = Get_iv();// valley[i] = Get_iv(); 

       count = count +1;
       ArrayParam[count] = Get_subb();// subband[i] = Get_subb();
       //Nhu vay count se di tu 1 den 6*n_used

    Label_End:
    } // End of for(i=1; i<=n_used; i++){ // Calculate for each particle

    //Sau khi simulate cac hat xong thi o tat ca cac processor ta co cac mang ArrayParam: Du lieu phan bo RAI RAC.
    // Can TAP HOP lai va de du lieu cau ArrayParam la GIONG NHAU o TAT CA CAC Processors
    Can them va tiep


  return;  
} // End of void emcd()
