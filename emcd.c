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
        
Starting date: March 17, 2010
Latest update: March 19, 2010 
***************************************************************************** */
#include <stdio.h>
#include <math.h>

void emcd(){ //Ensemble Monte Carlo for the Device 
    // Goi cac ham
    double **Get_p(),*Get_energy(); // THAM SO cua HAT
    int *Get_valley(),*Get_subband();// THAM SO cua HAT
    int Get_n_used(); // so hat dang su dung
    double Get_mesh_size_x();
    long Get_idum();
    float random2(long *idum); 
    double Get_dt(), *Get_max_gm(); // max_gm: Gamma_max cho tung SECTION
    int count_used_particles(); // Tinh SO HAT dang SU DUNG
    int Get_flag_ballistic_transport(); 
    
    void Set_iss_out(int out_p),Set_iss_eli(int eli_p), Set_iss_cre(int cre_p);
    void Set_idd_out(int out_p),Set_idd_eli(int eli_p), Set_idd_cre(int cre_p);
    void Set_kx(double k_momen_x),Set_dtau(double flight_time),Set_x_position(double position_in_x);
    void Set_electron_energy(double ener),Set_iv(int valley_pair_index),Set_subb(int subband_index);
    double Get_kx(),Get_dtau(),Get_x_position(),Get_electron_energy();
    int Get_iv(),Get_subb();
    
    //void Set_particle_i_th( int i); //Set (Lay) thu tu hat dang simulate For checking
      
    void drift(double tau);
    void scattering();
    
    // Cac bien local
    double **p     = Get_p();
    double *energy = Get_energy();
    int *valley = Get_valley();
    int *subband = Get_subband();
    double mesh_size_x = Get_mesh_size_x();
    long idum  = Get_idum();
    double dt = Get_dt();// Time step: observation time Ta TU DINH NGHIA doc tu Input.d
    double *max_gm = Get_max_gm();
    int flag_ballistic_transport;
    
// Thuc hien
    int i, iv_check;
    double dte=0.0,dt2=0.0,dte2=0.0,rr=0.0,dt3=0.0,dtp=0.0,tau=0.0;
    double x_update=0.0; // position update
    int s_update=0;// section update
    
    flag_ballistic_transport = Get_flag_ballistic_transport();//=1 la BALLISTIC
    
    // Reset the electron number
    Set_iss_out(0);// Nghia la iss_out = 0;
    Set_iss_eli(0);
    Set_iss_cre(0);
    Set_idd_out(0);
    Set_idd_eli(0);
    Set_idd_cre(0);//printf("\n idd_out=%d, idd_eli=%d, idd_cre=%d",Get_idd_out(),Get_idd_eli(),Get_idd_cre()); getchar();
    
    int n_used = Get_n_used(); // int ne = count_used_particles(); // Tim so hat dang su dung la bao nhieu
    //printf("\n Truoc khi chay emcd thuc su n_used=%d, ne= %d",n_used,ne);getchar();// Confirm 2 kieu la giong nhau
     
    for(i=1; i<=n_used; i++){ // Calculate for each particle (loop for all used particles);//for(i=1; i<=10; i++){  // Thu cho 10 hat thoi
                
        // Inverse mapping of particle atributes
        Set_kx(p[i][1]);         // kx = p[i][1]
        Set_x_position(p[i][2]); // x_position = p[i][2]
        Set_dtau(p[i][3]);       // dtau = p[i][3]
        // i_region = p[n][4]: Hien tai khong dung tham so nay
        Set_iv(valley[i]);       // iv = valley[i]
        Set_electron_energy(energy[i]);// electron_energy = energy[i]
        Set_subb(subband[i]);    // subb= subband[i]
        
        // Checking iv Neu ==9 thi HAT nay bi LOAI khong can thuc hien nua
        iv_check = Get_iv();
        if(iv_check==9) { goto L403;}// Bat dau 1 vong moi cua chu trinh ma khong thuc hien tiep cac lenh dang sau
        
        // BALLISTIC TRANSPORT or DIFFUSIVE TRANSPORT ?
        if( flag_ballistic_transport==1)//ballistic
           {
             drift(dt);// ballistic nen trong khoang thoi gian dt no chi drift thoi
             iv_check = Get_iv();// Do trong qua trinh Drift no CO THE THAY DOI iv
             if(iv_check==9) { goto L403;}
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
      
        while(dte <= dt){
             // Free-flight and scatter part
             dte2 = dte; // Da xet drift o tren roi

             // Lay thu tu hat dang simulate ma se co scattering
             //Set_particle_i_th(i); // Hat thu i dang simulate se co scattering

             scattering(); // Goi ham scattering
             
             do {
                  rr=random2(&idum);
                }
             while ((rr<=1.0e-6)||(rr>=1.0));
             
             // Do no da drift 1 doan nen co the vi tri cua HAT da sang SECTION KHAC 
             x_update = Get_x_position();
             s_update = (int)(round(x_update/mesh_size_x));
             dt3 = -log(rr)/max_gm[s_update]; // [s]
             dtp = dt - dte2;	// remaining time to scatter in dt-interval
         
             if(dt3<=dtp) {dt2 = dt3;}
             else {dt2 = dtp;}
      
             drift(dt2);
             //printf("\n Drift SAU KHI goi scattering() of particle %d with drift time Tau= %le",i,dt2);
             
             iv_check = Get_iv();// Do trong qua trinh Drift no CO THE THAY DOI iv
             if(iv_check==9) { goto L403;}
             	   
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
     p[i][1] = Get_kx();
     p[i][2] = Get_x_position();
     p[i][3] = Get_dtau();
     valley[i] = Get_iv(); 
     energy[i] = Get_electron_energy();   
     subband[i] = Get_subb();
   
    } // End of for(i=1; i<=n_used; i++){ // Calculate for each particle
  return;  
} // End of void emcd()
