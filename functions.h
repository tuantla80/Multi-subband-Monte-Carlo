
#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_ 
     // Material parameters
     void material_param();
	 double Get_ml();
	 double Get_mt();
	 double Get_constant_acoustic();
	 double Get_constant_optical_g_e();
	 double Get_constant_optical_g_a();
	 double Get_constant_optical_f_e();
	 double Get_constant_optical_f_a(); 
	 double Get_Vt();
	 double Get_intrinsic_carrier_density();
	 double Get_affinity_silicon();
         double Get_gateWF();
         double Get_Eg();
         double Get_BandgapOxide();
         double Get_nonparabolicity_factor();
         double Get_emax();
         double Get_hw0f_phonon();
         double Get_hw0g_phonon();
    
    //Device structure
    void device_structure();
         double Get_mesh_size_x();
         double Get_mesh_size_y();
         double Get_mesh_size_z();
         int Get_Tox_mesh();
         int Get_Tbox_mesh();
         int Get_Wox_mesh();
         int Get_nx0();
         int Get_nxa();
         int Get_nxb();
         int Get_nx1();
         int Get_ny0();
         int Get_nya();
         int Get_nyb();
         int Get_ny1();
         int Get_nz0();
         int Get_nza();
         int Get_nzb();
         int Get_nz1();
         int Get_nxgate0();
         int Get_nxgate1();
         double *Get_doping_density();
    
    // Find region 
    int find_region(int i,int j,int k);
    
    // Initialization for doping and potential
    void initial_potential_doping();
         void save_initial_potential_doping();
         void save_initial_charge_density_for_Poisson(); // Cai lay tu initial potentential
         void save_charge_density_for_Poisson();// Save charge density cai dua vao Poisson tai cac Interim
         void save_potential();// Save potential tai bat cu diem nao dat ham nay
 
    // Khoi tao so hat tai dau Mut cua Source va Drain contact
    void source_drain_carrier_number();
         int Get_nsource_side_carriers();
         int Get_ndrain_side_carriers();
         
    // Read voltages from an input file
    void read_voltages_input();
         double Get_V_source_input();
         double Get_Vd_start();
         double Get_Vd_end();
         double Get_Vd_step();
         double Get_Vg_start();
         double Get_Vg_end();
         double Get_Vg_step();
         void SetDrainVoltage(double v);
         double GetDrainVoltage();
         void SetGateVoltage(double v);
         double GetGateVoltage();

    // Read simulation list from Input.dat 
    void read_simulation_list();
         double Get_dt();
         double Get_tot_time();
         double Get_transient_time();
         double Get_length_plotted();
         double Get_depth_plotted();
         double Get_width_plotted();
         int Get_NSELECT();
         int Get_NTimeSteps();
         double Get_artificial_factor_cell_volume();
         
    // Read scattering and save list     
    void read_scattering_save_list();
        int Get_num_scat_used();
        int Get_flag_ballistic_transport(); 
        int Get_flag_acoustic();
        int Get_flag_zero_order_optical();
        int Get_flag_first_order_optical();
        int Get_flag_Coulomb();
        int Get_flag_Surface_roughness();
        int Get_save_eig_wavefuntion();
        int Get_save_doping_potential_init();
        int Get_save_electron_initlization();
        int Get_save_form_factor();
        int Get_save_scattering_rate();
        int Get_save_electron_population();
        int Get_save_electron_density();
        int Get_save_VelEnerCurr_all_time_steps();
        int Get_save_VelEnerCurr_after_transient();
        int Get_save_init_free_flight();
     
    // Applied voltage
    void apply_voltage();

    // Khoi tai so hat electrons trong device
    void electrons_initialization();
         void save_electron_parameters(char *fn);// Dat o dau thi save o day
         void save_electron_distribution(char *fnSubband,char *fnEnergy); // Dat o dau thi save o day
    // Khoi tao momentum, energy, valley, subband
    void init_kspace(int ne,int i);
    // Khoi tao position
    void init_realspace(int ne,int i);
    // Tinh so hat dang duoc su dung trong device voi valley KHAC 9
    int count_used_particles();
       
    // Ham random
    float random2(long *idum);
        
    // Define functions
    void initial_variables();
    void Set_n_used( int n);
    int Get_n_used();
    double **Get_p();
    double *Get_energy();
    int *Get_valley();
    int *Get_subband();
    long Get_idum();
    double *****Get_wave();
    double ***Get_eig();
    double ***Get_eig_jthSchro();
    double ***Get_doping();
    double ***Get_fai();
    double ***Get_fai_jthSchro();
    double *Get_pot();
    double **Get_trap_weights();
    double ****Get_form_factor();
    double *****Get_scat_table();
    double *Get_max_gm();
    void Set_iss_out(int out_p);
    int Get_iss_out();
    void Set_iss_eli(int eli_p);
    int Get_iss_eli();
    void Set_iss_cre(int cre_p);
    int Get_iss_cre();
    void Set_idd_out(int out_p);
    int Get_idd_out();
    void Set_idd_eli(int eli_p);
    int Get_idd_eli();
    void Set_idd_cre(int cre_p);
    int Get_idd_cre();
    void Set_kx(double k_momen_x);
    double Get_kx();
    void Set_dtau(double flight_time);
    double Get_dtau();
    void Set_x_position(double position_in_x);
    double Get_x_position();
    void Set_electron_energy(double ener);
    double Get_electron_energy();
    void Set_iv(int valley_pair_index);
    int Get_iv();
    void Set_subb(int subband_index);
    int Get_subb();
    void Set_particle_i_th( int i); //26/03/10 15:16
    int Get_particle_i_th();
    double ***Get_electron_density();
    double *Get_velocity_x_sum();
    double *Get_energy_sum();
    double *Get_current_sum();
    double *Get_velocity_x_sum_after_transient();
    double *Get_energy_sum_after_transient();
    double *Get_current_sum_after_transient();
    double *Get_pot_x_avg();
    double **Get_pot_yz_avg();


    // Lay Potential de dua vao 2D Schrodinger
    void get_potential(int s, int ny, int nz);
    
    // Using First order correction for eigen energy 
    void eigen_energy_correction(int StepsUsedCorrection);
    
    //Solved_2D_Schro_for_MSMC()
    void Solved_2D_Schro_Parallel_for_MSMC(int *argc, char ***argv);
    void Solved_2D_Schro_for_MSMC();
        void Schrodinger2D(int s, int v, double m1, double m2);
        void makeh(int ny,int nz,double ly,double lz,double *potential,double **ham,double m1, double m2);
        void diasym(int s, int v,int ny, int nz,double ly,double lz,int NSELECT, double **ham, double ***eig, double *****wave);
        void normalize_wave(int s,int v,int ny,int nz,int NSELECT, double *****wave);
        void save_eig_wave(int s,int v, char *fneigen, char *fnwave);
        double energy0(int ky,int kz,double ly,double lz,double m1, double m2);
        void multiply_matrix(int INDEX,double alpha, double **matrixA, double **matrixB, double **matrixC);
        void convert_mtx_to_array(int M, int N, double **mtx, double *a);
        void convert_array_to_mtx(int M, int N, double *a, double **mt);
        void print_matrix( char* desc, int m, int n, double* a, int lda );
        double psirealspace(double y,double z,int ny,int nz,double ly,double lz,double **vec, int index); 
        void save_subband_energy(char *fnsubband);
	 void save_potential_used_for_2D_Schrodinger(char *fn_x, char *fn_yz);
              
    // Form factor
    void form_factor_calculation();
       void trapezoidal_weights();
       void printf_form_factor_calculation(int nx, int NSELECT);
	void save_form_factor_calculation(int s_th, char *fn);
        
    // Making scattering table and normalize
    void scattering_table();
         void normalize_table();
         void save_scattering_table(int s_th, char *fn);
         void save_normalized_table();//For checking
    // Save results
    void save_results();
    
    // initial free flight time
    void init_free_flight();
    
    // Thuc hien Ensemble Monte Carlo (EMC) cho 1 step dt
    void emcd();
         void drift(double tau);
              void check_boundary();
              
    // Thuc hien scattering() sau khi hat da drift
    void scattering() ;
         void Acoustic_mechanism(int s_th,int valley,int e_step,double ei,double eff_mass,double rr);
         void ZeroOrder_Optical_f_Absorption(int s_th,int val_befo,int val_aft,int e_step,double ei,double eff_mass,double rr,int index);
         void ZeroOrder_Optical_f_Emission(int s_th,int val_befo,int val_aft,int e_step,double ei,double eff_mass,double rr,int index);
         void ZeroOrder_Optical_g_Absorption(int s_th,int val,int e_step,double ei,double eff_mass,double rr,int index); 
         void ZeroOrder_Optical_g_Emission(int s_th,int val,int e_step,double ei,double eff_mass,double rr,int index);
         void isotropic(int valley_final,int subband_final,double energy_final,double effective_mass);
    
    // Kiem tra so hat o S va D theo charge neutrality de neu can thi eliminate hoac create
    void check_source_drain_contacts();
    
    // Hat co subband[] =9 thi can delete
    void delete_particles();
    
    // De tinh electron density cho dau vao cua Poisson
    void electron_density_caculation();
         void save_electron_population(double ***electron_population, int nx_max, int v_max, int NSELECT);
	  void save_electron_density(char *fn);
	  void save_electron_density_all_sections(char *fn);
         
    // Tinh va Save velocity, energy and current cumulative
    void velocity_energy_cumulative(FILE *f1, FILE *f2, double time,int iteration_reference);
         void save_VeloEnerCurrent_all_time_steps(FILE *f1, double *velo,double *ener,double *curr,int steps,int nx_max,double mesh_size_x);
         void save_VeloEnerCurrent_after_transient(FILE *f2, double *velo,double *ener,double *curr,int n_ti,int nx_max,double mesh_size_x);
    
    
    // Thuc hien MSMC
    void multi_subbands_MC(FILE *logfile, int *argc, char ***argv);
    
    // Tinh current
    void current_calculation(double *cur_av);
         void linreg(double x[],double y[],double sigma[],int N
              ,double a_fit[],double sig_a[],double yy[],double *chisqr );
    // Tinh eigen correction
    void eigen_energy_correction(int StepsUsedCorrection);
         void save_eig_correction(int StepsUsedCorrection); 

    // Save potential
    void save_potential(int flag);
    void save_potential_average(FILE *f1, FILE *f2, double time, int iter_reference);
    

// For Poisson 3D
    /*
    // Tai file init.c co 2 ham duoi day
    void Initialize(int *argc, char ***argv);
    void Finalize(); 
    // Ket thuc Tai file init.c
    
    // Tai init_local.c co cac ham sau
    void initialize_general_constants();
    void set_position_x_direction();
    void set_position_y_direction();
    void set_position_z_direction();
    void set_regions();
    void divide_region(double h,double xlimit,double *Pos,int *K,double *xv);
    // Ket thuc Tai init_local.c
    
    // Tai init_pot.c co cac ham
    void set_initial_potential();
    void assign_built_in_potential(); // Referenced by SetPhi_Contact().
    void guess_potential(); // Referenced by set_initial_potential().
    void make_flat_zero_potential(int zflag); // KHONG Thay dung ham nay
    void mpi_distribute_potential_from_node0(int node, int np); // Referenced by set_initial_potential().
    // Ket thuc tai init_pot.c  
    
    // Tai doping.c co cac ham sau
    void SetDopingDensityFor(int convflag, char region, char doping_type, double doping_density);
    double GetDopingDensityFor(char region); 
    double GetDoping(int p);
    int signfunc(double x);
    // Ket thuc tai doping.c
    
    // Tai nq.c co cac ham
    double Func_Nq(int i,int j,int k);
    double Func_DNq(int i,int j,int k);
    double nq(int i,int j,int k);
    double Dnq(int i,int j,int k);
    double nq_linearized(int p);
    double nq_linearized2(int p);
    void read_charge_density();
    void mpi_distribute_nq_from_node0(double ***nq, int node, int np); 
    // Ket thuc tai nq.c

    // Tai plain_outs.c co cac ham
    void print_potential_plain(char *fn);
    void print_midline_potential_plain(char *fn);
    void mpi_gather_potential_at_node0();
    // Ke thuc tai plain_outs.c
    
    // Tai outs.c co cac ham
    void SetPRINT_MidPotFlag();
    void SetPRINT_Pot3dFlag();
    int GetPRINT_MidPotFlag();
    int GetPRINT_Pot3dFlag();
    void print_output();
    // Ket thuc tai outs.c
    
    // Tai define_poi.c xem poisson3d.h  
       
    void solve_3d_poisson() ;
    void solve_3d_poisson_linearized() ;
    void print_output() ;
    void read_charge_density() ;
    void (*SolvePoisson)() ;
    
    */
    // De chay duoc ham mai cua GS
    void Initialize() ;
    void Finalize() ;
    void solve_3d_poisson() ;
    void solve_3d_poisson_linearized() ;
    void print_output() ;
    void read_charge_density() ;
    void (*SolvePoisson)() ;
    
    // Solve Poisson at equilibrium
    void solve_Poisson_at_equilibrium(int *argc, char ***argv);
   
    void logfile_PoissonInput();
         void logging(FILE *log,const char *fmt,...);// Ham nay de tao ta logfile

    
#endif // _FUNCTIONS_H_
