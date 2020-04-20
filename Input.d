#DEVICE_PARAMETER
Device_type    = NanowireFET
L_Source/D[nm] = 4.0
L_Channel[nm]  = 4.0
L_gate[nm]     = 4.0
T_ox[nm]       = 1.0
T_box[nm]      = 1.0
W_ox[nm]       = 1.0
T_Si[nm]       = 2.0
W_Si[nm]       = 2.0
T_gate[T_ox+T_Si+T_box]    = 4.000000
W_gate[(W_ox+W_Si+W_ox)/2] = 2.00000
      
mesh_size_x[nm] = 0.2
mesh_size_y[nm] = 0.1
mesh_size_z[nm] = 0.1

Tox_mesh(discrete_nodes_for_Tox)   = 10
Tbox_mesh(discrete_nodes_for_Tbox) = 10
Wox_mesh(discrete_nodes_for_Wox)   = 10
	
n_region(TotalNumberActiveRegions) = 3
Source_Junction_Doping_type	       = n
doping_density(1)_Source[/cm3] 	= 1.0e+20
Drain_Junction_Doping_Type	       = n
doping_density(2)_Drain[/cm3]  	= 1.0e+20
Channel_Doping_Type		       = i
doping_density(3)_Channel[/cm3]	= -1.45e+10
	
#VOLTAGES_INPUT
V_Source_input[V]		      	   = 0.0	
V_Drain_input[V]...initial_final_step = 0.5,  0.5,  0.1 
V_Gate_input[V]...initial_final_step  = 0.5,  0.5,  0.1
        
#MATERIAL_PARAMETER
Temperature(K)			  = 300.0
Electron_energy_max[eV]_2x2_2.5eV	  = 4.0
Band_Gap_Si(eV)			  = 1.12
Band_Gap_of_Ox(eV)		     	  = 9.0
Longitudinal_effective_mass_(ml)     = 0.98
Tranverse_effective_mass_(mt)	  = 0.196
nonparabolicity_factor_[1/eV]	  = 0.5
Constant_dielectric_Silicon(eps_sc)  = 11.9
Constant_dielectric_Oxide(eps_oxide) = 3.9
intrinsic_density_Silicon(/cm3)	  = 1.45e10
gateTYPE_(metal_or_nPolysilicon)     = metal
Affinity_of_Silicon_(eV)_(*)	  = 4.00
Metal_gate_workfunction_(eV)_(**)    = 4.56
Gate_Workfunction_Offset(eV)(**-*)   = 0.56
Crystal_Density[kg/m3]		  = 2329.0
sound_velocity[m/s]		         = 9037.0	
Acoustic_Deformation_Pot[eV]	  = 12.0
f_phonon_deformation_pot[ev/m]	  = 1.1e11
f_phonon_energy[eV]		         = 0.059
g_phonon_deformation_pot[ev/m]	  = 8.0e10
g_phonon_energy[eV]		         = 0.063

#SCATTERING_LIST
	
(=1_Ballistic;_=0_Diffusive)           = 0
Acoustic(=1_choosing;_=0_not_choosing) = 1	
Zero_Order_Non_Polar_Optical_Phonon    = 1
First_Order_Non_Polar_Optical_Phonon   = 0				
Coulomb_scattering	                  = 0										
Surface_roughness	                  = 0

#SAVE_LIST
Save_eig_wavefunct(_1_Save_0_NotSave)   = 0
Save_doping_potential_initialization    = 0
Save_electron_initlization_distribution = 1
Save_form_factor_calculation	     = 1
Save_scattering_rate			     = 0
Save_electron_population_THUC_SU_chon_0 = 0
Save_electron_Density		     = 1
Save_velo_ener_curr_cumu_all_time_steps = 1
Save_velo_ener_curr_cumu_after_transie  = 1
Save_init_free_flight		     = 1

#SIMULATION_LIST 
dt...MC_time_step[fs]_(1fs=1.0e-15s)			= 0.1
Total_time[ps]_(1ps=1.0e-12s)				= 0.0004
Transient_time_after_which_results_are_calculated[ps]   = 0.0000
Length_where_3D_plots_in_cross_section_are_plotted[nm]  = 6.0
Depth_where_the_3D_plots_are_plotted[nm]			= 2.0
Width_where_the_3D_plots_are_plotted[nm]			= 2.0
Number_Subbands_for_2D_Schrodinger_(NSELECT)	       = 4
N_time_steps_to_Solve_2D_Schrodinger_Equation		= 1
artificial_factor_of_cell_volume_(default_1.0e-6)	= 3.0e-6	

#POISSON_INPUT
Max_Iteration_of_Newton_Method_GiaoSu_200 = 200
Convergence_Eps = 1.000000e-06
KSPType(Richardson,Chebychev,CG,BiCG,GMRES,BiCGSTAB,CGS,QMR1,QMR2,CR,LSQR,PreOnly) = BiCGSTAB
PCType(Jacobi,Block_Jacobi,SOR,SOR_Eis,ICC,ILU,ASM,LinearSolv,Combinations) = Block_Jacobi
KSP_Rtol = 1.000000e-09
GMRES_Restart = 50
Poisson_Linearized(n/y) = n
Transmission_Data_Size = 1000
Print_Midline_Potential(y/n) = n
Print_3D_Potential_at_Each_Bias_Point(y/n) = n

The end


