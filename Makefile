LOCAL_MKL_LIB = -L/opt/intel/mkl/10.1.1.019/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread
CC = /opt/intel/Compiler/11.0/081/bin/intel64/icc
CFLAGS           = -g 
MANSEC           = KSP

MAIN = main.o init.o init_local.o definefn.o material_param.o device_structure.o read_poisson_input.o read_voltages_input.o read_simulation_list.o read_scattering_save_list.o\
	initial_variable.o find_region.o source_drain_carrier_number.o  \
	init_kspace.o init_realspace.o random2.o electrons_initialization.o count_used_particles.o trapezoidal_weights.o Solved_2D_Schro_for_MSMC.o form_factor_calculation.o \
	scattering_table.o init_free_flight.o multi_subbands_MC.o  current_calculation.o emcd.o check_source_drain_contacts.o delete_particles.o electron_density_caculation.o \
	velocity_energy_cumulative.o eigen_energy_correction.o save_results.o drift.o scattering.o save_potential_parameters_from_input_file.o save_potential_average.o \
       initial_potential_doping.o logfile_PoissonInput.o

POISSON = init_pot.o poisson3d.o define_poi.o doping.o
SELFCONS = nq.o
AUX = outs.o plain_outs.o
MATH = my_util.o nrutil.o
LINKLIB = $(LOCAL_MKL_LIB) -lm
LIB = -lf2c -lm

include ${PETSC_DIR}/bmake/common/base

msmc: $(MAIN) $(POISSON) $(MTXH) $(NEGF) $(SELFCONS) $(AUX) $(MATH) $(DEBUG) chkopts
	-${CLINKER} -o  $@ $(MAIN) $(POISSON) $(MTXH) $(NEGF) $(SELFCONS) $(AUX) $(MATH) $(DEBUG) $(LIB) ${PETSC_KSP_LIB} ${RP_LIB}  $(LINKLIB)

remove:
	rm -f *.o
	rm -f *.o*
	rm -f *.e*
	rm -f *.*~
	rm -f msmc

clean:
	rm -f *.r*
	rm -f *.dat
	rm -f *.txt
	rm -f *.log
	rm -f *.o
	rm -f *.o*
	rm -f *.e*
	rm -f *.*~
	rm -f msmc

   	

