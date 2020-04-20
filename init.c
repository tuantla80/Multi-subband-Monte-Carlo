#include <stdio.h>
#include "petscksp.h"
#include "nanowire.h"

static char help[] = "Solves Poisson equation in parallel with KSP" ;

void Initialize(int *argc, char ***argv)
{
  // initialize MPI
  void read_poisson_input();
  void initialize_general_constants();
  void set_position_z_direction();
  void set_position_y_direction();
  void set_position_x_direction();
  void set_regions() ;	
  void initialize_poisson() ; 	

  PetscInitialize(argc,argv,(char *)0,help);// KHONG Lay ra o ham main luon
  
  read_poisson_input(); // Doc tham so #POISSON_INPUT tu Input file
  initialize_general_constants(); // Chuyen 1 so constant ve dang DON VI VAT LY day du
  //printf("\n Vuot qua ham  initialize_general_constants()");

  set_position_z_direction();
  //printf("\n Vuot qua ham  set_position_z_direction()");

  set_position_y_direction() ;
  //printf("\n Vuot qua ham  set_position_y_direction()");

  set_position_x_direction() ;
  //printf("\n Vuot qua ham  set_position_x_direction()");

  set_regions() ;
  //printf("\n Vuot qua ham  set_regions()");

  initialize_poisson() ;
  //printf("\n Vuot qua ham initialize_poisson()");

  //printf("\n Ket thuc Initialize() o file init.c");
}

void Finalize() 
{
  PetscFinalize() ;
}

