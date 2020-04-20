#include "petscksp.h"
#include "nanowire.h"

void print_potential_plain(char *fn)
{
  FILE *fout ;
  int i,j,k ;

  PoiNum N = GetPoiNum() ;
  double *X = GetPX() ;
  double *Y = GetPY() ;
  double *Z = GetPZ() ;
  double ***Phi = GetPot() ;

  void mpi_gather_potential_at_node0() ;
  mpi_gather_potential_at_node0() ;

  int node ;
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;

  if ( node==0 ) {

    fout = fopen(fn,"w") ;
    fprintf(fout,"# %d %d %d\n",N.x,N.y,N.z) ;

    for ( i=N.x0 ; i<=N.x1 ; i++ ) fprintf(fout,"%le\n",X[i]) ;
    for ( j=N.y0 ; j<=N.y1 ; j++ ) fprintf(fout,"%le\n",Y[j]) ;
    for ( k=N.z0 ; k<=N.z1 ; k++ ) fprintf(fout,"%le\n",Z[k]) ;

    for ( i=0 ; i<N.x ; i++ ) 
      for ( j=0 ; j<N.y ; j++ )
	for ( k=0 ; k<N.z ; k++ )
	  fprintf(fout,"%le\n",-Phi[i][j][k]) ;
    
    fclose(fout) ;
  }
}

// Midline Potential

void print_midline_potential_plain(char *fn)
{
  int i ;
  FILE *fout ;

  PoiNum N = GetPoiNum() ;
  double *X = GetPX() ;
  double ***Phi = GetPot() ;

  int node ;
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;

  void mpi_gather_potential_at_node0() ;
  mpi_gather_potential_at_node0() ;

  int yc = (N.y0+N.y1)/2 ;
  int zc = (N.z0+N.z1)/2 ;

  if ( node==0 ) {

    fout = fopen(fn,"w") ;
    fprintf(fout,"# %d %d %d\n",N.x,N.y,N.z) ;
    fprintf(fout,"# Vd = %lf Vg = %lf\n",GetDrainVoltage(),GetGateVoltage()) ;

    for ( i=0 ; i<N.x ; i++ ) 
      fprintf(fout,"%le %le\n",X[i],-Phi[i][yc][zc]) ;

    fprintf(fout,"\n") ;
    fflush(fout) ;
  }
}

// MPI

void mpi_gather_potential_at_node0()
{
  int node,np ;
  int trans_size = GetTransSize() ;
  MPI_Status status ;

  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  double ***Pot = GetPot() ;
  PoiNum N = GetPoiNum() ;

  int x0 = GetPx0() ;
  int y0 = N.y0-1 ;
  int y1 = N.y1+1 ;
  int z0 = N.z0-1 ;
  int z1 = N.z1+1 ;

  int n, n0, n1, data_size, ntrans, sum_transmitted, source, p0, tag ;
  void GetLocalN() ;

  if ( node>0 ) {
    GetLocalN(node,&n0,&n1) ;
    n0 += x0 - 1 ;
    n1 += x0 - 1 ;

    data_size = (n1-n0+1)*(y1-y0+1)*(z1-z0+1) ;
    ntrans = data_size/trans_size ;

    sum_transmitted = 0 ;
    p0 = 0 ;
    tag = 0 ;

    for ( n=1 ; n<=ntrans ; n++ ) {
      MPI_Send(&Pot[n0][y0][z0]+p0, trans_size, MPI_DOUBLE, 0, ++tag, PETSC_COMM_WORLD) ;
      p0 += trans_size ;
      sum_transmitted += trans_size ;
    }
    if ( sum_transmitted < data_size ) 
      MPI_Send(&Pot[n0][y0][z0]+p0, data_size-sum_transmitted, MPI_DOUBLE, 0, ++tag, PETSC_COMM_WORLD) ;

  }
  else {
    for ( source=1 ; source<np ; source++ ) {
      GetLocalN(source,&n0,&n1) ;
      n0 += x0 - 1 ;
      n1 += x0 - 1 ;

      data_size = (n1-n0+1)*(y1-y0+1)*(z1-z0+1) ;
      ntrans = data_size/trans_size ;
      sum_transmitted = 0 ;
      p0 = 0 ;
      tag = 0 ;

      for ( n=1 ; n<=ntrans ; n++ ) {
	MPI_Recv(&Pot[n0][y0][z0]+p0, trans_size, MPI_DOUBLE, source, ++tag, PETSC_COMM_WORLD, &status) ;
	p0 += trans_size ;
	sum_transmitted += trans_size ;
      }
      if ( sum_transmitted < data_size )
	MPI_Recv(&Pot[n0][y0][z0]+p0, data_size-sum_transmitted, MPI_DOUBLE, source, ++tag, PETSC_COMM_WORLD, &status) ;
    }
    
  } 
}
