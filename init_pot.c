
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"

static PoiNum N;
static int Ntot;
static int entry;
/******************************************************************************
 Thay doi 1 vai cai de mapping code cua ta voi 3D Poisson
 INPUT:  ***Phi = GetPot() cua Poisson
 
 Starting date: May 04, 2010
 Latest update: May 04, 2010 
******************************************************************************/


void set_initial_potential()
{
  if ( ++entry > 1 ) return ;

  void guess_potential() ;
  void make_flat_zero_potential() ;
  void mpi_distribute_potential_from_node0() ;

  // MPI

  int node, np ;

  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  // begin

  FILE *fphi = fopen("pot.r","r") ;
  double ***Phi = GetPot() ;
  N = GetPoiNum() ;
  Ntot = N.x*N.y*N.z ;

  int i,j,k,nx,ny,nz ;
  double v ;

  if ( fphi != (FILE *)NULL ) {
    if ( node==0 ) { 
      PetscPrintf(PETSC_COMM_WORLD,"\n*******\nREADING PHI FROM EXISTING FILE\n*******\n") ;
      fscanf(fphi,"%*s %d %d %d\n",&nx,&ny,&nz) ;
      if ( nx!=N.x || ny!=N.y ||  nz!=N.z ) report_error("Problems reading Phi's") ;
      
      for ( i=1 ; i<=nx ; i++ ) fscanf(fphi,"%*s") ;
      for ( j=1 ; j<=ny ; j++ ) fscanf(fphi,"%*s") ;
      for ( k=1 ; k<=nz ; k++ ) fscanf(fphi,"%*s") ;

      for ( i=0 ; i<N.x ; i++ )
        for ( j=0 ; j<N.y ; j++ )
          for ( k=0 ; k<N.z ; k++ ) {
            fscanf(fphi,"%le",&v) ;
	    Phi[i][j][k] = -v ;
	  }
    } 
    mpi_distribute_potential_from_node0(node,np) ;
  }
  else 
    guess_potential() ;

  // Print Midline Potential

  if ( GetPRINT_MidPotFlag() ) {
    void print_midline_potential_plain() ;
    void reset_print_midline_potential() ;
    
    print_midline_potential_plain("mpot.r0") ;
  }

}

void guess_potential() // Referenced by set_initial_potential().
{ 
  //****************Ta them vao*******************
  void mpi_distribute_potential_from_node0(int node, int np);


  // MPI
  int node, np ;
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;
  // begin
  N = GetPoiNum() ; // Neu khong du so du lieu thi lam sao ma ham guess_potential() chay dung duoc
  Ntot = N.x*N.y*N.z ;
  int nx,ny,nz ;
  double v ;
  //**********Het phan them vao*********************** 

  int p,i,j,k,region ;
  Region R ;

  double Vd = GetDrainVoltage() ;
  double phi_source = GetPhi_Source() ;
  double phi_drain = GetPhi_Drain() ;
  double phi_gate = GetPhi_Gate() ;
  double ***Phi = GetPot() ;

  for ( p=1 ; p<=Ntot ; p++ ) {
    get_ijk(p,&i,&j,&k) ;
    region = GetRegion(&R,p) ;
    if ( i<N.xa ) {
        Phi[i][j][k] = phi_source;
        //printf("\n phi_source tai guess_potential = %f", phi_source);
        }
    else if ( i>N.xb ){
         Phi[i][j][k] = phi_drain;
         //printf("\n phi_drain tai guess_potential = %f", phi_drain);
         }
    else {
         Phi[i][j][k] = (phi_drain-phi_source)/(N.xb-N.xa)*(i-N.xa)+phi_source;
         }

    if ( GATE_CONTACT ){
         Phi[i][j][k] = phi_gate;
         //printf("\n phi_gate tai guess_potential = %f", phi_gate);
         }
  }// End of for ( p=1 ; p<=Ntot ; p++ ) {

  //************Ta them vao*********************
   // Print Midline Potential
  if ( GetPRINT_MidPotFlag() ) {
    void print_midline_potential_plain() ;
    void reset_print_midline_potential() ;
    print_midline_potential_plain("mpot.r0") ;
  }

}// End of void guess_potential() 

/* Ham nay duoc dinh o definefn.c nhung de o day cho tien theo doi
void SetPhi_Contact() {

  void assign_built_in_potential() ;
  assign_built_in_potential() ;

  printf("\n ********** From SetPhi_Contact() *********");

  double Vd = GetDrainVoltage() ;
  printf("\n Drain Voltage Vd = %f", Vd);
  
  double Vg = GetGateVoltage() ;
  printf("\n Gate Voltage Vg  = %f", Vg);

  Phi_Source = E.biS ;
  printf("\n Phi_Source       = %f",Phi_Source);
  
  Phi_Drain  = E.biD + Vd ;
  printf("\n Phi_Drain        = %f",Phi_Drain);
  
  Phi_Gate   = Vg - P.offset_gate ; // for n-mos;  //Phi_Gate   = -Vg + P.offset_gate + Vd ; // for p-mos with CB_PEMT mode
  printf("\n Phi_Gate         = %f",Phi_Gate);
  
  //printf("\n ********** Ket thuc SetPhi_Contact() *********");//  getchar();
}

double GetPhi_Source() { 
  return(Phi_Source) ; 
}
double GetPhi_Drain() { 
  return(Phi_Drain) ; 
}
double GetPhi_Gate() {
  return(Phi_Gate) ; 
}

***************************************************************************************************************************/
void assign_built_in_potential() // Referenced by SetPhi_Contact().Ta da sua lai 1 it roi ( June 10, 2010)

{ 

  Energy *E = PGetEnergy() ;
  Param  *P = PGetParam() ;
  double kT = GetKT() ;
  
  //printf("\n ******From assign_built_in_potential() ******");

  double EF       = GetEF() ;   
  //printf("\n EF   = %f",EF);
  
  double EFnS     = GetEFnS() ; 
  //printf("\n EFnS = %f",EFnS);
  
  double EFnD     = GetEFnD() ; 
  //printf("\n EFnD = %f",EFnD);
 
  // Built-in Potentials
  /* Sua lai xem duoi. Cai nay la cua GS
  E->biS = EFnS - EF ;
  printf("\n Built-in potential E->biS = %f",E->biS);
  
  E->biD = EFnD - EF ;
  printf("\n Built-in potential E->biD = %f",E->biD);
  */
//*******************************************************************
  // Minh them vao de tinh built-in potential
  double *Get_doping_density();
  double *doping_density = Get_doping_density();//1,2 or 3 ->S,D or Channel

  double Get_Vt();
  double Vt = Get_Vt();

  double Get_intrinsic_carrier_density();
  double intrinsic_carrier_density = Get_intrinsic_carrier_density();
  
  double factor = 0.0, dop_temp =0.0, built_in_S =0.0, built_in_D=0.0 ; 

  // Tinh built-in potential tai Source
  dop_temp = doping_density[1];// Doping tai Source
  dop_temp = 0.5*dop_temp/intrinsic_carrier_density; // (May 26, 2010) Lam viec thieu can than la the day !
  if(dop_temp > 0.0){ // Hien nhien la the do ta dope la n-type
     factor = dop_temp + sqrt(dop_temp*dop_temp + 1.0);
     built_in_S =  Vt*log(factor);
  }
  else { // dop_term <= 0.0. Trong truong hop nay khong xay ra nhung van lam
        dop_temp = - dop_temp;
        factor = dop_temp + sqrt(dop_temp*dop_temp + 1.0);
        built_in_S = - Vt*log(factor);
  }

  // Tinh built-in potential tai Drain
  dop_temp = doping_density[2]; // Doping tai Drain
  dop_temp = 0.5*dop_temp/intrinsic_carrier_density; //(May 26, 2010)
  if(dop_temp > 0.0){ // Hien nhien la the do ta dope la n-type
     factor = dop_temp + sqrt(dop_temp*dop_temp + 1.0);
     built_in_D = Vt*log(factor);
  }
  else { // dop_temp <= 0.0. Trong truong hop nay khong xay ra nhung van lam
       dop_temp = - dop_temp;
       factor = dop_temp + sqrt(dop_temp*dop_temp + 1.0);
       built_in_D = - Vt*log(factor);
  }
  //printf("\n built_in_S =%f, built_in_D =%f",  built_in_S, built_in_D );  getchar();
  // Ket thuc phan minh them vao
//**************************************************************************************
  
 // Built-in Potentials. Can sua lai the nay

  E->biS = built_in_S ;//printf("\n Built-in potential E->biS = %f",E->biS);
  
  E->biD = built_in_D ;//  printf("\n Built-in potential E->biD = %f",E->biD);  getchar();
  //*******************************************************************
  
  // Gate Offset
  
  //printf("\n Before P->offset_gate = %f, E->c_si = %f ",P->offset_gate, E->c_si);
  //P->offset_gate -= (E->c_si - EF) ;
  // Do xay ra loi khi chay Poisson lan thu 2, 3, 4 vi gia tri P->offset_gate cu bi tru di. Nen ta set luon no la 0 cho tien

  P->offset_gate = 0.0;

  // printf("\n After P->offset_gate = %f",P->offset_gate);
  
  //printf("\n ******Ket thuc assign_built_in_potential() ******");
}

//********************************************************************************************

void make_flat_zero_potential(int zflag) // KHONG Thay dung ham nay
{
  int p,i,j,k,region ;
  Region R ;
  double ***Phi = GetPot() ;

  if ( zflag == 0 )
    PetscPrintf(PETSC_COMM_WORLD,"\n******\n FLAT ZERO POTENTIAL \n******\n") ;

  int Nyz = N.y*N.z ;

  for ( p=1 ; p<=Nyz ; p++ ) {
    get_ijk(p,&i,&j,&k) ;
    region = GetRegion(&R,p) ;

    Phi[i][j][k] = 0.0 ;
  }

}

void mpi_distribute_potential_from_node0(int node, int np) // Referenced by set_initial_potential().
{
  int x0,y0,z0,x1,y1,z1 ;
  int n,dest,datasize,ntrans,p0,sum_transmitted,tag ;
  int trans_size = GetTransSize(); //printf("\n Trans size = %d",trans_size);
  MPI_Status status ;
  double ***Pot ;

  Pot = GetPot() ;

  x0 = N.x0-1 ; x1 = N.x1+1 ;
  y0 = N.y0-1 ; y1 = N.y1+1 ;
  z0 = N.z0-1 ; z1 = N.z1+1 ;

  datasize = (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ;
  ntrans = datasize/trans_size ;

  if ( node==0 ) {
    for ( dest=1 ; dest<np ; dest++ ) {
      p0 = 0 ;
      sum_transmitted = 0 ;
      tag = 0 ;
      for ( n=1 ; n<=ntrans ; n++ ) {
        MPI_Send(&Pot[x0][y0][z0]+p0,trans_size,MPI_DOUBLE,dest,++tag,PETSC_COMM_WORLD) ;
        p0 += trans_size ;
        sum_transmitted += trans_size ;
      }
      if ( sum_transmitted < datasize )
        MPI_Send(&Pot[x0][y0][z0]+p0,datasize-sum_transmitted,MPI_DOUBLE,dest,++tag,PETSC_COMM_WORLD) ;
    }
  }
  else {
    p0 = 0 ;
    sum_transmitted = 0 ;
    tag = 0 ;
    
    for ( n=1 ; n<=ntrans ; n++ ) {
      MPI_Recv(&Pot[x0][y0][z0]+p0,trans_size,MPI_DOUBLE,0,++tag,PETSC_COMM_WORLD,&status) ;
      p0 += trans_size ;
      sum_transmitted += trans_size ;
    }
    if ( sum_transmitted < datasize )
      MPI_Recv(&Pot[x0][y0][z0]+p0,datasize-sum_transmitted,MPI_DOUBLE,0,++tag,PETSC_COMM_WORLD,&status) ;
  }
}
// End of void mpi_distribute_potential_from_node0(int node, int np) 
