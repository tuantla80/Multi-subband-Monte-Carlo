#include "nrutil.h"
#include "petscksp.h"
#include "nanowire.h"
#include "region.h"


static Mat MatA ;
static Vec Vecb,Vecx ;
static KSPType ksptype ;
static KSP ksp ;
static PC pc ;
static PCType pctype ;
static double ksprtol ;
static int gmres_restart ;

static PoiNum N ;
static int Ntot ;
static int *diag_entry ;
static double *alpha_xplus,*alpha_yplus,*alpha_zplus,*alpha_xminus,*alpha_yminus,*alpha_zminus ;
static double *alpha_center,*alpha_const,*alpha_q0e ;

double falpha_xplus(),falpha_yplus(),falpha_zplus() ;
double falpha_xminus(),falpha_yminus(),falpha_zminus() ;

#undef __FUNCT__
#define __FUNCT__ "initialize_poisson"

int initialize_poisson()
{
  N = GetPoiNum() ;
  Ntot = N.x*N.y*N.z ;

  PoiMem *M = PGetPoiMem() ;

  M->Phi = d3matrix_mpi(N.x0-1,N.x1+1,N.y0-1,N.y1+1,N.z0-1,N.z1+1) ;
  M->phi = d3matrix_mpi(N.x0-1,N.x1+1,N.y0-1,N.y1+1,N.z0-1,N.z1+1) ;
  M->n3d = d3matrix_mpi(N.x0,N.x1,N.y0,N.y1,N.z0,N.z1) ;
  M->p3d = d3matrix_mpi(N.x0,N.x1,N.y0,N.y1,N.z0,N.z1) ;
  
  void set_null_out_of_bound() ;
  set_null_out_of_bound() ;

  // Matrix Creation

  int np ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  int Nlocal = Ntot/np ;
  int d_nz,o_nz ;
  int *d_nnz = PETSC_NULL ;
  int *o_nnz = PETSC_NULL ;

  if ( Nlocal < N.z ) d_nz = 3  ;
  else if ( Nlocal < N.y*N.z ) d_nz = 5 ;
  else d_nz = 7 ;

  o_nz = 7-d_nz ; 

  PetscErrorCode ierr ;

  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,d_nz,d_nnz,o_nz,o_nnz,&MatA) ; CHKERRQ(ierr);
  ierr = MatSetFromOptions(MatA) ; CHKERRQ(ierr) ;

  // Assign Alphas

  void poisson_memory_allocation() ;
  poisson_memory_allocation() ;

  PoiParam CP = GetPoiParam() ;
  void set_alphas_1() ;

  set_alphas_1() ;

  // Assign Off-diagonal Elements

  int I,J,Istart,Iend,p,i,j,k ;
  int Ns = N.y*N.z ;
  double v ;

  ierr = MatGetOwnershipRange(MatA,&Istart,&Iend) ; CHKERRQ(ierr) ;

  for ( I=Istart ; I<Iend ; I++ ) {
    p = I+1 ;
    get_ijk(p,&i,&j,&k) ;
    if ( i-1>=N.x0 ) {
      J = I-Ns ;
      v = alpha_xminus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
    if ( j-1>=N.y0 ) {
      J = I-N.z ;
      v = alpha_yminus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
    if ( k-1>=N.z0 ) {
      J = I-1 ;
      v = alpha_zminus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
    if ( k+1<=N.z1 ) {
      J = I+1 ;
      v = alpha_zplus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
    if ( j+1<=N.y1 ) {
      J = I+N.z ;
      v = alpha_yplus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
    if ( i+1<=N.x1 ) {
      J = I+Ns ;
      v = alpha_xplus[p] ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }
  }

  // Create Vectors

  ierr = VecCreate(PETSC_COMM_WORLD,&Vecb) ; CHKERRQ(ierr) ;
  ierr = VecSetSizes(Vecb,PETSC_DECIDE,Ntot) ; CHKERRQ(ierr) ;
  ierr = VecSetFromOptions(Vecb) ; CHKERRQ(ierr) ;
  ierr = VecDuplicate(Vecb,&Vecx) ; CHKERRQ(ierr) ;

  // KSP Setup

  PCType GetPCType() ;
  KSPType GetKSPType() ;

  pctype = GetPCType() ;
  ksptype = GetKSPType() ;
  ksprtol = CP.ksprtol ;
  gmres_restart = CP.gmres_restart ;

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp) ; CHKERRQ(ierr) ;
  ierr = KSPGetPC(ksp,&pc) ; CHKERRQ(ierr) ;
  ierr = PCSetType(pc,pctype) ; CHKERRQ(ierr) ;
  ierr = KSPSetType(ksp,ksptype) ; CHKERRQ(ierr) ;
  ierr = KSPGMRESSetRestart(ksp,gmres_restart) ; CHKERRQ(ierr) ;
  ierr = KSPSetTolerances(ksp,ksprtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT) ; CHKERRQ(ierr) ;
  ierr = KSPSetFromOptions(ksp) ; CHKERRQ(ierr) ;

  // Set n0 and n1 for schrodinger and negf
  
  int I0,I1 ;
  int *I0s = ivector(0,np-1) ;
  int *I1s = ivector(0,np-1) ;

  get_ijk(Istart+1,&I0,&j,&k) ;
  if ( !(j==N.y0 && k==N.z0) ) I0++ ;

  get_ijk(Iend,&I1,&j,&k) ;

  MPI_Allgather(&I0,1,MPI_INT,&I0s[0],1,MPI_INT,PETSC_COMM_WORLD) ;
  MPI_Allgather(&I1,1,MPI_INT,&I1s[0],1,MPI_INT,PETSC_COMM_WORLD) ;

  void SetLocalN() ;
  SetLocalN(np,I0s,I1s) ;

  free_ivector(I0s,0,np-1) ;
  free_ivector(I1s,0,np-1) ;

  return 0 ;
}

#undef __FUNCT__
#define __FUNCT__ "solve_3d_poisson"

int solve_3d_poisson_linearized()
{
  // Skip If Read from Existing File

  int node,np ;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  // Set-Up

  int Istart,Iend,low,high,local_x0,local_x1,j,k ;
  double ***Phi = GetPot() ;
  void set_alphas_2(),fix_corners() ;

  PoiParam CP = GetPoiParam() ;
  int MaxIter = CP.MaxIter ;
  
  PetscErrorCode ierr ;
  PetscScalar *Vecx_local ;

  ierr = MatGetOwnershipRange(MatA,&Istart,&Iend) ; CHKERRQ(ierr) ;
  ierr = VecGetOwnershipRange(Vecb,&low,&high) ; CHKERRQ(ierr) ;

  get_ijk(Istart+1,&local_x0,&j,&k) ;
  get_ijk(Iend,&local_x1,&j,&k) ;

  fix_corners(Phi,local_x0,local_x1) ;
  
  set_alphas_2() ;
  
  // Iteration

  int iter,I,J,p,i,no_ksp_its ;
  double v,norm_local,norm,norm_old=-1.0e30,nq_linearized(),nq_linearized2(),fabs(),sqrt() ;
  void mpi_transfer_phi() ;

  // Assigning Diagonal Elements

  for ( I=Istart ; I<Iend ; I++ ) {
    p = I+1 ;
    J = I ;
    v = alpha_center[p]-alpha_q0e[p]*nq_linearized(p) ;
    ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
  }

  // Matrix Assemble
  
  ierr = MatAssemblyBegin(MatA,MAT_FINAL_ASSEMBLY) ; CHKERRQ(ierr) ;
  ierr = MatAssemblyEnd(MatA,MAT_FINAL_ASSEMBLY) ; CHKERRQ(ierr) ;
  ierr = MatSetOption(MatA,MAT_NO_NEW_NONZERO_LOCATIONS);CHKERRQ(ierr);
  
  // Assinging Right Hand Side Vector
  
  for ( I=low ; I<high ;  I++ ) {
    p = I+1 ;
    v = -alpha_const[p] + alpha_q0e[p]*nq_linearized2(p) ;
    ierr = VecSetValues(Vecb,1,&I,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
  }

  // KSP Procedure
  
  ierr = KSPSetOperators(ksp,MatA,MatA,SAME_NONZERO_PATTERN) ; CHKERRQ(ierr) ;
  ierr = KSPSolve(ksp,Vecb,Vecx) ; CHKERRQ(ierr) ;
  ierr = KSPGetIterationNumber(ksp,&no_ksp_its) ;
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP_Iter_No = %d\n",no_ksp_its) ; CHKERRQ(ierr) ;

  if ( no_ksp_its==0 ) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP_Iter_No = 0 Ocurred") ; 
    Finalize() ;
    exit(0) ;
  }

  // Update phi locallly

  ierr = VecGetArray(Vecx,&Vecx_local) ; CHKERRQ(ierr) ;
    
  for ( I=Istart ; I<Iend ; I++ ) {
    p = I+1 ;
    get_ijk(p,&i,&j,&k) ;
    Phi[i][j][k] = Vecx_local[I-Istart] ;
  }
  ierr = VecRestoreArray(Vecx,&Vecx_local) ; CHKERRQ(ierr) ;

  // Transfer Phi

  mpi_transfer_phi(Phi,node,np,Istart,Iend) ;
    
  // Fix Corners

  fix_corners(Phi,local_x0,local_x1) ;

  // Print Midline Potential

  if ( GetPRINT_MidPotFlag() ) {
    char fn[100] ;
    void print_midline_potential_plain() ;

    sprintf(fn,"mpot.r-Vd%.2lf-Vg%.2lf",GetDrainVoltage(),GetGateVoltage()) ;
    print_midline_potential_plain(fn,0) ;
  }

  return 0 ;
}

int solve_3d_poisson()
{
  // Skip If Read from Existing File

  int node,np ;
  double Func_Nq(), Func_DNq() ;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  // Set-Up

  int Istart,Iend,low,high,local_x0,local_x1,j,k ;
  double ***Phi = GetPot() ;
  double ***phi = GetPhiKth() ;
  void set_alphas_2(),fix_corners() ;

  PoiParam CP = GetPoiParam() ;
  int MaxIter = CP.MaxIter ;
  
  PetscErrorCode ierr ;
  PetscScalar *Vecx_local ;

  ierr = MatGetOwnershipRange(MatA,&Istart,&Iend) ; CHKERRQ(ierr) ;
  ierr = VecGetOwnershipRange(Vecb,&low,&high) ; CHKERRQ(ierr) ;

  get_ijk(Istart+1,&local_x0,&j,&k) ;
  get_ijk(Iend,&local_x1,&j,&k) ;

  fix_corners(Phi,local_x0,local_x1) ; 
  copy_R3mtx(phi,Phi,local_x0-1,local_x1+1,N.y0-1,N.y1+1,N.z0-1,N.z1+1) ;
  
  set_alphas_2() ;
  
  // Iteration

  int iter,I,J,p,i,no_ksp_its ;
  double v,norm_local,norm,norm_old=-1.0e30,Fnn(),DFnn(),fabs(),sqrt() ;
  void mpi_transfer_phi() ;

  PetscPrintf(PETSC_COMM_WORLD,"\n - Poisson Solver -\n") ; 

  for ( iter=1 ; iter<=MaxIter ; iter++ ) {

    // Assigning Diagonal Elements

    for ( I=Istart ; I<Iend ; I++ ) {
      p = I+1 ;
      J = I ;
      v = DFnn(p,Func_DNq) ;
      ierr = MatSetValues(MatA,1,&I,1,&J,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }

    // Matrix Assemble

    ierr = MatAssemblyBegin(MatA,MAT_FINAL_ASSEMBLY) ; CHKERRQ(ierr) ;
    ierr = MatAssemblyEnd(MatA,MAT_FINAL_ASSEMBLY) ; CHKERRQ(ierr) ;
    ierr = MatSetOption(MatA,MAT_NO_NEW_NONZERO_LOCATIONS);CHKERRQ(ierr);

    // Assinging Right Hand Side Vector

    for ( I=low ; I<high ;  I++ ) {
      p = I+1 ;
      v = Fnn(p,Func_Nq) ;
      ierr = VecSetValues(Vecb,1,&I,&v,INSERT_VALUES) ; CHKERRQ(ierr) ;
    }

    // KSP Procedure

    ierr = KSPSetOperators(ksp,MatA,MatA,SAME_NONZERO_PATTERN) ; CHKERRQ(ierr) ;
    ierr = KSPSolve(ksp,Vecb,Vecx) ; CHKERRQ(ierr) ;
    ierr = KSPGetIterationNumber(ksp,&no_ksp_its) ;

    if ( no_ksp_its==0 ) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP_Iter_No = 0 Ocurred") ; 
      Finalize() ;
      exit(0) ;
    }

    // Update phi locallly

    ierr = VecGetArray(Vecx,&Vecx_local) ; CHKERRQ(ierr) ;

    norm_local = 0.0 ;
    for ( I=Istart ; I<Iend ; I++ ) {
      p = I+1 ;
      get_ijk(p,&i,&j,&k) ;
      phi[i][j][k] += Vecx_local[I-Istart] ;
      norm_local += phi[i][j][k]*phi[i][j][k] ;
    }
    ierr = VecRestoreArray(Vecx,&Vecx_local) ; CHKERRQ(ierr) ;

    // Check Convergence 

    MPI_Allreduce(&norm_local,&norm,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD) ;
    norm = sqrt(norm) ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Poisson Iteration %d : Norm of phi = %le, KSP_Iter_No = %d\n",iter,norm,no_ksp_its) ; CHKERRQ(ierr) ;

    if ( fabs(norm-norm_old) < norm*CP.ConvEps ) break ;
    else norm_old = norm ;

    // Transfer intermediate phi

    mpi_transfer_phi(phi,node,np,Istart,Iend) ;
  }

  // Finalize 

  int i0,i1 ;
  get_ijk(Istart+1,&i0,&j,&k) ;
  get_ijk(Iend,&i1,&j,&k) ;

  copy_R3mtx(Phi,phi,i0-1,i1+1,N.y0,N.y1,N.z0,N.z1) ; 
  fix_corners(Phi,local_x0,local_x1) ; 

  // Print Midline Potential

  if ( GetPRINT_MidPotFlag() ) {
    char fn[100] ;
    void print_midline_potential_plain() ;

    sprintf(fn,"mpot.r-Vd%.2lf-Vg%.2lf",GetDrainVoltage(),GetGateVoltage()) ;
    print_midline_potential_plain(fn,0) ;
  }

  return 0 ;
}

void mpi_transfer_phi(double ***phi,int node, int np, int Istart, int Iend)
{
  int tag1,tag2,p0,n,ntrans,data_size,sum_transmitted ;
  int trans_size = GetTransSize() ;
  MPI_Status status ;

  int y0 = N.y0-1 ;
  int z0 = N.z0-1 ;
  int y1 = N.y1+1 ;
  int z1 = N.z1+1 ;
  int Ns = (y1-y0+1)*(z1-z0+1) ;

  int i0,j0,k0,i1,j1,k1 ;

  get_ijk(Istart+1,&i0,&j0,&k0) ;
  get_ijk(Iend+1,&i1,&j1,&k1) ;

  data_size = Ns ;
  ntrans = data_size/trans_size ;

  if ( node-1>=0 ) {
    p0 = 0 ;
    sum_transmitted = 0 ;
    tag1 = 1 ;
    tag2 = 2 ;
    for ( n=1 ; n<=ntrans ; n++ ) {
      MPI_Send(&phi[i0][j0][k0]+p0,trans_size,MPI_DOUBLE,node-1,tag1,PETSC_COMM_WORLD) ;
      MPI_Recv(&phi[i0][j0][k0]-Ns+p0,trans_size,MPI_DOUBLE,node-1,tag2,PETSC_COMM_WORLD,&status) ;
      p0 += trans_size ;
      sum_transmitted += trans_size ;
      tag1 += 2 ;
      tag2 += 2 ;
    }
    MPI_Send(&phi[i0][j0][k0]+p0,data_size-sum_transmitted,MPI_DOUBLE,node-1,tag1,PETSC_COMM_WORLD) ;
    MPI_Recv(&phi[i0][j0][k0]-Ns+p0,data_size-sum_transmitted,MPI_DOUBLE,node-1,tag2,PETSC_COMM_WORLD,&status) ;
  }

  if ( node+1<np ) {
    p0 = 0 ;
    sum_transmitted = 0 ;
    tag1 = 1 ;
    tag2 = 2 ; 
    for ( n=1 ; n<=ntrans ; n++ ) {
      MPI_Send(&phi[i1][j1][k1]-Ns+p0,trans_size,MPI_DOUBLE,node+1,tag2,PETSC_COMM_WORLD) ;
      MPI_Recv(&phi[i1][j1][k1]+p0,trans_size,MPI_DOUBLE,node+1,tag1,PETSC_COMM_WORLD,&status) ;
      p0 += trans_size ;
      sum_transmitted += trans_size ;
      tag1 += 2 ;
      tag2 += 2 ;
    }
    MPI_Send(&phi[i1][j1][k1]-Ns+p0,data_size-sum_transmitted,MPI_DOUBLE,node+1,tag2,PETSC_COMM_WORLD) ;
    MPI_Recv(&phi[i1][j1][k1]+p0,data_size-sum_transmitted,MPI_DOUBLE,node+1,tag1,PETSC_COMM_WORLD,&status) ;
  } 
}

void set_alphas_1()
{
  int p,region ;
  double sum ;
  Region R ;

  double e_si = GetESi() ;

  for ( p=1 ; p<=Ntot ; p++ ) {

    region = GetRegion(&R,p) ;

    sum = 0.0 ;
    sum += (alpha_xplus[p] = falpha_xplus(p)) ;
    sum += (alpha_yplus[p] = falpha_yplus(p)) ;
    sum += (alpha_zplus[p] = falpha_zplus(p)) ;
    sum += (alpha_xminus[p] = falpha_xminus(p)) ;
    sum += (alpha_yminus[p] = falpha_yminus(p)) ;
    sum += (alpha_zminus[p] = falpha_zminus(p)) ;
    alpha_center[p] = -sum ;

    if ( GATE_CONTACT || CORNERS ) alpha_center[p] = -1.0 ;

    if ( SILICON ) alpha_q0e[p] = NANOM*NANOM*q0/e_si ;
    else alpha_q0e[p] = 0.0 ;
  }

}

void set_alphas_2()
{
  int p,i,j,k,region ;
  Region R ;

  double ***Phi = GetPot() ;
  double phi_gate = GetPhi_Gate() ;

  for ( p=1 ; p<=Ntot ; p++ ) {

    region = GetRegion(&R,p) ;

    if ( GATE_CONTACT )  alpha_const[p] = phi_gate ;
    else if ( CORNERS ) {
      get_ijk(p,&i,&j,&k) ;
      alpha_const[p] = Phi[i][j][k] ;
    }
    else  alpha_const[p] = 0.0 ;
  }
}

double falpha_xplus(int p)
{
  int i,j,k,region ;
  double *hx,*hy,*hz,e_ox,e_si ;
  Region R ;

  hx = GetHx() ;
  hy = GetHy() ;
  hz = GetHz() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;
  
  if ( BULK ) return(2.0/(hx[i]+hx[i-1])/hx[i]) ;
  if ( SOURCE_WALL ) return(1.0) ;
  if ( SOURCE_CONTACT ) return(1.0) ;
  
  return(0.0) ;
}

double falpha_xminus(int p)
{
  int i,j,k,region ;
  double *hx,*hy,*hz,e_ox,e_si ;
  Region R ;

  hx = GetHx() ;
  hy = GetHy() ;
  hz = GetHz() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;
  
  if ( BULK ) return(2.0/(hx[i]+hx[i-1])/hx[i-1]) ;
  if ( DRAIN_WALL ) return(1.0) ;
  if ( DRAIN_CONTACT ) return(1.0) ;
  
  return(0.0) ;
}

double falpha_yplus(int p)
{
  int i,j,k,region ;
  double *PY,*hy,e_ox,e_si ;
  double hj,hjm1,y,yc,zc ;
  Region R ;

  PY = GetPY() ;
  hy = GetHy() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;

  hj = hy[j] ;
  hjm1 = hy[j-1] ;

  if ( BULK ) return(2.0/(hj+hjm1)/hj) ;
  if ( LEFT_INTERFACE  ) return(e_si/(e_ox*hj)) ;
  if ( RIGHT_INTERFACE ) return(1.0/hj) ;
  if ( LEFT_WALL ) return(1.0) ;
  
  return(0.0) ;
}

double falpha_yminus(int p)
{
  int i,j,k,region ;
  double *PY,*hy,e_ox,e_si ;
  double hj,hjm1,y,yc,zc ;
  Region R ;

  hy = GetHy() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;

  hj = hy[j] ;
  hjm1 = hy[j-1] ;

  if ( BULK ) return(2.0/(hj+hjm1)/hjm1) ;
  if ( LEFT_INTERFACE  ) return(1.0/hjm1) ;
  if ( RIGHT_INTERFACE ) return(e_si/(e_ox*hjm1)) ;
  if ( RIGHT_WALL ) return(1.0) ;
  
  return(0.0) ;
}

double falpha_zplus(int p)
{
  int i,j,k,region ;
  double *PZ,*hz,e_ox,e_si ;
  double hk,hkm1,z,yc,zc ;
  Region R ;

  hz = GetHz() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;

  hk = hz[k] ;
  hkm1 = hz[k-1] ;

  if ( BULK ) return(2.0/(hk+hkm1)/hk) ;
  if ( TOP_INTERFACE  ) return(e_si/(e_ox*hk)) ;
  if ( BOTT_INTERFACE ) return(1.0/hk) ;
  if ( TOP_WALL ) return(1.0) ;

  return(0.0) ;
}

double falpha_zminus(int p)
{
  int i,j,k,region ;
  double *PZ,*hz,e_ox,e_si ;
  double hk,hkm1,z,yc,zc ;
  Region R ;

  hz = GetHz() ;
  e_ox = GetEOx() ;
  e_si = GetESi() ;

  get_ijk(p,&i,&j,&k) ;
  region = GetRegion(&R,p) ;
  
  hk = hz[k] ;
  hkm1 = hz[k-1] ;

  if ( BULK ) return(2.0/(hk+hkm1)/hkm1) ;
  if ( TOP_INTERFACE  ) return(1.0/hkm1) ;
  if ( BOTT_INTERFACE ) return(e_si/(e_ox*hkm1)) ;

  if ( BOTT_WALL ) return(1.0) ;
  
  return(0.0) ;
}

double Fnn(p,nfunc)
     int p ;
     double (*nfunc)() ;
{
  int i,j,k ;
  double fv,***phi,exp(),GetDoping() ;

  phi = GetPhiKth() ;
  get_ijk(p,&i,&j,&k) ;

  fv = alpha_xminus[p]*phi[i-1][j][k] + alpha_xplus[p]*phi[i+1][j][k] \
    + alpha_yminus[p]*phi[i][j-1][k] + alpha_yplus[p]*phi[i][j+1][k] \
    + alpha_zminus[p]*phi[i][j][k-1] + alpha_zplus[p]*phi[i][j][k+1] \
    + alpha_center[p]*phi[i][j][k] + alpha_const[p] \
    + alpha_q0e[p]*(-(*nfunc)(i,j,k)-GetDoping(p)) ;

  return(-fv) ;
}

double DFnn(p,dnfunc)
     int p ;
     double (*dnfunc)() ;
{
  int i,j,k ;
  double ***phi,exp() ;

  phi = GetPhiKth() ;
  get_ijk(p,&i,&j,&k) ;

  return(alpha_center[p]+alpha_q0e[p]*(-(*dnfunc)(i,j,k))) ;
}

void poisson_memory_allocation()
{
  alpha_xplus  = dvector(1,Ntot) ;
  alpha_yplus  = dvector(1,Ntot) ;
  alpha_zplus  = dvector(1,Ntot) ;
  alpha_xminus = dvector(1,Ntot) ;
  alpha_yminus = dvector(1,Ntot) ;
  alpha_zminus = dvector(1,Ntot) ;
  alpha_center = dvector(1,Ntot) ;
  alpha_const  = dvector(1,Ntot) ;
  alpha_q0e    = dvector(1,Ntot) ;
}

void poisson_free_memory()
{
  free_dvector(alpha_xplus, 1,Ntot) ;
  free_dvector(alpha_yplus, 1,Ntot) ;
  free_dvector(alpha_zplus, 1,Ntot) ;
  free_dvector(alpha_xminus,1,Ntot) ;
  free_dvector(alpha_yminus,1,Ntot) ;
  free_dvector(alpha_zminus,1,Ntot) ;
  free_dvector(alpha_center,1,Ntot) ;
  free_dvector(alpha_const, 1,Ntot) ;
  free_dvector(alpha_q0e,   1,Ntot) ;
}

#define tiny 1.0e-06

int GetRegion(Region *PR, int p )
{
  int i,j,k,l,jj,kk,region ;
  double y,z,yc,zc,radius_cr ;
  double *PY = GetPY() ;
  double *PZ = GetPZ() ;
  Dimen D = GetDimen() ;
  PoiParam P = GetPoiParam() ;
  Region R ;
  int GetRegion_Bulk() ;

  *PR = GetR_() ;
  R = *PR ;
  get_ijk(p,&i,&j,&k) ;

  if ( (i==N.x0 && j==N.y0) || (i==N.x0 && j==N.y1) || \
       (i==N.x0 && k==N.z0) || (i==N.x0 && k==N.z1) || \
       (i==N.x1 && j==N.y0) || (i==N.x1 && j==N.y1) || \
       (i==N.x1 && k==N.z0) || (i==N.x1 && k==N.z1) || \
       (j==N.y0 && k==N.z0) || (j==N.y0 && k==N.z1) || \
       (j==N.y1 && k==N.z0) || (j==N.y1 && k==N.z1) ) {

    region = R._CORNERS ;
    return(region) ;
  }

  if ( i==N.x0 ) { // SOURCE
    if ( j<N.ya || j>N.yb || k<N.za || k>N.zb ) {
      region = R._SOURCE_WALL ;
      return(region) ;
    }
    else region = R._SOURCE_CONTACT ;
  }

  if ( i==N.x1 ) { // DRAIN
    if ( j<N.ya || j>N.yb || k<N.za || k>N.zb ) {
      region = R._DRAIN_WALL ;
      return(region) ;
    }
    else region = R._DRAIN_CONTACT ;
  }

  if ( i<N.xgate0 || i>N.xgate1 ) {
    if ( j==N.y0 ) { region = R._LEFT_WALL ; return(region) ; }
    if ( j==N.y1 ) { region = R._RIGHT_WALL ; return(region) ; }
    if ( k==N.z0 ) { region = R._TOP_WALL ; return(region) ; }
    if ( k==N.z1 ) { region = R._BOTT_WALL ; return(region) ; }
  }
  
  if ( i>=N.xgate0 && i<=N.xgate1 ) { // GATE
    if ( k==N.z0 ) { 
      region = R._GATE_CONTACT ; 
      return(region) ; 
    }
    if ( k==N.z1 ) {
      if ( j<=N.ygate0 || j>=N.ygate1 ) region = R._GATE_CONTACT ;
      else region = R._BOTT_WALL ;
      return(region) ;
    }
    if ( j==N.y0 ) {
      if ( k<=N.zgate ) region = R._GATE_CONTACT ;
      else region = R._LEFT_WALL ;
      return(region) ;
    }
    if ( j==N.y1 ) {
      if ( k<=N.zgate ) region = R._GATE_CONTACT ;
      else region = R._RIGHT_WALL ;
      return(region) ;
    } 
  }

  // BULK

  region = GetRegion_Bulk(i,j,k) ;

  if ( i==N.x0 ) {
    if ( OXIDE ) region = R._SOURCE_WALL ;
    else         region = R._SOURCE_CONTACT ;
  }
  else if ( i<N.xa ) {
    if ( CHANNEL_SILICON ) region = R._JUNCTION_SILICON ;
  }
  else if ( i>N.xb && i<N.x1 ) {
    if ( CHANNEL_SILICON ) region = R._JUNCTION_SILICON ;
  }
  else if ( i==N.x1 ) {
    if ( OXIDE ) region = R._DRAIN_WALL ;
    else         region = R._DRAIN_CONTACT ;
  }

  return(region) ;
}

int GetRegion_Bulk(int i, int j, int k)
{
  int region,jj,kk,l ;
  Region R = GetR_() ;

  if ( k==N.za && (j>N.ya && j<N.yb) ) { region = R._TOP_INTERFACE ; return(region) ; } 
  else if ( k==N.zb && (j>N.ya && j<N.yb) ) { region = R._BOTT_INTERFACE ; return(region) ; }
  else if ( j==N.ya && (k>N.za && k<N.zb) ) { region = R._LEFT_INTERFACE ; return(region) ; }
  else if ( j==N.yb && (k>N.za && k<N.zb) ) { region = R._RIGHT_INTERFACE ; return(region) ; }
  else if ( k<=N.za || k>=N.zb || j<=N.ya || j>=N.yb ) { region = R._OXIDE ; return(region) ; }
  else {
    region = R._CHANNEL_SILICON ;
    return(region) ;
  }
}

#undef tiny

void get_ijk(p,iv,jv,kv)
     int p,*iv,*jv,*kv ;
{
  int pp ;
  
  *iv = (p-1)/(N.y*N.z) ;
  pp = p-(*iv)*N.y*N.z ;
  *jv = (pp-1)/N.z ;
  *kv = pp-(*jv)*N.z-1 ;
}

int get_p(int i, int j, int k)
{
  return(i*N.y*N.z+j*N.z+k+1) ;
}

void set_null_out_of_bound()
{
  int i,j,k ;
  double ***Phi = GetPot() ;

  for ( j=N.y0-1 ; j<=N.y1+1 ; j++ ) 
    for ( k=N.z0-1 ; k<=N.z1+1 ; k++ ) {
      Phi[N.x0-1][j][k] = 0.0 ;
      Phi[N.x1+1][j][k] = 0.0 ;
    }

  for ( i=N.x0-1 ; i<=N.x1+1 ; i++ )
    for ( k=N.z0-1 ; k<=N.z1+1 ; k++ ) {
      Phi[i][N.y0-1][k] = 0.0 ;
      Phi[i][N.y1+1][k] = 0.0 ;
    }

  for ( i=N.x0-1 ; i<=N.x1+1 ; i++ )
    for ( j=N.y0-1 ; j<=N.y1+1 ; j++ ) {
      Phi[i][j][N.z0-1] = 0.0 ;
      Phi[i][j][N.z1+1] = 0.0 ;
    }
}

void fix_corners(double ***Phi, int local_x0, int local_x1)
{
  int i,j,k ;
  
  // Taking care of Corner Points

  if ( local_x0 == N.x0 ) {
    for ( j=N.y0 ; j<=N.y1 ; j++ ) {
      Phi[local_x0][j][N.z0] = Phi[local_x0][j][N.z0+1] ;
      Phi[local_x0][j][N.z1] = Phi[local_x0][j][N.z1-1] ;
    }
    for ( k=N.z0 ; k<=N.z1 ; k++ ) {
      Phi[local_x0][N.y0][k] = Phi[local_x0][N.y0+1][k] ;
      Phi[local_x0][N.y1][k] = Phi[local_x0][N.y1-1][k] ;
    }
  }
  
  if ( local_x1 == N.x1 ) { 
    for ( j=N.y0 ; j<=N.y1 ; j++ ) {
      Phi[local_x1][j][N.z0] = Phi[local_x1][j][N.z0+1] ;
      Phi[local_x1][j][N.z1] = Phi[local_x1][j][N.z1-1] ;
    }
    for ( k=N.z0 ; k<=N.z1 ; k++ ) {
      Phi[local_x1][N.y0][k] = Phi[local_x1][N.y0+1][k] ;
      Phi[local_x1][N.y1][k] = Phi[local_x1][N.y1-1][k] ;
    }
  }

  for ( i=local_x0 ; i<=local_x1 ; i++ ) {
    Phi[i][N.y0][N.z0] = 0.5*(Phi[i][N.y0+1][N.z0]+Phi[i][N.y0][N.z0+1]) ;
    Phi[i][N.y0][N.z1] = 0.5*(Phi[i][N.y0+1][N.z1]+Phi[i][N.y0][N.z1-1]) ;
    Phi[i][N.y1][N.z0] = 0.5*(Phi[i][N.y1-1][N.z0]+Phi[i][N.y1][N.z0+1]) ;
    Phi[i][N.y1][N.z1] = 0.5*(Phi[i][N.y1-1][N.z1]+Phi[i][N.y1][N.z1-1]) ;
  }
}
