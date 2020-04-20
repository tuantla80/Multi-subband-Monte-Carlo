#include "petscksp.h"
#include "nanowire.h"
#include "region.h"
#include "nrutil.h"

double Func_Nq(int i,int j,int k)
{
  double nq() ;
  return(nq(i,j,k)) ;
}

double Func_DNq(int i,int j,int k)
{
  double Dnq() ;
  return(Dnq(i,j,k)) ;
}

double nq(int i,int j,int k)
{
  double exp() ;

  double kT = GetKT() ;
  double ***phi = GetPhiKth() ;
  double ***phi_old = GetPot() ;
  double ***n3d = GetN3d() ;
  double ***p3d = GetP3d() ;
  double exp_beta = exp((phi[i][j][k]-phi_old[i][j][k])/kT) ;

  double x = n3d[i][j][k]*exp_beta - p3d[i][j][k]/exp_beta ;

  return(x) ;
}

double Dnq(int i,int j,int k)
{
  double exp() ;

  double kT = GetKT() ;
  double ***phi = GetPhiKth() ;
  double ***phi_old = GetPot() ;
  double ***n3d = GetN3d() ;
  double ***p3d = GetP3d() ;
  double exp_beta = exp((phi[i][j][k]-phi_old[i][j][k])/kT) ;

  double x= (n3d[i][j][k]*exp_beta + p3d[i][j][k]/exp_beta)/kT ;

  return(x) ;
}

double nq_linearized(int p)
{
  double kT = GetKT() ;
  double ***n3d = GetN3d() ;
  double ***p3d = GetP3d() ;

  int i,j,k ;
  get_ijk(p,&i,&j,&k) ;

  return((n3d[i][j][k]+p3d[i][j][k])/kT) ;
}

double nq_linearized2(int p)
{
  double GetDoping() ;

  double kT = GetKT() ;
  double ***n3d = GetN3d() ;
  double ***p3d = GetP3d() ;
  double ***phi = GetPot() ;

  int i,j,k ;
  get_ijk(p,&i,&j,&k) ;

  return(-(n3d[i][j][k]+p3d[i][j][k])*phi[i][j][k]/kT+(n3d[i][j][k]-p3d[i][j][k])+GetDoping(p)) ;
}
/******************************************************************************
 Thay doi 1 vai cai de mapping code cua ta voi 3D Poisson
 INPUT: double ***electron_density = Get_electron_density()cua ta dua vao  
        double ***n3d = GetN3d() cua Poisson
 De tu day trong ham Poisson thi electron_density la n3d. Co thanh phan Nd va Na
 la doping o S, D va Channel thi da duoc add Getdoping va tu do goi ham Fnn() cho solve_3d_poisson()
 hoac nq_linearized2() cho solve_3d_poisson_linearized()
 INPUT: int flag = 0 : Equilibrium Poisson
                 !=0 : Giai cho Non equilibrium Poisson tuc la qua trinh selfconsistently
 NOTE: Tuy vay ta TACH rieng ra 2 ham cho do lang nhang
 
 Starting date: May 04, 2010
 Latest update: May 04, 2010 
******************************************************************************/
void get_charge_density_for_Poisson()// Non equilibrium Poisson
{
  // Goi ham
  double ***Get_electron_density(); // ham cua ta
  //int find_region(int i,int j,int k);// neu trong vung S,D,Channel thi tinh duoc charge con neu o Oxide thi minh cho charge=0 chu
  // Bien local
  double ***electron_density = Get_electron_density();
  
  int Get_nya(),Get_nza();
  int nya = Get_nya();int nza = Get_nza();

  int np,node;
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  double ***n3d = GetN3d() ;
  double ***p3d = GetP3d();

  if( node==0 ) {
    int i,j,k; //n=0;
    PoiNum N = GetPoiNum();
    for ( i=N.x0 ; i<=N.x1 ; i++ ){
      for ( j=N.y0 ; j<=N.y1 ; j++ ){
         for ( k=N.z0 ; k<=N.z1 ; k++ ){
	   //n = find_region(i,j,k);
	   // if((n==1)||(n==2)||(n==3)){// S,D,Channel;
	   //  n3d[i][j][k] = electron_density[i][j-nya][k-nza]; // lay tu ham electron_density_calculation()
	      // CHU Y toa de cua mang n3d va electron_densuty la khac nhau. Can mapping voi nhau
	      // chi so i thi TRUNG NHAU, nhung chi so j va k thi KHAC NHAU do electron_density chi dinh nghia trong loi Silicon
	   //  }
	   //  else{ //   n3d[i][j][k] =0.0; }

	   n3d[i][j][k] = electron_density[i][j][k];
	   p3d[i][j][k] = 0.0; // do ta KHONG simulate hole
         }
      }
    }
  }// End of if(node==0)

  /*
   // Phan save de check o day
   if( node==0 ) {
    int i,j,k;
    PoiNum N = GetPoiNum();
    FILE *f; f = fopen("charge_density_for_Poisson.dat","w"); // "w" neu tep ton tai no se bi xoa
     if(f==NULL){ printf("\n Cannot open file charge_density_for_Poisson.dat"); return ; }
     fprintf(f,"\n #i j k  n3d[1/m3] p3d[1/m3] \n");
     for ( i=N.x0 ; i<=N.x1 ; i++ ){
      for ( j=N.y0 ; j<=N.y1 ; j++ ){
         for ( k=N.z0 ; k<=N.z1 ; k++ ){
	   fprintf(f," %d %d %d %le %le \n", i,j,k,n3d[i][j][k], p3d[i][j][k]);
	 }
      }
     }
     fclose(f);
   }

   */
 
  void mpi_distribute_nq_from_node0();// goi ham nay thoi
  mpi_distribute_nq_from_node0(n3d,node,np) ;
  mpi_distribute_nq_from_node0(p3d,node,np) ;

}// End of void get_charge_density_for_Poisson()// Non equilibrium Poisson
//**************************************************************************************
// Starting date: May  04, 2010
// Latest update: June 10, 2010 

void get_initial_charge_density_for_Poisson()// charge density tinh tu initial potential
{
   
  int find_region(int i,int j,int k);// neu trong vung S,D,Channel thi tinh duoc charge con neu o Oxide thi minh cho charge=0 chu

  double ***Phi = GetPot();
 
  double Get_Vt();
  double Vt = Get_Vt();

  // SAI double Get_intrinsic_carrier_density();//  double Ni = Get_intrinsic_carrier_density();
  double Nc_Si = 3.22*1.0e+25;//[m3] // Tai room temperature Nc cua Silicon la 3.22e+19/cm3.Semiconductor Device Fundamental - F.Pierret (page 51)
  double Ni = Nc_Si;
  
  int np,node;
  MPI_Comm_rank(PETSC_COMM_WORLD,&node) ;
  MPI_Comm_size(PETSC_COMM_WORLD,&np) ;

  double ***n3d = GetN3d();
  double ***p3d = GetP3d();

  double ***Get_electron_density(); // electron density lan dau tien la tinh tu first guess potential
  double ***electron_density = Get_electron_density();

  double GetDrainVoltage();
  double Vd = GetDrainVoltage();
  Energy *E = PGetEnergy() ;
  double Fermihalf (double x);

  if( node==0 ) {// Tinh tai 1 node roi sau do distribute
    int i,j,k,n=0;
    PoiNum N = GetPoiNum();
    for ( i=N.x0 ; i<=N.x1 ; i++ ){
      for ( j=N.y0 ; j<=N.y1 ; j++ ){
         for ( k=N.z0 ; k<=N.z1 ; k++ ){
	     n = find_region(i,j,k);
	     
	     if(n==1){ //S
	       // Boltzmann statistics
	       n3d[i][j][k] = Ni*exp((-Phi[i][j][k]+E->biS)/Vt);// vi tri toa de cua Phi va n3d la TUONG DUONG NHAU du Phi co them 2 index 2 dau
	       
	       //Fermi-Dirac statistics
	       //n3d[i][j][k] = Ni*Fermihalf((-Phi[i][j][k]+E->biS)/Vt);

	       electron_density[i][j][k] = n3d[i][j][k];// de ve ra
	     }

	     else if (n==2){//D
	       //Boltzmann statistics
	       n3d[i][j][k] = Ni*exp((-Phi[i][j][k]+Vd+ E->biD)/Vt);

	       //Fermi-Dirac statistics
	       //n3d[i][j][k] = Ni*Fermihalf((-Phi[i][j][k]+Vd+ E->biD)/Vt);

	       electron_density[i][j][k] =n3d[i][j][k];// de ve ra
	     }

	     else if (n==3){//Ch
	       //Boltzmann statistics
	       n3d[i][j][k] = Ni*exp((-Phi[i][j][k])/Vt);

	       //Fermi-Dirac statistics
	       //n3d[i][j][k] = Ni*Fermihalf((-Phi[i][j][k])/Vt);

	       electron_density[i][j][k] = n3d[i][j][k] ;// de ve ra
	      
	     }
	     else {// o Oxide
	        n3d[i][j][k] = 0.0;// Khong the gop lai vi e mu 0 bang 1
		electron_density[i][j][k]=0.0;
	     }

             p3d[i][j][k] = 0.0; // do ta KHONG simulate hole

	 }
      }
    }
  
  }// End of if( node==0 )

  void mpi_distribute_nq_from_node0();// goi ham nay thoi
  mpi_distribute_nq_from_node0(n3d,node,np) ;
  mpi_distribute_nq_from_node0(p3d,node,np) ;
}
// End of void get_initial_charge_density_for_Poisson()


void mpi_distribute_nq_from_node0(double ***nq, int node, int np) 
{
  int x0,y0,z0,x1,y1,z1 ;
  int n,dest,datasize,ntrans,p0,sum_transmitted,tag ;
  int trans_size = GetTransSize() ;
  PoiNum N = GetPoiNum() ;
  MPI_Status status ;

  x0 = N.x0 ; x1 = N.x1 ;
  y0 = N.y0 ; y1 = N.y1 ;
  z0 = N.z0 ; z1 = N.z1 ;

  datasize = (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ;
  ntrans = datasize/trans_size ;

  if ( node==0 ) {
    for ( dest=1 ; dest<np ; dest++ ) {
      p0 = 0 ;
      sum_transmitted = 0 ;
      tag = 0 ;
      for ( n=1 ; n<=ntrans ; n++ ) {
        MPI_Send(&nq[x0][y0][z0]+p0,trans_size,MPI_DOUBLE,dest,++tag,PETSC_COMM_WORLD) ;
        p0 += trans_size ;
        sum_transmitted += trans_size ;
      }
      if ( sum_transmitted < datasize )
        MPI_Send(&nq[x0][y0][z0]+p0,datasize-sum_transmitted,MPI_DOUBLE,dest,++tag,PETSC_COMM_WORLD) ;
    }
  }
  else {
    p0 = 0 ;
    sum_transmitted = 0 ;
    tag = 0 ;
    
    for ( n=1 ; n<=ntrans ; n++ ) {
      MPI_Recv(&nq[x0][y0][z0]+p0,trans_size,MPI_DOUBLE,0,++tag,PETSC_COMM_WORLD,&status) ;
      p0 += trans_size ;
      sum_transmitted += trans_size ;
    }
    if ( sum_transmitted < datasize )
      MPI_Recv(&nq[x0][y0][z0]+p0,datasize-sum_transmitted,MPI_DOUBLE,0,++tag,PETSC_COMM_WORLD,&status) ;
  }
}

