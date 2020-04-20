/* ****************************************************************************
This function is to INITIALIZE DEVICE STRUCTURE
	(1) defines source, drain, gate, etc. regions
	(2) defines different doping regions
	(3) defines mesh size
	(4) defines initial doping 
	       
Starting date: Feb 09, 2010
Update:        Feb 11, 2010
Latest update: April 29, 2010 (Them tham so de mapping voi 3D Poisson)
***************************************************************************** */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "constants.h"
#include "petscksp.h" // for Poisson 3D
#include "nanowire.h" // for Poisson 3D

static double mesh_size_x,mesh_size_y,mesh_size_z;
static int Tox_mesh,Tbox_mesh,Wox_mesh; // do la number of points
static int nx0, nxa, nxb, nx1;// so diem chia o tung vung cho truc x
static int ny0, nya, nyb, ny1;// truc y
static int nz0, nza, nzb, nz1;// truc z
static int nxgate0, nxgate1;  // so diem chia cho vi tri cua gate
static double *doping_density;

/* ****************** Reading device structure ********************************** */
void device_structure(){
    // Goi ham cho Poisson 3D
    void SetDopingDensityFor(int convflag, char region, char doping_type, double doping_density);
    double GetDopingDensityFor(char region);
     
    char device_type[40];// NanowireFET
    char cvar; // for Poisson 3D. doping dang n, p hay i
    double LsdInput,LchInput,LgInput,ToxInput,TboxInput,WoxInput,TsiInput,WsiInput;// tu Input file la nm
    double Lsd,Lch,Lg,Tox,Tbox,Wox,Tsi,Wsi; // chuyen tu [nm] -> [m] 
    double TgateInput, WgateInput; // for Poisson 3D
    int n_region;
    doping_density=dvector(1,3); //[1] S, [2] D; [3] Channel
    double *doping_densityInput = dvector(1,3);
    /* ******** READING DEVICE PARAMETER ********* */
    FILE *f;
    char s[180];
    f = fopen("Input.d","r");
    if(f == NULL){printf("Can't open file Input.d\n"); return;}
    //printf("\n\n Reading #DEVICE_PARAMETER");
    fscanf(f,"%s",s); // read the current row
	if(strcmp(s,"#DEVICE_PARAMETER")==0){
		
      	fscanf(f,"%*s %*s %s",device_type);// Xem kieu device_type la loai gi !
      	//printf("\n Device type is           %s",device_type); // Hien tai chi la NanowireFET

      	fscanf(f,"%*s %*s %lf",&LsdInput);
      	//printf("\n Source/Drain length [nm] = %5.2lf", LsdInput);
      	Lsd = LsdInput*1.0e-9; // [nm] -> [m]

        fscanf(f,"%*s %*s %lf",&LchInput);
        //printf("\n Channel length [nm]      = %5.2lf", LchInput);
        Lch = LchInput*1.0e-9; // [m]

      	fscanf(f,"%*s %*s %lf",&LgInput);
      	//printf("\n Gate length[nm]          = %5.2lf", LgInput);
      	Lg = LgInput*1.0e-9; // [m]

      	fscanf(f,"%*s %*s %lf",&ToxInput);
      	//printf("\n Tox [nm]                 = %5.2lf", ToxInput);
      	Tox = ToxInput*1.0e-9; // [m]

      	fscanf(f,"%*s %*s %lf",&TboxInput);
      	//printf("\n Tbox [nm]                = %5.2lf", TboxInput);
      	Tbox = TboxInput*1.0e-9; // [m]

      	fscanf(f,"%*s %*s %lf",&WoxInput); 
      	//printf("\n Wox [nm]                 = %5.2lf",WoxInput);
      	Wox = WoxInput*1.0e-9; // [m]
      	
      	fscanf(f,"%*s %*s %lf",&TsiInput);
      	//printf("\n Tsi [nm]                 = %5.2lf",TsiInput);
      	Tsi = TsiInput*1.0e-9; // [m]
      	
      	fscanf(f,"%*s %*s %lf",&WsiInput);
      	//printf("\n Wsi [nm]                 = %5.2lf",WsiInput);
      	Wsi = WsiInput*1.0e-9; // [m]
      	
      	fscanf(f,"%*s %*s %lf",&TgateInput);
      	//printf("\n T_gate [nm]              = %5.2lf", TgateInput);
      	// Khong chuyen sang [m] vi chi dung de mapping cho Poisson 3D ma thoi
      	
      	fscanf(f,"%*s %*s %lf",&WgateInput);
      	//printf("\n W_gate [nm]              = %5.2lf", WgateInput);
      	// Khong chuyen sang [m] vi chi dung de mapping cho Poisson 3D ma thoi
      	
      	fscanf(f,"%*s %*s %lf",&mesh_size_x); // dinh nghia do rong mesh size theo phuong x
      	//printf("\n mesh_size_x [nm]         = %5.2lf",mesh_size_x);
      	mesh_size_x = mesh_size_x*1.0e-9; // [m]

      	fscanf(f,"%*s %*s %lf",&mesh_size_y);
      	//printf("\n mesh_size_y [nm]         = %5.2lf",mesh_size_y);
      	mesh_size_y = mesh_size_y*1.0e-9; // [m]
      	
      	fscanf(f,"%*s %*s %lf",&mesh_size_z);
      	//printf("\n mesh_size_z [nm]         = %5.2lf",mesh_size_z);
      	mesh_size_z = mesh_size_z*1.0e-9; // [m]
      	
      	fscanf(f,"%*s %*s %d",&Tox_mesh);
      	//printf("\n Number of point Tox_mesh = %d",Tox_mesh);
      
      	fscanf(f,"%*s %*s %d",&Tbox_mesh);
      	//printf("\n Number of point Tbox_mesh= %d",Tbox_mesh);
      	
      	fscanf(f,"%*s %*s %d",&Wox_mesh);
      	//printf("\n Number of point Wox_mesh = %d",Wox_mesh);

      	fscanf(f,"%*s %*s %d",&n_region); //total_number_of_active_region
      	//printf("\n n_region                 = %d",n_region);
      	
      	// Doping at Source
        fscanf(f,"%*s %*s %c", &cvar);
      	//printf("\n Source doping type       = %c",cvar);
      	
        fscanf(f,"%*s %*s %lf",&doping_densityInput[1]); //Source
      	//printf("\n Source doping </cm3>     = %5.1le",doping_densityInput[1]);
      	doping_density[1] = doping_densityInput[1]*1.0e+6; // [/cm3] -> [1/m3]
      	SetDopingDensityFor(1,'s',cvar,doping_densityInput[1]);// cho Poisson 3D

      	// Doping at Drain
      	fscanf(f,"%*s %*s %c", &cvar);
      	//printf("\n Drain doping type        = %c",cvar);
      	
        fscanf(f,"%*s %*s %lf",&doping_densityInput[2]); // Drain
      	//printf("\n Drain doping </cm3>      = %5.1le",doping_densityInput[2]);
      	doping_density[2] = doping_densityInput[2]*1.0e+6; // [/cm3] -> [1/m3]
      	SetDopingDensityFor(1,'d',cvar,doping_densityInput[2]);// Cho Poisson 3D
       
        // Doping at Channel
        fscanf(f,"%*s %*s %c", &cvar);
      	//printf("\n Channel doping type      = %c",cvar);
      	
      	fscanf(f,"%*s %*s %lf",&doping_densityInput[3]); // Channel
      	//printf("\n Channel doping </cm3>    = %5.2le",doping_densityInput[3]);
      	doping_density[3] = doping_densityInput[3]*1.0e+6;// [/cm3] -> [1/m3]
      	SetDopingDensityFor(1,'c',cvar,doping_densityInput[3]);// Cho 3D Poisson
      	// Cua GS thi doping gia tri o Input file la + hay - thi luc vao deu duoc fabs()
 
    	fclose(f);
    } // End of if(strcmp(s,"#DEVICE_PARAMETER")==0){
      /* ******** END OF Reading #DEVICE_PARAMETER********* */
  	/* Specify reference points for source, drain and channel regions
   		and junction depth and oxide region (bellow gate), ect. */

    // Theo phuong x
    nx0 = 0; // goc TAO DO
    nxa = nx0 + (int)(Lsd/mesh_size_x + 0.5);// (int)(50.8)=50
    nxb = nxa + (int)(Lch/mesh_size_x + 0.5);
    nx1 = nxb + (int)(Lsd/mesh_size_x + 0.5);
    //printf("\n nx0=%d, nxa=%d, nxb=%d, nx1=%d",nx0, nxa, nxb, nx1);
    
    // Theo phuong y
    ny0 = 0; // goc toa do
    nya = Wox_mesh;
    nyb = nya + (int)(Wsi/mesh_size_y + 0.5);
    ny1 = nyb + Wox_mesh;
    //printf("\n ny0=%d, nya=%d, nyb=%d, ny1=%d",ny0, nya, nyb, ny1);
    
    // Theo phuong z
    nz0 = 0; // goc toa do
    nza = Tox_mesh;
    nzb = nza + (int)(Tsi/mesh_size_z + 0.5);
    nz1 = nzb + Tbox_mesh;
    //printf("\n nz0=%d, nza=%d, nzb=%d, nz1=%d",nz0, nza, nzb, nz1);

    // Toa do cho gate
    double x3 = Lsd+Lch+Lsd; // do dai device
    double x2 = x3 - Lg;
    nxgate0 = (int)(0.5*x2/mesh_size_x + 0.5);
    double x4 = x3 + Lg;
    nxgate1 = (int)(0.5*x4/mesh_size_x + 0.5);
    //printf("\n nxgate0=%d, nxgate1=%d",nxgate0, nxgate1);
    // Khi Lg=Lch thi nxgate0=nxa, nxgate1=nxb: Cach Kiem tra cong thuc tinh dung hay sai
    
    // Mapping for 3D Poisson
    Dimen *D;
    PoiNum *N;
    D = PGetDimen();
    N = PGetPoiNum();
    
    D->Lsrc       = LsdInput; // do dai Source
    N->Mx_src     = (int)(Lsd/mesh_size_x + 0.5);     // so diem chia o Source
    
    D->Lchannel   = LchInput; // do dai channel length
    N->Mx_channel = (int)(Lch/mesh_size_x + 0.5); // so diem chia cho Lch

    D->Tox   = ToxInput;
    N->Mz_ox = Tox_mesh;
    
    D->Tsi   = TsiInput;
    N->Mz_si = (int)(Tsi/mesh_size_z + 0.5);
    
    D->Tbox   = TboxInput;
    N->Mz_box = Tbox_mesh;

    D->Wox   = WoxInput;
    N->My_ox = Wox_mesh;
    
    D->Wsi   = WsiInput;
    N->My_si = (int)(Wsi/mesh_size_y + 0.5);

    D->Lgate = LgInput;
    D->Tgate = TgateInput;
    D->Wgate = WgateInput;
    
    // Kiem tra mapping for 3D Poisson
    //printf("\n Checking mapping for 3D Poisson at device_structure()");
    //printf("\n Doping at Source[1/m3] =%le",GetDopingDensityFor('s'));
    //printf("\n Doping at Channel[1/m3]=%le",GetDopingDensityFor('c'));
    //printf("\n Doping at Drain[1/m3]  =%le",GetDopingDensityFor('d'));
    //printf("\n Lsource[NoUnit] D->Lsrc     =%f,  Number of points N->Mx_src    =%d",D->Lsrc,N->Mx_src);
    //printf("\n Lchannel[NoUnit]D->Lchannel =%f, Number of points N->Mx_channel =%d",D->Lchannel,N->Mx_channel);
    //printf("\n Tox[NoUnit] D->Tox          =%f,   Number of points N->Mz_ox    =%d",D->Tox,N->Mz_ox);
    //printf("\n Tsi[NoUnit] D->Tsi          =%f,   Number of points N->Mz_si    =%d",D->Tsi,N->Mz_si);
    //printf("\n Tbox[NoUnit] D->Tbox        =%f,   Number of points N->Mz_box   =%d",D->Tbox,N->Mz_box);
    //printf("\n Wox[NoUnit] D->Wox          =%f,   Number of points N->My_ox    =%d",D->Wox,N->My_ox);
    //printf("\n Wsi[NoUnit] D->Wsi          =%f,   Number of points N->My_si    =%d",D->Wsi,N->My_si);
    //printf("\n Lgate[NoUnit] D->Lgate      =%f",D->Lgate);
    //printf("\n Tgate[NoUnit] D->Tgate      =%f",D->Tgate);
    //printf("\n Wgate[NoUnit] D->Wgate      =%f",D->Wgate);

    // Phan in ra o monitor khi rank=0
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    
    if(myrank ==0){ 
        printf("\n\n Reading DEVICE_PARAMETER");
    	printf("\n Device type is           %s",device_type); // Hien tai chi la NanowireFET
      	printf("\n Source/Drain length [nm] = %5.2lf", LsdInput);
        printf("\n Channel length [nm]      = %5.2lf", LchInput);
       	printf("\n Gate length[nm]          = %5.2lf", LgInput);
       	printf("\n Tox [nm]                 = %5.2lf", ToxInput);
       	printf("\n Tbox [nm]                = %5.2lf", TboxInput);
       	printf("\n Wox [nm]                 = %5.2lf",WoxInput);
       	printf("\n Tsi [nm]                 = %5.2lf",TsiInput);
       	printf("\n Wsi [nm]                 = %5.2lf",WsiInput);
      	printf("\n T_gate [nm]              = %5.2lf", TgateInput);
      	printf("\n W_gate [nm]              = %5.2lf", WgateInput);
       	printf("\n mesh_size_x [nm]         = %5.2lf",mesh_size_x/1.0e-9);
       	printf("\n mesh_size_y [nm]         = %5.2lf",mesh_size_y/1.0e-9);
      	printf("\n mesh_size_z [nm]         = %5.2lf",mesh_size_z/1.0e-9);
      	printf("\n Number of point Tox_mesh = %d",Tox_mesh);
       	printf("\n Number of point Tbox_mesh= %d",Tbox_mesh);
       	printf("\n Number of point Wox_mesh = %d",Wox_mesh);
      	printf("\n n_region                 = %d",n_region);
       	printf("\n Source doping </cm3>     = %5.1le",doping_densityInput[1]);
      	printf("\n Drain doping </cm3>      = %5.1le",doping_densityInput[2]);
       	printf("\n Channel doping </cm3>    = %5.2le",doping_densityInput[3]);
	printf("\n nx0=%d, nxa=%d, nxb=%d, nx1=%d",nx0, nxa, nxb, nx1);
	printf("\n ny0=%d, nya=%d, nyb=%d, ny1=%d",ny0, nya, nyb, ny1);
	printf("\n nz0=%d, nza=%d, nzb=%d, nz1=%d",nz0, nza, nzb, nz1);
	printf("\n nxgate0=%d, nxgate1=%d",nxgate0, nxgate1);
       
    }// End of if(myrank ==0)
    
   //Free local variable
   free_dvector(doping_densityInput,1,3);
    
return;
}
// End of void device_structure(){

double Get_mesh_size_x(){
       return mesh_size_x;
}

double Get_mesh_size_y(){
       return mesh_size_y;
}

double Get_mesh_size_z(){
       return mesh_size_z;
}

int Get_Tox_mesh(){
    return Tox_mesh;
}

int Get_Tbox_mesh(){
    return Tbox_mesh;
}

int Get_Wox_mesh(){
    return Wox_mesh;
}

int Get_nx0(){
    return nx0;
}

int Get_nxa(){
    return nxa;
}

int Get_nxb(){
    return nxb ;
}

int Get_nx1(){
    return nx1 ;
}

int Get_ny0(){
    return ny0 ;
}

int Get_nya(){
    return nya ;
}

int Get_nyb(){
    return nyb ;
}

int Get_ny1(){
    return ny1 ;
}

int Get_nz0(){
    return nz0 ;
}

int Get_nza(){
    return nza ;
}

int Get_nzb(){
    return nzb ;
}

int Get_nz1(){
    return nz1 ;
}

int Get_nxgate0(){
    return nxgate0 ;
}

int Get_nxgate1(){
    return nxgate1 ;
}

// Chu y cho doping
double *Get_doping_density(){
       return (doping_density);
}
       
   
       


