#ifndef _CONSTANT_ 
#define _CONSTANT_
	/* Physical Constant*/
	#define pi    3.14159265358979 /* pi number viet 3.14 co nghia la kieu du lieu la double */
	#define hbar  1.05459e-34      /* Planck's constant/2PI in J-sec */
	#define m0    9.10953e-31      /* Electron rest mass in kg */
	#define q     1.60219e-19      /* Electron charge */
	#define kb    1.3806505e-23    /* Boltzmann constant J/Kelvin */
	#define eps_0 8.85419e-12      /* Permittivity of free space. unit F/m */

	/* For Monte Carlo */
	#define max_electron_number  200000 /* maximum number of simulated electrons */
	#define max_subband_number   30    // tinh trong 1 valley 
	#define averaging_time       5.0e-13   /* [s] time to average/print data*/
	
	/* For 2D Schrodinger */
	#define beta   13.123288895  /* unit [1/eVnm**2], beta=m0/(hbar*hbar)*q*(1.0e-9*1.0e-9);

	/* Parameters of scattering table  */
	#define n_lev  1000 // # of energy levels in the scattering table
	#define n_scatt_max  10  // maximum # of scattering mechanisms
	
	/* For device geometry */
	#define nx_size  400    // number of points in x-direction 
	#define ny_size  200  // number of points in y-direction
	#define nz_size  200  // number of points in y-direction
	// Mac du doc tu Input file gia tri n_gate nhung ta chon max cua no la ny_size_negative
	
	/* De tao cac logfile: 1 Save logfie, KHAC 1 la Khong save */
	#define ENABLE_LOG 1

#endif //_CONSTANT_	







