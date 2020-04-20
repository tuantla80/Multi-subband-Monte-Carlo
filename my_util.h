#include <stdlib.h>

int nint() ;
int imax() ;
int imin() ;
double dmax() ;
double dmin() ;
void i_interchange() ;
void d_interchange() ;
char *sword() ;
char **svector() ;
double **dmatrix_mpi() ;
double ***d3matrix_mpi() ;
double ****d4matrix() ;
double ****d4matrix_mpi() ;
double *****d5matrix_mpi() ;
double maxV_of_Rvec() ;
int stop_when() ;
void report_error() ;
void get_time() ;
void reset_time() ;
void free_sword() ;
void free_svector() ;
void free_dvector() ;
void free_d4matrix_mpi() ;
void free_dmatrix_mpi() ;
//void nrerror() ;

void Rmtx_mul() ;
void rRmtx_mul() ;
void Rmtx_mul_ABt() ;
void Rmtx_mul_AtB() ;
void Rmtx_3mul() ;
void Rmtx_3mul_rABCt() ;
void Rmtx_3mul_AtBC() ;
void Rmtx_mul_dgemm() ;
void Rmtx_add() ;
void Rmtx_sub() ;
void r_times_Rmtx() ;
void r_times_Rmtx2() ;

void tdRmtx_mul() ;
void set_null_Rvec() ;
void set_null_Rmtx() ;
void set_null_R3mtx() ;
void Transpose_Rmtx() ;
void copy_Rvec() ;
void copy_Rmtx() ;
void copy_R3mtx() ;
void print_Rvec() ;
void print_Rmtx() ;
void print_Rmtx_condensed() ;
void Identity_Rmtx() ;
int ConvertBooleanCharToInt() ;
int ConvertBooleanToInt() ;

#define MaxChar 100

#define tiny_03 1.0e-03
#define tiny_04 1.0e-04
#define tiny_05 1.0e-05
#define tiny_06 1.0e-06
#define tiny_07 1.0e-07
#define tiny_08 1.0e-08
#define tiny_09 1.0e-09
#define tiny_10 1.0e-10
#define tiny_12 1.0e-12
#define tiny_15 1.0e-15
#define tiny_20 1.0e-20



