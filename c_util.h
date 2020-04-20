typedef struct {
   double r ;
   double i ;
} doublecomplex ;

doublecomplex *complex_dvector() ;
doublecomplex **complex_dmatrix() ;
doublecomplex ***complex_d3matrix() ;
doublecomplex ***complex_d3matrix_mpi() ;

void free_complex_dvector() ;
void free_complex_dmatrix() ;
void free_complex_d3matrix() ;

doublecomplex Complex() ;
doublecomplex Czero_() ;
doublecomplex Cneg() ;
doublecomplex Conj() ;
doublecomplex Cadd() ;
doublecomplex Csub() ;
doublecomplex Cmul() ;
doublecomplex Cdiv() ;
doublecomplex RCmul() ;
doublecomplex CRmul() ;
doublecomplex CRdiv() ;
doublecomplex RCadd() ;
doublecomplex CRadd() ;
doublecomplex Csqrt() ;

double Re() ;
double Cabs() ;
double Cabs2() ;

void Cmtx_mul() ;
void CmtxCj_mul() ;
void Cmtx_3mul() ;
void CmtxCj_3mul() ;
void RCmtx_mul() ;
void CRmtx_mul() ;
void RtimesCmtx() ;
void set_null_Cmtx() ;
void Cmtx_add() ;
void RCmtx_add() ;
void Cmtx_sub() ;
void Adjoint_mtx() ;
void copy_Cmtx() ;
void print_Cmtx() ;
void print_Cmtx_condensed() ;
void Identity_Cmtx() ;
void set_null_Cvec() ;

#define ComplexZero Complex(0.0,0.0)
