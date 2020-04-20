#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

double ***d3matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch);
int ***i3matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch);
double ****d4matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh);// Made by Tuan
//double ****d4matrix(int th, int zh,  int yh, int xh); // Ta dung cai nay (8 December 2009)
   // d4matrix(int th, int zh,  int yh, int xh) thuc chat d4matrix[0,th, 0, zh,0,yh,0, xh)
//double *****d5matrix(int dem1,int dem2,int dem3,int dem4,int dem5);// Made by Tuan
double *****d5matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh);
double ******d6matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh,int nql, int nqh);
double ****n4matrix(int nl,int nr,int ml,int mr, int kl,int kr,int jl,int jr);

void free_n4matrix(double ****m, int nl,int nr,int ml,int mr,int kl,int kr,int jl,int jr);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
void free_d3matrix( double ***m, int nl,int nh,int nrl,int nrh,int ncl,int nch) ;
void free_i3matrix(int ***m ,int nl,int nh,int nrl,int nrh,int ncl,int nch);
void free_d4matrix(double ****m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh);
void free_d5matrix(double *****m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh);
void free_d6matrix(double ******m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh,int nql, int nqh);
void free_n4matrix(double ****m, int nl,int nr,int ml,int mr,int kl,int kr,int jl,int jr);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */

