#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "my_util.h"

int imax(int a, int b)
{
  if ( a > b ) return(a) ;
  else return(b) ;
}

int imin(int a, int b)
{
  if ( a < b ) return(a) ;
  else return(b) ;
}

double dmax(double a, double b)
{
  if ( a > b ) return(a) ;
  else return(b) ;
}

double dmin(double a, double b)
{
  if ( a < b ) return(a) ;
  else return(b) ;
}

char *sword(len)
     int len ;
/* Allocate an char vector with range [nl..nh].*/
{
        char *v;
 
        v=(char *) malloc((unsigned) len*sizeof(char));
        if (!v) nrerror("allocation failure in svector()");
        return v ;
}

char **svector(nrl,nrh,len)
     int nrl,nrh,len ;
/*Allocate a char matrix with range [nrl..nrh][ncl..nch].*/
{
        int i;
        char **m;
                                                                                
        /*allocate pointers to rows.*/
        m=(char **) malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
        if (!m) nrerror("allocation failure 1 in smatrix()");
        m -= nrl;
                                                                                
        /*allocate rows and set pointers to them.*/
        for (i=nrl;i<=nrh;i++) {
                m[i]=(char *) malloc((unsigned) len*sizeof(char));
                if (!m[i]) nrerror("allocation failure 2 in smatrix()");
        }
                                                                                
        /*return pointer to array of pointers to rows.*/
        return m;
}

double **dmatrix_mpi(m0,m1,n0,n1)
     int m0,m1,n0,n1 ;
{
  int i,M,N ;
  double *v1,**v2 ;

  M = m1-m0+1 ;
  N = n1-n0+1 ;

  v1 = (double *) malloc((unsigned)(M*N)*sizeof(double)) ;
  v1 -= n0 ;

  v2 = (double **) malloc((unsigned)M*sizeof(double*)) ;
  v2 -= m0 ;

  for ( i=m0 ; i<M+m0 ; i++ )
    v2[i] = v1+(i-m0)*N ;

  return(v2) ;
}

double ***d3matrix_mpi(l0,l1,m0,m1,n0,n1) 
     int l0,l1,m0,m1,n0,n1 ;
{
  int i,L,M,N ;
  double *v1,**v2,***v3 ;

  L = (l1-l0+1) ;
  M = (m1-m0+1) ;
  N = (n1-n0+1) ;
  
  v1 = (double *) malloc((unsigned)(L*M*N)*sizeof(double)) ;
  v1 -= n0 ;

  v2 = (double **) malloc((unsigned)(L*M)*sizeof(double *)) ;
  v2 -= m0 ;

  for ( i=m0 ; i<L*M+m0 ; i++ )
    v2[i] = v1+(i-m0)*N ;
  
  v3 = (double ***) malloc((unsigned)L*sizeof(double **)) ;
  v3 -= l0 ;
  for ( i=l0 ; i<L+l0 ; i++ )
    v3[i] = v2+(i-l0)*M ;

  return(v3) ;
}

void free_sword(v)
     char *v ;
/*Frees a int vector allocated by ivector().*/
{
        free((char*) v);
}

void free_svector(m,nrl,nrh)
     char **m;
     int nrl,nrh ;
/*Frees a matrix allocated with svector.*/
{
        int i;
                                                                                
        for (i=nrh;i>=nrl;i--) free((char*) m[i]);
        free((char*) (m+nrl));
}


void free_dmatrix_mpi(double **v, int m0, int m1, int n0, int n1)
{
  free((char*)(v[m0]+n0)) ;
  free((char*)(v+m0)) ;
}

void copy_R3mtx(b,a,nx0,nx1,ny0,ny1,nz0,nz1)
     double ***a,***b ;
     int nx0,nx1,ny0,ny1,nz0,nz1 ;
{
  int i,j,k ;

  for ( i=nx0 ; i<=nx1 ; i++ )
    for ( j=ny0 ; j<=ny1 ; j++ )
      for ( k=nz0 ; k<=nz1 ; k++ )
        b[i][j][k] = a[i][j][k] ;
}

void report_error(error_text)
char error_text[];
{
        void exit();
 
        printf("%s\n",error_text);
        printf("...now exiting to system...\n");
        exit(1);
}

int nint(x)
double x ;
{
  int nx ;

  nx = (int) x ;
  if ( x > 0.0 ) {
    if ( (x - nx) < 0.5 ) return(nx) ;
    else return(nx+1) ;
  }
  else {
    if ( (nx - x) < 0.5 ) return(nx) ;
    else return(nx-1) ;
  }

}

int ConvertBooleanCharToInt(char c)
{
  if ( c=='y' ) return(1) ;
  else if ( c=='n' ) return(0) ;
  else report_error("Should be y or n") ;

  return 0 ;
}

int ConvertBooleanToInt(char *data)
{
  if ( !strcmp(data,"yes") ) return(1) ;
  else return(0) ;
}

static int c1,c2 ;

void get_time(s)
     char *s ;
{
  clock_t clock() ;
  static int entry=0 ;

  c2 = clock() ;
  printf("%s = %lf\n",s,(c2-c1)/(double)CLOCKS_PER_SEC) ;
  fflush(stdout) ;

  c1 = c2 ;
}

void reset_time()
{
  c1 = clock() ;
}
