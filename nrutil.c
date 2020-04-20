#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

#define FSLIOERR(x) { fprintf(stderr,"Error:: %s\n",(x)); fflush(stderr); exit(EXIT_FAILURE); }//8-Dec.2009 for d4matrix
   
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr," run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ***d3matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch)
     
{
  int i ;
  double ***m ;

  m = (double ***) malloc((unsigned)(nh-nl+1)*sizeof(double**)) ;
  if (!m) 
    {
      printf("allocation failure 1 in d3matrix()");
      exit(1);
    }
  m -= nl ;

  for ( i=nl ; i<=nh ; i++ )
    m[i] = dmatrix(nrl,nrh,ncl,nch) ;

  return m ;
}

//30/03/10 16:41 Tuan made i3matrix
int ***i3matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch)
{
  int i ;
  int ***m ;

  m = (int ***) malloc((unsigned)(nh-nl+1)*sizeof(int**)) ;
  if (!m) 
    {
      printf("allocation failure 1 in i3matrix()");
      exit(1);
    }
  m -= nl ;

  for ( i=nl ; i<=nh ; i++ )
    m[i] = imatrix(nrl,nrh,ncl,nch) ;

  return m ;
}
/*n4matrix allocates space for a four dimensional matrix*/
double ****n4matrix(int nl,int nr,int ml,int mr, int kl,int kr,int jl,int jr)
{
	int i,j,k;
	double ****m;
	
	m=(double ****) calloc((unsigned) (nr-nl+1),(unsigned) sizeof(double***));
	if (!m)
		nrerror("allocation failure 1 in n4matrix");
	m -= nl;
	for (i=nl; i<=nr; i++)
		{
		m[i]=(double ***) calloc((unsigned) (mr-ml+1),(unsigned) sizeof(double**));
		if (!m[i])
			nrerror("allocation failure 2 in n4matrix");
		m[i] -= ml;
		for (j=ml; j<=mr;j++)
			{
			m[i][j]=(double **) calloc((unsigned) (kr-kl+1),(unsigned) sizeof(double*));
			if (!m[i][j])
				nrerror("allocation failure 3 in n4matrix");
			m[i][j] -= kl;
			for(k=jl;k<=jr;k++)
				{
				m[i][j][k]=(double *) calloc((unsigned) (jr-jl+1),(unsigned) sizeof(double));
				if (!m[i][j])
					nrerror("allocation failure 4 in n4matrix");
				m[i][j][k] -= jl;
				}
			}
		}
	return m;
}

/*free_n4matrix free allocated space for n4matrix*/
void free_n4matrix(double ****m, int nl,int nr,int ml,int mr,int kl,int kr,int jl,int jr)
{
	int i,j,k;
	
	for (i=nl;i>=nr;i--)
		{
		for (j=ml;j>=mr;j--)
			{
			for (k=jl;k>=jr;k--)
				free((char *) (m[i][j][k]+jl));
			free((char *) (m[i][j]+kl));
			}
		free((char*) (m[i]+ml));
		 }
	free((char *) (m+nl));
}



void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d3matrix(double ***m ,int nl,int nh,int nrl,int nrh,int ncl,int nch)  
{
  int i ;
  for ( i=nh ; i>=nl ; i-- )
    free_dmatrix(m[i],nrl,nrh,ncl,nch) ;
  free((char*)(m+nl)) ;
}

//30/03/10 16:44 Tuan made free_i3matrix
void free_i3matrix(int ***m ,int nl,int nh,int nrl,int nrh,int ncl,int nch)  
{
  int i ;
  for ( i=nh ; i>=nl ; i-- )
    free_imatrix(m[i],nrl,nrh,ncl,nch) ;
  free((char*)(m+nl)) ;
}

/* Cac ham minh tu lam (Feb 23 and 24, 2010 // Tuan made */
///* Khong dung cach nay 
    //d4matrix allocates space for a four dimensional matrix
double ****d4matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh)
{
	int i;
	double ****m;
	
	m=(double ****) malloc((unsigned) (nh-nl+1)* sizeof(double***));
	if (!m){
		printf("allocation failure 1 in d4matrix()");
        exit(1);
      }
	m -= nl;
	for (i=nl; i<=nh; i++)
	    {
		  m[i] = d3matrix(nrl,nrh,ncl,nch,nkl,nkh);
        }
	return m;
}

/*free_d4matrix free allocated space for ndmatrix*/

void free_d4matrix(double ****m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh)
{
	int i;
	for ( i=nh ; i>=nl ; i-- ){
       free_d3matrix(m[i],nrl,nrh,ncl,nch,nkl,nkh);
      }
  free((char*)(m+nl)) ;
	
}

double *****d5matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh)
{
	int i;
	double *****m;
	
	m=(double *****) malloc((unsigned) (nh-nl+1)* sizeof(double****));
	if (!m){
		printf("allocation failure 1 in d5matrix()");
        exit(1);
      }
	m -= nl;
	for (i=nl; i<=nh; i++)
	    {
		  m[i] = d4matrix(nrl,nrh,ncl,nch,nkl,nkh,nml,nmh);
        }
	return m;
}

void free_d5matrix(double *****m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh)
{
	int i;
	for ( i=nh ; i>=nl ; i-- ){
       free_d4matrix(m[i],nrl,nrh,ncl,nch,nkl,nkh,nml,nmh);
      }
  free((char*)(m+nl)) ;
	
}
// Tuan included: March 08, 2010 
double ******d6matrix(int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh,int nql, int nqh)
{
	int i;
	double ******m;
	
	m=(double ******) malloc((unsigned) (nh-nl+1)* sizeof(double*****));
	if (!m){
		printf("allocation failure 1 in d6matrix()");
        exit(1);
      }
	m -= nl;
	for (i=nl; i<=nh; i++)
	    {
		  m[i] = d5matrix(nrl,nrh,ncl,nch,nkl,nkh,nml,nmh,nql,nqh);
        }
	return m;
}

void free_d6matrix(double ******m,int nl,int nh,int nrl,int nrh,int ncl,int nch,int nkl,int nkh,int nml,int nmh,int nql, int nqh)
{
	int i;
	for ( i=nh ; i>=nl ; i-- ){
       free_d5matrix(m[i],nrl,nrh,ncl,nch,nkl,nkh,nml,nmh,nql,nqh);
      }
  free((char*)(m+nl)) ;
	
}

/****************************************************************
 *  FUNCTION:  d5matrix -- dynamically allocate space in memory *
 *                         for a 5-dimensional matrix           *
 *  INPUTS: dem1 = the first dimension
            dem2 = the second dimension                            
            dem3 = the third dimension
            dem4 = the fourth dimension
            dem5 = the fifth dimension                          
 *  OUTPUT: A = dem1 x dem2 x dem3 x dem4 x dem5 matrix
 
 Made by Tuan - 11-Jun-2008 
 *****************************************************************/
/*
double *****d5matrix(int dem1,int dem2,int dem3,int dem4,int dem5)
{
int i, j, k,l;
double *****A;

A = (double *****) calloc(dem1, sizeof(double ****));
for (i =1; i <=dem1; i++)
    {
    A[i] = (double ****) calloc(dem2, sizeof(double ***));
    for (j = 1; j <= dem2; j++){
        A[i][j] = (double ***) calloc(dem3, sizeof(double **));
        for(k=1;k<=dem3;k++){
            A[i][j][k] = (double **) calloc(dem4, sizeof(double *));  
            for(l=1;l<=dem4;l++){
                A[i][j][k][l] = (double *) calloc(dem5, sizeof(double ));                                        
            }// End of for(l=1
         }// End of for(k=1
       } // End of for(j=1
     } // End of for(i=1  
return(A);
}  ///* end of function d5matrix 
*/

//*/
// Dung cach nay (8 December 2009)
/***************************************************************
 * d4matrix
 ***************************************************************/
/*! \fn double ****d4matrix(int th, int zh,  int yh, int xh)
    \brief allocate a 4D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..th][0..zh][0..yh][0..xh].  
        <br>To access all elements as a vector, use buf[0][0][0][i] where
        i can range from 0 to th*zh*yh*xh - 1.
        Adaptation of Numerical Recipes in C nrutil.c allocation routines. 

    \param th slowest changing dimension
    \param zh 2nd slowest changing dimension
    \param yh 2nd fastest changing dimension
    \param xh fastest changing dimension
    \return Pointer to 4D double array
 */
 /*
double ****d4matrix(int th, int zh,  int yh, int xh)
{

        int j;
        int nvol = th+1;
        int nslice = zh+1;
        int nrow = yh+1;
        int ncol = xh+1;
        double ****t;


        ///** allocate pointers to vols 
        t=(double ****) malloc((size_t)((nvol)*sizeof(double***)));
        if (!t) FSLIOERR("d4matrix: allocation failure");

        //** allocate pointers to slices 
        t[0]=(double ***) malloc((size_t)((nvol*nslice)*sizeof(double**)));
        if (!t[0]) FSLIOERR("d4matrix: allocation failure");

       // /** allocate pointers for ydim 
        t[0][0]=(double **) malloc((size_t)((nvol*nslice*nrow)*sizeof(double*)));
        if (!t[0][0]) FSLIOERR("d4matrix: allocation failure");


        ///** allocate the data blob
        t[0][0][0]=(double *) malloc((size_t)((nvol*nslice*nrow*ncol)*sizeof(double)));
        if (!t[0][0][0]) FSLIOERR("d4matrix: allocation failure");


        ///** point everything to the data blob
        for(j=1;j<nrow*nslice*nvol;j++) t[0][0][j]=t[0][0][j-1]+ncol;
        for(j=1;j<nslice*nvol;j++) t[0][j]=t[0][j-1]+nrow;
        for(j=1;j<nvol;j++) t[j]=t[j-1]+nslice;

        return t;
}
*/


#else /* ANSI */
/* traditional - K&R */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
long nch,ncl,nrh,nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m)	nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((unsigned int)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
long nch,ncl,nrh,nrl;
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
float ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}



#endif /* ANSI */

