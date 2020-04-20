// Ham ran2 o day, nhung ta thu su dung ran4 xem sao -> Ham random o duoi la ran4 day 
 //  Note #undef's at end of file 

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values)
float random2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {               //Initialize
		if (-(*idum) < 1) *idum=1;  // Be sure to prevent idum=0.
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {   //Load the shuffle table (after 8 warm-ups)
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;                      // Start here when not initilizing
	*idum=IA1*(*idum-k*IQ1)-k*IR1;      //Compute idum=(IA1*idum) % IM1 without
	if (*idum < 0) *idum += IM1;        //overflow by Schrage's method
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;      //Compute idum2=IA2*idum %IM2 likewise
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;                          //Will be in the range 0...NTAB-1.
	iy=iv[j]-idum2;                     //Here idum is shuffled, idum and idum2
	iv[j] = *idum;                      //are combine to generate output
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;// Because users don't expect endpoint values
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


// Ket thuc ham ran2. Ham ran4 below

/*

#define NITER 4

void psdes(unsigned long *lword, unsigned long *irword)
{
 unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
 static unsigned long c1[NITER]={
  0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
 static unsigned long c2[NITER]={
  0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

 for (i=0;i<NITER;i++) {
  ia=(iswap=(*irword)) ^ c1[i];
  itmpl = ia & 0xffff;
  itmph = ia >> 16;
  ib=itmpl*itmpl+ ~(itmph*itmph);
  *irword=(*lword) ^ (((ia = (ib >> 16) |
   ((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
  *lword=iswap;
 }
}
#undef NITER

// (C) Copr. 1986-92 Numerical Recipes Software z1+91. 
// Ham ran4 - ham tot nhat theo y tac gia nhung no chay hoi cham
float random(long *idum)
{
 //void psdes(unsigned long *lword, unsigned long *irword);
 unsigned long irword,itemp,lword;
 static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
 static unsigned long jflone = 0x00004080;
 static unsigned long jflmsk = 0xffff007f;
#else
 static unsigned long jflone = 0x3f800000;
 static unsigned long jflmsk = 0x007fffff;
#endif

 if (*idum < 0) {
  idums = -(*idum);
  *idum=1;
 }
 irword=(*idum);
 lword=idums;
 psdes(&lword,&irword);
 itemp=jflone | (jflmsk & irword);
 ++(*idum);
 return (*(float *)&itemp)-1.0;
}
// (C) Copr. 1986-92 Numerical Recipes Software z1+91. 

// */


