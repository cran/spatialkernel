#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <R.h>
#include "risk.h"

#ifdef LINUX
double Pi = M_PI; /*3.14159265358979323846 in i386-glibc21-linux*/
#else
double const Pi = 3.141592653589793238462643;
#endif    

extern void F77_NAME(kfun)(double*, double*, double*, double*, 
    double*, int*, double*);

/*standardized kernel functions*/
double kernelf(double *xpt, int *kernel)
{
    double x0=0, y0=0, bh=1, z=0;
    F77_CALL(kfun)(xpt, xpt+1, &x0, &y0, &bh, kernel, &z);
    return(z);
}

/*y[]==j elements of y should  be coincident to j.*/
double hatpij(int *i, int *j, double *xpts, int *y, int *n, double *h, int *kernel, 
	      double *c)
{
  double sum1=0, sum2=0, xpti[2], tmp;
  int k;
  for(k=0; k<*n; k++){
    if(k==*i) continue;
    xpti[0] = (xpts[*i]-xpts[k])/(*h); 
    xpti[1] = (xpts[*n+*i]-xpts[*n+k])/(*h);
    tmp = kernelf(xpti, kernel)/(h[0]*h[0]*(*(c+k)));
    if(y[k]==*j) sum1 += tmp;
    sum2 += tmp;
  }
  return (sum1/sum2);
}

void lc(double *xpts, int *y, int *n, double *h, int *kernel, double *c, 
	double *lc)
{
  int i;
  *lc = 0;
  for(i=0; i<*n; i++) *lc += log(hatpij(&i, &y[i], xpts, y, n, h, kernel, c));
}

/* c dependent to h c[]: (*n)*(*nh) matrix */
void lcn(double *xpts, int *y, int *n, double *h, int *nh, int *kernel, double *c, 
	double *lc)
{
  int i, j;
  for(j=0; j<*nh; j++) {
    lc[j] = 0;
    for(i=0; i<*n; i++) lc[j] += log(hatpij(&i, &y[i], xpts, y, n, &h[j], kernel, 
					   &c[*n*j]));
  }
}


/*IMPORTANT!!! make sure that {y[i] | i=1,2,...,n-1}=={0,1,2,...,n-1} */
void hatp(double *dpt, double *xpts, int *y, int *n, double *h, int *kernel, 
	   double *c, int *m, double *hatp)
{
  double *sum1, sum2=0, xpti[2], tmp;
  int i;
  sum1 = (double *)malloc((*m)*sizeof(double));
  if(sum1==NULL) {
    Rprintf("\nCannot allocate enough memory.\n");
    return; /*exit(1); FIXED*/
  } 
  for(i=0; i<*m; i++) sum1[i]=0;
  for(i=0; i<*n; i++){
    xpti[0] = (dpt[0] - xpts[i])/(*h); 
    xpti[1] = (dpt[1] - xpts[*n+i])/(*h);
    tmp = (kernelf(xpti, kernel)/(h[0]*h[0])); /*(*(c+i));*/
    sum1[y[i]] += tmp;
    sum2 += tmp;
  }
  for(i=0; i<m[0]; i++) hatp[i] = sum1[i]/sum2;
  free(sum1);
}

/* calculate hatp at many points once; hatp[npts*m], (xpts, y) is the data */
void hatpn(double *dpts, int *ndpts, double *xpts, int *y, int *n, double *h, 
	   int *kernel, double *c, int *m, double *hatp)
{
  double *sum1, sum2, xpti[2], tmp;
  int i, j;
  sum1 = (double *)malloc((*m)*sizeof(double));
  if(sum1==NULL) {
    Rprintf("\nCannot allocate enough memory.\n");
    return; /*exit(1); FIXED*/
  } 

  for(j=0; j<*ndpts; j++) {
    for(i=0; i<*m; i++) sum1[i]=0;
    sum2 = 0;
    for(i=0; i<*n; i++){
      xpti[0] = (dpts[j] - xpts[i])/(*h); 
      xpti[1] = (dpts[*ndpts+j] - xpts[*n+i])/(*h);
      tmp = (kernelf(xpti, kernel)/(h[0]*h[0])); /*/(*(c+i));*/
      sum1[y[i]] += tmp;
      sum2 += tmp;
    }
    for(i=0; i<*m; i++) hatp[j+i*(*ndpts)] = sum1[i]/sum2;
  }
  free(sum1);
}

/* known phat at xpts, wrsp[n]--workspace */
void varphat(double *dpts, int *ndpts, double *xpts, int *y, double *phat, int *n, double *h, 
	   int *kernel, double *c, int *m, double *wrsp, double *pvar)
{
  double sum, xpti[2];
  int i, j, k;
  for(j=0; j<*ndpts; j++) {
    sum = 0;
    for(i=0; i<*n; i++) {
      xpti[0] = (dpts[j] - xpts[i])/(*h); 
      xpti[1] = (dpts[*ndpts+j] - xpts[*n+i])/(*h);
      wrsp[i] = (kernelf(xpti, kernel)/(h[0]*h[0])); /*/(*(c+i));*/
      sum += wrsp[i];
    }
    for(i=0; i<*n; i++) wrsp[i] = wrsp[i]*wrsp[i]/(sum*sum);    
    for(k=0; k<*m; k++) {
      for(i=0; i<*n; i++) {
        pvar[j+k*(*ndpts)] += wrsp[i]*phat[i+k*(*n)]*
          (1-phat[i+k*(*n)]);
      }
    }
  }
}

void hat_lambda_c(double *pts, int *npts, double *dpts, int *ndpts, 
		  double *h, int *kernel, double *c, double *lmb)
{
  int i , j;
  double x[2];
  for(i=0; i<*npts; i++){
    lmb[i]=0;
    for(j=0; j<*ndpts; j++){
      x[0] = (pts[i] - dpts[j])/(*h);
      x[1] = (pts[i+*npts] - dpts[j+*ndpts])/(*h);	
      /*lmb[i] += kernelf(x, kernel)/(h[0]*h[0]*(*(c+j))); */
      lmb[i] += kernelf(x, kernel);
    }
    lmb[i] = lmb[i]/(h[0]*h[0]*c[i]);
  }
  return;
} 

/*Baddeley leave-one-out version of lambdathat specifically for Kinhat*/
void hat_lambda_b(double *pts, int *npts, double *h, int *kernel, double *c, double *lmb)
{
  int i , j;
  double x[2];
  for(i=0; i<*npts; i++){
    lmb[i]=0;
    for(j=0; j<*npts; j++){
      if(j==i) continue;
      x[0] = (pts[i] - pts[j])/(*h);
      x[1] = (pts[i+*npts] - pts[j+*npts])/(*h);	
      lmb[i] += kernelf(x, kernel);
    }
    lmb[i] = lmb[i]/(h[0]*h[0]*c[i]);
  }
  return;
} 

void area_poly(double *pts, int *npts, double *area)
{
  double *x, *y;
  int i;
  *area=0;
  x=pts; 
  y=pts+*npts;
   *area += x[0]*(y[1]-y[*npts-1]);
  for(i=1; i<*npts-1; i++) *area += x[i]*(y[i+1]-y[i-1]);
  *area += x[*npts-1]*(y[0]-y[*npts-2]);
  *area /= 2.0;   /* 2*A(P) */
  return;
}

/* intensity estimation for source point process */
void lambda_source_pp(double *pts, int *npts, double *c, double *spts, int *nspts,
		      double *h, int *kernel, double *lmb)
{
  int i, j;
  double x[2];
  for(i=0; i<*npts; i++) {
    lmb[i] = 0;
    for(j=0; j<*nspts; j++) {
      x[0] = (pts[i] - spts[j])/(*h); 
      x[1] = (pts[*npts+i] - spts[*nspts+j])/(*h);
      lmb[i] += c[j]*kernelf(x, kernel);
    }
    lmb[i] = lmb[i]/(h[0]*h[0]);
  }
  return;
}
