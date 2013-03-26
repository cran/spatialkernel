#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include "risk.h"
#include "gpc.h"
#include"R.h"

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include "triangle.h"

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
#endif /* _STDLIB_H_ */

void realfree(REAL *p)
{
  if(p!=(REAL *)NULL) free(p);
}

void intfree(int *p){
  if(p!=(int *) NULL) free(p);
}

#define abs(x) ((x) < 0 ? -(x) : (x) )
#define min(x, y) ((x)<(y) ? (x) : (y))
#define max(x, y) ((x)<(y) ? (y) : (x))

/*circumcircle through (xtri[i], ytri[i]), i=1,2,3;
return centre (p[0], p[1]) and radius r; 
return r2=r*r for convience
 r2 = -1 means the three points are collinear*/
void circumcircle(double *xtri, double *ytri, double *p, double *r2)
{
  double A,B,C,D,E,F,G;

  A = xtri[1] - xtri[0];
  B = ytri[1] - ytri[0];
  C = xtri[2] - xtri[0];
  D = ytri[2] - ytri[0];

  G = 2.0*(A*(ytri[2] - ytri[1])-B*(xtri[2] - xtri[1]));
  if(G==0) {
    r2[0] = -1;
    return;
  }

  E = A*(xtri[0] + xtri[1]) + B*(ytri[0] + ytri[1]);
  F = C*(xtri[0] + xtri[2]) + D*(ytri[0] + ytri[2]);

  p[0] = (D*E - B*F) / G;
  p[1] = (A*F - C*E) / G;
  r2[0] = (xtri[0] - p[0])*(xtri[0] - p[0]) + (ytri[0] - p[1])*(ytri[0] - p[1]);
}

/*rangexy: [0]-[1] is the x-range; [2]-[3] is the y-range;
passing in value of addindx[0] is the length of xgrid and ygrid;
for checking if xgrid and ygrid is long enough to hold;
returning addindx[] indicate the locations: number of added
between 0 to 1, 1 to 2, ..., npoly-1 to 0; the last is the number added inside;
addindx[npoly+1]; (xgrid, ygird) not include polygon points;
*/
void add_grids(double* poly, int* npoly, double *rangxy, double *h, 
	       double *xgrid, double *ygrid, int* addindx)
{
  double xmin, xmax, ymin, ymax;
  double x1, x2, y1, y2;
  int i, j, inout, ngrids=0, checksum;
  double tmp, tmpx, tmpy, deltax, deltay;
  
  checksum = addindx[0];
  xmin=rangxy[0];
  xmax=rangxy[1]; 
  ymin=rangxy[2]; 
  ymax=rangxy[3];
  
  /*add steiner points on the boundary of polygon;*/
  for(i=0; i<npoly[0]; i++){
    addindx[i] = 0;
    x1 = poly[i]; 
    y1 = poly[npoly[0]+i];
    if(i<npoly[0]-1){
      x2 = poly[i+1]; 
      y2 = poly[npoly[0]+i+1];    
    } else {
      x2 = poly[0]; 
      y2 = poly[npoly[0]];    
    }
    tmp = ceil(max(abs(x2-x1), abs(y2-y1))/h[0]);
    if(tmp > 1){
      deltax = (x2-x1)/tmp;
      deltay = (y2-y1)/tmp;
      for(j=1; j<tmp; j++) {
	xgrid[ngrids] = x1+((double)j)*deltax;
	ygrid[ngrids] = y1+((double)j)*deltay;
	addindx[i]++;
	ngrids++;
      }
    }
  }
  
  /*add points inside */
  addindx[npoly[0]] = 0;
  for(tmpx=xmin; tmpx<xmax; tmpx+=h[0]){
    for(tmpy=ymin; tmpy<ymax; tmpy+=h[0]){
      pnpoly_(&tmpx, &tmpy, poly, &poly[npoly[0]], npoly, &inout);
      if(inout==1){
	xgrid[ngrids] = tmpx;
	ygrid[ngrids] = tmpy;
	addindx[npoly[0]]++;
	ngrids++;
      }
    }
  }
  
  if(ngrids>checksum){
    error("\nMemory for added points in add_grids too small.\n");
    return; /* exit(1) FIXED*/;
  }
  return;
}

/*pts[2*npts], tri[3*ntri];
pts[0], pts[1] are the first point coordinates
tri[0], tri[1], tri[2] are three integer index for the first triangle, ...
*/
void constrained(double* poly, int* npoly, double *rangxy, double *h, 
  double* pts, int* npts, int *tri, int* ntri)
{
  struct triangulateio in, out;
  int i, j, nc1, nc2;
  /*double xgrid[npts[0]-npoly[0]], ygrid[npts[0]-npoly[0]];
  int addindx[npoly[0]+1], nedge, ninside;*/
  int *addindx, nedge, ninside;
  REAL *xgrid, *ygrid;
  addindx=(int *)malloc((npoly[0]+1)*sizeof(int));
  xgrid=(REAL *)malloc((npts[0]-npoly[0])*sizeof(REAL));
  ygrid=(REAL *)malloc((npts[0]-npoly[0])*sizeof(REAL));
  if(addindx==NULL || xgrid==NULL || ygrid==NULL) {
    error("\nCannot allocate enough memory.\n");
    return; /* exit(1); FIXED */
  } 
  addindx[0] = npts[0]-npoly[0];  

  if(npoly[0]>3000){
    error("\nNumber of polygon points exceed 3000.\n");
    return;/* exit(1); FIXED */
  }

  add_grids(poly, npoly, rangxy, h, xgrid, ygrid, addindx);
  /* for(i=0; i<npoly[0]+1; i++) Rprintf("\naddindx[%d]=%d", i, addindx[i]);*/

  nedge=0;
  for(i=0; i<npoly[0]; i++) nedge += addindx[i];
  ninside=addindx[npoly[0]];

  in.numberofpoints = npoly[0]+nedge+ninside;
  in.numberofpointattributes = 0;

  npts[0] = in.numberofpoints;

  in.pointlist = (REAL *) malloc(in.numberofpoints * 2*sizeof(REAL));
  nc1=0; /*index for in.pointlist;*/
  nc2=0; /*index for xgrid and ygrid;*/
  for(i=0; i<npoly[0]; i++) {
    /*polygon point;*/
    in.pointlist[2*nc1]=poly[i];
    in.pointlist[2*nc1+1]=poly[npoly[0]+i];
    nc1++;
    /*added points on segment between point i to i+1;*/
    for(j=0; j<addindx[i]; j++){
      in.pointlist[2*nc1]=xgrid[nc2];
      in.pointlist[2*nc1+1]=ygrid[nc2];
      nc1++;
      nc2++;
    }
  }

  /*added points inside the polygon;*/
  for(i=0; i<ninside; i++){
    in.pointlist[2*nc1]=xgrid[nc2];
    in.pointlist[2*nc1+1]=ygrid[nc2];
    nc1++;
    nc2++;
  }

  in.pointmarkerlist = (int *) NULL;
  
  in.numberofsegments = npoly[0]+nedge;
  in.segmentlist = (int *) malloc(2*in.numberofsegments*sizeof(int));
  in.segmentlist[0]=in.numberofsegments-1;
  in.segmentlist[1]=0;
  for(i=1; i<in.numberofsegments; i++){
    in.segmentlist[2*i] = i-1;
    in.segmentlist[2*i+1] = i;
  }
  in.segmentmarkerlist = (int *) NULL;
  
  in.numberofholes = 0;  
  in.numberofregions = 0;
  
  out.pointlist = (REAL *) NULL; 
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL; 
  
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;

  triangulate("pzBNQ", &in, &out, (struct triangulateio *) NULL);
 
  if(ntri[0]<out.numberoftriangles) {
    error("\nMemory for generated triangles too small.\nGenerated %d triangles.\n", out.numberoftriangles);
    return;/* exit(1); FIXED */
  } else ntri[0]=out.numberoftriangles;

  for(i=0; i<2*npts[0]; i++) pts[i] = in.pointlist[i];

  for(i=0; i<3*ntri[0]; i++) tri[i]=out.trianglelist[i];

  realfree(in.pointlist);
  intfree(in.segmentlist);
  intfree(in.pointmarkerlist);
  intfree(in.segmentmarkerlist);
  realfree(out.pointlist);
  intfree(out.trianglelist);
  
  intfree(out.pointmarkerlist);
  intfree(out.segmentlist);
  intfree(out.segmentmarkerlist);
  return;
}

/* n vertices polygon to triangles (index) by connecting
 vertex[0] to the others */
void adaptconvex(double *xpoly, double *ypoly, int *npoly, 
		 double *x, double *y, double *h, 
		 int *kernel, double *eps, double *err, int *mcalls, 
		 int *ncalls, int *ier, int *nw, double *w,
		 double *ans)
{
  int i;
  int ier0, ncalls0;
  double tri[6], eps0, err0, res; /*cubtri_ ask tri in special order*/
  tri[0] = xpoly[0]; tri[1] = ypoly[0]; /* vertex[0] fixed */
  ier[0]=ier[1]=ier[2]=ier[3]=ier[4]=ier[5]=0;
  eps0 = *eps/(*npoly-2);
  *err = 0;
  *ncalls = 0;
  *ans = 0;
  for(i=1; i<*npoly-1; i++) {
    tri[2] = xpoly[i]; tri[3] = ypoly[i]; /*vertex[1] */
    tri[4] = xpoly[i+1]; tri[5] = ypoly[i+1]; /*vertex[2]*/
    ncalls0 = 0;
    cubtri_(tri, &eps0, mcalls, &res, &err0, &ncalls0, w, nw,
	    x, y, h, &ier0, kernel);
    ier[ier0]++;
    *err += err0;
    if(*ncalls < ncalls0) *ncalls = ncalls0;
    *ans += res;
  }
}

/*calculate c(x_i, h), for all x_i, i=1,...,n and one h together
xpts are coordinates of x_i; poly and xpts: first x coordinates, then y...
cxh[n] as long as xpts[n*2];
c1*h to grid the polygon area; c2*h away from xpt integrand is 0;*/
void adaptpoly(double *poly, int *npoly, double *xpts, int *nxpts, 
	       double *h, int *kernel, double *c1, double *c2, double *rngxy,
	       double *eps, double *err, int *mcalls, int *ncalls, int *ier,
	       double *cxh)
{
  double *tripts, c1h, c2h, neps, err0;
  int ntripts, ntri, *tri, ncalls0; /* , ier0;*/
  int nxgrid, nygrid;
  double tri0[6], tri1[6], *w, ans; /*, p[2], r2;*/
  int nw;
  double sq[8], clip[22];
  int three=3, four=4, nclip, ier1[6];
  int i, j, k;
  c1h = c1[0]*h[0];
  c2h = c2[0]*h[0];
  nxgrid = (int) ceil((rngxy[1]-rngxy[0])/(c1h));
  nygrid = (int) ceil((rngxy[3]-rngxy[2])/(c1h));
  ntripts = 2*(nxgrid*nygrid+npoly[0]);
  tripts = (double *) malloc(2*ntripts*sizeof(double));
  ntri = 2*ntripts;
  tri = (int *) malloc(3*ntri*sizeof(int));
  constrained(poly, npoly, rngxy, &c1h, tripts, &ntripts, tri, &ntri);
  neps = eps[0]/((double) ntri);
  nw = (int) ceil(3.0*(19.0+3.0*mcalls[0])/38.0);
  w = (double *) malloc(nw*sizeof(double));
  for(i=0; i<nxpts[0]; i++) {
    ncalls[i] = 0;
    cxh[i] = 0;
    err[i] = 0;
  }
  for(j=0; j<6; j++) ier[j] = 0;
  for(j=0; j<ntri; j++){ /*each triangle*/
    /* vertice index tri[3*j], tri[3*j+1], tri[3*j+2]
     point-xy: tripts[2*tri[3*j]], tripts[2*tri[3*j]+1], ... */
    for(k=0; k<3; k++) {
      tri0[2*k] = tripts[2*tri[3*j+k]]; /*x*/
      tri0[2*k+1] = tripts[2*tri[3*j+k]+1]; /*y*/
      tri1[k] = tri0[2*k]; 
      tri1[3+k] = tri0[2*k+1]; 
    }
    /* circumcircle(tri1, &tri1[3], p, &r2); 
    p circumcircle centre, r2=r*r, r radius
    d=distance of circumcircle centre to integration centre point
    if d > c2h+sqrt(r2) then integration=0
    r2=-1 means tri1 collinear*/
    for(i=0; i<nxpts[0]; i++) { 
      /*if( (r2!=-1)&&( (p[0]-xpts[i])*(p[0]-xpts[i])+(p[1]-xpts[nxpts[0]+i])
      	      *(p[1]-xpts[nxpts[0]+i]) < 
      	      (c2h+sqrt(r2))*(c2h+sqrt(r2)) ) ) {
      	if(*kernel!=1) { // non-Gaussian */
      sq[0] = xpts[i] - c2h; sq[1] = xpts[i] + c2h; 
      sq[2] = xpts[i] + c2h; sq[3] = xpts[i] - c2h; 
      sq[4] = xpts[nxpts[0]+i] - c2h; sq[5] = xpts[nxpts[0]+i] - c2h; 
      sq[6] = xpts[nxpts[0]+i] + c2h; sq[7] = xpts[nxpts[0]+i] + c2h; 
      simple_polygon_clip(tri1, &tri1[3], &three, sq, &sq[4], &four, 
			  clip, &clip[11], &nclip);
      if(nclip > 0){
	adaptconvex(clip, &clip[11], &nclip, &xpts[i], &xpts[nxpts[0]+i],
		    h, kernel, eps, &err0, mcalls, &ncalls0, ier1, &nw, w, 
		    &ans);
	ier[0] += ier1[0]; ier[1] += ier1[1]; ier[2] += ier1[2]; 
	ier[3] += ier1[3]; ier[4] += ier1[4]; ier[5] += ier1[5];
	if(ncalls[i]<ncalls0) ncalls[i]=ncalls0;
	cxh[i] += ans;
      }
    }
  }
  free(tripts); 
  free(tri);
  free(w);
}
 
/* int main()
{
  REAL poly[]={0,1,1,0,0,0,1,1}, rngxy[]={0,1,0,1};
  REAL xpts[]={0.25, 0.5, 0.5, 0.25}, cxh[]={0, 0};
  double h=20, c1=10.1, c2=20, eps=1.0e-6, err[2];
  int npoly=4, nxpts=2, kernel=1; 
  int mcalls=1000, ncalls[2], ier[6];
  int i;
  for(i=0; i<1000; i++){
  adaptpoly(poly, &npoly, xpts, &nxpts, &h, &kernel, &c1, &c2, rngxy, 
	     &eps, err, &mcalls, ncalls, ier, cxh);
  }
  Rprintf("\nPoint      Cubature        Ncalls        Est abs error\n");
  for(i=0; i<nxpts; i++) {
    Rprintf("%5d%15.8g%10d%22.6g\n", i+1, cxh[i], ncalls[i], err[i]);
  }
  Rprintf("\nFrequency of integration indicators:\n");
  for(i=0; i<6; i++) Rprintf("%6d", i);
  Rprintf("\n");
  for(i=0; i<6; i++) Rprintf("%6d", ier[i]);
  Rprintf("\n\nPI/2=%11.10g\n", asin(1.0));
  return 0;
}
*/
