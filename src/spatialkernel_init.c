#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void adaptpoly(double *, int *, double *, int *, double *, int *, 
  double *, double *, double *, double *, double *, int *, int *, int *, 
  double *);
extern void area_poly(double *, int *, double *);
extern void hat_lambda_b(double *, int *, double *, int *, double *, double *);
extern void hat_lambda_c(double *, int *, double *, int *, double *, int *, 
  double *, double *);
extern void hatpn(double *, int *, double *, int *, int *, double *, int *, 
  double *, int *, double *);
extern void lcn(double *, int *, int *, double *, int *, int *, double *, 
  double *);
extern void varphat(double *, int *, double *, int *, double *, int *, double *, 
  int *, double *, int *, double *, double *);

/* .Fortran calls */
extern void F77_NAME(dokinhat)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pnpoly)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(psnpoly)(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"adaptpoly",    (DL_FUNC) &adaptpoly,    15},
  {"area_poly",    (DL_FUNC) &area_poly,     3},
  {"hat_lambda_b", (DL_FUNC) &hat_lambda_b,  6},
  {"hat_lambda_c", (DL_FUNC) &hat_lambda_c,  8},
  {"hatpn",        (DL_FUNC) &hatpn,        10},
  {"lcn",          (DL_FUNC) &lcn,           8},
  {"varphat",      (DL_FUNC) &varphat,      12},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"dokinhat", (DL_FUNC) &F77_NAME(dokinhat), 10},
  {"pnpoly",   (DL_FUNC) &F77_NAME(pnpoly),    6},
  {"psnpoly",  (DL_FUNC) &F77_NAME(psnpoly),   7},
  {NULL, NULL, 0}
};

void R_init_spatialkernel(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
