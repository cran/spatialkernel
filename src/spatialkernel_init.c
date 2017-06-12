#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* References for what this means
https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols
https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Registering-native-routines
*/

/* .C calls */
extern void adapt_poly(double *, int *, double *, int *, double *, int *, 
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
extern void var_phat(double *, int *, double *, int *, double *, int *, double *, 
  int *, double *, int *, double *, double *);

/* .Fortran calls */
extern void F77_NAME(dokinhat)(double *, double *, int *, double *, double *, double *, int *, double *, int *, double *);
extern void F77_NAME(pnpoly)(double *, double *, double *, double *, int *, int *);
extern void F77_NAME(psnpoly)(double *, double *, int *, double *, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
  {"adapt_poly",    (DL_FUNC) &adapt_poly,    15},
  {"area_poly",    (DL_FUNC) &area_poly,     3},
  {"hat_lambda_b", (DL_FUNC) &hat_lambda_b,  6},
  {"hat_lambda_c", (DL_FUNC) &hat_lambda_c,  8},
  {"hatpn",        (DL_FUNC) &hatpn,        10},
  {"lcn",          (DL_FUNC) &lcn,           8},
  {"var_phat",      (DL_FUNC) &var_phat,      12},
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
