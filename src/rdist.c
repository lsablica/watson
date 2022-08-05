
#include <R.h>
#include <Rdefines.h>
#include "rdist.h"

/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

#ifdef __GNUC__
#  define ATTRIBUTE__UNUSED        __attribute__ ((unused))
#else
#  define ATTRIBUTE__UNUSED
#endif


/*---------------------------------------------------------------------------*/
/* Marginal Distribution                                                     */
/*---------------------------------------------------------------------------*/

double watson_lpdf (double x, const void *params) {
  /* watson_lpdf <- function(x, kappa, d) { kappa*x*x+(d-3)*log(1-x^2)/2 } */
  double kappa = ((const double*) params)[0];
  double d = ((const double*) params)[1];

  if (d == 3.) {
    return (kappa*x*x);
  }
  else {
    return (kappa*x*x + (d-3.)*log(1-x*x)/2.);
  }
}

double watson_dlpdf (double x, const void *params) {
  /* watson_dlpdf <- function(x, kappa, d) { 2*kappa*x -(d-3)*x/(1-x^2) } */
  double kappa = ((const double*) params)[0];
  double d = ((const double*) params)[1];

  if (d == 3.) {
    return (2.*kappa*x);
  }
  else {
    return (2.*kappa*x - (d-3.)*x/(1.-x*x));
  }
}

double watson_d2lpdf (double x, const void *params) {
  /* watson_d2lpdf <- function(x, kappa, d) { 2*kappa -(d-3)*(1+x^2)/(1-x^2)^2 } */
  double kappa = ((const double*) params)[0];
  double d = ((const double*) params)[1];
  double xsq = x*x;

  if (d == 3.) {
    return (2.*kappa);
  }
  else {
    return (2.*kappa - (d-3.)*(1.+xsq)/((1-xsq)*(1-xsq)));
  }
}


/*---------------------------------------------------------------------------*/
/* Create TinflexC object                                                    */
/* Wrapper for Tinflex_lib_setup()                                           */
/*---------------------------------------------------------------------------*/


typedef TINFLEX_GEN * TINFLEX_SETUP_FUNC (TINFLEX_FUNCT *lpdf, TINFLEX_FUNCT *dlpdf, TINFLEX_FUNCT *d2lpdf,
					  const void *params,
					  int n_ib, const double *ib,
					  int n_c, const double *c,
					  double rho, int max_intervals);

typedef SEXP TINFLEX_SAMPLE_FUNC (TINFLEX_GEN *gen, int n);
typedef double TINFLEX_SAMPLE_DOUBLE_FUNC (TINFLEX_GEN *gen);


/*---------------------------------------------------------------------------*/
/* Create TinflexC object                                                    */
/*---------------------------------------------------------------------------*/

SEXP Tinflexsampler_tag(void) {
  static SEXP tag = NULL; 
  if (!tag) tag = install("R_TINFLEX_C_TAG");
  return tag;
} /* end Tinflexsampler_tag() */

/*...........................................................................*/

void Tinflexsampler_free (SEXP sexp_gen)
{
  TINFLEX_GEN *gen;

  static void * (*free_func)(TINFLEX_GEN *) = NULL;
  if (free_func == NULL)
    free_func = R_GetCCallable("Tinflex", "Tinflex_lib_free");
  
  gen = R_ExternalPtrAddr(sexp_gen);
  free_func (gen);
  R_ClearExternalPtr(sexp_gen);
} /* end of Tinflexsampler_free() */ 

/*...........................................................................*/

SEXP Tinflexsampler_setup (SEXP sexp_obj,
			   SEXP sexp_params,
			   SEXP sexp_ib, SEXP sexp_c,
			   SEXP sexp_rho, SEXP sexp_max_intervals)
{
  const double *params;
  const double *ib;
  int n_ib;
  const double *c;
  int n_c;
  double rho;
  int max_intervals;

  TINFLEX_GEN *gen;
  SEXP sexp_gen = R_NilValue;

  /* set PDF */
  TINFLEX_FUNCT *lpdf = watson_lpdf;
  TINFLEX_FUNCT *dlpdf = watson_dlpdf;
  TINFLEX_FUNCT *d2lpdf = watson_d2lpdf;

  /* get Tinflex setup function */
  static TINFLEX_SETUP_FUNC *Tinflex_setup = NULL;
  if (Tinflex_setup == NULL)
    Tinflex_setup = (TINFLEX_SETUP_FUNC*) R_GetCCallable("Tinflex", "Tinflex_lib_setup");

  /* extract arguments */
  params = REAL(sexp_params);
  ib = REAL(sexp_ib);
  n_ib = length(sexp_ib);
  c = REAL(sexp_c);
  n_c = length(sexp_c);
  rho = NUMERIC_VALUE(sexp_rho);
  max_intervals = INTEGER_VALUE(sexp_max_intervals);

  /* run setup */
  gen = Tinflex_setup (lpdf, dlpdf, d2lpdf, params,
		       n_ib, ib, n_c, c, rho, max_intervals);

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, Tinflexsampler_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, Tinflexsampler_free);

  /* return pointer to R */
  UNPROTECT(1);
    
  return (sexp_gen);

} /* end of Tinflexsampler_setup() */

/*---------------------------------------------------------------------------*/

SEXP Tinflexsampler_sampler (SEXP sexp_n,
			     SEXP sexp_params,
			     SEXP sexp_ib, SEXP sexp_c,
			     SEXP sexp_rho, SEXP sexp_max_intervals)
{
  int n;
  const double *params;
  const double *ib;
  int n_ib;
  const double *c;
  int n_c;
  double rho;
  int max_intervals;

  TINFLEX_GEN *gen;

  SEXP sexp_res;

  /* set PDF */
  TINFLEX_FUNCT *lpdf = watson_lpdf;
  TINFLEX_FUNCT *dlpdf = watson_dlpdf;
  TINFLEX_FUNCT *d2lpdf = watson_d2lpdf;

  /* get Tinflex setup function */
  static TINFLEX_SETUP_FUNC *Tinflex_setup = NULL;
  if (Tinflex_setup == NULL)
    Tinflex_setup = (TINFLEX_SETUP_FUNC*) R_GetCCallable("Tinflex", "Tinflex_lib_setup");

  /* get Tinflex sampler function */
  static TINFLEX_SAMPLE_FUNC *Tinflex_sample = NULL;
  if (Tinflex_sample == NULL)
    Tinflex_sample = (TINFLEX_SAMPLE_FUNC*) R_GetCCallable("Tinflex", "Tinflex_lib_sample");

  /* get Tinflex free function */
  static void * (*Tinflex_free)(TINFLEX_GEN *) = NULL;
  if (Tinflex_free == NULL)
    Tinflex_free = R_GetCCallable("Tinflex", "Tinflex_lib_free");

  /* extract arguments */
  n = INTEGER_VALUE(sexp_n);
  params = REAL(sexp_params);
  ib = REAL(sexp_ib);
  n_ib = length(sexp_ib);
  c = REAL(sexp_c);
  n_c = length(sexp_c);
  rho = NUMERIC_VALUE(sexp_rho);
  max_intervals = INTEGER_VALUE(sexp_max_intervals);

  /* run setup */
  gen = Tinflex_setup (lpdf, dlpdf, d2lpdf, params,
		       n_ib, ib, n_c, c, rho, max_intervals);

  /* generate random sample */
  PROTECT(sexp_res = Tinflex_sample (gen, n));

  /* destroy generator object */
  Tinflex_free (gen);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_res);

} /* end of Tinflexsampler_sampler() */

/*---------------------------------------------------------------------------*/
