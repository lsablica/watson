#include "Rcpp.h"
// Define the method signature

#ifdef __cplusplus
extern "C" {
#endif
   
#include "rdist.h"
   
#ifdef __cplusplus
}
#endif


// [[Rcpp::export]]
Rcpp::NumericVector Tinflexsampler_sampler_from_c(Rcpp::NumericVector& sexp_n,
                                                  Rcpp::NumericVector& sexp_params,
                                                  Rcpp::NumericVector& sexp_ib,
                                                  Rcpp::NumericVector& sexp_c,
                                                  Rcpp::NumericVector& sexp_rho,
                                                  Rcpp::NumericVector& sexp_max_intervals) {
   
   // Compute the result in _C_
   // Import it into R
   SEXP ab = Tinflexsampler_sampler(sexp_n,
                                    sexp_params,
                                    sexp_ib, sexp_c,
                                    sexp_rho, sexp_max_intervals);
   
   // Cast as a NumericVector 
   Rcpp::NumericVector result( ab );
   
   
   // Return result
   return result;
}
