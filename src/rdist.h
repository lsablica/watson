#include "Tinflex_lib.h"

SEXP Tinflexsampler_setup (SEXP sexp_obj,
			   SEXP sexp_params,
			   SEXP sexp_ib, SEXP sexp_c,
			   SEXP sexp_rho, SEXP sexp_max_intervals);

SEXP Tinflexsampler_sampler (SEXP sexp_n,
			     SEXP sexp_params,
			     SEXP sexp_ib, SEXP sexp_c,
			     SEXP sexp_rho, SEXP sexp_max_intervals);
