#' @useDynLib watson
#' @importFrom Rcpp sourceCpp
#' @import Tinflex
NULL
#' @export
print.watfit <- function(x, ...){
  cat("Fitted " , x$K ,"-components Watson mixture: \n\n", sep = "")
  cat("Weights:", x$weights, "\n")
  cat("Kappa:", x$kappa_vector, "\n\n")
  cat("Mu: \n")
  print(x$mu_matrix)
  cat("\nLog-likelihood: ", x$L, ",  Average log-likelihood: ", x$L/x$details$n  , "\n", sep = "")
}
#' @export
coef.watfit <- function(object, ...) object[c("weights", "kappa_vector", "mu_matrix")]

#' @export
logLik.watfit <- function(object, newdata, ...){
  if (missing(newdata)) {
    return(object$ll)
  } else {
    if(inherits(newdata, "Matrix")){
      log_like2(newdata, object$kappa_vector, object$mu_matrix, object$weights, object$K, object$details$beta, nrow(newdata))
    } else{
      log_like1(newdata, object$kappa_vector, object$mu_matrix, object$weights, object$K, object$details$beta, nrow(newdata))
    }
  }
}
#' @export
predict.watfit <- function(object, newdata = NULL, type = c("class_ids", "memberships"), ...){
    type <- match.arg(type)
    M <- if(is.null(newdata)){
      object$a_posterior
    } else {
      if(inherits(newdata, "Matrix")){
        predictC2(newdata, object$kappa_vector, object$mu_matrix, object$weights, object$details$E, object$K)
      } else{
        predictC1(newdata, object$kappa_vector, object$mu_matrix, object$weights, object$details$E, object$K)
      }
    }
    v <- if(type == "class_ids") max.col(M) else M
    attr(v, "loglik") <- attr(M, "loglik")
    v
}

#' @title Diametrical clustering
#' @description  \code{diam_clus} clusters axial data on sphere using  the algorithm proposed in Dhillon et al. (2003). 
#' @param x a numeric data matrix, with rows corresponding to observations. Can be a dense matrix,
#'          or any of the supported sparse matrices from \code{\link[Matrix]{Matrix}} package by \code{\link[Rcpp]{Rcpp}}.
#' @param k an integer giving the number of mixture components.
#' @param niter integer indicating the number of iterations of the diametrical clustering algorithm, default: 100.
#' @return a matrix with the concentration directions with an attribute "id" defining the classified categories.
#' @examples
#' ## Generate a sample 
#' a <- rmwat(n = 200, weights = c(0.5,0.5), kappa = c(20, 20),
#'                         mu = matrix(c(1,1,-1,1),nrow = 2))
#' ## Fit basic model
#' q <- diam_clus(a, 2)
#' @rdname diam_clus
#' @references Inderjit S Dhillon, Edward M Marcotte, and Usman Roshan. Diametrical clustering 
#' for identifying anti-correlated gene clusters. Bioinformatics, 19(13):1612-1619, 2003.
#' @export
diam_clus <- function(x, k, niter = 100){
  if(inherits(x, "Matrix")){
    a = diam_clus2(x,k,niter)
  } else{
    a = diam_clus1(x,k,niter)
  }
  a
}


#' @title Fit Mixtures of Watson Distributions
#' @description  \code{watson} fits a finite mixture of multivariate Watson distributions.
#' @param x a numeric data matrix, with rows corresponding to observations. Can be a dense matrix,
#'          or any of the supported sparse matrices from  \code{\link[Matrix]{Matrix}} package by \code{\link[Rcpp]{Rcpp}}.
#' @param k an integer giving the number of mixture components.
#' @param control a list of control parameters. See Details.
#' @param ... a list of control parameters (overriding those specified in control).
#' @return An object of class "watfit" representing the fitted mixture, which is a list containing the fitted
#'         weights, concentrations parameters (kappa), concentrations directions (mu) and further metadata.
#' @details watson returns an object of class "watfit" representing the fitted mixture of Watson
#' distributions model. Available methods for such objects include \code{\link[stats]{coef}}, \code{\link[stats]{logLik}}, \code{\link[base]{print}} and \code{\link[stats]{predict}}.
#' \code{\link[stats]{predict}} has an extra type argument with possible values \code{"class_ids"} (default) and \code{"memberships"}
#' for indicating hard or soft prediction, respectively.
#'
#' The mixture of Watson distributions is fitted using EM variants as specified by control
#' option E (specifying the E-step employed), with possible values "softmax" (default), "hardmax" or
#' "stochmax". For "stochmax", class assignments are drawn from the posteriors for each observation in the E-step as
#'  outlined as SEM in Celeux and Govaert (1992). The stopping criterion for this algorithm is by default
#'  changed to not check for convergence (logical control option converge), but to return the parameters with
#'  the maximum likelihood encountered.
#'
#'  In the M-step, the parameters of the respective component distributions
#'  are estimated via maximum likelihood, which is accomplished by solving the equation
#'   \deqn{g(\alpha,\beta, \kappa)=r,}
#'   where
#'   \deqn{0<\alpha<\beta,  \  0\leq r\leq 1}
#'   and
#'    \deqn{g(\alpha, \beta, \kappa) =  (\alpha/\beta)M(\alpha+1, \beta+1,  \kappa)/M(\alpha, \beta,  \kappa),} with M being the Kummer's function.
#'   Via control argument M, one can specify how to (approximately) solve these equations.
#'   The possible methods are:
#'
#'   "Sra_2007"
#'   uses the approximation of Sra (2007).
#'
#'   "BBG"
#'   uses the approximation of Bijral et al. (2007), without the correction term.
#'
#'   "BBG_c"
#'  uses the approximation of Bijral et al. (2007), with the correction term.
#'
#'   "Sra_Karp_2013"
#'   uses the bounds derived in Sra and Karp (2013), with the decision rule .
#'
#'   "bisection"
#'   uses a bisection to solve the problem using evaluation proposed in Writeup1 (2018).
#'
#'   "newton"
#'   uses a bracketet type of Neton algorithm to solve the problem using evaluation proposed in Writeup1 (2018). (Default.)
#'
#'   "lognewton"
#'   uses a bracketet type of Neton algorithm to solve the problem \eqn{log(g((\alpha, \beta, \kappa)) = log(r)}
#'   using evaluation proposed in Writeup1 (2018).
#'
#'   Additional control parameters are as follows.
#'
#'   maxiter
#'   an integer giving the maximal number of EM iterations to be performed, default: 100.
#'
#'   reltol
#'   the minimum relative improvement of the objective function per iteration. If improvement is less, the EM
#'   algorithm will stop under the assumption that no further significant improvement can be made, defaults
#'   to sqrt(.Machine$double.eps).
#'
#'   ids
#'   either a vector of class memberships or TRUE which implies that the class memberships are obtained from
#'   the attribute named "id" of the input data; these class memberships are used for initializing the EM
#'   algorithm and the algorithm is stopped after the first iteration.
#'
#'   init_iter
#'   a numeric vector setting the number of diametrical clustering iterations to do, before the EM starts, default: 0.
#'
#'   start
#'   a specification of the starting values to be employed. Can be a list of matrices giving the memberships
#'    of objects to components. This information is combined with the \code{init_iter} parameter and together form
#'    the initialization procedure. If nothing is given, the starting values are drwan randomly.
#'
#'   If several starting values are specified, the EM algorithm is performed individually to each starting
#'   value, and the best solution found is returned.
#'
#'   nruns
#'   an integer giving the number of EM runs to be performed. Default: 1. Only used if start is not given.
#'
#'   minweight
#'   a numeric indicating the minimum prior probability. Components falling below this threshold are removed
#'   during the iteration. If is greater than 1, the value is taken as the minimal number of observations in a component, default: 0 if
#'   E = "softmax" and 2 if other type of E-method is used .
#'
#'   converge
#'   a logical, if TRUE the EM algorithm is stopped if the reltol criterion is met and the current parameter
#'   estimate is returned. If FALSE the EM algorithm is run for maxiter iterations and the parametrizations
#'   with the maximum likelihood encountered during the EM algorithm is returned. Default: TRUE, changed to
#'   FALSE if E="stochmax".
#'
#'   N
#'   an integer indicating number of iteration used when the Kummer function is approximate, default: 30.
#'   
#'   verbose
#'   a logical indicating whether to provide some output on algorithmic progress, default: FALSE.
#'
#' @examples
#' \donttest{
#' ## Generate a sample with two orthogonal circles (negative kappas)
#' a <- rmwat(n = 200, weights = c(0.5,0.5), kappa = c(-200,-200),
#'                         mu = matrix(c(1,1,1,-1,1,1),nrow = 3))
#' ## Fit basic model
#' q <- watson(a, 2)
#' ## Fit the models, giving the true categories
#' q <- watson(a, 2, ids=TRUE)
#' ## Fit a model with hard-assignment, and 50 runs
#' q <- watson(a, 2, E = "hard", nruns = 50)
#' ## Print details
#' q
#' ## Extract coefficients
#' coef(q)
#' ## Calculate likelihood for new data
#' a2 <- rmwat(n = 20, weights = c(0.5,0.5), kappa = c(-200,-200),
#'                         mu = matrix(c(1,1,1,-1,1,1),nrow = 3))
#' logLik(q, a2)
#' ## Compare the fitted classes to the true ones:
#' table(True = attr(a, "id"), Fitted = predict(q))
#' }
#' @rdname watson
#' @export

watson <- function(x, k, control = list(), ...){
  control <- c(control, list(...))
  
  idd <- attr(x, "id")
  x <- as.matrix(x)
  attr(x, "id") <- idd
  n <- nrow(x)
  p <- ncol(x)
  sparse <- inherits(x, "Matrix")

  maxiter <- control$maxiter
  if(is.null(maxiter)) maxiter <- 100L

  reltol <- control$reltol
  if(is.null(reltol)) reltol <- sqrt(.Machine$double.eps)


  E_methods <- c("softmax", "hardmax", "stochmax")
  E_type <- control$E
  if(is.null(E_type)){
    E_type <- "softmax"
  } else {
    pos <- pmatch(tolower(E_type), tolower(E_methods))
    if(is.na(pos))
      stop("Invalid E-step method.")
    E_type <- E_methods[pos]
  }

  M_methods <- c("newton","lognewton", "bisection", "BBG", "Sra_2007", "BBG_c", "Sra_Karp_2013")
  M_type <- control$M
  if(is.null(M_type)){
    M_type <- "newton"
  } else {
    pos <- pmatch(tolower(M_type), tolower(M_methods))
    if(is.na(pos))
      stop("Invalid M-step method.")
    M_type <- M_methods[pos]
  }

  minweight <- control$minweight
  if(is.null(minweight)){
    minweight <- if(E_type == "softmax") 0 else 2/n
  } else if(minweight >= 1 && n>minweight){
    minweight <- minweight/n
  } else if(n <= minweight || minweight<0){
    stop("minweight must be in [0,1) or [1,n)")
  }

  converge <- control$converge
  if(is.null(converge)) {
    converge <- if(E_type == "stochmax") FALSE else TRUE
  }

  N <- control$N
  if(is.null(N)) {
    N <- 30L
  }
  
  verbose <- control$verbose
  if (is.null(verbose)) verbose <- getOption("verbose")

  ids <- control$ids
  start <- control$start
  init_iter <- control$init_iter
  nruns <- control$nruns
  if(!is.null(ids)) {
    if(identical(ids, TRUE)){
      ids <- attr(x, "id")
      ids <- as.integer(as.factor(ids))
    } else {
      ids <- as.integer(as.factor(ids))
      if(length(ids)!=n) stop("length of 'ids' needs to match the number of observations")
    }
    k <- max(ids) 
    maxiter <- 0L
    init_iter <- if(is.null(init_iter)) 0L else init_iter[1]
    if(!is.null(nruns) || !is.null(start))
      warning("control arguments 'nruns' and 'start' ignored because 'ids' are specified")
    nruns <- 1L
    beta_mat <- matrix(0, n, k)
    beta_mat[cbind(seq_len(n), ids)] <- 1
    start <- list(list(given = TRUE, init_iter = init_iter, matrix = beta_mat))
  } else{
    if(is.null(start)) {
      if(is.null(nruns)){
        nruns <- if(is.null(init_iter)) 1L else length(init_iter)
      }else{
        if(nruns >= 1L){
          nruns <- floor(nruns)
        } else{
          stop("nruns must be greater or equal than 1")
        }
      }
      init_iter <- if(is.null(init_iter)) rep_len(as.list(0L),nruns) else rep_len(as.list(init_iter),nruns)
      given <- rep(list(FALSE), nruns)
      matrix <- rep(list(matrix(0,1,1)),nruns)
      start <- lapply(1:nruns, function(i) list(given = given[[i]], init_iter = init_iter[[i]], matrix = matrix[[i]]))
    }
    else {
      matrix <- if(!is.list(start)) rep(list(start),nruns) else start
      nruns <- length(matrix)
      init_iter <- if(is.null(init_iter)) rep_len(as.list(0L),nruns) else rep_len(as.list(init_iter),nruns)
      given <- rep(list(TRUE), nruns)
      start <- lapply(1:nruns, function(i) list(given = given[[i]], init_iter = init_iter[[i]], matrix = matrix[[i]]))
    }
  }
  obj <- if(sparse){
    EM2(data = x, K = k, E_type = E_type, M_type = M_type, minalpha = minweight, convergence = converge,
       maxiter = maxiter, N = N, reltol = reltol, start = start, verbose = verbose)
  } else{
    EM1(data = x, K = k, E_type = E_type, M_type = M_type, minalpha = minweight, convergence = converge,
       maxiter = maxiter, N = N, reltol = reltol, start = start, verbose = verbose)
  }
  details <- list(p = p, n = n, beta = p/2, reltol = reltol, maxiter = maxiter, E = E_type, M = M_type, minweight = minweight, N = N, converge = converge, sparse = sparse)
  K <- dim(obj[[2]])[2]
  ll <- structure(obj[[5]], class = "logLik", df = (p+1)*K - 1, nobs = n)
  obj <- c(list(x), list(K), obj, list(ll), list(details))
  names(obj) <- c("data", "K" ,"a_posterior", "kappa_vector", "mu_matrix", "weights", "L", "ll", "details")
  nam <- paste("clus",1:K, sep = "_")
  colnames(obj$a_posterior) <- nam
  colnames(obj$kappa_vector) <- nam
  colnames(obj$mu_matrix) <- nam
  colnames(obj$weights) <- nam
  class(obj) <- "watfit"
  obj
}
