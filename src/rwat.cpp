#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "Rcpp.h"
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif
   
#include "rdist.h"
   
#ifdef __cplusplus
}
#endif

double ACGvsTinflex(int n, double kappa, double d){
   if(d == 2){
      return 10; //acg always
   }
   Environment pkg = Environment::namespace_env("watson");
   arma::cube tinflexdata = pkg["resultTinflex"];
   arma::cube acgdata = pkg["resultACG"];
   arma::vec ns = {1, 3, 5, 10, 20, 50, 100, 500, 1000, 10000};
   arma::vec kappas = {-100, -50, -10, -1, 1, 10, 50, 100};
   arma::vec ds = {3, 5, 10, 20, 50, 100, 200, 1000};
   double npos = sum(ns <= n);
   if(npos == ns.n_elem) npos-- ;
   double kappapos = sum(kappas <= kappa);
   if(kappapos == kappas.n_elem) kappapos-- ;
   if(kappapos == 0) kappapos++ ;
   double dpos = sum(ds <= d);
   if(dpos == ds.n_elem) dpos-- ;
   double npart = (n - ns(npos-1))/(ns(npos)-ns(npos-1));
   double kappaart = (kappa - kappas(kappapos-1))/(kappas(kappapos)-kappas(kappapos-1));
   double dpart = (d - ds(dpos-1))/(ds(dpos)-ds(dpos-1));

   double Tinflexnmkm = tinflexdata(npos-1,kappapos-1,dpos-1) + dpart*(tinflexdata(npos-1,kappapos-1,dpos) - tinflexdata(npos-1,kappapos-1,dpos-1));
   double Tinflexnmk = tinflexdata(npos-1,kappapos,dpos-1) + dpart*(tinflexdata(npos-1,kappapos,dpos) - tinflexdata(npos-1,kappapos,dpos-1));
   double Tinflexnkm = tinflexdata(npos,kappapos-1,dpos-1) + dpart*(tinflexdata(npos,kappapos-1,dpos) - tinflexdata(npos,kappapos-1,dpos-1));
   double Tinflexnk = tinflexdata(npos,kappapos,dpos-1) + dpart*(tinflexdata(npos,kappapos,dpos) - tinflexdata(npos,kappapos,dpos-1));
   
   double Tinflexnm = Tinflexnmkm + kappaart*(Tinflexnmk - Tinflexnmkm);
   double Tinflexn = Tinflexnkm + kappaart*(Tinflexnk - Tinflexnkm);
   
   double Tinflex = Tinflexnm + npart*(Tinflexn - Tinflexnm);
   
   
   double ACGnmkm = acgdata(npos-1,kappapos-1,dpos-1) + dpart*(acgdata(npos-1,kappapos-1,dpos) - acgdata(npos-1,kappapos-1,dpos-1));
   double ACGnmk = acgdata(npos-1,kappapos,dpos-1) + dpart*(acgdata(npos-1,kappapos,dpos) - acgdata(npos-1,kappapos,dpos-1));
   double ACGnkm = acgdata(npos,kappapos-1,dpos-1) + dpart*(acgdata(npos,kappapos-1,dpos) - acgdata(npos,kappapos-1,dpos-1));
   double ACGnk = acgdata(npos,kappapos,dpos-1) + dpart*(acgdata(npos,kappapos,dpos) - acgdata(npos,kappapos,dpos-1));
   
   double ACGnm = ACGnmkm + kappaart*(ACGnmk - ACGnmkm);
   double ACGn = ACGnkm + kappaart*(ACGnk - ACGnkm);
   
   double ACG = ACGnm + npart*(ACGn - ACGnm);
   
   //Rcout << "ACG:" << ACG << std::endl;
   //Rcout << "Tinflex:" << Tinflex << std::endl;
   //Rcout << "Tinflex/ACG:" << Tinflex/ACG << std::endl;
   
   return Tinflex/ACG; // 1> Tinflex better, >1 ACG better
}

NumericVector Tinflexsampler_sampler_from_c(int n,
                                            double kappa,
                                            double d,
                                            double cT,
                                            double rho) {
   NumericVector sexp_params = {kappa,d};
   NumericVector sexp_c = {cT};
   NumericVector sexp_rho = {rho};
   IntegerVector sexp_n = {n};
   NumericVector sexp_ib = {0,1};
   NumericVector sexp_max_intervals = {1001};
   SEXP ab = Tinflexsampler_sampler(sexp_n, sexp_params, sexp_ib, sexp_c,
                                    sexp_rho, sexp_max_intervals);
   Rcpp::NumericVector result( ab );
   return result;
}

// [[Rcpp::export]]
arma::mat rwatTinflex(int n, double kappa, arma::vec &mu, double cT, double rho){
   double norm = as_scalar(sum(pow(mu,2)));
   int p = mu.n_elem;
   arma::mat A(n, p, arma::fill::randn);
   if(kappa == 0 || norm == 0){/*uniform*/
      return normalise(A,2,1);
   }
   mu = mu/sqrt(norm);
   NumericVector v = Tinflexsampler_sampler_from_c(n, kappa, p, cT, rho);
   arma::vec w = as<arma::vec>(wrap(v));
   arma::vec choice = {-1, 1};
   arma::vec index = RcppArmadillo::sample(choice, n, true); 
   w = w % index;
   A = A - A*mu*mu.t() ;
   A = arma::normalise(A, 2, 1);
   A = A.each_col() % sqrt(1 - w%w);
   A = w*mu.t() + A;
   return A;
} 

// [[Rcpp::export]]
arma::mat rwatACG(int n, double kappa, arma::vec &mu, double b = -10){
   double norm = as_scalar(sum(pow(mu,2)));
   int p = mu.n_elem;
   arma::mat A(n, p);
   if(kappa == 0 || norm == 0){/*uniform*/
      return normalise(A.randn(),2,1);
   }
   mu = mu/sqrt(norm);
   int count = 0;
   int Nt = 0;
   double beta1bar, beta2bar, mutx, kappacross2, unif, known, norm2, b1mutx2, ratio;
   arma::vec candidate;
   if(b<=0) b = p*0.5 + kappa + sqrt(kappa*kappa + p*p*0.25 - kappa*p + 2*kappa) ;
   // Rcout << "b:" << b << std::endl;
   known = -p*0.5*log(p) + 0.5*(p-b);
   beta2bar = -1 + sqrt(b/(b-2*kappa));
   beta1bar = beta2bar*(beta2bar + 2);
   while(count<n){
      candidate = arma::randn<arma::vec>(p);
      mutx = arma::dot(mu, candidate) ;
      b1mutx2 = beta1bar*mutx*mutx;
      norm2 = arma::dot(candidate,candidate) + b1mutx2;
      kappacross2 = kappa*(mutx*mutx + b1mutx2)/norm2;
      unif = arma::randu<double>();
      ratio = kappacross2 + p*0.5*log(b-2*kappacross2) + known;
      if(log(unif)<ratio){
         candidate = (candidate + beta2bar*mutx*mu)/sqrt(norm2);
         A.row(count) = arma::trans(candidate);
         count += 1;
      }
      Nt += 1;
      if(Nt % 1000000 == 0) Rcpp::checkUserInterrupt();
   }  
   return A;
} 

//' @title Random Sampling from a Mixture of Watson Distributions
//' @description \code{rmwat} generates a random sample from a mixture of multivariate Watson distributions.
//' @param n an integer giving the number of samples to draw.
//' @param weights a numeric vector with non-negative elements giving the mixture probabilities.
//' @param kappa a numeric vector giving the kappa parameters of the mixture components.
//' @param mu a numeric matrix with columns giving the mu parameters of the mixture components.
//' @param method a string indicating whether ACG sampler (\code{method = "acg"}), Tinflex sampler (\code{method = "tinflex"}) or automatic selection (\code{method = "auto"}) of the sampler should be used, default: "acg".  
//' @param b a positive numeric hyper-parameter used in the sampling. If not a positive value is given, optimal choice of b is used, default: -10.
//' @param rho performance parameter: requested upper bound for ratio of area below hat to area below squeeze (numeric). See \code{\link[Tinflex]{Tinflex.setup}}, default: 1.1.
//' @return  A matrix with rows equal to the generated values.
//' @details The function generates samples from finite mixtures of Watson distributions,
//'          using methods from Sablica, Hornik and Leydold (2022) \url{https://research.wu.ac.at/en/publications/random-sampling-from-the-watson-distribution}.
//' @examples
//'
//' ## simulate from Watson distribution
//' sample1 <- rmwat(n = 20, weights = 1, kappa = 20, mu = matrix(c(1,1,1),nrow = 3))
//'
//' ## simulate from a mixture of Watson distributions
//' sample2 <- rmwat(n = 20, weights = c(0.5,0.5), kappa = c(-200,-200),
//'                             mu = matrix(c(1,1,1,-1,1,1),nrow = 3))
//' @rdname rmwat
//' @references Sablica, Hornik and Leydold (2022). Random Sampling from the Watson Distribution \url{https://research.wu.ac.at/en/publications/random-sampling-from-the-watson-distribution}.
//' @export
// [[Rcpp::export]]
NumericMatrix rmwat(int n, arma::vec &weights, arma::vec kappa, arma::mat &mu, String method = "acg",
                    double b = -10, double rho=1.1){
  weights = arma::normalise(weights, 1);
  int p = mu.n_rows;
  int K = mu.n_cols;
  arma::mat A(n, p);
  arma::uvec sample = RcppArmadillo::sample(arma::regspace<arma::uvec>(0, K-1), n, true, weights);
  int size;
  String method2;
  arma::uvec which;
  arma::vec mus;
  for(int i = 0; i < K; i++) {
    which = arma::find(sample==i);
    size = which.n_elem;
    mus = mu.col(i);
    if(size>0){
       method2 = (method == "auto") ? ((ACGvsTinflex(size, kappa(i), p) < 1) ? "tinflex" : "acg") : method;
       if(method2 == "acg"){
          A.rows(which) = rwatACG(size, kappa(i), mus, b);
       } else {
          A.rows(which) = rwatTinflex(size, kappa(i), mus, 0, rho);
       } 
    }
  }
  sample = sample + 1;
  IntegerVector fac = wrap(sample);
  IntegerVector lev = seq(1,K);
  fac.attr("dim") = R_NilValue;
  fac.attr("levels") = as<CharacterVector>(lev);
  fac.attr("class") = "factor";
  
  NumericMatrix res = wrap(A);
  res.attr("id") = fac;
  res.attr("class") = "rmwat";
  return res;
}
