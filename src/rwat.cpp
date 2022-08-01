#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rwat(int n, double kappa, arma::vec &mu, double b = -10){
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
//' @param b a positive numeric hyper-parameter used in the sampling. If not a positive value is given, optimal choice of b is used, default: -10.
//' @return  A matrix with rows equal to the generated values.
//' @details The function generates samples from finite mixtures of Watson distributions,
//'          using adjusted BACG algorithm of Kent (2013) for the case of Watson distribution. The algorithm is of the
//'          rejection-sampling form.
//' @examples
//'
//' ## simulate from Watson distribution
//' sample1 <- rmwat(n = 20, weights = 1, kappa = 20, mu = matrix(c(1,1,1),nrow = 3))
//'
//' ## simulate from a mixture of Watson distributions
//' sample2 <- rmwat(n = 20, weights = c(0.5,0.5), kappa = c(-200,-200),
//'                             mu = matrix(c(1,1,1,-1,1,1),nrow = 3))
//' @rdname rmwat
//' @references Kent J.T., Ganeiber A.M. and Mardia K.V. (2013). A new method to simulate the Bingham and related distributions
//'   in directional data analysis with applications \url{http://arxiv.org/pdf/1310.8110v1.pdf}
//' @export
// [[Rcpp::export]]
NumericMatrix rmwat(int n, arma::vec &weights, arma::vec kappa, arma::mat &mu, double b = -10){
  weights = arma::normalise(weights, 1);
  int p = mu.n_rows;
  int K = mu.n_cols;
  arma::mat A(n, p);
  arma::uvec sample = RcppArmadillo::sample(arma::regspace<arma::uvec>(0, K-1), n, true, weights);
  int size;
  arma::uvec which;
  arma::vec mus;
  for(int i = 0; i < K; i++) {
    which = arma::find(sample==i);
    size = which.n_elem;
    // Rcout << size << std::endl;
    mus = mu.col(i);
    if(size>0){
      A.rows(which) = rwat(size, kappa(i), mus, b);
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
