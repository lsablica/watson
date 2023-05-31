#include <RcppArmadillo.h>
#include "hyper.h"
using namespace Rcpp;
using namespace std;

double gneg(double x, double alpha, double beta, int N){
  double start1 = 2*(alpha+N)/((beta+N)-1-x+sqrt(pow(x-(beta+1+N), 2.0) + 4*x*(alpha+1+N)));
  double start2 = 1-2*(beta-alpha)/((beta+N)-1 + x + sqrt(-4*((beta-alpha)+1)*x + pow(x + (beta+N+1),2)) );
  
  for(int i = N-1; i >= 0; --i) {
    start1 = (alpha+i)/((beta+i) - x + x * start1);
    start2 = (alpha+i)/((beta+i) - x + x * start2);
  }
  return (start2+start1)/2;
}

// [[Rcpp::export]]
double g(double alpha, double beta, double x, int N = 30) {
  if(x==0){
    return alpha/beta ;
  }
  else if(x>0){
    alpha = beta - alpha;
    return 1 - gneg(-x, alpha, beta, N);
  }
  else{
    return gneg(x, alpha, beta, N);
  }
}

// [[Rcpp::export]]
double kummerM(double alpha, double beta, double r) {
  double res; 
  gsl_sf_hyperg_1F1_e(alpha, beta, r, res);
  return res;
}

double B(double x,double c, double d, double alpha, double beta){
   return ((alpha/beta)*(c+d))/(sqrt(x*x+c*c+2*x*((alpha/beta)*(c+d)-d))-x+d);
}

double integal(double r, double alpha, double beta, double c, double d){
   return (beta*beta*(-c + d)*r + (alpha - beta)*(beta*(c - d) + alpha*(c + d))*std::log(1 - r) - alpha*alpha*(c + d)*std::log(r))/(2*alpha*beta);
} 

int log_hyperg_1F1_bounds(double alpha, double beta, double r, double &res){
   double w, rr, c, d, c2, d2, d3, g1, g4, BU, BL, z = 0;
   if(r<0){
      alpha = beta - alpha;
      r = z = -r;
   } 
   w = (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
   rr = alpha*(w*(beta+2)-1-alpha)/(w*alpha*(2+beta)+(1+beta)*(beta-2*alpha));
   c = beta + 1;
   d = (c*(alpha+beta)-2*alpha*beta)/(beta-alpha);
   c2 = beta + w;
   d2 =  (c2*(alpha+beta)-2*alpha*beta)/(beta-alpha);
   d3 = (beta+1)*(beta*beta-alpha*(beta-2))/((beta-alpha)*(beta+2));
   g1 = B(r, c, d, alpha, beta);
   g4 = std::max(B(r, c2, d2, alpha, beta), B(r, c, d3, alpha, beta));
   BU = r*g1 - z - (integal(g1, alpha, beta, c, d) - integal(alpha/beta, alpha, beta, c, d));
   if(g4 > rr){
      BL = r*g4 - z - (integal(g4, alpha, beta, c2, d2) - integal(rr, alpha, beta, c2, d2) + integal(rr, alpha, beta, c, d3) - integal(alpha/beta, alpha, beta, c, d3));
   } else{
      BL = r*g4 - z - (integal(g4, alpha, beta, c, d3) - integal(alpha/beta, alpha, beta, c, d3));
   }
   
   res = (BL+BU)/2; 
   return (BU - BL < 0.03960525)? 0 : 1; // log(1.02)*2 = 0.03960525
}

double log_hyperg_1F1_iter(double alpha, double beta, double r, int N = 10){
  double alpha2 = beta - alpha, summ = 0;
  int flor = floor(alpha2);
  
  if (flor == alpha2){ 
    flor -= 1;
  }
  
  for(int i = 1; i <= flor; i++){
    summ += std::log(beta - i) - std::log(alpha2 - i) + std::log(g(alpha2 - i, beta - i, - r, N));
  }
  if(r>0){
    return summ + r + std::log(kummerM(alpha2-flor, beta-flor, -r));
  } else{
    return summ + std::log(kummerM(alpha , beta-flor, r));
  }
}

// [[Rcpp::export]]
double log_hyperg_1F1(double alpha, double beta, double r, int N = 10){
  double res; 
  int success = gsl_sf_hyperg_1F1_e(alpha, beta, r, res);
  
  if(success==0){
    return std::log(res);
  } else{
    success = gsl_sf_hyperg_1F1_e(beta - alpha, beta, -r, res); 
  }
  
  if(success==0){
     return r+std::log(res);
  } else{
     success = log_hyperg_1F1_bounds(alpha, beta, r, res); 
  }
  
  if(success==0){
    return res;
  } else{
    return log_hyperg_1F1_iter(alpha, beta, r, N);
  }
}

arma::sp_mat extract_rows(const arma::sp_mat &a, arma::uvec &x, double i){
  arma::sp_mat b = a.t(); /*because no .rows for sp_mat*/
  arma::uvec e = find(x==i);
  arma::sp_mat B(b.n_rows, e.n_elem);
  arma::uvec::iterator it     = e.begin();
  arma::uvec::iterator it_end = e.end();
  int j = 0;
  for(; it != it_end; ++it && ++j){
    B.col(j) = b.col(*it);
  }
  return B.t();
}
arma::mat extract_rows(const arma::mat &a, arma::uvec &x, double i){
  return a.rows(find(x==i));
}

template <class T>
void diamclus_internal(const T &data, arma::mat &beta_matrix, arma::mat &mu_matrix, int K, int n, int maxiter){
  K = K-1;
  arma::uvec part = arma::randi<arma::uvec>(n, arma::distr_param(0, K));
  T A;
  for(int i = 0; i <= K; i++) {
    A = extract_rows(data, part, i);
    int nn = A.n_rows;
    arma::rowvec r(sum(A,0)/nn);
    mu_matrix.col(i) = normalise(trans(r));
  }
  int iter = 0;
  arma::uvec oldpart;
  while(iter < maxiter){
    oldpart = part;
    part = index_max(pow(data*mu_matrix, 2), 1); // E-step
    if(arma::all(oldpart == part) && iter != 0) break;
    for(int i = 0; i < K; i++) {     // M-step
      A = extract_rows(data, part, i);
      mu_matrix.col(i) = normalise((trans(A)*A)*mu_matrix.col(i));
    }
    iter +=1;
  }
  arma::uvec j(1);
  for(int i = 0; i <= K; i++) {
    j(0) = i;
    beta_matrix.submat(arma::find(part==i),j).ones();
  }
}

template <class T>
NumericMatrix diam_clus( T &data, int K, int maxiter = 100) {
  data = normalise(data, 2, 1);
  int p = data.n_cols;
  int n = data.n_rows;
  arma::mat beta_matrix(n, K, arma::fill::zeros), mu_matrix(p, K);
  diamclus_internal(data, beta_matrix, mu_matrix, K, n, maxiter);
  arma::uvec categ = index_max(beta_matrix, 1);
  NumericVector V = wrap(categ);
  NumericMatrix M = wrap(mu_matrix);
  V.attr("dim") = R_NilValue;
  M.attr("id") = V+1;
  return M;
}
//[[Rcpp::export]]
NumericMatrix diam_clus1(arma::mat &data, int K, int maxiter = 100){
  return diam_clus(data, K, maxiter);
}
// [[Rcpp::export]]
NumericMatrix diam_clus2(arma::sp_mat &data, int K, int maxiter = 100){
  return diam_clus(data, K, maxiter);
}

void hard(arma::mat &beta_matrix, int K, int n){
  arma::uvec j(1), maxindex = index_max( beta_matrix, 1);
  beta_matrix.zeros();
  for(int i = 0; i < K; i++) {
    j(0) = i;
    beta_matrix.submat(arma::find(maxindex==i),j).ones();
  }
  return;
}
void soft(arma::mat &beta_matrix, int K, int n){
  return;
}
void stoch(arma::mat &beta_matrix, int K, int n){
  arma::uvec j(1), maxindex = arma::sum(arma::repelem(arma::randu(n),1,K)>arma::cumsum(beta_matrix, 1), 1);
  beta_matrix.zeros();
  for(int i = 0; i < K; i++) {
    j(0) = i;
    beta_matrix.submat(arma::find(maxindex==i),j).ones();
  }
  return;
}
template <class T>
bool E_step(const T &data, arma::mat &beta_matrix, arma::vec &kappa_vector, arma::mat &mu_matrix,
            const arma::rowvec &pi_vector, void (*E_method)(arma::mat&, int, int), int &K, bool convergence, double minalpha,
            double beta, int n, double p, double &lik, double reltol, arma::mat &max_beta_matrix, 
            arma::vec &max_kappa_vector, arma::mat &max_mu_matrix, arma::rowvec &max_pi_vector, double &max_log_lik){
  arma::mat A = pow(data*mu_matrix,2);
  A.each_row() %= kappa_vector.t();
  arma::rowvec kummer_vector(K);
  for(int i = 0; i < K; i++) {
    kummer_vector(i) = log_hyperg_1F1(0.5, beta, kappa_vector(i), 10);
  }
  A += arma::repelem(log(pi_vector),n,1) - arma::repelem(kummer_vector,n,1);
  
  arma::vec maxx = max( A, 1);
  maxx += log(sum(exp(A.each_col() - maxx),1));
  double lik_new = sum(maxx); 
  if(convergence && std::abs(lik - lik_new) < reltol * (std::abs(lik) + reltol)){
    lik = lik_new;
    return true;
  } else{
    if(!convergence && lik_new > max_log_lik){
      max_beta_matrix  = beta_matrix;
      max_mu_matrix    = mu_matrix;
      max_pi_vector    = pi_vector;
      max_kappa_vector = kappa_vector;
      max_log_lik      = lik_new;
    }
    lik = lik_new;
    A.each_col() -= maxx;
    beta_matrix = exp(A);
    E_method(beta_matrix, K, n);
    if(minalpha>0){
      beta_matrix = normalise(beta_matrix.cols(find((sum(beta_matrix)/n) >minalpha)),1,1);
      K = beta_matrix.n_cols;
      mu_matrix.resize(p,K);
      kappa_vector.resize(K);
    }
    return false;
  }
}

template <class T>
NumericMatrix predictC(T &data, arma::vec &kappa_vector, arma::mat &mu_matrix, arma::rowvec &pi_vector, String E_type, int K){
  data = normalise(data, 2, 1);
  double p = data.n_cols;
  int n = data.n_rows;
  double beta  = p/2;
  void (*E_method)(arma::mat&, int, int);
  
  if(E_type == "softmax"){
    E_method = soft;
  } else if(E_type == "hardmax"){
    E_method = hard;
  } else{ // stochmax
    E_method = stoch;
  }
  arma::mat beta_matrix(n, K, arma::fill::zeros);
  double log_lik = -1e11, max_log_lik = 1e16;
  E_step(data, beta_matrix, kappa_vector, mu_matrix, pi_vector, E_method, K, false, 0, beta, n, p, log_lik, 0, 
         beta_matrix, kappa_vector, mu_matrix, pi_vector, max_log_lik);
  NumericMatrix betamat = wrap(beta_matrix);
  betamat.attr("loglik") = log_lik;
  return betamat;
}


// [[Rcpp::export]]
NumericMatrix predictC1(arma::mat &data, arma::vec &kappa_vector, arma::mat &mu_matrix, arma::rowvec &pi_vector, String E_type, int K){
  return predictC(data, kappa_vector, mu_matrix, pi_vector, E_type, K);
}
// [[Rcpp::export]]
NumericMatrix predictC2(arma::sp_mat &data, arma::vec &kappa_vector, arma::mat &mu_matrix, arma::rowvec &pi_vector, String E_type, int K){
  return predictC(data, kappa_vector, mu_matrix, pi_vector, E_type, K);
}

double hybridnewton(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  double x,a,b,rr,w;
  bool change = true;
  if(r < alpha/beta){
     alpha = beta-alpha; 
     r = 1-r;
    change = false;
  }
   
  w =  (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
  rr =  (w*w*alpha+w*(beta-alpha)*alpha)/((beta+1)*(beta-alpha)+w*w*alpha); 
   
  if(r > rr){
    double c3,c4;
    c3 =  beta + (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
    c4 =  (c3*(alpha+beta)-2*alpha*beta)/(beta-alpha);
    a = -(r*beta-alpha)/(r*(1-r))*((c3*(beta*r+alpha)-c4*(beta*r-alpha))/(2*alpha*beta));
    alpha = beta-alpha; 
    r = 1-r;
  } else{
    r = 1-r;
    alpha = beta-alpha; 
    a = (r*beta-alpha)*(1+sqrt(1+(4*(beta+1)*r*(1-r))/(alpha*(beta-alpha))))/(2*r*(1-r));
  }
  b = (r*beta-alpha)*(1+r/alpha)/(r*(1-r));
  x = (a+b)/2;
  double z, oldx = x + 20;
  int niter = 0;
  
  while(tol < std::abs(x-oldx) && niter < maxiter){
    z = g(alpha, beta, x, N);
    oldx = x; 
    x -= (z-r)/(((1-beta/x)*z + alpha/x - pow(z,2)));
    if(x < a || x > b){
      if(z>r){
        b = oldx;
      }
      else{
        a = oldx;
      }
      x= (b+a)/2;
    }
    niter +=1;
  }
  if(change){
    x = -x;
  } 
  return x;
}

double hybridlognewton(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  double x,a,b,rr,w;
  bool change = true;
  if(r < alpha/beta){
    alpha = beta-alpha; 
    r = 1-r;
    change = false;
  }
   
  w =  (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
  rr =  (w*w*alpha+w*(beta-alpha)*alpha)/((beta+1)*(beta-alpha)+w*w*alpha); 
   
  if(r > rr){
    double c3,c4;
    c3 =  beta + (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
    c4 =  (c3*(alpha+beta)-2*alpha*beta)/(beta-alpha);
    a = -(r*beta-alpha)/(r*(1-r))*((c3*(beta*r+alpha)-c4*(beta*r-alpha))/(2*alpha*beta));
    alpha = beta-alpha; 
    r = 1-r;
  } else{
    r = 1-r;
    alpha = beta-alpha; 
    a = (r*beta-alpha)*(1+sqrt(1+(4*(beta+1)*r*(1-r))/(alpha*(beta-alpha))))/(2*r*(1-r));
  }
  b = (r*beta-alpha)*(1+r/alpha)/(r*(1-r));
  x = (a+b)/2;
  r = log(r);
  double z, lz, oldx = x + 20;
  int niter = 0;
  
  while(tol < std::abs(x-oldx) && niter < maxiter){
    z = g(alpha, beta, x, N);
    oldx = x; 
    lz = log(z);
    x -=  (lz-r)/( ((1-beta/x) + alpha/(x*z) - z) ) ;
    if(x < a || x > b){
      if(lz>r){
        b = oldx;
      }
      else{
        a = oldx;
      }
      x= (b+a)/2;
    }
    niter +=1;
  }
  if(change){
    x = -x;
  } 
  return x;
}

double bisection(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  double x,a,b,rr,w;
  bool change = true;
  if(r < alpha/beta){
    alpha = beta-alpha; 
    r = 1-r;
    change = false;
  }
   
  w =  (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
  rr =  (w*w*alpha+w*(beta-alpha)*alpha)/((beta+1)*(beta-alpha)+w*w*alpha); 
   
  if(r > rr){
    double c3,c4;
    c3 =  beta + (sqrt(16*alpha*beta+8*alpha+1)+4*alpha+1)/(8*alpha);
    c4 =  (c3*(alpha+beta)-2*alpha*beta)/(beta-alpha);
    a = -(r*beta-alpha)/(r*(1-r))*((c3*(beta*r+alpha)-c4*(beta*r-alpha))/(2*alpha*beta));
    alpha = beta-alpha; 
    r = 1-r;
  } else{
    r = 1-r;
    alpha = beta-alpha; 
    a = (r*beta-alpha)*(1+sqrt(1+(4*(beta+1)*r*(1-r))/(alpha*(beta-alpha))))/(2*r*(1-r));
  }
  b = (r*beta-alpha)*(1+r/alpha)/(r*(1-r));
  x = (a+b)/2;
  double z, oldx = x + 20;
  int niter = 0;
  
  while(tol < std::abs(x-oldx) && niter < maxiter){
    z = g(alpha, beta, x, N);
    niter +=1;
    if(z<r){
      a = x;
    } else if(z>r){
      b = x;
    }
    else{
      break;
    }
    oldx = x;
    x = (a+b)/2;
  }
  if(change){
    x = -x;
  } 
  return x;
}

double BBG(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  return (beta*r-alpha)/(r*(1-r));
}

double BBG_c(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  return (beta*r-alpha)/(r*(1-r))+r/(2*beta*(1-r));
}

double Sra_2007(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  return (alpha+beta-1)*(1/(1-r) - alpha/(r*(beta-1)));
}
double Sra_2013(double r, double alpha, double beta, int N = 30, double tol = 1e-9, int maxiter = 100) {
  if(r < alpha/(2*beta)){
    return (r*beta-alpha)*(1+r/alpha)/(r*(1-r));
  } else if(r > 2*alpha/sqrt(beta)){
    return (r*beta-alpha)*(1+sqrt(1+(4*(beta+1)*r*(1-r))/(alpha*(beta-alpha))))/(2*r*(1-r));
  } else{
    return (r*beta-alpha)*(1+(r-1)/(beta-alpha))/(r*(1-r));
  }
}
template <class T>
void M_step(const T &data, double (*M_method)(double, double, double, int, double, int),
            const arma::mat &beta_matrix, arma::vec &kappa_vector, arma::mat &mu_matrix,
            arma::rowvec &pi_vector, int K, int N, double reltol, double p, int n, double beta){ 
  K = K - 1;
  arma::rowvec beta_sum = sum(beta_matrix);
  pi_vector = beta_sum/n;
  
  arma::vec eigva, mu1, mup;
  arma::mat eigve, S;
  double kappa1, kappap, v1, v2, lam1, lamp;
  for(int i = 0; i <= K; i++){
    S = trans(data)*diagmat(beta_matrix.col(i))*data/beta_sum(i);
    if(S.has_nan()) {
      kappa_vector(i) = 0;
      return;
    }
    arma::eig_sym(eigva, eigve, S);
    lam1 = eigva(p-1); lamp = eigva(0);
    if(lam1 <= 0){
      stop("group defined by null-matrix in data occurred, please remove rows containing only zeros");
    } else if(lam1 >= 1){
      warning("cluster with only one observation occurred, results can be unprecise");
      v1 =   -std::numeric_limits<double>::infinity();
      kappa1 = 1e16;
      mu1 = eigve.col(p-1);
    } else{
      kappa1 = M_method(lam1, 0.5, beta, N, reltol, 100);
      if(kappa1>=0){
        mu1 = eigve.col(p-1);
      } else{
        mu1 = eigve.col(0);
        lam1 = lamp;
      }
      v1 = kappa1*lam1 - log_hyperg_1F1(0.5, beta, kappa1);
    }
    if(lamp <= 0){
      mu_matrix.col(i) = mu1;
      kappa_vector(i) = kappa1;
      continue; 
    } else if(lamp >= 1){
      stop("group defined by pure idenity-matrix in data occurs, likelihood is infinity, try different methods");
    } else{
      kappap = M_method(lamp, 0.5, beta, N, reltol, 100);
      if(kappap>=0){
        mup = eigve.col(p-1);
        lamp = eigva(p-1);
      } else{
        mup = eigve.col(0);
      }
      v2 = kappap*lamp - log_hyperg_1F1(0.5, beta, kappap);
    }
    if(v1>v2){
      mu_matrix.col(i) = mu1;
      kappa_vector(i) = kappa1;
    } else{
      mu_matrix.col(i) = mup;
      kappa_vector(i) = kappap;
    }
  }
}

template <class T>
double log_like(const T &data, const arma::vec &kappa_vector, const arma::mat &mu_matrix,
                const arma::rowvec &pi_vector, int K, double beta, int n){
  arma::mat A = pow(data*mu_matrix,2);
  A.each_row() %= kappa_vector.t();
  arma::rowvec kummer_vector(K);
  for(int i = 0; i < K; i++) {
    kummer_vector(i) = log_hyperg_1F1(0.5, beta, kappa_vector(i), 10);
  }
  A = arma::repelem(log(pi_vector),n,1) + A - arma::repelem(kummer_vector,n,1);
  arma::vec maxx = max(A, 1);
  maxx += log(sum(exp(A.each_col() - maxx),1));
  return sum(maxx); 
}
// [[Rcpp::export]]
double log_like1(arma::mat &data, const arma::vec &kappa_vector, const arma::mat &mu_matrix,
                 const arma::rowvec &pi_vector, int K, double beta, int n){
  data = normalise(data, 2, 1);
  return log_like(normalise(data, 2, 1), kappa_vector, mu_matrix, pi_vector, K, beta, n); 
}
// [[Rcpp::export]]
double log_like2(arma::sp_mat &data, const arma::vec &kappa_vector, const arma::mat &mu_matrix,
                 const arma::rowvec &pi_vector, int K, double beta, int n){
  data = normalise(data, 2, 1);
  return log_like(data, kappa_vector, mu_matrix, pi_vector, K, beta, n); ;
}

template <class T>
double log_like_hardinit(const T &data, const arma::vec &kappa_vector, 
                         const arma::mat &mu_matrix, int K, double beta, int n, const arma::mat &beta_matrix){
   arma::uvec clus(n);
   arma::mat A;
   arma::vec like;
   double overall = 0;
   for(int i = 0; i < K; i++) {
      clus = arma::conv_to<arma::uvec>::from(beta_matrix.col(i));
      A = extract_rows(data, clus, 1);
      like = kappa_vector(i)*pow(A*mu_matrix.col(i),2) - log_hyperg_1F1(0.5, beta, kappa_vector(i), 10);
      overall += sum(like);
   }
   return overall; 
}

template <class T>
double init(const T &data, arma::mat &beta_matrix, arma::vec &kappa_vector, arma::mat &mu_matrix, arma::rowvec &pi_vector,
            double (*M_method)(double, double, double, int, double, int), List start,
            int K, int N, double reltol, double p, int n, double beta, int maxiter){
  beta_matrix.set_size(n, K);
  mu_matrix.set_size(p, K);
  kappa_vector.set_size(K);
  pi_vector.set_size(K);
  bool given = start["given"];
  int init_iter = start["init_iter"];
  if(given){
    beta_matrix = as<arma::mat>(start["matrix"]);
    if(init_iter>0){
      diamclus_internal(/*c-r*/data, /*r*/beta_matrix,/*r*/mu_matrix, K, n, init_iter);
    }
    M_step(/*c-r*/data, M_method,/*c-r*/beta_matrix,/*r*/kappa_vector,/*r*/mu_matrix,/*r*/pi_vector,
           /*r*/K, N, reltol, p, n, beta);
  } else{
    if(init_iter>0){
      beta_matrix.zeros();
      diamclus_internal(/*c-r*/data, /*r*/beta_matrix,/*r*/mu_matrix, K, n, init_iter);
      M_step(/*c-r*/data, M_method,/*c-r*/beta_matrix,/*r*/kappa_vector,/*r*/mu_matrix,/*r*/pi_vector,
             /*r*/K, N, reltol, p, n, beta);
    } else{
      beta_matrix = normalise(beta_matrix.randu(),1,1);
      mu_matrix = normalise(mu_matrix.randn(),2,0);
      kappa_vector.randn();
      pi_vector = sum(beta_matrix)/n;
    }
  }
  if(maxiter == 0){
    return log_like_hardinit(data,/*c-r*/kappa_vector,/*c-r*/mu_matrix, K, beta, n, beta_matrix);
  }else{
    return -1e10;
  }
  
}
template <class T>
void watson(const T &data, int K, void (*E_method)(arma::mat&, int, int), double (*M_method)(double, double, double, int, double, int),
            arma::mat &beta_matrix, arma::vec &kappa_vector, arma::mat &mu_matrix, arma::rowvec &pi_vector,
            double minalpha, bool convergence, int maxiter, int N, double reltol, double p, int n, double beta, double& log_lik,
            arma::mat &max_beta_matrix, arma::vec &max_kappa_vector, arma::mat &max_mu_matrix, arma::rowvec &max_pi_vector){
  
  int i = 0;
  double max_log_lik = -1e11;
  bool stop;
  while(i<maxiter){
    stop = E_step(/*c-r*/data,/*r*/beta_matrix,/*c-r*/kappa_vector,/*c-r*/mu_matrix,/*c-r*/ pi_vector, E_method, K, convergence,
                  minalpha, beta, n, p, log_lik, reltol, max_beta_matrix, max_kappa_vector, max_mu_matrix, max_pi_vector, max_log_lik);
    if(stop) break;
    M_step(/*c-r*/data, M_method,/*c-r*/beta_matrix,/*r*/kappa_vector,/*r*/mu_matrix,/*r*/pi_vector,/*r*/K, N, reltol, p, n, beta);
    i += 1;
  }
  if(!convergence){
    beta_matrix  = max_beta_matrix;
    mu_matrix    = max_mu_matrix;
    pi_vector    = max_pi_vector;
    kappa_vector = max_kappa_vector;
    log_lik      = max_log_lik;
  }
  return;/*List::create(beta_matrix, kappa_vector, mu_matrix, pi_vector, log_lik)*/
}


template <class T>
List EM(T &data, int K, String E_type, String M_type,  double minalpha,
        bool convergence, int maxiter, int N, double reltol, List &start, bool verbose = false){
  data = normalise(data, 2, 1); // normalize in case
  double p = data.n_cols;
  int n = data.n_rows;
  double beta  = p/2;
  int N_runs = start.size();
  List best_result = List(1);
  
  double (*M_method)(double, double, double, int, double, int);
  void (*E_method)(arma::mat&, int, int);
  
  if(E_type == "softmax"){
    E_method = soft;
  } else if(E_type == "hardmax"){
    E_method = hard;
  } else{ // stochmax
    E_method = stoch;
  }
  
  if(M_type == "newton"){
    M_method = hybridnewton;
  } else if(M_type == "lognewton"){
    M_method = hybridlognewton;
  } else if(M_type == "bisection"){
    M_method = bisection;
  } else if(M_type == "BBG"){
    M_method = BBG;
  } else if(M_type == "Sra_2007"){
    M_method = Sra_2007;
  } else if(M_type == "Sra_Karp_2013"){
    M_method = Sra_2013;
  } else{ // "BBG-c"
    M_method = BBG_c;
  }
  
  arma::mat beta_matrix(n, K), max_beta_matrix(n, K);
  arma::mat mu_matrix(p, K), max_mu_matrix(p, K);
  arma::vec kappa_vector(K), max_kappa_vector(K);
  arma::rowvec pi_vector(K), max_pi_vector(K);
  double log_lik, best_log_lik = -1e11;
  for(int i=0; i<N_runs; i++){
    if(verbose) Rcout << "Run: " << i+1 << "/" << N_runs << std::endl; 
    log_lik = init(data, beta_matrix, kappa_vector, mu_matrix, pi_vector, M_method, start[i], K, N, reltol, p, n, beta, maxiter);
    watson(data, K, E_method, M_method, beta_matrix, kappa_vector, mu_matrix, pi_vector, minalpha, convergence, 
           maxiter, N, reltol, p, n, beta, log_lik, max_beta_matrix, max_kappa_vector, max_mu_matrix, max_pi_vector);
    if(log_lik>best_log_lik){
      best_log_lik = log_lik;
      best_result = List::create(beta_matrix, kappa_vector.t(), mu_matrix, pi_vector, log_lik);
    }
    if(i % 5 == 0) Rcpp::checkUserInterrupt();
  }
  
  return best_result;
}
// [[Rcpp::export]]
List EM1(arma::mat &data, int K, String E_type, String M_type,  double minalpha = 0,
         bool convergence = true, int maxiter = 100, int N = 30, double reltol = 1e-9, List start = R_NilValue, bool verbose = false){
  return EM(data, K, E_type, M_type,  minalpha, convergence , maxiter,  N, reltol, start, verbose);
}
// [[Rcpp::export]]
List EM2(arma::sp_mat &data, int K, String E_type, String M_type,  double minalpha = 0,
         bool convergence = true, int maxiter = 100, int N = 30, double reltol = 1e-9, List start = R_NilValue, bool verbose = false){
  return EM(data, K, E_type, M_type,  minalpha, convergence , maxiter,  N, reltol, start, verbose);
}