/* The following code is based on the source code of GNU GSL */

/* specfunc/gamma.c, specfunc/hyperg_1F1.c, specfunc/hyperg.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2010 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#include <RcppArmadillo.h>
#include "hyper.h"

const double GSL_DBL_EPSILON = 2.2204460492503131e-16;
const double GSL_DBL_MAX = 1.7976931348623157e+308;
const double SUM_LARGE = (1.0e-5*GSL_DBL_MAX);
const double GSL_LOG_DBL_MAX = 7.0978271289338397e+02;
const double GSL_LOG_DBL_MIN = (-7.0839641853226408e+02);
const double GSL_SQRT_DBL_MAX = 1.3407807929942596e+154;
const double GSL_SQRT_DBL_MIN = 1.4916681462400413e-154;
const double LogRootTwoPi_ = 0.9189385332046727418;
double MPI = 3.14159265358979323846264338328;
double ME = 2.71828182845904523536028747135;
const double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
  -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
  -176.61502916214059906584551354,
  12.507343278686904814458936853,
  -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};


int gsl_sf_hyperg_1F1_series_e(const double a, const double b, const double x, double &result){
  double an  = a;
  double bn  = b;
  double n   = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double sum_val = 1.0;
  
  while(abs_del/std::abs(sum_val) > 0.25*GSL_DBL_EPSILON) { 
    double u, abs_u;
    
    if (n > 10000.0) {
      result = sum_val;
      return 1; /*not cvg*/
    }
    u = x * (an/(bn*n));
    abs_u = std::abs(u);
    if(abs_u > 1.0 && max_abs_del > GSL_DBL_MAX/abs_u) {
      result = sum_val;
      return 1; /*overflow*/
    }
    del *= u;
    sum_val += del;
    if(std::abs(sum_val) > SUM_LARGE) {
      result = sum_val;
      return 1; /*overflow*/
    }
    
    abs_del = std::abs(del);
    if(abs_del > max_abs_del) max_abs_del = abs_del;
    an += 1.0;
    bn += 1.0;
    n  += 1.0;
  }
  
  result = sum_val;
  return 0;
}



int lngamma_lanczos(double x, double &result){
    int k;
    double Ag, term1, term2;
    
    x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */
    
    Ag = lanczos_7_c[0];
    for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }
    
    /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
    term1 = (x+0.5)*std::log((x+7.5)/ME);
    term2 = LogRootTwoPi_ + std::log(Ag);
    result = term1 + (term2 - 7.0);
    return 0;
  }

int gsl_sf_lngamma_e(double x, double &result_lg){
  return lngamma_lanczos(x, result_lg);
}

int gsl_sf_hyperg_2F0_series_e(const double a, const double b, const double x, int n_trunc, double &result){
  const int maxiter = 2000;
  double an = a;
  double bn = b;
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double last_abs_del = 1.0;
  
  while(abs_del/std::abs(sum) > GSL_DBL_EPSILON && n < maxiter) {
    
    double u = an * (bn/n * x);
    double abs_u = std::abs(u);
    
    if(abs_u > 1.0 && (max_abs_del > GSL_DBL_MAX/abs_u)) {
      result = sum;
      return 1; /*overflow*/
    }
    del *= u;
    sum += del;
    abs_del = std::abs(del);
    if(abs_del > last_abs_del) break; /* series is probably starting to grow */
    
    last_abs_del = abs_del;
    if(abs_del > max_abs_del) max_abs_del = abs_del;
    
    an += 1.0;
    bn += 1.0;
    n  += 1.0;
    
    if(an == 0.0 || bn == 0.0) break;        /* series terminated */
    
    if(n_trunc >= 0 && n >= n_trunc) break;  /* reached requested timeout */
  }
  
  result = sum;
  if(n >= maxiter){
    return 1; /*not cvg*/;
  } else{
    return 0;
  }
} 

int gsl_sf_exp_mult_err_e(const double x, const double y, double &result){
  const double ay  = std::abs(y);
  
  if(y == 0.0) {
    result = 0.0;
    return 0;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)&& (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)) {
    double ex = std::exp(x);
    result = y * ex;
    return 0;
  }
  else {
    const double ly  = std::log(ay);
    const double lnr = x + ly;
    
    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      return 1; /*overflow*/
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      return 1; /*underflow*/
    }
    else {
      const double sy  = y >= 0 ? 1 : -1;
      const double M   = std::floor(x);
      const double N   = std::floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      const double eMN = std::exp(M+N);
      const double eab = std::exp(a+b);
      result = sy * eMN * eab;
      return 0;
    }
  }
}

int hyperg_1F1_asymp_posx(const double a, const double b, const double x, double &result){
  double lg_b, lg_a;
  
  int stat_b = gsl_sf_lngamma_e(b, lg_b);
  int stat_a = gsl_sf_lngamma_e(a, lg_a);
  
  if(stat_a + stat_b == 0) {
    double F;
    int stat_F = gsl_sf_hyperg_2F0_series_e(b-a, 1.0-a, 1.0/x, -1, F);
    if(stat_F == 0 && F != 0) {
      double lnx = std::log(x);
      double ln_term_val = (a-b)*lnx;
      double ln_pre_val = lg_b - lg_a + ln_term_val + x;
      int stat_e = gsl_sf_exp_mult_err_e(ln_pre_val, F, result);
      return stat_e;
    }
    else {
      result = 0;
      return stat_F;
    }
  }
  else {
    return 1;
  }
}

int hyperg_1F1_largebx(const double a, const double b, const double x, double &result){
  double y = x/b;
  double f = std::exp(-a*std::log1p(-y));
  double t1 = -((a*(a+1.0))/(2*b))*std::pow((y/(1.0-y)),2.0);
  double t2 = (1/(24*b*b))*((a*(a+1)*y*y)/std::pow(1-y,4))*(12+8*(2*a+1)*y+(3*a*a-a-2)*y*y);
  double t3 = (-1/(48*b*b*b*std::pow(1-y,6)))*a*((a + 1)*((y*((a + 1)*(a*(y*(y*((y*(a - 2) + 16)*(a - 1)) + 72)) + 96)) + 24)*std::pow(y, 2)));
  result = f * (1 + t1 + t2 + t3);
  return 0;
}


int hyperg_1F1_large2bm4a(const double a, const double b, const double x, double &result){
  double eta    = 2.0*b - 4.0*a;
  double cos2th = x/eta;
  double sin2th = 1.0 - cos2th;
  double th = std::acos(std::sqrt(cos2th));
  double pre_h  = 0.25*MPI*MPI*eta*eta*cos2th*sin2th;
  double lg_b;
  int stat_lg = gsl_sf_lngamma_e(b, lg_b);
  double t1 = 0.5*(1.0-b)*std::log(0.25*x*eta);
  double t2 = 0.25*std::log(pre_h);
  double lnpre_val = lg_b + 0.5*x + t1 - t2;
  double s1 = std::sin(a*MPI);
  double s2 = std::sin(0.25*eta*(2.0*th - std::sin(2.0*th)) + 0.25*MPI);
  double ser_val = s1 + s2;
  int stat_e = gsl_sf_exp_mult_err_e(lnpre_val, ser_val,  result);
  return stat_e + stat_lg;
}

int hyperg_1F1_asymp_negx(const double a, const double b, const double x, double &result){
  double lg_b, lg_bma;
  
  int stat_b   = gsl_sf_lngamma_e(b, lg_b);
  int stat_bma = gsl_sf_lngamma_e(b-a, lg_bma);
  
  if(stat_b + stat_bma == 0) {
    double F;
    int stat_F = gsl_sf_hyperg_2F0_series_e(a, 1.0+a-b, -1.0/x, -1, F);
    if(F != 0) {
      double ln_term_val = a*std::log(-x);
      double ln_pre_val = lg_b - lg_bma- ln_term_val;
      int stat_e = gsl_sf_exp_mult_err_e(ln_pre_val, F, result);
      return stat_e;
    }
    else {
      result = 0.0;
      return stat_F;
    }
  }
  else {
    return 1;
  }
}

int hyperg_1F1_luke(const double a, const double c, const double xin, double &result){
    const double RECUR_BIG = 1.0e+50;
    const int nmax = 5000;
    int n = 3;
    const double x  = -xin;
    const double x3 = x*x*x;
    const double t0 = a/c;
    const double t1 = (a+1.0)/(2.0*c);
    const double t2 = (a+2.0)/(2.0*(c+1.0));
    double F = 1.0;
    double prec;
    
    double Bnm3 = 1.0;                                  /* B0 */
    double Bnm2 = 1.0 + t1 * x;                         /* B1 */
    double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
    
    double Anm3 = 1.0;                                                      /* A0 */
    double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
    double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */
    
    while(1) {
      double npam1 = n + a - 1;
      double npcm1 = n + c - 1;
      double npam2 = n + a - 2;
      double npcm2 = n + c - 2;
      double tnm1  = 2*n - 1;
      double tnm3  = 2*n - 3;
      double tnm5  = 2*n - 5;
      double F1 =  (n-a-2) / (2*tnm3*npcm1);
      double F2 =  (n+a)*npam1 / (4*tnm1*tnm3*npcm2*npcm1);
      double F3 = -npam2*npam1*(n-a-2) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
      double E  = -npam1*(n-c-1) / (2*tnm3*npcm2*npcm1);
      
      double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
      double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
      double r = An/Bn;
      
      prec = std::abs((F - r)/F);
      F = r;
      
      if(prec < GSL_DBL_EPSILON || n > nmax) break;
      
      if(std::abs(An) > RECUR_BIG || std::abs(Bn) > RECUR_BIG) {
        An   /= RECUR_BIG;
        Bn   /= RECUR_BIG;
        Anm1 /= RECUR_BIG;
        Bnm1 /= RECUR_BIG;
        Anm2 /= RECUR_BIG;
        Bnm2 /= RECUR_BIG;
        Anm3 /= RECUR_BIG;
        Bnm3 /= RECUR_BIG;
      }
      else if(std::abs(An) < 1.0/RECUR_BIG || std::abs(Bn) < 1.0/RECUR_BIG) {
        An   *= RECUR_BIG;
        Bn   *= RECUR_BIG;
        Anm1 *= RECUR_BIG;
        Bnm1 *= RECUR_BIG;
        Anm2 *= RECUR_BIG;
        Bnm2 *= RECUR_BIG;
        Anm3 *= RECUR_BIG;
        Bnm3 *= RECUR_BIG;
      }
      
      n++;
      Bnm3 = Bnm2;
      Bnm2 = Bnm1;
      Bnm1 = Bn;
      Anm3 = Anm2;
      Anm2 = Anm1;
      Anm1 = An;
    }
    
    result = F;
    return 0;
  }

int hyperg_1F1_small_a_bgt0(const double a, const double b, const double x, double &result){
  const double bma = b-a;
  const double oma = 1.0-a;
  const double ap1mb = 1.0+a-b;
  const double abs_bma = std::abs(bma);
  const double abs_oma = std::abs(oma);
  const double abs_ap1mb = std::abs(ap1mb);
  const double ax = std::abs(x);
  
  if(b >= 1.4*ax) {
    return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
  }
  else if(x > 0.0) {
    if(x > 100.0 && abs_bma*abs_oma < 0.5*x) {
      return hyperg_1F1_asymp_posx(a, b, x, result);
    }
    else if(b < 5.0e+06) {
      /* Recurse backward on b from
      * a suitably high point.
      */
      const double b_del = std::ceil(1.4*x-b) + 1.0;
      double bp = b + b_del;
      double r_Mbp1;
      double r_Mb;
      double Mbp1;
      double Mb;
      double Mbm1;
      int stat_0 = gsl_sf_hyperg_1F1_series_e(a, bp+1.0, x, r_Mbp1);
      int stat_1 = gsl_sf_hyperg_1F1_series_e(a, bp,     x, r_Mb);
      Mbp1 = r_Mbp1;
      Mb   = r_Mb;
      while(bp > b+0.1) {
        /* Do backward recursion. */
        Mbm1 = ((x+bp-1.0)*Mb - x*(bp-a)/bp*Mbp1)/(bp-1.0);
        bp -= 1.0;
        Mbp1 = Mb;
        Mb   = Mbm1;
      }
      result = Mb;
      return stat_0+stat_1;
    }
    else if (std::abs(x) < std::abs(b) && std::abs(a*x) < sqrt(std::abs(b)) * std::abs(b-x)) {
      return hyperg_1F1_largebx(a, b, x, result);
    } else {
      return hyperg_1F1_large2bm4a(a, b, x, result);
    }
  }
  else {
    /* x < 0 and b not large compared to |x|
    */
    if(ax < 10.0 && b < 10.0) {
      return gsl_sf_hyperg_1F1_series_e(a, b, x, result);
    }
    else if(ax >= 100.0 && (abs_ap1mb > 1.0 ? abs_ap1mb : 1.0) < 0.99*ax) {
      return hyperg_1F1_asymp_negx(a, b, x, result);
    }
    else {
      return hyperg_1F1_luke(a, b, x, result);
    }
  }
}

int gsl_sf_exp_e(const double x, double &result){
  if(x > GSL_LOG_DBL_MAX) {
    return 1;
  } else if(x < GSL_LOG_DBL_MIN) {
    return 1;
  } else {
    result = std::exp(x);
    return 0;
  }
}


int gsl_sf_hyperg_1F1_e(const double a, const double b, const double x, double &result){
  if(x == 0.0) {
    result = 1.0;
    return 0;
  } else if(a == b) {
    return gsl_sf_exp_e(x, result);
  } else {
    return hyperg_1F1_small_a_bgt0(a, b, x, result);
  }
}

// // [[Rcpp::export]]
// Rcpp::List h(const double a, const double b, const double x){
//   double result;
//   int r = gsl_sf_hyperg_1F1_e(a, b, x, result);
//   return Rcpp::List::create(r, result);
// }
