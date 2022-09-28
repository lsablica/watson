## Test environments
* local Ubuntu 18.04 LTS, R 4.1.2
* win-builder 
* R-hub

## R CMD check results
There were no ERRORs or WARNINGs.

It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is above the 1MB threshold.
My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp and RcppArmadillo. Without the speed up gained from those C++ functions, this package would become impractical.

## Downstream dependencies
none



