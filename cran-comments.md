## Test environments
* local Ubuntu 24.04 LTS, R 4.5.2
* win-builder 
* GitHub Actions (ubuntu-24.04, macOS-latest, windows-latest) https://github.com/lsablica/watson/actions

## R CMD check results
There were no ERRORs or WARNINGs.

The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.

It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is above the 1MB threshold.
My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp and RcppArmadillo. Without the speed up gained from those C++ functions, this package would become impractical.

## Downstream dependencies
none


