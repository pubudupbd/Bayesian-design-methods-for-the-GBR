############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Calculate trace for KLD calculation
## Author: Pubudu Thialn
trace <- function(M){
  tr <- sum(diag(M))
  tr
}
