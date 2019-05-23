############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Calculate KLD utility
## Author: Pubudu Thialn
#############################################################################

kld.post <- function(X,Y,inits){
  Sol<-LP_approx(X,Y,inits) 
  mu_post <- Sol$par
  Hes_mat=Sol$hessian #calculate Hesian matrix
  
  if(is.non.singular.matrix(Hes_mat)==TRUE){
    if(is.positive.definite(Hes_mat)==FALSE){
      Hes_mat=make.positive.definite(Hes_mat,tol=1e-10)
    }
    
    Sigma_post<-solve(Hes_mat,tol=1e-40)
    det.out <- 0.5*(trace(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log(dSigma_prior/det(Sigma_post)))#Kld utility value
  }
  else
  {
    det.out=0.00000001
  }
  if(is.finite(det.out)==TRUE)
  {
    det.out=det.out
  }
  else
  {
    det.out=0.00000001
  }
  return(det.out)
}