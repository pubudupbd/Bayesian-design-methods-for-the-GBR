############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Calculate KLD utility
## Author: Pubudu Thialn
#############################################################################

kld.post <- function(X,Y,inits,mu_prior,Sigma_prior,Sigma_priorb,iSigma_prior,dSigma_prior,v0,leny,euclidean_nextY2,num_par,L,r_rand2){
  Sol<-LP_approx(X,Y,inits,mu_prior,Sigma_prior,Sigma_priorb,v0,leny,euclidean_nextY2,num_par,L,r_rand2)
  mu_post <- Sol$par
  Hes_mat=Sol$hessian #hessian matrix
    
  Sigma_post<-solve(Hes_mat,tol=1e-40)
  det.out <- 0.5*(trace(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log(dSigma_prior/det(Sigma_post))) #formular for KLD
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