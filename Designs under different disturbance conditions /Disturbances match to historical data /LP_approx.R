############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Use optim function to find minimum of log posterior function
## Author: Pubudu Thialn

LP_approx <- function(X,Y,inits,mu_prior,Sigma_prior,Sigma_priorb,v0,leny,euclidean_nextY2,num_par,L,r_rand2){
  lb=c(mu_prior[1]-3*sqrt(Sigma_priorb[1]),mu_prior[2]-3*sqrt(Sigma_priorb[2]),mu_prior[3]-3*sqrt(Sigma_priorb[3]),mu_prior[4]-3*sqrt(Sigma_priorb[4]),mu_prior[5]-3*sqrt(Sigma_priorb[5]),mu_prior[6]-3*sqrt(Sigma_priorb[6]),mu_prior[7]-3*sqrt(Sigma_priorb[7]),mu_prior[8]-3*sqrt(Sigma_priorb[8]),mu_prior[9]-3*sqrt(Sigma_priorb[9]),mu_prior[10]-3*sqrt(Sigma_priorb[10]),mu_prior[11]-3*sqrt(Sigma_priorb[11]),mu_prior[12]-3*sqrt(Sigma_priorb[12]),mu_prior[13]-3*sqrt(Sigma_priorb[13]),mu_prior[14]-3*sqrt(Sigma_priorb[14])) #lower bound for optim
  ub=c(mu_prior[1]+3*sqrt(Sigma_priorb[1]),mu_prior[2]+3*sqrt(Sigma_priorb[2]),mu_prior[3]+3*sqrt(Sigma_priorb[3]),mu_prior[4]+3*sqrt(Sigma_priorb[4]),mu_prior[5]+3*sqrt(Sigma_priorb[5]),mu_prior[6]+3*sqrt(Sigma_priorb[6]),mu_prior[7]+3*sqrt(Sigma_priorb[7]),mu_prior[8]+3*sqrt(Sigma_priorb[8]),mu_prior[9]+3*sqrt(Sigma_priorb[9]),mu_prior[10]+3*sqrt(Sigma_priorb[10]),mu_prior[11]+3*sqrt(Sigma_priorb[11]),mu_prior[12]+3*sqrt(Sigma_priorb[12]),mu_prior[13]+3*sqrt(Sigma_priorb[13]),mu_prior[14]+3*sqrt(Sigma_priorb[14])) # upper bound for optim
  lp.approx <- optim(par=inits,Y=Y,X=X,L=L,num_par=num_par,v0=v0,leny=leny,mu_prior=mu_prior,Sigma_priorb=Sigma_priorb,euclidean_nextY2=euclidean_nextY2,r_rand2=r_rand2,fn=log_posterior, hessian=TRUE,method="L-BFGS-B",lower=lb,upper=ub)
  lp.approx
}