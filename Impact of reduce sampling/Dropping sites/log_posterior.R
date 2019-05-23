############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Calculate log posterior
## Author: Pubudu Thialn
#############################################################################
log_posterior <- function(par,Y,X,mu_prior,Sigma_priorb,v0,leny,euclidean_nextY2,num_par,L,r_rand2) {
  tcov2<-v0*diag(rep(1,leny))+exp(par[2])*exp(-euclidean_nextY2/(exp(par[3]))^2) #Guassian covariance 
  zb<-(chol(tcov2))%*%r_rand2
  
  lp_tem<-X%*%par[4:num_par]
  lp<-rep(lp_tem,L)+zb

    m1=exp(lp)/(1+exp(lp))
    m2=1/exp(par[1])  #change
    log_like=mean(colSums(dbeta(Y, m1*m2,(1-m1)*m2, log = T))) #log likelihood
    log_prior<-sum(dnorm(x=par,mean=mu_prior,sd=sqrt(Sigma_priorb),log=TRUE)) #log prior
    log_post<- log_like+log_prior #log posterior

  if(is.finite(log_post)==TRUE) {log_post=log_post}
  else{log_post=10000000}
 -1*log_post
}
