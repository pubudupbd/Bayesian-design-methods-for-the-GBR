############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Scenario 1, environmental disturbances were simulated to match the historical data in the Whitsunday region  
## Author: Pubudu Thialn
#############################################################################

library("foreach",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("inline",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("iterators",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("parallel")
library("doParallel",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("matrixcalc",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")#is.singular
library("Matrix",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("mvtnorm",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("sp",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")
library("corpcor",lib.loc="/home/n9592385/Year_2/paper1_code/intel_pack")


set.seed(20)

source('log_posterior.R')
source('LP_approx.R')
source('kld.post.R')
source('trace.R')
source('covariate_generate.R')
source('euc_dist_cal.R')

B<<-100#T posterior distributions need to be sampled from or approximated
leny <<- 135 # number of data points for next sampling year
L<<-100 #number of samples for Monte carlo intergration
v0<<-0.00001 #nugect

#Prior for the design
mu_prior <<- c(-2.5232,-5.9825,-1.1189,-1.2726,0.153,0.9066,0.2782,-0.1062,-0.7962,-0.23,-0.0093,-0.449,-0.2152,-0.0431)
Sigma_prior <<- diag(c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2))
Sigma_priorb <<- c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2)

#simulate parameters
iSigma_prior <<- solve(Sigma_prior)
dSigma_prior <<- det(Sigma_prior)
num_par <<- length(mu_prior)

theta_sim <- matrix(0,nrow=B,ncol=length(mu_prior))
for (j in 1:length(mu_prior)){
  theta_sim[,j] <-rnorm(B,mu_prior[j],sqrt(Sigma_priorb[j]))
}

theta_sim <- data.frame(theta_sim)
theta_sim0 <- as.data.frame(sweep(theta_sim,2,colMeans(theta_sim),"-"))  # shift the mean 
C <- chol(nearPD(cov(theta_sim0))$mat)
SN <- matrix(rnorm(B * ncol(theta_sim)), ncol(theta_sim))#standard normals
X <- t(C) %*% SN
theta_mul <- data.frame(as.matrix(t(X)))
names(theta_mul) <- names(theta_sim)
theta_mul <- as.data.frame(sweep(theta_mul,2,colMeans(theta_sim),"+")) 
theta_sim_all<-as.matrix(theta_mul)
r_rand1<- matrix(rnorm(leny*B,0,1),leny,B)

# KLD utility claculation
kld.utility<-function(d,B,theta_sim,r_rand1){
  kld.post.val=matrix(0,B)
  for (j in 1:B){
    r_rand2<- matrix(rnorm(leny*L,0,1),leny,L)
    X<-cbind(1,Get_Covariate(d,20)) # covariates
    theta = theta_sim[,(4:num_par)][j,]
    euclidean_nextY2=euc_dist_cal(d)
    cov<-v0*diag(rep(1,leny))+exp(theta_sim[,2][j])*exp(-euclidean_nextY2/(exp(theta_sim[,3][j]))^2)
    zb <-(chol(cov))%*%r_rand1[,j] # spatial random effects
    
    lp = X%*%theta+zb
    mu=exp(lp)/(1+exp(lp))
    phi=1/exp(theta_sim[,1][j])
    Y =rbeta(leny,mu*phi,(1-mu)*phi)
    Y=ifelse(1-Y< 1e-5,0.9999990,Y)
    Y=ifelse(Y< 1e-5,0.000001,Y) #coral cover data
    
    inits=theta_sim[j,]
    kld_out =kld.post(X,Y,inits,mu_prior,Sigma_prior,Sigma_priorb,iSigma_prior,dSigma_prior,v0,leny,euclidean_nextY2,num_par,L,r_rand2)
    kld.post.val[j]=kld_out
  }
  return(kld.post.val)
}

#cycle through sites
max_site_cal<-function(n_ind,i,curr_x,curr_u,site,kld.utility,B,theta_sim,r_rand1)
{
  prop_x<-curr_x
  max_site<-foreach(j=1:27,.combine='cbind',.export = c("B","n_ind","i","curr_x","curr_u","prop_x","site","leny","L","v0","theta_sim","r_rand1","num_par","mu_prior","Sigma_prior","Sigma_priorb","iSigma_prior","dSigma_prior"),.packages = c("matrixcalc","Matrix","mvtnorm","sp"))%dopar%
  {
    tem_out=c(mean(kld.utility(replace(prop_x,i,site[j]),B,theta_sim,r_rand1)),site[j])
    tem_out
  }
  max_site=data.frame(max_site)
  print(t(max_site))
  MS=max(data.frame(max_site)[1,])
  if(MS>curr_u){
    prop_x[i]=which(abs((max_site[1,]-MS))<1e-5)
    curr_x <- prop_x
    curr_u <- MS
  }
  return(list(curr_x,curr_u))
}
#######################################################################
# Coordinate exchange algorithm
optimal_des<-function(n_ind,kld.utility,curr_x,site,B,theta_sim,r_rand1)
{
  curr_u <- mean(kld.utility(curr_x,B,theta_sim,r_rand1))
  for(k in 1:5){
    for(i in 1:27){
      print(c(k,i))
      out<-max_site_cal(n_ind,i,curr_x,curr_u,site,kld.utility,B,theta_sim,r_rand1)
      curr_x<-out[[1]]
      curr_u<-out[[2]]
      print(paste0("Current ud: ", curr_u))
      print(matrix(c(round(curr_x,0)),byrow=T,nrow=1,ncol=27))
    }
    curr_x <- curr_x
    curr_u <- curr_u
  }
  return(list(curr_x,curr_u))
}
#######################################################################
site<-seq(1,27,1)

E=5#independent runs
rsite<-sample(site, 27*E, replace = T, prob = NULL)
index <- Sys.getenv("PBS_ARRAY_INDEX")
n_ind <- as.numeric(index)

cl <- makeCluster(27,outfile=paste("Independent_Run",n_ind,".txt",sep=""))
clusterEvalQ(cl,.libPaths("/home/n9555293/pubudu/intel_pack"))
registerDoParallel(cl)
getDoParWorkers()

clusterCall(cl, function() { source("log_posterior.R") })
clusterCall(cl, function() { source("LP_approx.R") })
clusterCall(cl, function() { source("kld.post.R") })
clusterCall(cl, function() { source("trace.R") })
clusterCall(cl, function() { source("covariate_generate.R") })
clusterCall(cl, function() { source("euc_dist_cal.R") })

#   # starting from 1 , n_ind <- 1
index <- (1+(n_ind-1)*27):(n_ind*27)
index <- index[index<=E*27]
rstart.d<-rsite[index]

v_ut <- "B100"
v_ut <- paste(v_ut,n_ind,sep="")

optimal_design<-optimal_des(n_ind,kld.utility,rstart.d,site,B,theta_sim_all,r_rand1)

v_fileName <- paste(v_ut,"RData",sep=".")
save(list= c("optimal_design"),file = v_fileName)

stopCluster(cl)

