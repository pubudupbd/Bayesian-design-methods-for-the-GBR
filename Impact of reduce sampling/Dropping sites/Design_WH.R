############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: The effect of visiting fewer LTMP sites by dropping reefs  in a given region  
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

B<<-100 #T posterior distributions need to be sampled from or approximated
leny <<- 90 # number of data points for next sampling year
L<<-100 #number of samples for Monte carlo intergration
v0<<-0.00001 #nugect

#prior for the design
mu_prior <<- c(-2.5232,-5.9825,-1.1189,-1.2726,0.153,0.9066,0.2782,-0.1062,-0.7962,-0.23,-0.0093,-0.449,-0.2152,-0.0431)
Sigma_prior <<- diag(c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2))
Sigma_priorb <<- c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2)
iSigma_prior <<- solve(Sigma_prior)
dSigma_prior <<- det(Sigma_prior)
num_par <<- length(mu_prior)

#Simulate parameters
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

kld.utility<-function(d,B,theta_sim,r_rand1){
  kld.post.val=matrix(0,B)
  for (j in 1:B){
    r_rand2<- matrix(rnorm(leny*L,0,1),leny,L)
    X<-cbind(1,Get_Covariate(d,20)) # simulate covariates
    theta = theta_sim[,(4:num_par)][j,]
    euclidean_nextY2=euc_dist_cal(d)
    cov<-v0*diag(rep(1,leny))+exp(theta_sim[,2][j])*exp(-euclidean_nextY2/(exp(theta_sim[,3][j]))^2)
    zb <-(chol(cov))%*%r_rand1[,j] #spatial random effects
    
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
  return(kld.post.val) #kld utility value
}

#circle through sites
max_site_cal<-function(new_design,i,curr_x,curr_u,site,kld.utility,B,theta_sim,r_rand1)
{
  prop_x<-curr_x
  max_site<-foreach(jf=seq(1,5,2),.combine='cbind',.export = c("new_design","B","i","curr_x","curr_u","prop_x","site","leny","L","v0","theta_sim","r_rand1","num_par","mu_prior","Sigma_prior","Sigma_priorb","iSigma_prior","dSigma_prior"),.packages = c("matrixcalc","Matrix","mvtnorm","sp"))%dopar%
  {
    tem_out=c(mean(kld.utility(replace(prop_x,c(2*i-1,2*i),as.numeric(as.vector(paste(new_design[jf:(jf+1),i])))),B,theta_sim,r_rand1)),jf)
    tem_out
  }
  
  max_site=data.frame(max_site)
  MS=max(data.frame(max_site)[1,])
  if(MS>curr_u){
    sID<-seq(1,5,2)
    keep_site=which(abs((max_site[1,]-MS))<1e-5)
    print(keep_site)
    curr_x <- replace(prop_x,c(2*i-1,2*i),as.numeric(as.vector(paste(new_design[sID[keep_site]:(sID[keep_site]+1),i]))))
    curr_u <- MS
  }
  return(list(curr_x,curr_u))
}
#######################################################################

# Coordinate exchange algorithm
optimal_des<-function(new_design,rstart.d_id,kld.utility,site,B,theta_sim,r_rand1)
{
  random_dis<-list()
  for(j in 1:9)
  {
    random_dis[[j]]<-new_design[rstart.d_id[j]:(rstart.d_id[j]+1),j]
  }
  rstart.d<-matrix(do.call(cbind,random_dis),ncol = 1,nrow=18)
  curr_x<-as.numeric(as.vector(paste(rstart.d)))

  curr_u <- mean(kld.utility(curr_x,B,theta_sim,r_rand1))
  for(k in 1:30){
    for(i in 1:9){
      print(c(k,i))
      out<-max_site_cal(new_design,i,curr_x,curr_u,site,kld.utility,B,theta_sim,r_rand1)
      curr_x<-out[[1]]
      curr_u<-out[[2]]
      print(paste0("Current ud: ", curr_u))
      print(matrix(c(round(curr_x,0)),byrow=T,nrow=1,ncol=18))
    }
    curr_x <- curr_x
    curr_u <- curr_u
  }
  return(list(curr_x,curr_u))
}
#######################################################################
site<-seq(1,27,1)
E=10#independent runs
indexv<-seq(1,5,2)
rsite<-sample(indexv,9*E, prob = NULL,replace = T)
index <- Sys.getenv("PBS_ARRAY_INDEX")
n_ind <- as.numeric(index)

#   # starting from 1 , n_ind <- 1
index <- (1+(n_ind-1)*9):(n_ind*9)
index <- index[index<=E*9]
rstart.d_id<-rsite[index]

v_ut <- "twosite"
v_ut <- paste(v_ut,n_ind,sep="")

cl <- makeCluster(3,outfile=paste("Independent_Run",n_ind,".txt",sep=""))
clusterEvalQ(cl,.libPaths("/home/n9592385/Year_2/paper1_code/intel_pack"))
registerDoParallel(cl)
getDoParWorkers()

clusterCall(cl, function() { source("log_posterior.R") })
clusterCall(cl, function() { source("LP_approx.R") })
clusterCall(cl, function() { source("kld.post.R") })
clusterCall(cl, function() { source("trace.R") })
clusterCall(cl, function() { source("covariate_generate.R") })
clusterCall(cl, function() { source("euc_dist_cal.R") })

WH_reef_no<- read.csv(file="WH_reef_no.csv")

design_com=site #Change design here

reef_sites<-list()
reef_sites_com<-list()

#create possible designs
for (i in 1:9)
{
  t=design_com[i]
  unique(WH_reef_no$Site[WH_reef_no$reef_no==t])
  reef_sites[[i]]<-unique(WH_reef_no$Site[WH_reef_no$reef_no==t])
  reef_sites_com[[i]]<-as.list(combn(reef_sites[[i]], 2))
  
}
new_design=do.call(cbind,reef_sites_com)
new_design=as.matrix(new_design,ncol = 9,nrow=6)


optimal_design<-optimal_des(new_design,rstart.d_id,kld.utility,site,B,theta_sim_all,r_rand1)

v_fileName <- paste(v_ut,"RData",sep=".")
save(list= c("optimal_design"),file = v_fileName)

stopCluster(cl)

