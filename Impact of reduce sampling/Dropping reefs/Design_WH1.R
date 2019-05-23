############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: The effect of visiting fewer LTMP sites by dropping reefs  in a given region  
## Author: Pubudu Thialn
#############################################################################

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


E=20# independent evaluations
B=1000
leny <<- 120 # number of data points for next sampling year
L<<-100
v0<<-0.00001 #nugect

mu_prior <<- c(-2.5232,-5.9825,-1.1189,-1.2726,0.153,0.9066,0.2782,-0.1062,-0.7962,-0.23,-0.0093,-0.449,-0.2152,-0.0431)
Sigma_prior <<- diag(c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2))
Sigma_priorb <<- c(0.0421^2,0.4846^2,0.0596^2,0.0779^2,0.0833^2,0.2077^2,0.0864^2,0.0227^2,0.0969^2,0.0534^2,0.0063^2,0.0528^2,0.067^2,0.0289^2)
iSigma_prior <<- solve(Sigma_prior)
dSigma_prior <<- det(Sigma_prior)
num_par <<- length(mu_prior)

# Simulate from prior 
theta_sim <- matrix(0,nrow=B*E,ncol=length(mu_prior))
for (j in 1:length(mu_prior)){
  theta_sim[,j] <-rnorm(B*E,mu_prior[j],sqrt(Sigma_priorb[j]))
}

theta_sim <- data.frame(theta_sim)
theta_sim0 <- as.data.frame(sweep(theta_sim,2,colMeans(theta_sim),"-"))  # shift the mean 
C <- chol(nearPD(cov(theta_sim0))$mat)
SN <- matrix(rnorm(B * ncol(theta_sim)*E), ncol(theta_sim))#standard normals
X <- t(C) %*% SN
theta_mul <- data.frame(as.matrix(t(X)))
names(theta_mul) <- names(theta_sim)
theta_mul <- as.data.frame(sweep(theta_mul,2,colMeans(theta_sim),"+")) 
theta_sim_all<-as.matrix(theta_mul)


strand_all=matrix(rnorm(leny*B*E,0,1),leny*E,B)


### Projecting coordinates into UTM coordinates system
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

# Calculate Euclidean distances
eu.dist<-function(X)
{
  tmp<-dist(X, diag = TRUE,upper = TRUE,method = "euclidean", p = 2)
  ans<-as.matrix(tmp)
  return(ans)
}

# #Euclidean distance 
# ###########################################################
WH_site<- read.csv(file="WH_site.csv",head=TRUE,sep=",")
euc_dis_sites<-function(WH_site,d)
{
  WH_site=WH_site[d,]
  WH_cord<-rbind(WH_site,WH_site,WH_site,WH_site,WH_site)
  x1<-WH_cord$longitude
  y1<-WH_cord$latitude
  
  UTM1<-LongLatToUTM(x1,y1,55)
  
  log1<-UTM1$X
  lat1<-UTM1$Y
  
  #standardization of distances
  mean_log1<-mean(log1)
  sd_log1<-sd(log1)
  standard_log1<-(log1-mean_log1)/sd_log1
  
  mean_lat1<-mean(lat1)
  sd_lat1<-sd(lat1)
  standard_lat1<-(lat1-mean_lat1)/sd_lat1
  
  standard_cord1<-t(rbind(standard_log1,standard_lat1))
  euclidean_nextY<-eu.dist(standard_cord1)
  euclidean_nextY2<-euclidean_nextY^2
  
  return(euclidean_nextY2)
}

kld.utility<-function(d,B,theta_sim,r_rand1){
  for (j in 1:B){
    r_rand2<<- matrix(rnorm(leny*L,0,1),leny,L)
    X<-cbind(1,Get_Covariate(d,20))
    theta = theta_sim[,(4:num_par)][j,]
    
    cov<-v0*diag(rep(1,leny))+exp(theta_sim[,2][j])*exp(-euclidean_nextY2/(exp(theta_sim[,3][j]))^2)
    zb <-(chol(cov))%*%r_rand1[,j]
    
    lp = X%*%theta+zb
    mu=exp(lp)/(1+exp(lp))
    phi=1/exp(theta_sim[,1][j])
    Y =rbeta(leny,mu*phi,(1-mu)*phi)
    Y=ifelse(1-Y< 1e-5,0.9999990,Y)
    Y=ifelse(Y< 1e-5,0.000001,Y)
    
    inits=theta_sim[j,]
    kld_out =kld.post(X,Y,inits)
    kld.post.val[j]=kld_out
  }
  return(kld.post.val)
}

all_possible_design<-read.csv("all_possible_reef.csv",header = T)
design_com=all_possible_design[,1] #Change design here

new_design_tem<-list()

for (i in 1:8)
{
  t=design_com[i]
  new_design_tem[[i]]<-unique(WH_reef_no$Site[WH_reef_no$reef_no==t])
}

new_design=do.call(rbind,new_design_tem)

new_design=matrix(new_design,ncol = 24,nrow=1)

euclidean_nextY2<<-euc_dis_sites(WH_site,new_design)

index <- Sys.getenv("PBS_ARRAY_INDEX")
n_ind <- as.numeric(index)

index <- (1+(n_ind-1)*B):(n_ind*B)
index <- index[index<=E*B]
theta_sim<-theta_sim_all[index,]

index1 <- (1+(n_ind-1)*leny):(n_ind*leny)
index1 <- index1[index1<=E*B*leny]
r_rand1<- strand_all[index1,]

kld.post.val=matrix(0,B)
KLD_udB1000=mean(kld.utility(new_design,B,theta_sim,r_rand1))

v_ut <- "design1_"
v_ut <- paste(v_ut,n_ind,sep="")

v_fileName <- paste(v_ut,"RData",sep=".")
save(list= c("KLD_udB1000"),file = v_fileName)




