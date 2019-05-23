############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Simulate covariates
## Author: Pubudu Thialn
#############################################################################

WH<- read.csv(file="all_data_final_site.csv",head=TRUE,sep=",") #standardised covariates
WH$bathy<-(WH$bathy-mean(WH$bathy))/sd(WH$bathy)
WH$chloro_water_column<-(WH$chloro_water_column-mean(WH$chloro_water_column))/sd(WH$chloro_water_column)
WH$CRS_T_AV<-(WH$CRS_T_AV-mean(WH$CRS_T_AV))/sd(WH$CRS_T_AV)

dis_mat<- read.csv(file="distribution_mat.csv",head=TRUE,sep=",")# The estimated parameters for the distributions of time-varying disturbance covariates at each site in the Whitsunday region
next_year_covariate<-list()

Get_Covariate=function(site,seed_val)
{
  set.seed(seed_val)
  for(k in 1:5){
  covariate_sub<-matrix(0,18,10)
  for(i in 1:18){
    t=site[i]
    covariate_sub[i,1]<-ifelse(WH$SHELF[WH$Site==t][1]=='M',1,0)# site specific
    covariate_sub[i,2]<-ifelse(WH$SHELF[WH$Site==t][1]=='O',1,0)
    covariate_sub[i,3]<-WH$notake[WH$Site==t][1]
    covariate_sub[i,4]<-WH$bathy[WH$Site==t][1]
    covariate_sub[i,5]<-WH$chloro_water_column[WH$Site==t][1]
    covariate_sub[i,6]<-WH$CRS_T_AV[WH$Site==t][1]
    r<-runif(1,0,1) #time varying covariates
    if(r<dis_mat[t,3]){
      covariate_sub[i,7]<-0
    }
    else{covariate_sub[i,7] <-rnorm(1,dis_mat[t,4],dis_mat[t,5])}
    
    covariate_sub[i,8] <- sample(c(0,1), 1, replace=TRUE, prob=c(dis_mat[t,1], (1-dis_mat[t,1])) )
    covariate_sub[i,9] <- sample(c(0,1), 1, replace=TRUE, prob=c(dis_mat[t,2], (1-dis_mat[t,2])) )
    covariate_sub[i,10]<-1.50537634 #standardized year 2017
  }
  next_year_covariate[[k]]<-covariate_sub
}
next_year_t=do.call(rbind,next_year_covariate)
return(next_year_t)
}




