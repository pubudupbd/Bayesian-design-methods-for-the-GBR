############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Calculate Euclidean distances
## Author: Pubudu Thialn
#############################################################################


WH_site<- read.csv(file="site_cord.csv",head=TRUE,sep=",") # site locations; long anf lat

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

euc_dist_cal<-function(cord){
WH_site_pick<-WH_site[cord,][,2:3]
WH_cord<-rbind(WH_site_pick,WH_site_pick,WH_site_pick,WH_site_pick,WH_site_pick)
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

