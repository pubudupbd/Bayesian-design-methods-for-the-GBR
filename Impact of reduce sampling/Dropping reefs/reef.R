############################################################################
## Project: Bayesian optimal design: improving the effectiveness of reef monitoring with time varying covariate information
## Script purpose: Create possible designs
## Author: Pubudu Thialn

# set number for reef names
WH<- read.csv(file="all_data_final_site.csv",head=TRUE,sep=",")
reef_no<-read.csv("reef_name_no.csv",head=TRUE,sep=",")
WH_reef_no<-merge(x = WH, y = reef_no, by = "Name", all = TRUE)
write.csv(WH_reef_no,"WH_reef_no.csv")

# All possible reef combinations
reef<-seq(1,9,1)
all_possible_reef=combn(reef, 8)

write.csv(all_possible_reef,"all_possible_reef.csv")


