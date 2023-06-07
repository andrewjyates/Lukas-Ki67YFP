
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

packages <- c("tidyverse", "ggplot2", "readr", "forcats")

tmp=lapply(packages, require, character.only = TRUE)

#----F test for adding p as a free parameter 

load("DataObjects/flow_sim_frac_pointest_p_free.Rdata")
load("DataObjects/flow_sim_frac_pointest_p_zero.Rdata")
load("DataObjects/flow_frac_pointest_p_zero.Rdata")
load("DataObjects/flow_frac_pointest_p_free.Rdata")

F<-function(SSR1,SSR2,p1,p2,n){
  f1<-(SSR1-SSR2)/(p2-p1)
  f2<-SSR2/(n-p2)
  f<-f1/f2
  p_val<-1-pf(f, p2-p1, n-p2, lower.tail=F)
  return(p_val)
}

# compare p free and p=0 fits for the flow frac model with fixed thymic functions

p1 <- rep(5,2)#Theta,mu, beta, lp, N
p2 <- rep(6,2)#Theta,mu, beta, p, lp, N

n <- c(43,42)

SSR1<-map(flow_frac_pointest_p_zero, ~.$ssr)  # enter SSR model with less parameters
SSR2<-map(flow_frac_pointest_p_free, ~.$ssr)  # enter SSR model with more parameters

F_flow_frac <- map(c(1,2), ~ F(SSR1[[.]],SSR2[[.]],p1[[.]],p2[[.]],n[[.]])) 

###############################################
# same for simultaneous fit: 

p1 <- rep(10,2)#f,a,b,c,a2,b2,c2,theta,mu,beta
p2 <- rep(11,2)#f,a,b,c,a2,b2,c2,theta,mu,beta,p

n <- c(42,42)

SSR1<-map(flow_sim_frac_pointest_p_zero, ~.$ssr)  # enter SSR model with less parameters
SSR2<-map(flow_sim_frac_pointest_p_free, ~.$ssr)  # enter SSR model with more parameters

F_flow_frac_sim <- map(c(1,2), ~ F(SSR1[[.]],SSR2[[.]],p1[[.]],p2[[.]],n[[.]]))  # give p-value for H0: hance that data generated under null hyp (=simpler model) , 
# if <0.05 then "model with less parameters" is better

cat("p-values (improvement by adding p nonzero?)\n")
cat("CD4 flow frac", F_flow_frac[[1]], "\n")
cat("CD8 flow frac", F_flow_frac[[2]], "\n")
cat("CD4 flow frac sim", F_flow_frac_sim[[1]], "\n")
cat("CD4 flow frac sim", F_flow_frac_sim[[2]], "\n")

