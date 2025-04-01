install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, repos = "http://cran.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

install_and_load(c(
  "data.table", "readxl", "dplyr", "tidyverse", "ggplot2", "usmap",
  "maps", "sf", "sp", "spdep", "SDPDmod", "spatialreg", "splm", "plm"
))

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

rm(list = ls())
setwd("/Users/kexinxie/Downloads/GitHub/SpatialSpilloverEffect_Measles/")

## prepare dataset
model="mean" #c(mean,median,Q3)
tau_factor=TRUE
vhi_factor=TRUE

source("code/01_datascript.R")
source("code/02_GSAR_impacts.R")

######################################################################################################
############################################### cost model ###########################################
######################################################################################################


####### gamma model without SAR #######

pdata.inla <- pdata
n <- nrow(pdata.inla)
pdata.inla$idx <- 1:n

### Model definition 
pdata.inla.scale <- pdata.inla
pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")] <- scale(pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")])


### Model definition 
f0.ols <- totalCost/pop ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment

### Zero-variance for error term
zero.variance =  list(prec=list(initial = 10,fixed=FALSE))


### fitting the gamma lag model without SAR
slmm0.ols<- inla( f0.ols,
                  data=pdata.inla.scale, family="gamma",
                  scale=100000,
                  #control.family = list(variant=0),
                  control.family = list(hyper=zero.variance),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,return.marginals=TRUE,return.marginals.predictor=TRUE)
                  #verbose = TRUE,
                  #offset=log(pop)
)

summary(slmm0.ols) 

### RMSE
print(sqrt(sum((slmm0.ols$summary.fitted.values$`0.5quant`-pdata.inla$totalCost/pdata.inla$pop)^2)/length(pdata.inla$totalCost))) #33.71079

### moran test
print(moran.test(pdata.inla$totalCost/pdata.inla$pop-slmm0.ols$summary.fitted.values$`0.5quant`,listw=W_listw)) 


####### gamma model with SAR #######

pdata.inla <- pdata
n <- nrow(pdata.inla)
pdata.inla$idx <- 1:n

pdata.inla.scale <- pdata.inla
pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")] <- scale(pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")])

### Model definition 
f0 <- totalCost/pop ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment

### Covariate matrix
mmatrix0 <- model.matrix(f0,pdata.inla.scale)

### Zero-variance for error term
zero.variance = list(prec=list(initial = 10, 
                               prior = "loggamma",
                               param= c(10,10)))


### Compute eigenvalues for slm model, used to obtain rho.min and rho.max
e = eigenw(W_listw)
re.idx = which(abs(Im(e)) < 1e-6)
rho.max = 1/max(Re(e[re.idx]))
rho.min = 1/min(Re(e[re.idx]))
rho = mean(c(rho.min, rho.max))

### Precision matrix for beta coeffients’ prior
betaprec <- .001
Q.beta = Diagonal(n=ncol(mmatrix0), betaprec)

### Priors on the hyperparameters
hyper0 = list(
  prec = list(
    prior = "loggamma",
    param = c(10, 10)),
  rho = list(
    initial=0,
    prior = "normal",
    param = c(0,100)))

### fitting the spatial lag gamma model:
slmm0 <- inla( totalCost/pop ~ -1 +
                 f(idx, model="slm",
                   args.slm=list(
                     rho.min = rho.min,
                     rho.max = rho.max,
                     W=W,
                     X=mmatrix0,
                     Q.beta=Q.beta),
                   hyper=hyper0),
               data=pdata.inla.scale, family="gamma",
               scale=100000,
               #control.family = list(variant=0),
               control.family = list(hyper=zero.variance),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,return.marginals=TRUE,return.marginals.predictor=TRUE)
               #verbose = TRUE,
               #offset=log(pop)
)

summary(slmm0) #DIC=77445.46 waic=76221.7

### Re-scale rho to real scale
rhomarg <- inla.tmarginal(function(x){rho.min+x*(rho.max-rho.min)},
                          slmm0$marginals.hyperpar[[3]])
inla.zmarginal(rhomarg)
plot(rhomarg, type = "l", main = "Spatial autocorrelation")

### Summary of the coefficients (at the end of the vector of random effects)
beta.summary <- round(slmm0$summary.random$idx[n+1:ncol(mmatrix0),],3)
print(beta.summary)
beta.marginal <- slmm0$marginals.random$idx[n+1:ncol(mmatrix0)]

### fitted value
fitted.value0 <- slmm0$summary.fitted.values

### RMSE
print(sqrt(sum((slmm0$summary.fitted.values$`0.5quant`-pdata.inla$totalCost/pdata.inla$pop)^2)/length(pdata.inla$totalCost))) #0.01901173

### moran test
print(moran.test(pdata.inla$totalCost/pdata.inla$pop-slmm0$summary.fitted.values$`0.5quant`,listw=W_listw))


#### Impacts SLM ####
library(parallel)
options(mC.cores = 3)

samp_slmm0 <- inla.posterior.sample(1000, slmm0)
imp_slmm0_all <- compute_impacts_slm_gamma_all(samp_slmm0,
                                               n.areas = n, W = W, n.var = 9, mmatrix = mmatrix0)


#### prediction ####
pred_51760_quantile<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix0)]
  rho <-  rho.min + theta[3] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix0[1:133,]
  pred0.dat[,2] <- 0
  # linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- solve(diag(133) - rho * w) %*% pred0.dat %*% coeff
  fitted.value0<-exp(linear.fitted0)
  
  pred.dat<-mmatrix0[1:133,]
  pred.dat[,2] <- 0
  pred.dat[which(grepl("51760",row.names(pred.dat))),2]=1
  #linear.fitted1 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted1 <- solve(diag(133) - rho * w) %*% pred.dat %*% coeff
  fitted.value1<-exp(linear.fitted1)
  
  inci_cost <- fitted.value1-fitted.value0
  
  return(inci_cost)
}),n=n,matrix0=mmatrix0,samp_slmm0)


### Prediction plot
pred_51760_quantile.plot <- as.data.frame(t(apply(pred_51760_quantile,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(inci_cost_median=`50%`*pop,
         inci_cost_per_pop_median=`50%`,
         inci_cost_UB=`97.5%`*pop,
         inci_cost_per_pop_UB=`97.5%`,
         inci_cost_LB=`2.5%`*pop,
         inci_cost_per_pop_LB=`2.5%`)

pred_51091_quantile<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix0)]
  rho <-  rho.min + theta[3] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix0[1:133,]
  pred0.dat[,2] <- 0
  #linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- solve(diag(133) - rho * w) %*% pred0.dat %*% coeff
  fitted.value0<-exp(linear.fitted0)
  
  pred.dat<-mmatrix0[1:133,]
  pred.dat[,2] <- 0
  pred.dat[which(grepl("51091",row.names(pred.dat))),2]=1
  #linear.fitted1 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted1 <- solve(diag(133) - rho * w) %*% pred.dat %*% coeff
  fitted.value1<-exp(linear.fitted1)
  
  inci_cost <- fitted.value1-fitted.value0
  
  return(inci_cost)
}),n=n,matrix0=mmatrix0,samp_slmm0)

pred_51091_quantile.plot <- as.data.frame(t(apply(pred_51091_quantile,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(inci_cost_median=`50%`*pop,
         inci_cost_per_pop_median=`50%`,
         inci_cost_UB=`97.5%`*pop,
         inci_cost_per_pop_UB=`97.5%`,
         inci_cost_LB=`2.5%`*pop,
         inci_cost_per_pop_LB=`2.5%`)

#### plot median together 
color<-viridisLite::plasma(12)
color0<-viridisLite::viridis(12)

p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_median, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)",x= "",y="Richmond City")+
  theme_bw()+theme(legend.position = "bottom")

p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_median, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)",x= "",y=" ")+ #title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p3<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_median, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((C))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)",x= "",y="Highland County")+#,title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario")+
  theme_bw()+theme(legend.position = "bottom")

p4<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_median, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>% 
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((D))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)",x= "",y=" ")+ #title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p13<-ggpubr::ggarrange(p1, p3, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
p24<-ggpubr::ggarrange(p2, p4, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

pdf("result/plots/cost_median_all_plots.pdf",width=18, height=9)
ggpubr::ggarrange(p13,p24, ncol = 2, nrow = 1)
dev.off()

#### plot UB together
color<-viridisLite::plasma(12)
color0<-viridisLite::viridis(12)
p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_UB, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)")+
  theme_bw()+theme(legend.position = "bottom")

p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_UB, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)")+ #title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p3<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_UB, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((C))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)")+#,title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario")+
  theme_bw()+theme(legend.position = "bottom")

p4<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_UB, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>% 
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((D))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)")+ #title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p13<-ggpubr::ggarrange(p1, p3, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
p24<-ggpubr::ggarrange(p2, p4, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

pdf("result/plots/cost_UB_all_plots.pdf",width=18, height=9)
ggpubr::ggarrange(p13,p24, ncol = 2, nrow = 1)
dev.off()


#### plot LB together 
color<-viridisLite::plasma(12)
color0<-viridisLite::viridis(12)

p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_LB, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)")+
  theme_bw()+theme(legend.position = "bottom")


p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_LB, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)")+ #title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")


p3<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_LB, breaks=c(0.001,0.01,0.05,0.1,0.5,1,10,100,500,1000,10000,50000,70000))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((C))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)")+#,title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario")+
  theme_bw()+theme(legend.position = "bottom")


p4<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_LB, breaks=c(1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5))) %>%
  filter(GEOID!=51515) %>% 
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color0,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((D))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)")+ #title="Predicted incremental cost when MMR level reduced by 1% at Highland county at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p13<-ggpubr::ggarrange(p1, p3, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
p24<-ggpubr::ggarrange(p2, p4, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

pdf("result/plots/cost_LB_all_plots.pdf",width=18, height=9)
ggpubr::ggarrange(p13,p24, ncol = 2, nrow = 1)
dev.off()

#### plot all spillover effects
county_list <- va.shp[which(va.shp$GEOID!=51515),]$GEOID
pred_spillover<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix0)]
  rho <-  rho.min + theta[3] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix0[1:133,]
  pred0.dat[,2] <- 0
  inver_mat <- solve(diag(133) - rho * w)
  # linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- inver_mat %*% pred0.dat %*% coeff
  fitted.value0<-exp(linear.fitted0)
  
  
  indirect_eff <- matrix(0,nrow=133,ncol=1)
  k <- 1
  for(ind in county_list){
    pred.dat<-mmatrix0[1:133,]
    pred.dat[,2] <- 0
    pred.dat[which(grepl(ind,row.names(pred.dat))),2]=1
    #linear.fitted1 <- matrix(0, nrow = 133, ncol = m)
    linear.fitted1 <- inver_mat %*% pred.dat %*% coeff
    fitted.value1<-exp(linear.fitted1)
    inci_cost <- fitted.value1-fitted.value0
    indirect_eff[k,] <- apply(inci_cost,2,sum)
    k <- k+1
  }
  return(indirect_eff)
}),n=n,matrix0=mmatrix0,samp_slmm0)

pred_spillover_quantile.plot <- as.data.frame(t(apply(pred_spillover,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(inci_cost_median=`50%`*pop,
         inci_cost_per_pop_median=`50%`,
         inci_cost_UB=`97.5%`*pop,
         inci_cost_per_pop_UB=`97.5%`,
         inci_cost_LB=`2.5%`*pop,
         inci_cost_per_pop_LB=`2.5%`)


p1_spillover<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_spillover_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_median, breaks=c(1000,5000,10000,25000,50000,75000,100000,500000,550000))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=viridisLite::plasma(8),drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Cumulative incremental cost (US$)",x= "",y="")+
  theme_bw()+theme(legend.position = "bottom")

p2_spillover<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_spillover_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(inci_cost_per_pop_median, breaks=c(0.01,0.05,0.1,0.25,0.5,0.75,1,1.5,2,2.5,3))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=viridisLite::viridis(10),drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Incremental cost per person (US$)",x= "",y="")+ #title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

pdf("result/plots/spillover_cost_median_all_plots.pdf",width=10, height=10)
ggpubr::ggarrange(p1_spillover,p2_spillover, ncol = 1, nrow = 2)
dev.off()