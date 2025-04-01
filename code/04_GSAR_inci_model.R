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
  "maps", "sf", "sp", "spdep", "SDPDmod", "spatialreg", "splm", "plm", "INLA"
))

rm(list = ls())
setwd("/Users/kexinxie/Downloads/GitHub/SpatialSpilloverEffect_Measles/")

## prepare dataset
model="mean" #c(mean,median,Q3)
tau_factor=TRUE
vhi_factor=TRUE

source("code/01_datascript.R")
source("code/02_GSAR_impacts.R")

######################################################################################################
############################################### incidence model ######################################
######################################################################################################

####### poisson model without SAR #######
pdata.inla <- pdata
n <- nrow(pdata.inla)
pdata.inla$idx <- 1:n

pdata.inla.scale <- pdata.inla
pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")] <- scale(pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")])


### Model definition 
f1.ols <- incidence ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment


### fitting the poisson model without SAR:
slmm1.ols <- inla( f1.ols,
                   data=pdata.inla.scale, 
                   E=pop,
                   family="xpoisson", #family="xpoisson", #family="zeroinflatedpoisson1", #
                   #control.family = list(hyper=theta.prior),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,
                                        return.marginals=TRUE,return.marginals.predictor=TRUE),
                   control.predictor=list(compute = T,link=1)
)

summary(slmm1.ols) 
plot(slmm1.ols)

### fitted value
fitted.value1ols.median <- slmm1.ols$summary.fitted.values$`0.5quant` * slmm1.ols$.args$E

### RMSE
print(sqrt(sum((fitted.value1ols.median-pdata.inla$incidence)^2)/length(pdata.inla$incidence))) 

### moran test
print(moran.test(pdata.inla$incidence-fitted.value1ols.median,listw=W_listw)) #0.1685


####### poisson lag model#######
library(INLA)
pdata.inla <- pdata
n <- nrow(pdata.inla)
pdata.inla$idx <- 1:n
pdata.inla$incidence_int <- round(pdata.inla$incidence,0)


pdata.inla.scale <- pdata.inla
pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")] <- scale(pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")])


### Model definition 
f1 <- incidence ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment

### Covariate matrix###
mmatrix1 <- model.matrix(f1,pdata.inla.scale)


#### prior for logit(prob) term
theta.prior = list(theta=list(
  initial = -10,
  prior = "gaussian",
  param = c(-10,100)))


e = eigenw(W_listw)
re.idx = which(abs(Im(e)) < 1e-6)
rho.max = 1/max(Re(e[re.idx]))
rho.min = 1/min(Re(e[re.idx]))
rho = mean(c(rho.min, rho.max))

### Precision matrix for beta coeffients’ prior
betaprec <- .001
Q.beta = Diagonal(n=ncol(mmatrix1), betaprec)

### Priors on the hyperparameters for zeroinflatedpoisson1
hyper = list(
  prec = list(
    prior = "loggamma",
    param = c(100,100)),
  rho = list(
    initial=0,
    prior = "normal",
    param = c(0,100)))

### R-code for fitting the spatial lag poisson model:
slmm1 <- inla( incidence ~ -1 +
                 f(idx, model="slm",
                   args.slm=list(
                     rho.min = rho.min,
                     rho.max = rho.max,
                     W=W,
                     X=mmatrix1,
                     Q.beta=Q.beta),
                   hyper=hyper),
               data=pdata.inla.scale,  
               E=pop,
               family="xpoisson", #family="xpoisson", #family="zeroinflatedpoisson1", #
               #control.family = list(hyper=theta.prior),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,
                                    return.marginals=TRUE,return.marginals.predictor=TRUE),
               control.predictor=list(compute = T,link=1)
)

summary(slmm1) 

### Summary of the coefficients (at the end of the vector of random effects)
beta.summary <- round(slmm1$summary.random$idx[n+1:ncol(mmatrix1),],3)
print(beta.summary)

beta.marginal <- slmm1$marginals.random$idx[n+1:ncol(mmatrix1)]

### Re-scale rho to real scale 
rhomarg <- inla.tmarginal(function(x){rho.min+x*(rho.max-rho.min)},
                          slmm1$marginals.hyperpar[[2]]) #for poisson
inla.zmarginal(rhomarg)
plot(rhomarg, type = "l", main = "Spatial autocorrelation")

### fitted value
fitted.value1.median <- slmm1$summary.fitted.values$`0.5quant`*slmm1$.args$E

### RMSE
print(sqrt(sum((fitted.value1.median-pdata.inla$incidence)^2)/length(pdata.inla$incidence))) #1.017716

### moran test
print(moran.test(pdata.inla$incidence-fitted.value1.median,listw=W_listw))


#### Impacts SLM ####
library(parallel)
options(mC.cores = 3)

samp_slmm1 <- inla.posterior.sample(1000, slmm1)

imp_slmm1_all <- compute_impacts_slm_poisson_all(samp_slmm1,
                                                 n.areas = n, W = W, n.var = 9, mmatrix = mmatrix1)


#### prediction ####
pred_51760_quantile<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix1)]
  rho <-  rho.min + theta[2] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix1[1:133,]
  pred0.dat[,2] <- 0
  linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- solve(diag(133) - rho * w) %*% pred0.dat %*% coeff
  fitted.value0<-exp(linear.fitted0)
  
  pred.dat<-mmatrix1[1:133,]
  pred.dat[,2] <- 0
  pred.dat[which(grepl("51760",row.names(pred.dat))),2]=1
  linear.fitted1 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted1 <- solve(diag(133) - rho * w) %*% pred.dat %*% coeff
  fitted.value1<-exp(linear.fitted1)
  
  inci_cost <- fitted.value1-fitted.value0
  
  return(inci_cost)
}),n=n,mmatrix1=mmatrix1,samp_slmm1)

### Prediction plot
pred_51760_inci_quantile.plot <- as.data.frame(t(apply(pred_51760_quantile,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(median=`50%`*1000000,
         UB=`97.5%`*1000000,
         LB=`2.5%`*1000000)

pred_51091_quantile<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix1)]
  rho <-  rho.min + theta[2] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix1[1:133,]
  pred0.dat[,2] <- 0
  linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- solve(diag(133) - rho * w) %*% pred0.dat %*% coeff#+matrix(rep(slmm1$offset.linear.predictor[1:133],m),nrow=133,ncol=m)
  fitted.value0<-exp(linear.fitted0)
  
  pred.dat<-mmatrix1[1:133,]
  pred.dat[,2] <- 0
  pred.dat[which(grepl("51091",row.names(pred.dat))),2]=1
  linear.fitted1 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted1 <- solve(diag(133) - rho * w) %*% pred.dat %*% coeff#+matrix(rep(slmm1$offset.linear.predictor[1:133],m),nrow=133,ncol=m)
  fitted.value1<-exp(linear.fitted1)
  
  inci_cost <- fitted.value1-fitted.value0
  
  return(inci_cost)
}),n=n,mmatrix1=mmatrix1,samp_slmm1)

pred_51091_inci_quantile.plot <- as.data.frame(t(apply(pred_51091_quantile,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(median=`50%`*1000000,
         UB=`97.5%`*1000000,
         LB=`2.5%`*1000000) 

### plot all median together
color<-terrain.colors(12)#viridisLite::cividis(12)

p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(median, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  #scale_fill_viridis(option = "C")+
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Richmond City")+#title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(median, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Highland County")+
  theme_bw()+theme(legend.position = "bottom")

pdf("result/plots/inci_median_all_plots.pdf",width=10, height=10)
ggpubr::ggarrange(p1,p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

### plot all LB together
color<-terrain.colors(12)#viridisLite::cividis(12)

p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(LB, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  #scale_fill_viridis(option = "C")+
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Richmond City")+#title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(LB, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Highland County")+
  theme_bw()+theme(legend.position = "bottom")

pdf("result/plots/inci_LB_all_plots.pdf",width=10, height=10)
ggpubr::ggarrange(p1,p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()


### plot all UB together
color<-terrain.colors(12)#viridisLite::cividis(12)

p1<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51760_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(UB, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  #scale_fill_viridis(option = "C")+
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((A))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Richmond City")+#title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

p2<-va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_51091_inci_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(UB, breaks=c(1e-8,1e-7,1e-6,5e-6,1e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.3))) %>% 
  filter(GEOID!=51515) %>% #10%:c(0,1e-6,1e-3,1e-2,1e-1,10,100,1000,10000,100000,1000000,1300000
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  scale_fill_manual(values=color,drop=F)+
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="Highland County")+
  theme_bw()+theme(legend.position = "bottom")

pdf("result/plots/inci_UB_all_plots.pdf",width=10, height=10)
ggpubr::ggarrange(p1,p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()


#### plot all spillover effects
county_list <- va.shp[which(va.shp$GEOID!=51515),]$GEOID
pred_spillover_incid<-inla.posterior.sample.eval((function(...) {
  coeff <- idx[n+1:ncol(mmatrix1)]
  rho <-  rho.min + theta[2] * (rho.max - rho.min)
  
  pred0.dat<-mmatrix1[1:133,]
  pred0.dat[,2] <- 0
  inver_mat <- solve(diag(133) - rho * w)
  # linear.fitted0 <- matrix(0, nrow = 133, ncol = m)
  linear.fitted0 <- inver_mat %*% pred0.dat %*% coeff
  fitted.value0<-exp(linear.fitted0)
  
  
  indirect_eff <- matrix(0,nrow=133,ncol=1)
  k <- 1
  for(ind in county_list){
    pred.dat<-mmatrix1[1:133,]
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
}),n=n,matrix1=mmatrix1,samp_slmm1)

pred_spillover_incid_quantile.plot <- as.data.frame(t(apply(pred_spillover_incid,1,function(x){quantile(x,c(0.025,0.5,0.975))}))) %>%
  mutate(county_fips=unique.county$county_fips,
         pop=unique.county$pop) %>%
  mutate(median=`50%`*1000000,
         UB=`97.5%`*1000000,
         LB=`2.5%`*1000000) 

p3_spillover <- va.shp %>% mutate(GEOID=as.integer(GEOID)) %>% 
  left_join(pred_spillover_incid_quantile.plot,by=c("GEOID"="county_fips")) %>%
  mutate(pred_group=cut(median, breaks=c(0.01,0.05,0.1,0.25,0.5,1,2.5,5,7.5,10))) %>%
  filter(GEOID!=51515) %>%
  ggplot() +
  geom_sf(aes(fill=pred_group),color="grey") +
  #scale_fill_viridis(option = "C")+
  scale_fill_manual(values=terrain.colors(9),drop=F)+  #viridisLite::viridis(9)
  annotate("text", x=-83.8, y=39.4, label= "bold((B))",size=5,parse = TRUE) +
  labs(fill="Incremental incidence case per 1,000,000 population",x= "",y="")+#title="Predicted incremental cost when MMR level reduced by 1% at Richmond City at benchmark scenario",
  theme_bw()+theme(legend.position = "bottom")

pdf("result/plots/spillover_incid_median_all_plots.pdf",width=10, height=7)
p3_spillover
dev.off()
