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
model="mean" #c(mean,median,Q3)
tau_factor=TRUE
vhi_factor=TRUE

source("code/01_datascript.R")
source("code/02_GSAR_impacts.R")

##################################### cross validation predictive model ####################################
rmse_cost_lm<- c()
rmse_cost_gsar<- c()
rmse_incidence_lm <- c()
rmse_incidence_gsar <- c()

k<-1

for (k in 1:133){
  
  pdata.cv <- load_cv_data(k=k)
  
  pdata.inla <- pdata.cv$pdata
  n <- nrow(pdata.inla)
  pdata.inla$idx <- 1:n
  
  pdata.inla.scale <- pdata.inla
  pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")] <- scale(pdata.inla.scale[,c("avg_income","prop_gender","prop_less_5","prop_employment")])
  
  ## Cost
  ### Model definition 
  f0 <- totalCost/pop ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment
  
  ### Covariate matrix
  mmatrix0 <- model.matrix(f0,pdata.inla.scale)
  
  ### Zero-variance for error term
  zero.variance = list(prec=list(initial = 10, 
                                 prior = "loggamma",
                                 param= c(10,10)))
  
  
  ### Compute eigenvalues for slm model, used to obtain rho.min and rho.max
  e = eigenw(pdata.cv$W_listw)
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
  
  slmm0.ols<- inla( f0,
                    data=pdata.inla.scale, family="gamma",
                    scale=100000,
                    #control.family = list(variant=0),
                    control.family = list(hyper=list(prec=list(initial = 10,fixed=FALSE))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,return.marginals=TRUE,return.marginals.predictor=TRUE)
                    #verbose = TRUE,
                    #offset=log(pop)
  )
  
  ### fitting the spatial lag gamma model:
  slmm0 <- inla( totalCost/pop ~ -1 +
                   f(idx, model="slm",
                     args.slm=list(
                       rho.min = rho.min,
                       rho.max = rho.max,
                       W=pdata.cv$W,
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
  
  
  ### fitted value
  fitted.value0 <- slmm0$summary.fitted.values
  
  # RMSE and store results
  rmse_cost_lm <- c(rmse_cost_lm,sqrt(sum((slmm0.ols$summary.fitted.values$`0.5quant`-pdata.inla$totalCost/pdata.inla$pop)^2)/length(pdata.inla$totalCost)))
  rmse_cost_gsar <- c(rmse_cost_gsar,sqrt(sum((slmm0$summary.fitted.values$`0.5quant`-pdata.inla$totalCost/pdata.inla$pop)^2)/length(pdata.inla$totalCost)))
  
  ## Incidence
  ### Model definition 
  f1 <- incidence ~ alpha+tau+vhi+avg_income+prop_gender+prop_less_5+prop_employment
  
  ### Covariate matrix###
  mmatrix1 <- model.matrix(f1,pdata.inla.scale)
  
  
  #### prior for logit(prob) term
  theta.prior = list(theta=list(
    initial = -10,
    prior = "gaussian",
    param = c(-10,100)))
  
  
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
  
  slmm1.ols <- inla( f1,
                     data=pdata.inla.scale, 
                     E=pop,
                     family="xpoisson", #family="xpoisson", #family="zeroinflatedpoisson1", #
                     #control.family = list(hyper=theta.prior),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,config = TRUE,
                                          return.marginals=TRUE,return.marginals.predictor=TRUE),
                     control.predictor=list(compute = T,link=1)
  )
  
  fitted.value1ols.median <- slmm1.ols$summary.fitted.values$`0.5quant` * slmm1.ols$.args$E
  
  ### fitting the spatial lag poisson model:
  slmm1 <- inla( incidence ~ -1 +
                   f(idx, model="slm",
                     args.slm=list(
                       rho.min = rho.min,
                       rho.max = rho.max,
                       W=pdata.cv$W,
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
  
  fitted.value1.median <- slmm1$summary.fitted.values$`0.5quant`*slmm1$.args$E
  
  ### RMSE
  rmse_incidence_lm <- c(rmse_incidence_lm, sqrt(sum((fitted.value1ols.median-pdata.inla$incidence)^2)/length(pdata.inla$incidence)))
  
  rmse_incidence_gsar <- c(rmse_incidence_gsar, sqrt(sum((fitted.value1.median-pdata.inla$incidence)^2)/length(pdata.inla$incidence)))
  
  
}

cost.sar.rmse <- sqrt(mean(rmse_cost_gsar)) 
cost.lm.rmse <- sqrt(mean(rmse_cost_lm)) 
incidence.sar.rmse <- sqrt(mean(rmse_incidence_gsar)) 
incidence.lm.rmse <- sqrt(mean(rmse_incidence_lm)) 