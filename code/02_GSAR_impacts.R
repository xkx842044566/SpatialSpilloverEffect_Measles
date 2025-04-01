#-----------------------------------------------------------
# Compute Impacts for a Single Covariate (Poisson Model)
#-----------------------------------------------------------
compute_impacts_slm_poisson <- function(..., n.areas, W, n.var, intercept = 1, mmatrix) {
  
  # Extract coefficient for the nth variable
  coeff <- idx[n.areas + intercept + n.var]
  
  # Linear predictor and its exponential transformation
  eta <- Predictor[1:n.areas]
  deta <- exp(eta)
  
  # Spatial autocorrelation parameter
  rho <- rho.min + theta[2] * (rho.max - rho.min)
  
  # Spatial multiplier matrix
  SW <- solve(Diagonal(n.areas, x = 1) - rho * W) * coeff * deta
  
  # Compute direct, indirect, and total impacts
  total.impact <- sum(SW) / n.areas
  dir.impact   <- sum(diag(SW)) / n.areas
  indir.impact <- total.impact - dir.impact
  
  return(c(dir.impact, indir.impact, total.impact))
}


#-----------------------------------------------------------
# Compute Impacts for All Covariates (Poisson Model)
#-----------------------------------------------------------
compute_impacts_slm_poisson_all <- function(inla_samples, n.areas, W, n.var, intercept = TRUE, mmatrix) {
  
  res <- sapply(1:n.var, function(X) {
    
    impacts <- inla.posterior.sample.eval(
      compute_impacts_slm_poisson,
      inla_samples,
      n.areas = n.areas,
      W = W,
      n.var = X,
      mmatrix = mmatrix
    )
    
    # Summary statistics
    stats <- c(
      apply(impacts, 1, mean),
      apply(impacts, 1, sd),
      apply(impacts, 1, function(x) quantile(x, 0.025)),
      apply(impacts, 1, function(x) quantile(x, 0.975))
    )
    
    return(stats)
  })
  
  res <- as.data.frame(res)
  colnames(res) <- colnames(mmatrix)[intercept + 1:n.var]
  rownames(res) <- c("direct (mean)", "indirect (mean)", "total (mean)",
                     "direct (s.d.)", "indirect (s.d.)", "total (s.d.)",
                     "direct (2.5%)", "indirect (2.5%)", "total (2.5%)",
                     "direct (97.5%)", "indirect (97.5%)", "total (97.5%)")
  
  return(res)
}


#-----------------------------------------------------------
# Compute Impacts for a Single Covariate (Gamma Model)
#-----------------------------------------------------------
compute_impacts_slm_gamma <- function(..., n.areas, W, n.var, intercept = 1, mmatrix) {
  
  # Extract coefficient
  coeff <- idx[n.areas + intercept + n.var]
  
  # Linear predictor and exponential transform
  eta  <- Predictor[1:n.areas]
  deta <- exp(eta)
  
  # Spatial autocorrelation parameter (for Gamma)
  rho <- rho.min + theta[3] * (rho.max - rho.min)
  
  # Spatial multiplier matrix
  SW <- solve(Diagonal(n.areas, x = 1) - rho * W) * coeff * deta
  
  # Compute impacts
  total.impact <- sum(SW) / n.areas
  dir.impact   <- sum(diag(SW)) / n.areas
  indir.impact <- total.impact - dir.impact
  
  return(c(dir.impact, indir.impact, total.impact))
}


#-----------------------------------------------------------
# Compute Impacts for All Covariates (Gamma Model)
#-----------------------------------------------------------
compute_impacts_slm_gamma_all <- function(inla_samples, n.areas, W, n.var, intercept = TRUE, mmatrix) {
  
  res <- sapply(1:n.var, function(X) {
    
    impacts <- inla.posterior.sample.eval(
      compute_impacts_slm_gamma,
      inla_samples,
      n.areas = n.areas,
      W = W,
      n.var = X,
      mmatrix = mmatrix
    )
    
    # Summary statistics
    stats <- c(
      apply(impacts, 1, mean),
      apply(impacts, 1, sd),
      apply(impacts, 1, function(x) quantile(x, 0.025)),
      apply(impacts, 1, function(x) quantile(x, 0.975))
    )
    
    return(stats)
  })
  
  res <- as.data.frame(res)
  colnames(res) <- colnames(mmatrix)[intercept + 1:n.var]
  rownames(res) <- c("direct (mean)", "indirect (mean)", "total (mean)",
                     "direct (s.d.)", "indirect (s.d.)", "total (s.d.)",
                     "direct (2.5%)", "indirect (2.5%)", "total (2.5%)",
                     "direct (97.5%)", "indirect (97.5%)", "total (97.5%)")
  
  return(res)
}
