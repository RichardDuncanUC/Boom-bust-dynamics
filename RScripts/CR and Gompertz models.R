
#-------------------------------------------------------------------------------
# Consumer resource model - fixed log_K_mean
#-------------------------------------------------------------------------------

CR <- nimbleCode({
  
  ## ------------------------------------------------------------
  ## Hierarchical spatial priors on CR parameters
  ## ------------------------------------------------------------

  ## r_d > 0  (log scale) - intrinsic growth rate
  log_r_d_mean ~ dnorm(0, sd = 1)  
  r_d_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.r_d, zero_mean = 1)
  for(i in 1:Nsite) {
    log_r_d[i] <- log_r_d_mean + r_d_raw[i]
    r_d[i]     <- exp(log_r_d[i])
  }
  
  ## K > 0  (log scale) - carrying capacity
  log_K_mean <- log(100)
  K_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.K, zero_mean = 1)
  for(i in 1:Nsite) {
    log_K[i] <- log_K_mean + K_raw[i]
    K[i]     <- exp(log_K[i])
  }
  
  ## alpha > 0  (log scale) - attack rate
  log_alpha_mean ~ dnorm(0, sd = 1)
  alpha_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.alpha, zero_mean = 1)
  for(i in 1:Nsite) {
    log_alpha[i] <- log_alpha_mean + alpha_raw[i]
    alpha[i]     <- exp(log_alpha[i])
  }
  
  ## eps in (0,1)  (logit scale) - conversion efficiency
  logit_eps_mean ~ dnorm(0, sd = 1)
  eps_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.eps, zero_mean = 1)
  for(i in 1:Nsite) {
    logit_eps[i] <- logit_eps_mean + eps_raw[i]
    eps[i]       <- ilogit(logit_eps[i])
  }
  
  ## mu in (0,1)  (logit scale) - mortality rate 
  logit_mu_mean ~ dnorm(0, sd = 1)
  mu_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.mu, zero_mean = 1)
  for(i in 1:Nsite) {
    logit_mu[i] <- logit_mu_mean + mu_raw[i]
    mu[i]       <- ilogit(logit_mu[i])
  }
  
  ## CAR precisions (parameterized via SD for better mixing)
  sd.r_d   ~ dexp(1)
  sd.K     ~ dexp(1)
  sd.alpha ~ dexp(1)
  sd.eps   ~ dexp(1)
  sd.mu    ~ dexp(1)
  tau.r_d   <- 1 / sd.r_d^2
  tau.K     <- 1 / sd.K^2
  tau.alpha <- 1 / sd.alpha^2
  tau.eps   <- 1 / sd.eps^2
  tau.mu    <- 1 / sd.mu^2

  ## ------------------------------------------------------------
  ## Site-specific consumer–resource trajectories starting at first_year
  ## Model density
  ## ------------------------------------------------------------
  for (i in 1:Nsite) {
    
    ## Initial states
    logC[i, first_year[i]] ~ dnorm(alpha0, tau.init.C)  
    log(lambda[i, first_year[i]]) <- logC[i, first_year[i]] + log_area_offset[i, first_year[i]]
    R[i, first_year[i]] <- K[i]
    
    for (j in (first_year[i] + 1):Nperiod) {
      
      ## Observation model: Negative Binomial 
      y[i, j] ~ dnbinom(size = size, prob = prob[i, j])
      prob[i, j] <- size / (size + lambda[i, j])
      log(lambda[i, j]) <- logC[i, j] + log_area_offset[i, j]
      
      # Process model
      ## Consumer update (log scale) with process error 
      mu.logC[i, j] <- logC[i, j-1] + 
                       log(1 - mu[i] + eps[i] * alpha[i] * R[i, j-1])
      logC[i, j] ~ dnorm(mu.logC[i, j], tau.process.C)
      
      ## Resource update (natural scale) with floor constraint
      R_temp[i, j] <- R[i, j-1] +
                      r_d[i] * R[i, j-1] * (1 - R[i, j-1] / K[i]) - 
                      alpha[i] * R[i, j-1] * exp(logC[i, j-1])
      ## Floor constraint to prevent negative resources
      R[i, j] <- max(R_temp[i, j], 0.001)
    }
    
    ## Link to density
    for (j in first_year[i]:Nperiod) {
      log_density[i, j] <- logC[i,j]
      density[i, j] <- exp(log_density[i, j])
    }
  }
  
  ## ------------------------------------------------------------
  ## Hyperpriors and observation priors
  ## ------------------------------------------------------------
  
  ## Initial conditions
  alpha0     ~ dnorm(-4, sd = 2)       # small initial consumer population
  tau.init.C ~ dgamma(2, 1)
  sd.init.C <- 1 / sqrt(tau.init.C)
  
  ## Process error for consumer dynamics
  tau.process.C ~ dgamma(2, 1)
  sd.process.C <- 1 / sqrt(tau.process.C)
  
  ## Overdispersion parameter
  size ~ dexp(0.1)  
})

################################################################################
#-------------------------------------------------------------------------------
# Gompertz model 
#-------------------------------------------------------------------------------

Gompertz <- nimbleCode({
  
  ## ------------------------------------------------------------
  ## Hierarchical spatial priors on Gompertz parameters
  ## ------------------------------------------------------------
  
  ## r > 0  (log scale) - intrinsic growth rate
  log_r_mean ~ dnorm(0, sd = 1)
  r_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.r, zero_mean = 1)
  for(i in 1:Nsite) {
    log_r[i] <- log_r_mean + r_raw[i]
    r[i]     <- exp(log_r[i])
  }
  
  ## K > 0  (log scale) - carrying capacity
  log_K_mean ~ dnorm(5, sd = 2)
  K_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.K, zero_mean = 1)
  for(i in 1:Nsite) {
    log_K[i] <- log_K_mean + K_raw[i]
    K[i]     <- exp(log_K[i])
  }
  
  ## CAR precisions (parameterized via SD for better mixing)
  sd.r ~ dexp(1)
  sd.K ~ dexp(1)

  tau.r <- 1 / sd.r^2
  tau.K <- 1 / sd.K^2

  ## ------------------------------------------------------------
  ## Site-specific Gompertz trajectories
  ## ------------------------------------------------------------
  for (i in 1:Nsite) {
    
    ## Initial state
    logC[i, first_year[i]] ~ dnorm(alpha0, tau.init.C)  

    for (j in (first_year[i] + 1):Nperiod) {
      
      ## Observation model: Negative Binomial 
      y[i, j] ~ dnbinom(size = size, prob = prob[i, j])
      prob[i, j] <- size / (size + lambda[i, j])
      log(lambda[i, j]) <- logC[i, j] + log_area_offset[i, j]

      # Process model
      ## log(N[t]) = log(N[t-1]) + r*(log(K) - log(N[t-1]))
      mu.logC[i, j] <- logC[i, j-1] + r[i] * (log_K[i] - logC[i, j-1])
      logC[i, j] ~ dnorm(mu.logC[i, j], tau.process.C)
    }
    
    ## Link to density
    for (j in first_year[i]:Nperiod) {
      log_density[i, j] <- logC[i,j]
      density[i, j] <- exp(log_density[i, j])
    }
  }
  
  ## ------------------------------------------------------------
  ## Hyperpriors and observation priors
  ## ------------------------------------------------------------
  
  ## Initial conditions
  alpha0     ~ dnorm(-4, sd = 2)       # small initial population
  tau.init.C ~ dgamma(2, 1)
  sd.init.C <- 1 / sqrt(tau.init.C)
  
  ## Process error for population dynamics
  tau.process.C ~ dgamma(2, 1)
  sd.process.C <- 1 / sqrt(tau.process.C)

})
