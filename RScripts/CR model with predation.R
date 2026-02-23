
# Nimble code for hare consumer-resource model with fox predation

CR_joint_fox <- nimbleCode({
  
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
  
  ## alpha > 0  (log scale) - attack rate on resource
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
  
  ## mu in (0,1)  (logit scale) - baseline mortality rate 
  logit_mu_mean ~ dnorm(0, sd = 1)
  mu_raw[1:Nsite] ~ dcar_normal(adj[], weights[], num[], tau.mu, zero_mean = 1)
  for(i in 1:Nsite) {
    logit_mu[i] <- logit_mu_mean + mu_raw[i]
    mu[i]       <- ilogit(logit_mu[i])
  }
  
  ## ------------------------------------------------------------
  ## Fox predation parameters (Type II functional response)
  ## ------------------------------------------------------------
  ## Attack rate: predation efficiency per fox per hare
  log_attack_rate ~ dnorm(0, sd = 1)
  attack_rate <- exp(log_attack_rate)
  
  ## Handling time: time to kill and consume one hare
  log_handling_time ~ dnorm(0, sd = 1)
  handling_time <- exp(log_handling_time)
  
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
      
      ## Type II functional response for fox predation
      ## Per-capita predation rate saturates with hare density
      ## Total predation = per-capita rate × fox density
      hare_density_prev[i, j-1] <- exp(logC[i, j-1])
      per_capita_kill_rate[i, j-1] <- (attack_rate * hare_density_prev[i, j-1]) / 
                                       (1 + attack_rate * handling_time * hare_density_prev[i, j-1])
      
      ## Total fox predation mortality
      ## = per-capita kills × fox density / hare density
      ## Automatically = 0 when fox_density = 0
      fox_mortality[i, j-1] <- (per_capita_kill_rate[i, j-1] * fox_density[i, j-1]) / 
                               hare_density_prev[i, j-1]
      
      ## Total mortality (baseline + predation)
      total_mortality[i, j-1] <- mu[i] + fox_mortality[i, j-1]
      mortality_rate[i, j-1] <- min(total_mortality[i, j-1], 0.99)
      
      ## Consumer update (log scale) with process error 
      mu.logC[i, j] <- logC[i, j-1] + 
                       log(1 - mortality_rate[i, j-1] + eps[i] * alpha[i] * R[i, j-1])
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


