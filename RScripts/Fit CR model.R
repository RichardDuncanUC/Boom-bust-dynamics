
library(tidyverse)
library(ozmaps)
library(sf)
library(spdep)
library(nimble)
library(coda)
library(Matrix)

#-------------------------------------------------------------------------------

rm(list = ls())
set.seed(1234)

# select either fox or hare 
mod.sp <- "fox"

# number of years prior to first sighting/bounty record to set as first_year
# prior to first_year, species are assumed absent from the district
buffer_years <- 2
if(mod.sp == "hare") buffer_years <- 7

# Model file
source("./RScripts/CR and Gompertz models.R")

#------------------------------------------------------------
# Read in hare and fox data
#------------------------------------------------------------

# read in data and filter to species
dat <- read.csv("./data/Fox hare data.csv") %>%
  filter(species == mod.sp) %>%
  glimpse()

# identify boards with at least one bounty record
all_board <- tapply(!is.na(dat$bounty_n), dat$board, sum)
in_board <- names(all_board)[all_board > 0]

# filter to those boards
dat <- dat |>
  filter(board %in% in_board) |>
  glimpse()

#-------------------------------------------------------------------------------
# read in location data for all boards (including those with no fox and hare data)
loc <- read.csv("./data/Board locations.csv") |>
  glimpse()

loc_all_boards <- loc # keep original for neighbourhood reconstruction

#------------------------------------------------------------
# Build Voronoi polygons for ALL boards
#------------------------------------------------------------

st_voronoi_point <- function(points) {
  if (!all(st_geometry_type(points) == "POINT"))
    stop("Input must be POINT geometries")
  g <- st_combine(st_geometry(points))
  v <- st_voronoi(g) %>%
    st_collection_extract()
  v[unlist(st_intersects(points, v))]
}

#------------------------------------------------------------
# Load NSW + ACT boundary using ozmaps and clip Voronoi
#------------------------------------------------------------

aus <- ozmaps::ozmap_states

nsw <- ozmaps::ozmap_states %>%
  filter(NAME %in% c("New South Wales", "Australian Capital Territory")) %>%
  st_union() %>%
  st_make_valid()

p <- st_as_sf(loc_all_boards, coords = c("x", "y"), crs = 4326) %>%
  st_transform(st_crs(nsw))

v <- st_voronoi_point(p)
v <- st_transform(v, st_crs(nsw))

vc <- st_intersection(v, nsw)
plot(vc)

# Calculate centroids
centroids <- st_centroid(vc)

# Extract coordinates from centroids
coords <- st_coordinates(centroids)

# Add text labels
text(coords[, "X"], coords[, "Y"], 
     labels = p$board,  
     cex = 0.6,         
     col = "black")

#------------------------------------------------------------
# Full neighbourhood structure (before filtering)
#------------------------------------------------------------
Wards_nb_full <- poly2nb(v)
nb_full <- nb2WB(Wards_nb_full)

#------------------------------------------------------------
# Filter loc, v to boards with at least one bounty record
#------------------------------------------------------------
keep_idx <- which(loc_all_boards$board %in% in_board)

loc <- loc_all_boards %>%
  filter(board %in% in_board)

v  <- v[keep_idx]
vc <- vc[keep_idx]
plot(vc)

# Add text labels
centroids <- st_centroid(vc)
coords <- st_coordinates(centroids)
text(coords[, "X"], coords[, "Y"], 
     labels = loc$board,  
     cex = 0.6,         
     col = "black")

#------------------------------------------------------------
# Build adjacency
#------------------------------------------------------------

old_to_new <- rep(NA_integer_, length(Wards_nb_full))
old_to_new[keep_idx] <- seq_along(keep_idx)

neighbors_list   <- Wards_nb_full[keep_idx]
mapped_neighbors <- lapply(neighbors_list, function(x) old_to_new[x])
mapped_neighbors <- lapply(mapped_neighbors, function(x) x[!is.na(x)])

new_adj     <- unlist(mapped_neighbors, use.names = FALSE)
new_num     <- lengths(mapped_neighbors)
new_weights <- rep(1, length(new_adj))

nb <- list(adj = new_adj, weights = new_weights, num = new_num)

#------------------------------------------------------------
# Species × year matrix for bounty counts
#------------------------------------------------------------
spp <- dat %>%
  dplyr::select(board, year, bounty_n) %>%
  pivot_wider(names_from = year,
              values_from = bounty_n,
              names_prefix = "y")

#------------------------------------------------------------
# Convert to matrix form
#------------------------------------------------------------
spp.mat <- spp %>%
  dplyr::select(-board) %>%
  as.matrix() %>%
  round()

Nsite   <- nrow(spp.mat)
Nperiod <- ncol(spp.mat)

Nsite
Nperiod

#------------------------------------------------------------
# Area matrix + log offset
#------------------------------------------------------------
area_vec <- loc$area
names(area_vec) <- loc$board

area_mat <- matrix(
  rep(area_vec[spp$board], Nperiod),
  nrow = Nsite, ncol = Nperiod
)

log_area_offset <- log(area_mat)

#-------------------------------------------------------------------------------
# first bounty record year in each board
first <- dat |>
  filter(!is.na(bounty_n)) |>
  group_by(board, first_newspaper) |>
  summarise(first_bounty = min(year)) |>
  ungroup() |>
  mutate(dif = first_newspaper - first_bounty,
         first_record = pmin(first_newspaper, first_bounty, na.rm = T)) |>
  glimpse()

# for foxes, difference between first newspaper sighting and first bounty record
mean(first$dif, na.rm = T)

# First year to start modelling arrival accounting for arrival date with buffer
first_year <- first$first_record - buffer_years

## Convert calendar years to period indices
first_year_index <- first_year - min(dat$year) + 1  
first_year_index <- ifelse(first_year_index < 1, 1, first_year_index)

#-------------------------------------------------------------------------------
# Nimble model inputs

myData <- list(
  y = spp.mat
)

myConstants <- list(
  Nsite   = Nsite,
  Nperiod = Nperiod,
  adj     = nb$adj,
  weights = nb$weights,
  num     = nb$num,
  first_year = first_year_index,
  log_area_offset = log_area_offset
)

## ------------------------------------------------------------
## Parameters to monitor
## ------------------------------------------------------------

myParam <- c(
  "alpha0",
  "size",
  "sd.process.C",
  "sd.init.C",

  ## hierarchical CR parameters
  "r_d", "K", "alpha", "eps", "mu", 

  ## hyperparameters
  "log_r_d_mean", "log_K_mean", "log_alpha_mean",
  "logit_eps_mean", "logit_mu_mean",

  ## CAR precisions
  "sd.r_d", "sd.K", "sd.alpha", "sd.eps", "sd.mu", 
  
  "density", "lambda", "R"
)

## ------------------------------------------------------------
## Initial values
## ------------------------------------------------------------
## Generate K_raw for calculating initial R values
K_raw <- rnorm(Nsite, 0, 0.1)

## Initialize R matrix
R_init <- matrix(NA, nrow = Nsite, ncol = Nperiod)
for(i in 1:Nsite) {
  ## Calculate K for this site based on K_raw
  K_i <- exp(K_raw[i])  # log_K_mean = 0
  
  ## Set R at colonization year and prior to K
  for(j in 1:first_year_index[i]) {
    R_init[i, j] <- K_i
  }

  ## Fill subsequent years with K (will be updated by model)
  for(j in (first_year_index[i] + 1):Nperiod) {
    R_init[i, j] <- K_i * 0.95  # Slightly below K
  }
}

## Initialize logC matrix
logC_init <- matrix(NA, nrow = Nsite, ncol = Nperiod)
for(i in 1:Nsite) {
  ## Set logC at colonization year and prior
  for(j in 1:first_year_index[i]) {
    logC_init[i, j] <- -6 
  }

  ## Fill subsequent years with slight increase
  for(j in (first_year_index[i] + 1):Nperiod) {
    logC_init[i, j] <- -6 + 0.1 * (j - first_year_index[i])
  }
}

## Complete initial values list
myInits <- list(
  ## CR hyper-means
  log_r_d_mean    = 1,
  log_alpha_mean  = -3,
  logit_eps_mean  = -1,
  logit_mu_mean    = -1.4,  
  
  ## CAR raw fields for parameters
  r_d_raw    = rnorm(Nsite, 0, 0.1),
  K_raw      = K_raw,
  alpha_raw  = rnorm(Nsite, 0, 0.1),
  eps_raw    = rnorm(Nsite, 0, 0.1),
  mu_raw     = rnorm(Nsite, 0, 0.1),
  
  ## observation model
  alpha0     = -8,
  tau.init.C = 1,
  size = 1.6,

  ## CAR precisions - parameterized as SDs
  sd.r_d   = 1,
  sd.K     = 1,
  sd.alpha = 1,
  sd.eps   = 1,
  sd.mu    = 1,
  
  ## Process error for consumer dynamics
  tau.process.C = 1,
  
  ## State variables
  R = R_init,
  logC = logC_init
)
#-------------------------------------------------------------------------------
# Fit model in Nimble

n.chains <- 3
n.iter   <- 410000
n.burn   <- 10000
n.thin   <- 20

myMod <- nimbleModel(code = CR,
                     constants = myConstants,
                     data = myData,
                     inits = myInits)

CmyMod <- compileNimble(myMod)

CmyModConf <- configureMCMC(CmyMod,
                            monitors = myParam,  
                            useConjugacy = FALSE,
                            enableWAIC = TRUE)

myMCMC <- buildMCMC(CmyModConf)
CmyMCMC <- compileNimble(myMCMC, project = myMod)

mod1 <- runMCMC(CmyMCMC,
                niter = n.iter,
                nburnin = n.burn,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE)

#===============================================================================
# Posterior summaries for ALL parameters
#===============================================================================

all_draws <- do.call(rbind, lapply(mod1$samples, as.matrix))

sum.out <- apply(all_draws, 2,
                 function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)) %>%
  t() %>%
  `colnames<-`(c("lcl", "lcl50", "med", "ucl50", "ucl"))

#===============================================================================
# Convergence diagnostics
#===============================================================================

params_to_check <- c(
  "alpha0",
  "size",
  "sd.process.C",
  "sd.init.C",

  ## hyperparameters
  "log_r_d_mean", "log_K_mean", "log_alpha_mean",
  "logit_eps_mean", "logit_mu_mean", 

  ## CAR precisions
  "sd.r_d", "sd.K", "sd.alpha", "sd.eps", "sd.mu")



mod1_subset <- lapply(mod1$samples, function(chain) chain[, params_to_check, drop = FALSE]) %>%
  as.mcmc.list()

rhat_vals <- gelman.diag(mod1_subset, multivariate = FALSE)$psrf[, 1]
ess_vals <- effectiveSize(mod1_subset)

outp <- sum.out[params_to_check, ] %>%
  as.data.frame() %>%
  mutate(var = rownames(.),
    rhat = rhat_vals[var],
    ess  = ess_vals[var])

outp

#===============================================================================
# Extract posterior median density (foxes/km²) and counts
#===============================================================================

density_medians <- sum.out[substr(rownames(sum.out), 1, 7) == "density", "med"]
count_medians <- sum.out[substr(rownames(sum.out), 1, 6) == "lambda", "med"]

results <- data.frame(board = rep(in_board, Nperiod),
                      year = rep(min(dat$year):max(dat$year), each = Nsite),
                      density_est = ifelse(density_medians > 0, density_medians, NA),
                      count_est = ifelse(count_medians > 0, count_medians, NA),
                      count_obs = as.vector(spp.mat)) |>
  glimpse()

# include bounty data
db <- dat %>%
  dplyr::select(board, year, bounty_price) |>
  arrange(year, board) |>
  glimpse()

results <- full_join(results, db) |>
  arrange(year, board) |>
  ungroup() |>
  glimpse()

dmat <- matrix(density_medians, nrow = Nsite, ncol = Nperiod, byrow = FALSE)
cmat <- matrix(count_medians, nrow = Nsite, ncol = Nperiod, byrow = FALSE)

#-------------------------------------------------------------------------------
# extract draws of the CR  parameters

# function to extract posterior draws
pd <- function(var, exact = FALSE) {
  lv <- nchar(var)
  temp <- all_draws[, substr(colnames(all_draws), 1, lv) == var]
  return(temp)
}

rd <- pd("r_d")
K <- pd("K")
alpha <- pd("alpha[")
eps <- pd("eps")
alpha0 <- pd("alpha0")
mu <- pd("mu[")

# Posterior array of fitted counts
lambda <- pd("lambda")
Ndraw <- dim(lambda)[1]
lambda_mat <- array(dim = c(Nsite, Nperiod, Ndraw))

# and fitted R
R_fitted <- pd("R[")  # Add this to your monitored params
R_mat <- array(dim = c(Nsite, Nperiod, Ndraw))

for(k in 1:Ndraw) {
  lambda_mat[, , k] <- matrix(lambda[k, ], nrow = Nsite, ncol = Nperiod, byrow = FALSE)
  lambda_mat[, , k][lambda_mat[, , k] == 0] <- NA
  R_mat[, , k] <- matrix(R_fitted[k, ], nrow = Nsite, ncol = Nperiod, byrow = FALSE)
}

#-------------------------------------------------------------------------------
# Calculate deterministic trajectory for each draw
# and one step ahead process residuals

# logC is the deterministic density
# logL is the deterministic count

R <- logC <- logL <- array(NA, c(Nsite, Nperiod, Ndraw))
# Process residuals: realized state vs deterministic expectation
process_residuals <- array(NA, dim = c(Nsite, Nperiod, Ndraw))

# loop through the posterior draws
for(k in 1:Ndraw) {

# Run deterministic dynamics
for (i in 1:Nsite) {
  
  logC[i, 1:(first_year_index[i]-1), k] <- NA
  logL[i, 1:(first_year_index[i]-1), k] <- NA
  R[i, 1:(first_year_index[i]-1), k] <- NA
  logC[i, first_year_index[i], k] <- log(dmat[i, first_year_index[i]])
  logL[i, first_year_index[i], k] <- log(cmat[i, first_year_index[i]])
  R[i, first_year_index[i], k] <- K[k, i]

  for (j in (first_year_index[i] + 1):Nperiod) {
    logC[i, j, k] <- logC[i, j-1, k] + 
                  log(1 - mu[k, i] + eps[k, i] * alpha[k, i] * R[i, j-1, k])
      
    R_increment <- rd[k, i] * R[i, j-1, k] * (1 - R[i, j-1, k] / K[k, i]) - 
                   alpha[k, i] * R[i, j-1, k] * exp(logC[i, j-1, k])
      
    R[i, j, k] <- max(R[i, j-1, k] + R_increment, 0.001)
    
    # Calculate predicted log count 
    logL[i, j, k] <- logC[i, j, k] + log_area_offset[i, j]
    
    # OSA residuals on counts
    # Realized consumer state
    logC_realized <- log(lambda_mat[i, j, k])
    # Deterministic expectation from process model
    logC_prev <- log(lambda_mat[i, j-1, k])
    mu_i <- mu[k, i]
    eps_i <- eps[k, i]
    alpha_i <- alpha[k, i]
    R_prev <- R_mat[i, j-1, k]
    mu_logC <- logC_prev + log(1 - mu_i + eps_i * alpha_i * R_prev)
    # OSA process residual 
    process_residuals[i, j, k] <- logC_realized - mu_logC
  }
}}

# Summarise
R_med <- apply(R, c(1, 2), median)
logL_med <- apply(logL, c(1, 2), median)
logC_med <- apply(logC, c(1, 2), median)
resid_med <- apply(process_residuals, c(1,2), median, na.rm=TRUE)

# add to results file
results$density_det <- as.vector(exp(logC_med))
results$count_det <- as.vector(exp(logL_med))
results$R_det <- as.vector(R_med)
results$osa_resid <- as.vector(resid_med)

results <- results |>
  full_join(loc) |>
  arrange(board, year) |>
  glimpse()

#===============================================================================
# Save results
#===============================================================================
# CR parameters
cr_param <- data.frame(rd = apply(rd, 2, median, na.rm = T),
                       K = apply(K, 2, median, na.rm = T),
                       alpha = apply(alpha, 2, median, na.rm = T),
                       eps = apply(eps, 2, median, na.rm = T),
                       mu = apply(mu, 2, median, na.rm = T))

out <- list()
out$results <- results
out$param <- outp
out$cr_param <- cr_param
out$mod <- c(n.chains = n.chains,
             n.iter   = n.iter,
             n.burn   = n.burn,
             n.thin   = n.thin)

assign(paste0(mod.sp, "_cr"), out)

save(list = c(paste0(mod.sp, "_cr")),
     file = "./model output/", mod.sp, " CR results.RData")

