
# ------------------------------------------------------------
# Consumer-resource reaction-diffusion model
# dR/dt = r*R*(1 - R/K) - a*R*C           (resource)
# dC/dt = e*a*R*C - m*C                   (consumer)
#   a = per-capita resource ingestion rate 
#   e = assimilation efficiency
# Solved with deSolve (time) + ReacTran (space, 2D diffusion)
# ------------------------------------------------------------

library(deSolve)
library(ReacTran)

rm(list = ls())

# -----------------------------
# Parameters 
# -----------------------------

# hare parameters
parms <- list(
  r = 0.10,
  K = 100,
  a = 0.54,
  e = 0.03,
  m = 0.447,
  DC = 767            # km^2/year  (gives ~60 km/year wave speed)
)

# ---------------
# Spatial domain
# ---------------
nx <- 80; ny <- 80
cell_km <- 30          # km per grid cell
Lx <- nx * cell_km; Ly <- ny * cell_km   # domain: 2000 x 2000 km

xgrid <- setup.grid.1D(x.up = 0, L = Lx, N = nx)
ygrid <- setup.grid.1D(x.up = 0, L = Ly, N = ny)
grid  <- setup.grid.2D(x.grid = xgrid, y.grid = ygrid)

parms$nx <- nx; parms$ny <- ny; parms$nn <- nx * ny
parms$grid <- grid

# ---------------------
# Initial conditions
# ---------------------
# Resource at carrying capacity everywhere
R0 <- matrix(parms$K, nrow = ny, ncol = nx)

# Consumer invasion: single central cell at low density
cx <- round(nx/2); cy <- round(ny/2)
C_seed <- 0.001
C0 <- matrix(0, nrow = ny, ncol = nx)
C0[cy, cx] <- C_seed

y0 <- c(as.vector(R0), as.vector(C0))

# -------------
# Model 
# -------------
CR_2D <- function(t, y, parms) {
  with(parms, {
    R <- matrix(y[1:nn], nrow = ny, ncol = nx)
    C <- matrix(y[(nn + 1):(2 * nn)], nrow = ny, ncol = nx)

    intake <- a * R * C
    reacR  <- r * R * (1 - R / K) - intake
    reacC  <- e * intake - m * C

    # Consumer dispersal only — resource is spatially fixed
    diffC <- tran.2D(
      C = C, D.x = DC, D.y = DC, grid = grid,
      flux.x.up = 0, flux.x.down = 0, flux.y.up = 0, flux.y.down = 0
    )

    dR <- reacR
    dC <- diffC$dC + reacC

    list(c(as.vector(dR), as.vector(dC)))
  })
}

# -------------
# Integrate
# -------------
tmax  <- 30
times <- seq(0, tmax, by = 0.5)
system.time({
  out <- ode(y = y0, times = times, func = CR_2D, parms = parms, method = "lsoda")
})

# -------------
# Post-processing & plotting
# -------------
get_state <- function(out, i, nx, ny) {
  vec <- out[i, -1]
  R <- matrix(vec[1:(nx*ny)], nrow = ny, ncol = nx)
  C <- matrix(vec[(nx*ny + 1):(2*nx*ny)], nrow = ny, ncol = nx)
  list(R = R, C = C)
}

snap_times <- c(10, 15, 20, 25)
snap_idx   <- sapply(snap_times, function(tt) which.min(abs(times - tt)))

pal_orange <- grDevices::colorRampPalette(c("#fff5eb","#fdae6b","#7f2704"))

# Half-axis coordinates in km: centre (0 km) to right edge (Lx/2 km)
half_idx <- cx:nx
x_coords <- seq(0, Lx / 2, length.out = length(half_idx))

# Consistent y-axis ceiling across all 1D panels
C_max_snap <- max(sapply(snap_idx, function(i) {
  st <- get_state(out, i, nx, ny)
  max(st$C[cy, half_idx])
}))

#-------------------------------------------------------------------------------
op <- par(mfcol = c(3, length(snap_times)),
          mar   = c(6, 5, 2, 1),
          oma   = c(0, 0, 2, 0))

for (j in seq_along(snap_idx)) {
  st <- get_state(out, snap_idx[j], nx, ny)

  # Row 1: 2D spatial image — lower half only (semi-circle spreading downward)
  # use 1:ny for a full circle
  image(st$C[cy:ny, cx:nx], col = pal_orange(120), axes = FALSE, main = "")
    axis(1, at = seq(0, 1200, 200)/1200, labels = rep("", 7))
    if(j == 1) axis(2, at = seq(0, 1200, 200)/1200, labels = seq(0, 1200, 200))
    if(j > 1) axis(2, at = seq(0, 1200, 200)/1200, labels = rep("", 7))
    mtext(paste0("t = ", snap_times[j], " yr"), side = 3, line = 0.5, cex = 1.2, font = 2)
    if(j == 1) mtext("Distance from introduction (km)", side = 2, outer = TRUE, cex = 1, adj = 0.97, line = -1.5)
  
  # Row 2: 1D profile along central row, centre to right edge
  if(j == 1) {
  plot(x_coords, st$C[cy, half_idx],
       type = "l", lwd = 2, col = "#e6550d",
       xlim = c(0, Lx / 2), ylim = c(0, C_max_snap * 1.05),
       xlab = "", ylab = "",
       main = "", bty = "l", xaxs = "i", yaxs = "i")
  }

  if(j > 1) {
  plot(x_coords, st$C[cy, half_idx],
       type = "l", lwd = 2, col = "#e6550d",
       xlim = c(0, Lx / 2), ylim = c(0, C_max_snap * 1.05),
       xlab = "", ylab = "", yaxt = "n",
       main = "", bty = "l", xaxs = "i", yaxs = "i")
    axis(2, at = seq(0, 1.2, 0.2), labels = rep("", 7))
  }

  # Row 3: empty panel
  plot.new()
}

mtext("Distance from introduction (km)", side = 1, outer = TRUE, cex = 1, font = 1, line = -25)
mtext("Consumer density", side = 2, outer = TRUE, cex = 1, font = 1, line = -1.5)
par(op)
