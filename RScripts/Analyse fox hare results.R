
library(tidyverse)
library(sf)
library(ozmaps)
library(fields)

rm(list = ls())

unlogit <- function(x) exp(x) / (1 + exp(x))

# colours
fox_color <- rgb(0.871, 0.176, 0.149, 0.7)
hare_color <- rgb(0.192, 0.510, 0.741, 0.7) 

# read in model outputs
load("./model output/fox CR results.RData")
load("./model output/fox gompertz results.RData")
load("./model output/hare CR results.RData")
load("./model output/hare gompertz results.RData")

# Hare with fox predation
load("./model output/hare CR joint fox results.RData")

# Number of boards
length(table(fox_cr$results$board))
length(table(hare_cr$results$board))

# first bounty records
min(fox_cr$results$year[fox_cr$results$count_obs > 0], na.rm = T)
min(hare_cr$results$year[hare_cr$results$count_obs > 0], na.rm = T)

#-------------------------------------------------------------------------------
# Returns per year by board (Figure S1)

fr <- fox_cr$results |>
  filter(!is.na(count_obs)) |>
  dplyr::select(board, year, y) |>
  unique() |>
  mutate(spp = "Fox") |>
  glimpse()

hr <- hare_cr$results |>
  filter(!is.na(count_obs)) |>
  dplyr::select(board, year, y) |>
  unique() |>
  mutate(spp = "Hare") |>
  glimpse()

fhr <- bind_rows(fr, hr) |>
  mutate(board= fct_reorder(board, y))

ggplot(fhr, aes(y = board, x = year)) +
  geom_point() +
  ylab("Board") +
  xlab("Year") +
  facet_wrap(~ spp, ncol = 2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12))

#-------------------------------------------------------------------------------
# Model comparison (Figure S2)
# one step ahead process residuals

fox_cr_resid <- fox_cr$results$osa_resid
fox_gompertz_resid <- fox_gompertz$results$osa_resid
hare_cr_resid <- hare_cr$results$osa_resid
hare_gompertz_resid <- hare_gompertz$results$osa_resid

# calculate mean squared OSA residuals
sfcr <- round(mean(fox_cr_resid^2, na.rm = T), 3)
sfgz <- round(mean(fox_gompertz_resid^2, na.rm = T), 1)
shcr <- round(mean(hare_cr_resid^2, na.rm = T), 3)
shgz <- round(mean(hare_gompertz_resid^2, na.rm = T), 1)

# plot OSA residuals
par(mfrow = c(3, 2), mar = c(4, 4, 2, 2))
hist(fox_gompertz_resid, col = fox_color, main = paste0("Fox Gompertz: MS = ", sfgz),
     xlab = "Residual", breaks = 15)
mtext("A", adj = 0, cex = 1.5)

hist(hare_gompertz_resid, col = hare_color, main = paste0("Hare Gompertz: MS = ", shgz),
     xlab = "Residual", breaks = 15)
mtext("B", adj = 0, cex = 1.5)

hist(fox_cr_resid, col = fox_color, main = paste0("Consumer-resource: MS = ", sfcr),
     xlab = "Residual", breaks = 15)
mtext("C", adj = 0, cex = 1.5)

hist(hare_cr_resid, col = hare_color, main = paste0("Consumer-resource: MS = ", shcr),
     xlab = "Residual", breaks = 15)
mtext("D", adj = 0, cex = 1.5)

# plot observed vs predicted
plot(fox_cr$results$count_obs ~ fox_cr$results$count_est, pch = 19, col = fox_color,
     log = "xy", xlim = c(1, 10000), bty = "l",
     ylab = "Observed", xlab = "Predicted")
  abline(0, 1)
mtext("E", adj = 0, cex = 1.5)

plot(hare_cr$results$count_obs ~ hare_cr$results$count_est, pch = 19, col = hare_color,
     log = "xy", xlim = c(10, 100000), bty = "l",
     ylab = "Observed", xlab = "Predicted")
abline(0, 1)
mtext("F", adj = 0, cex = 1.5)

#-------------------------------------------------------------------------------
# plot the fit of consumer-resource models to the data (Figures S3 & 4)

fox_cr$results$board = fct_reorder(fox_cr$results$board, fox_cr$results$y)

ggplot(fox_cr$results, aes(y = count_obs/area, x = year)) +
  geom_point() +
  geom_line(aes(y = density_est), col = fox_color, lwd = 1.1) +
  ylab("Harvest density (km-2)") +
  xlab("Year") +
  facet_wrap(~ board) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


hare_cr$results$board = fct_reorder(hare_cr$results$board, hare_cr$results$y)

ggplot(hare_cr$results, aes(y = count_obs/area, x = year)) +
  geom_point() +
  geom_line(aes(y = density_est), col = hare_color, lwd = 1.1) +
  ylab("Harvest density (km-2)") +
  xlab("Year") +
  facet_wrap(~ board) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-------------------------------------------------------------------------------
# compare hare CR with CR adding fox predation (Figure S5)

hare_cr_resid <- hare_cr$results$osa_resid
hare_cr_fox_resid <- hare_cr_joint_fox$results$osa_resid

# calculate mean squared residuals
shcr <- round(mean(hare_cr_resid^2, na.rm = T), 4)
shcrf <- round(mean(hare_cr_fox_resid^2, na.rm = T), 4)

par(mfrow = c(2, 2))
hist(hare_cr_resid, col = hare_color, main = paste0("Hare CR: MS = ", shcr),
     xlab = "Residual", breaks = 15)
mtext("A", adj = 0, cex = 1.5)

hist(hare_cr_fox_resid, col = hare_color, main = paste0("Hare CR with fox predation: MS = ", shcrf),
     xlab = "Residual", breaks = 15)
mtext("B", adj = 0, cex = 1.5)

#-------------------------------------------------------------------------------
# Table S1 summarising posterior draws of selected parameters

# names of parameters to include in table
sel_var <- c("log_r_d_mean", "log_alpha_mean", "logit_eps_mean", "logit_mu_mean",
             "sd.r_d", "sd.K", "sd.alpha", "sd.eps", "sd.mu",
             "size", "sd.process.C")

# function to produce table
tab <- function(fp) {
  fp <- fp |>
    dplyr::select(var, med, lcl, ucl, rhat, ess) |>
    filter(var %in% sel_var) |>
    arrange(match(var, sel_var)) 

  # transform to raw scale
  fp[fp$var == "log_r_d_mean", c("med", "lcl", "ucl")] <- exp(fp[fp$var == "log_r_d_mean", c("med", "lcl", "ucl")])
  fp[fp$var == "log_alpha_mean", c("med", "lcl", "ucl")] <- exp(fp[fp$var == "log_alpha_mean", c("med", "lcl", "ucl")])
  fp[fp$var == "logit_eps_mean", c("med", "lcl", "ucl")] <- unlogit(fp[fp$var == "logit_eps_mean", c("med", "lcl", "ucl")])
  fp[fp$var == "logit_mu_mean", c("med", "lcl", "ucl")] <- unlogit(fp[fp$var == "logit_mu_mean", c("med", "lcl", "ucl")])

  fp <- fp |>
    mutate(med = round(med, 3),
           lcl = round(lcl, 3), 
           ucl = round(ucl, 3),
           rhat = round(rhat, 4),
           ess = round(ess, 0))
  return(fp)
}

# create table and save to csv file
fout <- tab(fox_cr$param)
hout <- tab(hare_cr$param)

fout$spp <- "fox"
hout$spp <- "hare"
out <- bind_rows(fout, hout)
write.csv(out, "./model output/Output parameters.csv", row.names = F)

################################################################################
# Draw Figures 2 and 3 in main text
#-------------------------------------------------------------------------------
# Get NSW boundary and southern border
aus <- ozmaps::ozmap_states

nsw <- ozmaps::ozmap_states %>%
  filter(NAME %in% c("New South Wales", "Australian Capital Territory")) %>%
  st_union() %>%
  st_make_valid()

# get NSW-Vic border for plotting Boards ordered from south-north
# get NSW southern border
vic <- ozmaps::ozmap_states %>% 
  filter(NAME == "Victoria") %>% 
  st_union() %>% 
  st_make_valid()

nsw_border <- st_boundary(nsw)
vic_border <- st_boundary(vic)

# The NSW–Victoria border as a LINESTRING 
nsw_vic_border <- st_intersection(nsw_border, vic_border)
nsw_vic_border_proj <- st_transform(nsw_vic_border, 4283)

################################################################################
# Figure 2

# years of records for each species
fr <- fox_cr$results |>
  filter(!is.na(count_obs) & count_obs > 0) |>
  group_by(board, x, y) |>
  summarise(nyear = n()) |>
  glimpse()

hr <- hare_cr$results |>
  filter(!is.na(count_obs) & count_obs > 0) |>
  group_by(board, x, y) |>
  summarise(nyear = n()) |>
  glimpse()

# Save original par settings
op <- par(no.readonly = TRUE)

# Create plot layout
par(mfrow = c(2, 2))

# Top left: Fox
plot(nsw, main = "Fox")
points(fr$x, fr$y, pch = 19, col = fox_color, cex = fr$nyear/6)
mtext("B", adj = 0, cex = 1.5)

# Top right: Hare
plot(nsw, main = "Hare")
points(hr$x, hr$y, pch = 19, col = hare_color, cex = hr$nyear/6)
mtext("C", adj = 0, cex = 1.5)
ps <- c(1, 10, 20, 30) / 6
legend(x = 137, y = -29, 
       legend = c(1, 10, 20, 30), pch = 19, col = "grey", pt.cex = ps, xpd = NA, 
       y.intersp = 1.8, x.intersp = 1.5, bty = "n",
       title = "Years")

# Function to plot density offtake over time
dens.first <- function(dat, mt) {
  
  first.rec <- dat %>%
    mutate(obs_density = count_obs / area) |>
    arrange(board, year) %>%
    group_by(board) %>%
    mutate(first_positive_year = min(year[count_obs > 0], na.rm = TRUE),
      years_since_first   = if_else(year >= first_positive_year,
                                    year - first_positive_year,
                                    NA_real_)) %>%
    ungroup()
  
  dens.mean <- first.rec |>
    group_by(years_since_first) |>
    summarise(mean.dens = mean(obs_density, na.rm = T)) |>
    filter(!is.nan(mean.dens) & !is.na(years_since_first))

  dat_name <- deparse(substitute(dat))
  point_color <- ifelse(dat_name == "fox_cr$results", fox_color, hare_color)
  
  plot(obs_density ~ years_since_first, data = first.rec, pch = 19, bty = "l",
       col = rgb(0.7, 0.7, 0.7, 0.5), 
       ylab = "Harvest density (km-2)", xlab = "Years since first bounty",
       xlim = c(0, 35))
    points(mean.dens ~ years_since_first, data = dens.mean, pch = 19,
           col = point_color, cex = 1.5)
    fit <- smooth.spline(dens.mean$mean.dens ~ dens.mean$years_since_first,
                         spar = 0.5)
    lines(fit, col = point_color, lwd = 2)
    mtext(mt, adj = 0, cex = 1.5)
}

# Plot fox and hare density offtake over time
par(mar = c(5, 4, 4, 2))
dens.first(fox_cr$results, "D")
dens.first(hare_cr$results, "")

# Now overlay the Australia map in the center
par(fig = c(0.375, 0.625, 0.45, 0.70), new = TRUE)
par(mar = c(0, 0, 0, 0))
plot(st_geometry(aus), col = NA, border = "black", lwd = 1.5)
plot(st_geometry(nsw), col = "grey", border = "black", add = TRUE, lwd = 1.5)

# add point for were first releases occurred
points(144.96, -38.2, pch = 19, col = "red", cex = 1.5)

# Get bounding box of Australia
aus_mainland <- aus %>%
  st_crop(xmin = 113, xmax = 154, ymin = -44, ymax = -10)
bbox_aus <- st_bbox(aus_mainland)

# Extract coordinates
xmin <- bbox_aus["xmin"]
xmax <- bbox_aus["xmax"]
ymin <- bbox_aus["ymin"]
ymax <- bbox_aus["ymax"]

# Draw rectangle just outside the Australia bounds
# Add small buffer 
x_buffer <- 0.1 * (xmax - xmin)
y_buffer <- 0.1 * (ymax - ymin)

rect(xmin - x_buffer, ymin - y_buffer, 
     xmax + x_buffer, ymax + y_buffer, 
     border = "black", lwd = 1.5)
mtext("A", adj = 0.15, line = -5.5, cex = 1.5)

# NOW add the "D" after the overlay
par(fig = c(0.5, 1, 0, 0.5), new = TRUE)  # Return to bottom-right panel
par(mar = c(5, 4, 4, 2))
mtext("E", adj = 0, cex = 1.5)

# Restore original par settings
par(op)

################################################################################
# Figure 3

# function to plot a group of trajectories
plot_group <- function(dat, gp_num, xlm = "", ylm = "") {
  dat.gp <- filter(dat, gp == gp_num)
  br <- levels(factor(dat.gp$board))
  temp <- filter(dat.gp, board == br[1])
  plot(density_est ~ year, data = temp, type = "l", col = cols[gp_num], 
       ylim = yl, bty = "l", ylab = ylm, xlab = xlm, xlim = xl, lwd = 2)
  for(j in 2:length(br)) {
    temp <- filter(dat.gp, board == br[j])
    lines(density_est ~ year, data = temp, col = cols[gp_num], lwd = 2)
  }
}
  
#-------------------------------------------------------------------------------
# Board locations for fox
pts_proj <- st_as_sf(
  data.frame(lon = fox_cr$results$x, lat = fox_cr$results$y),
  coords = c("lon", "lat"),
  crs = 4283   
)

# distance of each Board from NSW-Vic border lumped into 5 groups
dist <- as.numeric(st_distance(pts_proj, nsw_vic_border_proj)[, 1])
fox_cr$results$gp <- as.numeric(cut_number(dist, 5))

# unique locations
fr <- fox_cr$results |>
  group_by(board, x, y, gp) |>
  summarise(max_pred = max(density_est, na.rm = TRUE),
            ny = n(), .groups = "drop")

# plot layout
op <- par(no.readonly = TRUE)
layout(matrix(c(1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15), ncol = 3, byrow = F),
       height = c(1, 1, 1, 1, 1, 0.1, 2, 0.1, 2, 0.1, 1, 1, 1, 1, 1))
par(mar = c(4, 4, 1, 0))

# fox colours
cols <- colorRampPalette(c("#fc9272", "#de2d26", "#a50f15"))(5)

yl <- c(0, max(fox_cr$results$density_est, na.rm = T))
xl <- c(1893, 1930)

# plot fox trajectories in the five groups
plot_group(fox_cr$results, 5) 
plot_group(fox_cr$results, 4) 
plot_group(fox_cr$results, 3, ylm = "Bounty offtake density (km-2)"); 
plot_group(fox_cr$results, 2) 
plot_group(fox_cr$results, 1, xlm = "Year")

plot.new()

# plot max density of each Board on map of NSW
par(mar = c(0, 0, 0, 0))
plot(nsw)
mtext("Fox", side = 3, cex = 1.2, font = 2)
points(fr$x, fr$y, pch = 19, col = cols[fr$gp], cex = log10(fr$max_pred + 1.3)*14)

legend_vals <- c(0.05, 0.2, 0.4)  
legend_sizes <- log10(legend_vals + 1.3) * 14

legend(x = 141, y = -37,
       legend = legend_vals,
       pch = 19, col = "grey",
       pt.cex = legend_sizes,
       title = "Maximum density (km-2)",
       cex = 1.2,
       bty = "n",  # no box
       y.intersp = 1.8,  # spacing between items
       x.intersp = 1.5, xpd = NA)

arrows(140.5, -32, 138, -32, lwd = 4, length = 0.25, angle = 20, xpd = NA)

plot.new()

#-------------------------------------------------------------------------------
# same again for hare
pts_proj <- st_as_sf(
  data.frame(lon = hare_cr$results$x, lat = hare_cr$results$y),
  coords = c("lon", "lat"),
  crs = 4283   
)

dist <- as.numeric(st_distance(pts_proj, nsw_vic_border_proj)[, 1])
hare_cr$results$gp <- as.numeric(cut_number(dist, 5))

# unique locations
hr <- hare_cr$results |>
  group_by(board, x, y, gp) |>
  summarise(max_pred = max(density_est, na.rm = TRUE),
            ny = n(), .groups = "drop")

cols <- colorRampPalette(c("#9ecae1", "#3182bd", "#08519c"))(5)

plot(nsw)
mtext("Hare", side = 3, cex = 1.2, font = 2)
points(hr$x, hr$y, pch = 19, col = cols[hr$gp], cex = log10(hr$max_pred + 1.3)*3)

legend_vals <- c(1, 5, 10)  # Adjust based on your actual range
legend_sizes <- log10(legend_vals + 1.3) * 3

legend(x = 141, y = -37,
       legend = legend_vals,
       pch = 19, col = "grey",
       pt.cex = legend_sizes,
       title = "Maximum density (km-2)",
       cex = 1.2,
       bty = "n",  # no box
       y.intersp = 1.8,  # spacing between items
       x.intersp = 1.5, xpd = NA)

arrows(153.5, -32, 156, -32, lwd = 4, length = 0.25, angle = 20, xpd = NA)

plot.new()

par(mar = c(4, 3, 1, 1))
yl <- c(0, max(hare_cr$results$density_est, na.rm = T))
xl <- c(1883, 1930)

plot_group(hare_cr$results, 5) 
plot_group(hare_cr$results, 4) 
plot_group(hare_cr$results, 3); 
plot_group(hare_cr$results, 2) 
plot_group(hare_cr$results, 1, xlm = "Year")

par(op)

#-------------------------------------------------------------------------------
# Figure S6
# Function to plot max density over NSW

fox_cols <- colorRampPalette(c("#fff5f0", "#fcbba1", "#fc9272", "#de2d26", "#a50f15", "#67000d"))(64)
hare_cols <- colorRampPalette(c("#f7fbff", "#c6dbef", "#9ecae1", "#3182bd", "#08519c", "#08306b"))(64)

# Function to plot maximum density surface
maxdens <- function(dat, cl, title) {
  # Fit thin plate spline
  fit <- Tps(cbind(dat$x, dat$y), dat$max_pred)

  # Create prediction grid within NSW bounding box
  bb <- st_bbox(nsw)
  grid_lon <- seq(bb["xmin"], bb["xmax"], length.out = 200)
  grid_lat <- seq(bb["ymin"], bb["ymax"], length.out = 200)
  grid_xy <- expand.grid(lon = grid_lon, lat = grid_lat)

  # Predict surface
  grid_xy$z <- predict(fit, x = as.matrix(grid_xy[, c("lon", "lat")]))

  # Mask to NSW: keep only grid cells inside the boundary
  grid_sf <- st_as_sf(grid_xy, coords = c("lon", "lat"), crs = 4283)
  inside <- st_intersects(grid_sf, nsw, sparse = FALSE)[, 1]
  grid_xy$z[!inside] <- NA

  # Reshape to matrix for image()
  z_mat <- matrix(grid_xy$z, nrow = length(grid_lon), ncol = length(grid_lat))

  # Plot
  image.plot(grid_lon, grid_lat, z_mat,
             col = cl,
             xlab = "Longitude", ylab = "Latitude",
             main = title)
  plot(st_geometry(nsw), add = TRUE, border = "black", lwd = 1.5)
#  points(fr$x, fr$y, pch = 16, cex = 0.4, col = "grey20")
}

par(mfrow = c(2, 2), mar = c(5, 4, 4, 4))
maxdens(fr, fox_cols, "Fox maximum density")
maxdens(hr, hare_cols, "Hare maximum density")

