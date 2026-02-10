##### Here we are extracting Sea Ice Coverage (SIC) from https://figshare.com/collections/_/5709767 paper - https://www.nature.com/articles/s41558-021-01254-9

# --- Seasonal SIC at sites using geodesic-nearest valid grid cell (per season) ---
library(ncdf4)
library(geosphere)
library(ggOceanMaps)
library(ggspatial)
library(sf)

sites <- read.csv("metadata/sitemetadata.csv")  # columns: Site, Lat, Lon

wrap180 <- function(x) ifelse(x > 180, x - 360, x)

season_files <- c(
  DJF = "rawdata/climate/SeaIce/conc_recon_ens_DJF.nc",
  MAM = "rawdata/climate/SeaIce/conc_recon_ens_MAM.nc",
  JJA = "rawdata/climate/SeaIce/conc_recon_ens_JJA.nc",
  SON = "rawdata/climate/SeaIce/conc_recon_ens_SON.nc"
)

# output holders
out_list <- list()
nearest_all <- list()

# plot settings
lims <- c(-70, -56, -70, -60)
shape_vals <- c(DJF=16, MAM=17, JJA=15, SON=3)

for (seas in names(season_files)) {
  
  f <- season_files[[seas]]
  
  nc <- nc_open(f)
  latg <- ncvar_get(nc, "latitutde")
  long <- wrap180(ncvar_get(nc, "longitude"))
  sic  <- ncvar_get(nc, "ensmean")                # [time, lon, lat]
  fill <- ncatt_get(nc, "ensmean", "_FillValue")$value
  nc_close(nc)
  
  # valid cells: any year has real value (zeros are valid)
  valid <- apply(sic, c(2,3), function(v) any(is.finite(v) & v != fill))
  valid_vec <- as.vector(valid)
  
  # flattened valid grid for fast geodesic nearest
  grid_lon <- as.vector(long)
  grid_lat <- as.vector(latg)
  grid_xy_valid <- cbind(grid_lon[valid_vec], grid_lat[valid_vec])
  valid_idx <- which(valid_vec)
  
  nlon <- dim(long)[1]
  flat_to_ij <- function(k) c(i_lon = ((k - 1) %% nlon) + 1, i_lat = ((k - 1) %/% nlon) + 1)
  
  nt <- dim(sic)[1]
  years <- 1905:(1905 + nt - 1)
  
  nearest_df <- data.frame(Site = sites$Site, Season = seas,
                           SiteLon = sites$Lon, SiteLat = sites$Lat,
                           NearestLon = NA_real_, NearestLat = NA_real_,
                           NearestDist_km = NA_real_)
  
  for (s in seq_len(nrow(sites))) {
    p <- c(sites$Lon[s], sites$Lat[s])
    
    d_m <- geosphere::distHaversine(grid_xy_valid, p)
    k_valid <- which.min(d_m)
    k_full <- valid_idx[k_valid]
    ij <- flat_to_ij(k_full)
    i_lon <- ij["i_lon"]; i_lat <- ij["i_lat"]
    
    nearest_df$NearestLon[s] <- long[i_lon, i_lat]
    nearest_df$NearestLat[s] <- latg[i_lon, i_lat]
    nearest_df$NearestDist_km[s] <- d_m[k_valid] / 1000
    
    # extract SIC time series for this season/site
    ts <- sic[, i_lon, i_lat]
    out_list[[length(out_list) + 1]] <- data.frame(
      Site = sites$Site[s],
      Season = seas,
      Year = years,
      SIC = as.numeric(ts),
      NearestLon = long[i_lon, i_lat],
      NearestLat = latg[i_lon, i_lat],
      NearestDist_km = d_m[k_valid] / 1000
    )
  }
  
  nearest_all[[seas]] <- nearest_df
}

sic_sites <- do.call(rbind, out_list)
nearest_points <- do.call(rbind, nearest_all)

# write output table

write.csv(sic_sites, "rawdata/climate/SeaIce/site_season_year_SIC_nearest.csv", row.names = FALSE)
# build sf lines: one per row of nearest_points
lines_sfc <- st_sfc(lapply(seq_len(nrow(nearest_points)), function(i) {
  st_linestring(matrix(c(nearest_points$SiteLon[i],    nearest_points$SiteLat[i],
                         nearest_points$NearestLon[i], nearest_points$NearestLat[i]),
                       ncol = 2, byrow = TRUE))
}), crs = 4326)

lines_sf <- st_sf(Season = nearest_points$Season, geometry = lines_sfc)

pdf("mapping/sites_and_nearest_by_season.pdf", width = 8, height = 12)

p <- basemap(limits = lims, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  
  geom_sf(data = lines_sf, aes(colour = Season),
          linewidth = 0.6, alpha = 0.8, inherit.aes = FALSE) +
  
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),
                                data = sites, color = "black",
                                cex = 4.2, crs = 4326) +
  
  ggspatial::geom_spatial_point(aes(x = NearestLon, y = NearestLat, shape = Season, colour = Season),
                                data = nearest_points, cex = 3.2, alpha = 0.9, crs = 4326) +
  
  scale_shape_manual(values = shape_vals) +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm"))

print(p)
dev.off()


## Here is a chat made plot


## Base R faceted time series: points (transparent) + season-wise loess per site
## Assumes data.frame sic_sites with columns: Site, Season, Year, SIC
stopifnot(all(c("Site","Season","Year","SIC") %in% names(sic_sites)))

d <- sic_sites
d <- d[is.finite(d$Year) & is.finite(d$SIC) & !is.na(d$Site) & !is.na(d$Season), ]

d$Site   <- factor(d$Site)
d$Season <- factor(d$Season, levels = c("DJF","MAM","JJA","SON"))

season_cols <- c(
  DJF = "#fee8c9",
  MAM = "#fdbb84",
  SON = "#74a9cf",
  JJA = "#045a8d"
)

pt_alpha <- 0.20

sites <- levels(d$Site)
ns <- length(sites)

nr <- floor(sqrt(ns + 1))
nc <- ceiling((ns + 1) / nr)

xlim_all <- range(d$Year, na.rm = TRUE)
ylim_all <- range(d$SIC,  na.rm = TRUE)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

pdf("rawdata/climate/SeaIce/PerSiteSeasons.pdf",width = 8,height = 6)

par(mfrow = c(nr, nc), mar = c(3.6, 3.8, 2.2, 0.8), xaxs = "i", yaxs = "i")

for (i in seq_len(ns + 1)) {
  
  if (i <= ns) {
    s  <- sites[i]
    ds <- d[d$Site == s, ]
    
    plot(NA, xlim = xlim_all, ylim = ylim_all,
         xlab = "", ylab = "", main = paste("Site", s), las = 1)
    
    for (sea in levels(d$Season)) {
      dss <- ds[ds$Season == sea, ]
      if (nrow(dss) == 0) next
      
      points(dss$Year, dss$SIC,
             pch = 16, cex = 0.7,
             col = adjustcolor(season_cols[sea], alpha.f = pt_alpha))
      
      ok <- is.finite(dss$Year) & is.finite(dss$SIC)
      if (sum(ok) >= 10) {
        lo <- loess(SIC ~ Year, data = dss[ok, ], span = 0.35, degree = 2,
                    control = loess.control(surface = "direct"))
        xx <- seq(min(dss$Year[ok]), max(dss$Year[ok]), length.out = 200)
        yy <- predict(lo, newdata = data.frame(Year = xx))
        lines(xx, yy, col = season_cols[sea], lwd = 2)
      }
    }
    
    box(col = "grey40")
    
  } else {
    plot.new()
    legend("center", bty = "n",
           legend = levels(d$Season),
           col = season_cols[levels(d$Season)],
           lwd = 2, pch = 16, pt.cex = 0.9)
  }
}
dev.off()

