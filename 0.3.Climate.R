########################################################
####=.     Analysis of WP5  Antarctic climate      =####
####==== Luke E. Holman====25.08.2025====################
########################################################

####====0.0 Packages====####

# For SIC 
library(ncdf4)
library(geosphere)
library(ggOceanMaps)
library(ggspatial)
library(sf)


####====0.1 Functions====####

# month name/abbr/number -> integer 1..12
month_num <- function(x){
  if(is.numeric(x)) return(pmin(pmax(as.integer(x),1),12))
  x  <- as.character(x)
  nx <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(nx), nx,
                match(toupper(substr(x,1,3)), toupper(month.abb)))
  as.integer(out)
}

# continuous year, centered on a month (default July=7)
make_continuous_year <- function(year, month, centre=7){
  m  <- month_num(month)
  c0 <- month_num(centre)[1]
  year + (m - c0)/12
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



####====0.2 Data in====####
#read in metadata 
metadat <- read.csv("metadata/RepMetaData.csv")
dates <- read.csv("metadata/dates.csv")
coredata <- read.csv("metadata/sitemetadata.csv")
coredata$col <- c("#D55E00","#cc79a7","#0072B2","#f0e442","#009e73")
climateout <- dates[,c(1,2,7)]
minval <- 1770




####====1.0 Southern Annular Mode data ====####

## read in data
model <- read.table("rawdata/climate/SAM_1000yr.txt",header = 1)
obs <- read.table("rawdata/climate/BES_SAMannual.txt",header = TRUE,  comment.char = "#",fill = TRUE,strip.white = TRUE)

## make shared dataframe of observational data 
shared <- data.frame("year"=intersect(model$age_CE,obs$year),
                     "modelObs"=model$Marshall_SAM_Index[match(intersect(model$age_CE,obs$year),model$age_CE)],
                     "Obs"=obs$ANN[match(intersect(model$age_CE,obs$year),obs$year)])

cor(shared$modelObs,shared$Obs)

shift <- mean(shared$modelObs) - mean(shared$Obs)
shift2 <- mean(shared$Obs) - mean(shared$modelObs)


obs$adjusted <- obs$ANN+shift




## Now we make a combined dataset to use - this uses observations adjusted to the normalisation used in Abrams from 1957 -> 2025

combinedSAM <- data.frame("year"=sort(unique(c(obs$year,model$age_CE)),decreasing = TRUE),
                          #"SAM_Index_OLD"=c(rev(obs$adjusted),model$SAM_Index[match(1956:1000,model$age_CE)]),
                          "SAM_Index"=c(rev(obs$ANN),model$SAM_Index[match(1956:1000,model$age_CE)]+shift2))

## now we make some rolling means & loess for plotting 
combinedSAM$SAM_Index_3yr <- stats::filter(combinedSAM$SAM_Index, rep(1/3, 3), sides = 2)
combinedSAM$SAM_Index_5yr <- stats::filter(combinedSAM$SAM_Index, rep(1/5, 5), sides = 2)
combinedSAM$SAM_Index_7yr <- stats::filter(combinedSAM$SAM_Index, rep(1/7, 7), sides = 2)
fit <- loess(SAM_Index ~ year, data = combinedSAM, span = 0.2)
combinedSAM$SAM_loess <- predict(fit)

#make a minin dataset for plotting
combinedSAM.s <- combinedSAM[combinedSAM$year>=minval,]
## SAM plot
pdf("figures/climate/clim.SAM.pdf",width=12,height=4)
par(xpd = NA) 
plot(combinedSAM.s$year,combinedSAM.s$SAM_Index,pch=16,col="grey90",cex=1.2,
     xlim=c(minval,2020),ylim=c(-4.5,4.5),bty="n",axes=FALSE,ylab="SAM Index",xlab="")
points(combinedSAM.s$year[combinedSAM.s$year>1956],combinedSAM.s$SAM_Index[combinedSAM.s$year>1956],pch=16,col="grey5",cex=0.9)
points(combinedSAM.s$year[combinedSAM.s$year<1957],combinedSAM.s$SAM_Index[combinedSAM.s$year<1957],pch=16,col="grey50",cex=0.9)
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_3yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.4))
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_5yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.6))
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_7yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.8))
points(combinedSAM.s$year,combinedSAM.s$SAM_loess,type="l",lwd=3,col=add.alpha("darkorchid3",0.8))
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1775,-25),labels = seq(2000,1775,-25) )  
axis(2)
dev.off()

### SAM output for statistics here

climateout$SAM_Index <- combinedSAM$SAM_Index[match(round(climateout$CalYear),combinedSAM$year)]
climateout$SAM_Index_3yr <- combinedSAM$SAM_Index_3yr[match(round(climateout$CalYear),combinedSAM$year)]
climateout$SAM_Index_5yr <- combinedSAM$SAM_Index_5yr[match(round(climateout$CalYear),combinedSAM$year)]
climateout$SAM_Index_7yr <- combinedSAM$SAM_Index_7yr[match(round(climateout$CalYear),combinedSAM$year)]


####====2.0 ENSO data ====####

nino34 <- read.csv("rawdata/climate/nino34.long.anom.csv",header=TRUE, stringsAsFactors = FALSE, check.names = FALSE,comment.char = "#",na.strings ="-9999")
#trim any white space
nino34$Date <- trimws(as.character(nino34$Date))
nino34$Date2 <- as.Date(trimws(as.character(nino34$Date)))
is_slash <- grepl("/", nino34$Date, fixed = TRUE)
is_dash  <- grepl("-", nino34$Date, fixed = TRUE)

## "1870-01-01"
nino34$Date2[is_dash]  <- as.Date(nino34$Date[is_dash], format = "%Y-%m-%d")
## "01/12/1900" = DD/MM/YYYY in your dump (day is always 01)
nino34$Date2[is_slash] <- as.Date(nino34$Date[is_slash], format = "%d/%m/%Y")


nino34$year <- as.integer(format(nino34$Date2, "%Y"))
nino34$month <- as.integer(format(nino34$Date2, "%m"))
## create ENSO year (Jan belongs to current year; Nov–Dec belong to next year)
enso_year <- nino34$year
enso_year[nino34$month %in% c(11, 12)] <- enso_year[nino34$month %in% c(11, 12)] + 1

## keep only NDJ months
ndj <- nino34[nino34$month %in% c(11, 12, 1), ]

## attach ENSO year
ndj$enso_year <- enso_year[nino34$month %in% c(11, 12, 1)]

## compute NDJ mean Niño3.4
nino34_NDJ <- aggregate(value ~ enso_year, data = ndj, FUN = mean, na.rm = TRUE)

## rename cleanly
names(nino34_NDJ) <- c("year", "Nino34_NDJ")

## now we make some rolling means & loess for plotting 
nino34_NDJ$value_3yr <- stats::filter(nino34_NDJ$Nino34_NDJ, rep(1/3, 3), sides = 2)
nino34_NDJ$value_5yr <- stats::filter(nino34_NDJ$Nino34_NDJ, rep(1/5, 5), sides = 2)
nino34_NDJ$value_7yr <- stats::filter(nino34_NDJ$Nino34_NDJ, rep(1/7, 7), sides = 2)
fit <- loess(Nino34_NDJ ~ year, data = nino34_NDJ, span = 0.2)
nino34_NDJ$ENSO_loess <- predict(fit)

### now lets plot 
pdf("figures/climate/clim.ENSO.pdf",width=12,height=4)
par(xpd = NA) 
plot(nino34_NDJ$year,nino34_NDJ$Nino34_NDJ,pch=16,col="grey90",cex=1.2,
     xlim=c(minval,2020),ylim=c(-3,3),bty="n",axes=FALSE,ylab="ENSO NDJ",xlab="")
points(nino34_NDJ$year,nino34_NDJ$Nino34_NDJ,pch=16,col="grey5",cex=0.9)
points(nino34_NDJ$year, nino34_NDJ$value_3yr, type="l", lwd=2, col=add.alpha("firebrick",0.4))
points(nino34_NDJ$year, nino34_NDJ$value_5yr, type="l", lwd=2, col=add.alpha("firebrick",0.6))
points(nino34_NDJ$year, nino34_NDJ$value_7yr, type="l", lwd=2, col=add.alpha("firebrick",0.8))
points(nino34_NDJ$year,nino34_NDJ$ENSO_loess,type="l",lwd=3,col=add.alpha("darkorchid3",0.8))
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1775,-25),labels = seq(2000,1775,-25) )  
axis(4)
dev.off()

climateout$ENSO <- nino34_NDJ$Nino34_NDJ[match(round(climateout$CalYear),nino34_NDJ$year)]
climateout$ENSO_3yr <- nino34_NDJ$value_3yr[match(round(climateout$CalYear),nino34_NDJ$year)]
climateout$ENSO_5yr <- nino34_NDJ$value_5yr[match(round(climateout$CalYear),nino34_NDJ$year)]
climateout$ENSO_7yr <- nino34_NDJ$value_7yr[match(round(climateout$CalYear),nino34_NDJ$year)]


### SAM and ENSO 

pdf("figures/climate/clim.SAM&ENSO.pdf",width=12,height=5)
par(xpd = NA,mfrow=c(2,1),mar=c(1, 4, 0.5, 3),xpd = NA) 
plot(combinedSAM.s$year,combinedSAM.s$SAM_Index,pch=16,col="grey90",cex=1.2,
     xlim=c(1770,2026),ylim=c(-4.6,5.14),bty="n",axes=FALSE,ylab="",xlab="", xaxs = "i", yaxs = "i")
points(combinedSAM.s$year[combinedSAM.s$year>1956],combinedSAM.s$SAM_Index[combinedSAM.s$year>1956],pch=16,col="grey5",cex=0.9)
points(combinedSAM.s$year[combinedSAM.s$year<1957],combinedSAM.s$SAM_Index[combinedSAM.s$year<1957],pch=16,col="grey50",cex=0.9)
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_3yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.4))
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_5yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.6))
points(combinedSAM.s$year,combinedSAM.s$SAM_Index_7yr,type="l",lwd=2,col=add.alpha("dodgerblue",0.8))
points(combinedSAM.s$year,combinedSAM.s$SAM_loess,type="l",lwd=3,col=add.alpha("darkorchid3",0.8))
#axis(1,lwd=1,lwd.ticks = 1,at = seq(2000,1775,-25),labels = seq(2000,1775,-25) )  
axis(2)


plot(nino34_NDJ$year,nino34_NDJ$Nino34_NDJ,pch=16,col="grey90",cex=1.2,
     xlim=c(1770,2026),ylim=c(-3,3),bty="n",axes=FALSE,ylab="",xlab="")
points(nino34_NDJ$year,nino34_NDJ$Nino34_NDJ,pch=16,col="grey5",cex=0.9)
points(nino34_NDJ$year, nino34_NDJ$value_3yr, type="l", lwd=2, col=add.alpha("firebrick",0.4))
points(nino34_NDJ$year, nino34_NDJ$value_5yr, type="l", lwd=2, col=add.alpha("firebrick",0.6))
points(nino34_NDJ$year, nino34_NDJ$value_7yr, type="l", lwd=2, col=add.alpha("firebrick",0.8))
points(nino34_NDJ$year,nino34_NDJ$ENSO_loess,type="l",lwd=3,col=add.alpha("darkorchid3",0.8))
#axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1775,-25),labels = seq(2000,1775,-25) )  
axis(4)
dev.off()




####====3.0 SIC data ====####
##### Here we are extracting Sea Ice Coverage (SIC) from https://figshare.com/collections/_/5709767 paper - https://www.nature.com/articles/s41558-021-01254-9

# --- Seasonal SIC at sites using geodesic-nearest valid grid cell (per season) ---

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

pdf("mapping/SIC.sites_and_nearest_by_season.pdf", width = 8, height = 12)

basemap(limits = lims, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  
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


dev.off()


### Now get rocking with the sea ice concentration (SIC) 

stopifnot(all(c("Site","Season","Year","SIC") %in% names(sic_sites)))

d <- subset(
  sic_sites,
  is.finite(Year) & is.finite(SIC) & !is.na(Site) & !is.na(Season)
)

d$Site   <- factor(d$Site)
d$Season <- factor(d$Season, levels = c("DJF","MAM","JJA","SON"))

season_cols <- c(DJF="#fee8c9", MAM="#fdbb84", SON="#74a9cf", JJA="#045a8d")
pt_alpha <- 0.50

sites    <- levels(d$Site)
xlim_all <- c(1900, 2020)
ylim_all <- range(d$SIC, na.rm = TRUE)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

pdf("figures/climate/PerSiteSeasons.pdf", width = 3.5, height = 5)
par(mfrow = c(5, 1), mar = c(2, 2, 0.2, 1), xaxs = "i", yaxs = "i", bty = "n")

for (s in sites) {
  ds <- d[d$Site == s, ]
  
  plot(NA, xlim = xlim_all, ylim = ylim_all, xlab = "", ylab = "", main = "", las = 1)
  
  for (sea in levels(d$Season)) {
    dss <- ds[ds$Season == sea, ]
    if (!nrow(dss)) next
    
    
    points(1903,80,pch=16,col=coredata$col[match(s,coredata$Site)],cex=2)
    
    points(
      dss$Year, dss$SIC,
      pch = 16, cex = 0.7,
      col = adjustcolor(season_cols[sea], alpha.f = pt_alpha)
    )
    
    ok <- is.finite(dss$Year) & is.finite(dss$SIC)
    if (sum(ok) >= 10) {
      lo <- loess(
        SIC ~ Year, data = dss[ok, ],
        span = 0.35, degree = 2,
        control = loess.control(surface = "direct")
      )
      xx <- seq(min(dss$Year[ok]), max(dss$Year[ok]), length.out = 200)
      yy <- predict(lo, newdata = data.frame(Year = xx))
      lines(xx, yy, col = season_cols[sea], lwd = 2)
    }
  }
  
  box(col = "grey40")
}

dev.off()

#### What about the average across seasons
out <- aggregate(SIC ~ Site + Year, data = d, FUN = mean, na.rm = TRUE)


pdf("figures/climate/PerSiteSIC.pdf", width = 12, height = 4)
par(mfrow = c(1, 1), mar = c(4, 4, 1, 3), xaxs = "i", yaxs = "i", bty = "n",xpd = NA)

plot(out$Year,out$SIC,col=add.alpha(coredata$col[match(out$Site,coredata$Site)],0.3),pch=16,xlim=c(1770,2026),ylab="",xlab="",axes=FALSE)
#axis(1,lwd=1,lwd.ticks = 1,at = seq(2025,1775,-25),labels = seq(2025,1775,-25) )  
axis(4)



for (site in coredata$Site){
  loopdat <- out[out$Site==site,]
  looploess <- loess(SIC~Year,data=loopdat,
                     span = 0.35, degree = 2,
                     control = loess.control(surface = "direct"))
  xx <- seq(1900, 2020,0.1)
  yy <- predict(looploess,newdata = data.frame(Year = xx))  
  
  lines(xx, yy, col = coredata$col[coredata$Site==site], lwd = 3)
  
}

dev.off()



### now lets pull the data out for analysis 

climateout$SIC_DJF <- NA
climateout$SIC_MAM <- NA
climateout$SIC_JJA <- NA
climateout$SIC_SON <- NA

for (line in 1:length(climateout$Site)){
  climateout$SIC_DJF[line] <- sic_sites$SIC[match(paste0(toupper(climateout$Site[line]),round(climateout$CalYear[line]),"DJF"),paste0(sic_sites$Site,sic_sites$Year,sic_sites$Season))]
  climateout$SIC_MAM[line] <- sic_sites$SIC[match(paste0(toupper(climateout$Site[line]),round(climateout$CalYear[line]),"MAM"),paste0(sic_sites$Site,sic_sites$Year,sic_sites$Season))]
  climateout$SIC_JJA[line] <- sic_sites$SIC[match(paste0(toupper(climateout$Site[line]),round(climateout$CalYear[line]),"JJA"),paste0(sic_sites$Site,sic_sites$Year,sic_sites$Season))]
  climateout$SIC_SON[line] <- sic_sites$SIC[match(paste0(toupper(climateout$Site[line]),round(climateout$CalYear[line]),"SON"),paste0(sic_sites$Site,sic_sites$Year,sic_sites$Season))]
}
climateout$SIC_mean <- rowMeans(
  climateout[, c("SIC_DJF","SIC_MAM","SIC_JJA","SIC_SON")],
  na.rm = FALSE
)

plot(climateout[,12:16])


####====4.0 PCA of climate variables ====####
  
## First the 3 year

vars <- c(
  "SAM_Index_3yr",
  "ENSO_3yr",
  "SIC_mean"
)
pca_dat <- climateout[, vars]
pca_dat <- pca_dat[complete.cases(pca_dat), ]

pca <- prcomp(pca_dat, center = TRUE, scale. = TRUE)



summary(pca)
pca$rotation
scores <- pca$x
climate_pca <- cbind(
  climateout[complete.cases(climateout[, vars]), ],
  scores
)
#### loadings plot
pdf("figures/climate/PCA_loadings.pdf",height = 4.5,width = 4)
load <- pca$rotation[, 1:2, drop = FALSE]
plot(load[,1], load[,2], xlim=c(-1,1), ylim=c(-1,1),
     xlab="Loading on PC1", ylab="Loading on PC2", type="n")
abline(h=0, v=0, lty=3)
arrows(0, 0, load[,1], load[,2], length=0.08)
text(load[,1], load[,2], labels=c("SAM","ENSO","Sea Ice"), pos=3, cex=0.8)
dev.off()

#var explained
ve <- (pca$sdev^2) / sum(pca$sdev^2)
pdf("figures/climate/PCA_var.pdf",height = 4,width = 2.8)
barplot(ve, names.arg=paste0("PC", seq_along(ve)),
        ylab="Proportion variance explained")
dev.off()

# plot over time 
year <- climate_pca$CalYear
plot(year, climate_pca$PC2, xlab="Year", ylab="Climate PC1", pch=16)
m1 <- lm(PC2 ~ splines::ns(CalYear, df=4), data=climate_pca)
yy <- seq(min(year, na.rm=TRUE), max(year, na.rm=TRUE), length.out=300)
pred <- predict(m1, newdata=data.frame(CalYear=yy))
lines(yy, pred, lwd=2)


pdf("figures/climate/PCA_PC1.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 1), mar = c(5, 5, 0.6, 0.6), xaxs = "i", yaxs = "i", bty = "n")
plot(climate_pca$CalYear,climate_pca$PC1,col=add.alpha(coredata$col[match(toupper(climate_pca$Site),coredata$Site)],0.8),pch=16,ylab="PC1",xlab="Year",xlim=c(1900,2026),ylim=c(-2.8,2.8))
dev.off()

pdf("figures/climate/PCA_PC2.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 1), mar = c(5, 5, 0.6, 0.6), xaxs = "i", yaxs = "i", bty = "n")
plot(climate_pca$CalYear,climate_pca$PC2,col=add.alpha(coredata$col[match(toupper(climate_pca$Site),coredata$Site)],0.8),pch=16,ylab="PC2",xlab="Year",xlim=c(1900,2026),ylim=c(-3.8,3.8))
dev.off()

write.csv(climateout,"cleandata/climate.data.all.csv")
write.csv(climate_pca,"cleandata/climate.data.PCA.csv")



########### BASEMENT CODE


#### What about the average across seasons
out <- aggregate(SIC ~ Site + Year, data = d, FUN = mean, na.rm = TRUE)


pdf("figures/climate/PerSiteSIC.pdf", width = 5, height = 4)
par(mfrow = c(1, 1), mar = c(5, 5, 0.2, 0.2), xaxs = "i", yaxs = "i", bty = "n")

plot(out$Year,out$SIC,col=add.alpha(coredata$col[match(out$Site,coredata$Site)],0.3),pch=16,xlim=c(1900,2026),ylab="Sea Ice Coverage (%)",xlab="Year")

for (site in coredata$Site){
  loopdat <- out[out$Site==site,]
  looploess <- loess(SIC~Year,data=loopdat,
                     span = 0.35, degree = 2,
                     control = loess.control(surface = "direct"))
  xx <- seq(1900, 2020,0.1)
  yy <- predict(looploess,newdata = data.frame(Year = xx))  
  
  lines(xx, yy, col = coredata$col[coredata$Site==site], lwd = 3)
  
}

dev.off()


### OLD ENSO

obs <- read.csv("rawdata/climate/nina34.anom.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,comment.char = "#",na.strings ="-9999")
model <- read.table("rawdata/climate/ENSO_700yr.txt",header = 1)

## if obs$date is like "1948-01-01"
obs$date <- as.Date(obs$date)

## if it is like "01/01/1948"
## obs$date <- as.Date(obs$date, format = "%d/%m/%Y")

obs$year  <- as.integer(format(obs$date, "%Y"))
obs$month <- as.integer(format(obs$date, "%m"))

## seasonal year: Nov/Dec counted into next year
obs$season_year <- ifelse(obs$month %in% c(11,12), obs$year + 1L, obs$year)

## NDJ = Nov, Dec, Jan (Li says NDJ)
ndj <- obs[obs$month %in% c(11,12,1), c("season_year","month","obs")]
ndj <- ndj[is.finite(ndj$obs), ]

## require all 3 months present (strict); otherwise NA
nino_ndj <- aggregate(obs ~ season_year, data = ndj, FUN = function(v) {
  if (length(v) == 3) mean(v) else NA_real_
})
names(nino_ndj) <- c("year","Nino34_NDJ")
nino_ndj <- nino_ndj[is.finite(nino_ndj$Nino34_NDJ), ]

## overlap with reconstruction
yrs <- intersect(model$age_AD, nino_ndj$year)

r <- model$sst.anom[match(yrs, model$age_AD)]
o <- nino_ndj$Nino34_NDJ[match(yrs, nino_ndj$year)]

ok <- complete.cases(r, o)
r <- r[ok]; o <- o[ok]; yrs <- yrs[ok]

## diagnostics
cor(r, o)
coef(lm(o ~ r))

fit <- lm(o ~ r)

a <- coef(fit)[1]   # intercept
b <- coef(fit)[2]   # slope

model$sst_anom_adj <- a + b * model$sst.anom



plot(o~r)



## SAM plot
pdf("figures/climate/clim.SAM.pdf",width=12,height=4)
plot(combinedSAM$year,combinedSAM$SAM_Index,pch=16,col="grey90",cex=1.2,
     xlim=c(minval,2020),ylim=c(-4,4),bty="n",axes=FALSE,ylab="SAM Index",xlab="")
points(combinedSAM$year,combinedSAM$SAM_Index,pch=16,col="grey10",cex=1.2)
points(combinedSAM$year,combinedSAM$SAM_Index,pch=16,col="grey40",cex=0.9)
points(SAM.long$age_CE,SAM.long$SAM_Index,xlim=c(2016,minval),pch=16,col=add.alpha("grey70",0.6))
points(SAM.long$age_CE,SAM.long$SAM_Index_5yr,type="l",lwd=2,col=add.alpha("black",0.3))
points(SAM.long$age_CE,SAM.long$SAM_loess,type="l",col="darkorchid3",lwd=2)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
dev.off()

#### Now match up the data 

climateout$CalYear


SAM.long$age_CE

match(climateout$CalYear,SAM.long$age_CE)



####====0.1 Functions====####
# Calculate the PZIPSO25 index from IPSO25 and HBI-III concentrations.
# ip25, hbi3: numeric vectors (ng/g dry sediment).
# c: optional balance factor; if NULL it is calculated as mean(ip25)/mean(hbi3).
# Returns a numeric vector of index values (0–1) with "c" stored as an attribute.
pzi_ipso25 <- function(ip25, hbi3, c = NULL, na.rm = TRUE) {
  if (is.null(c)) {
    c <- mean(ip25, na.rm = na.rm) / mean(hbi3, na.rm = na.rm)
  }
  out <- ip25 / (ip25 +(c * hbi3))
  attr(out, "c") <- c
  out
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


####====0.2 Data in====####
#read in metadata 
metadat <- read.csv("metadata/RepMetaData.csv")
dates <- read.csv("metadata/dates.csv")
coredata <- read.csv("metadata/sitemetadata.csv")
coredata$col <- c("#D55E00","#cc79a7","#0072B2","#f0e442","#009e73")


# data in
clim.p56 <- read.table("rawdata/climate/PS97_056-1_geochem_biomarker.txt",skip = 44,header=TRUE, sep="\t")
clim.p68 <- read.table("rawdata/climate/PS97_068-2_geochem_biomarker.txt",skip = 45,header=TRUE, sep="\t")
clim.SC <- read.csv("rawdata/climate/Pienkowski et al_data.csv")
clim.SC$PZIPSO25index <-pzi_ipso25(clim.SC$IPSO25..ng.g.,clim.SC$HBI_III..ng.g.)

minval <- 1770


### Southern Annular Mode data 
### now lets add in the climate data 
  
x <- readLines("rawdata/climate/BES_SAM.txt", warn = FALSE)   # read file lines
x <- trimws(x)                                               # trim whitespace
x <- x[nzchar(x) & !grepl("^#", x)]                          # drop blank and comment lines

hdr <- x[1]                                                  # header with month names
months <- scan(text = hdr, what = character(), quiet = TRUE) # split header to tokens
stopifnot(length(months) == 12)                              # require 12 months

dat <- gsub("\\s+", " ", trimws(x[-1]))                      # normalize spacing on data rows

sam <- read.table(text = dat, header = FALSE, fill = TRUE, stringsAsFactors = FALSE,
                  col.names = c("Year", months), check.names = FALSE)  # parse to table
yearlymean <- data.frame("Year"=sam[,1],
                         "SAMmean"= rowMeans(sam[,-1], na.rm = TRUE))
months <- colnames(sam)[-1] 
y <- c(t(as.matrix(sam[,-1]))) 
y <- y[!is.na(y)]  

sam.data <- data.frame("SAM"=c(t(as.matrix(sam[,-1]))),
                       "year"=rep(sam$Year, each = 12),
                       "month"=rep(months,length(sam[,1])))
sam.data <- sam.data[!is.na(sam.data$SAM),]

##make a continous month

# month name/abbr/number -> integer 1..12
month_num <- function(x){
  if(is.numeric(x)) return(pmin(pmax(as.integer(x),1),12))
  x  <- as.character(x)
  nx <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(nx), nx,
                match(toupper(substr(x,1,3)), toupper(month.abb)))
  as.integer(out)
}

# continuous year, centered on a month (default July=7)
make_continuous_year <- function(year, month, centre=7){
  m  <- month_num(month)
  c0 <- month_num(centre)[1]
  year + (m - c0)/12
}

sam.data$SAM.cont <- make_continuous_year(sam.data$year, sam.data$month, centre="JUL")

## now SAM modelled long term data 

SAM.long <- read.table("rawdata/climate/SAM_1000yr.txt",header = 1)
SAM.long$SAM_Index_5yr <- filter(SAM.long$SAM_Index, rep(1/5, 5), sides = 2)
fit <- loess(SAM_Index ~ age_CE, data = SAM.long, span = 0.2)
SAM.long$SAM_loess <- predict(fit)


## SAM plot
pdf("figures/climate/clim.SAM.pdf",width=12,height=4)
plot(sam.data$SAM.cont,sam.data$SAM,type="l",col="grey90",
     xlim=c(minval,2020),ylim=c(-4,4),bty="n",axes=FALSE,ylab="SAM Index",xlab="")
points(yearlymean$Year,yearlymean$SAMmean,pch=16,col="grey10",cex=1.2)
points(yearlymean$Year,yearlymean$SAMmean,pch=16,col="grey40",cex=0.9)
points(SAM.long$age_CE,SAM.long$SAM_Index,xlim=c(2016,minval),pch=16,col=add.alpha("grey70",0.6))
points(SAM.long$age_CE,SAM.long$SAM_Index_5yr,type="l",lwd=2,col=add.alpha("black",0.3))
points(SAM.long$age_CE,SAM.long$SAM_loess,type="l",col="darkorchid3",lwd=2)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
dev.off()

### now MEI - the ENSO index

# file with the pasted block, e.g. "rawdata/climate/MEI.txt"
months <- toupper(month.abb)

enso <- read.table("rawdata/climate/MEI.txt",
                   header = FALSE, fill = TRUE, comment.char = "#",
                   na.strings = "-999",
                   col.names = c("Year", months), check.names = FALSE)

# yearly mean
enso.yearly <- data.frame(Year = enso$Year,
                          ENSOmean = rowMeans(enso[,-1], na.rm = TRUE))

# long format
vals <- c(t(as.matrix(enso[,-1])))
enso.data <- data.frame(
  ENSO  = vals,
  year  = rep(enso$Year, each = 12),
  month = rep(months, times = nrow(enso))
)
enso.data <- enso.data[!is.na(enso.data$ENSO), ]

enso.data$ENSO.cont <- make_continuous_year(enso.data$year, enso.data$month)


## now the long data 

ENSO.long <- read.table("rawdata/climate/ENSO_700yr.txt",header = 1)
ENSO.long$ENSO_Index_5yr <- filter(ENSO.long$sst.anom, rep(1/5, 5), sides = 2)
fit <- loess(sst.anom ~ age_AD, data = ENSO.long, span = 0.2)
ENSO.long$ENSO_loess <- predict(fit)




## ENSO plot
pdf("figures/climate/clim.ENSO.pdf",width=12,height=4)
plot(enso.data$ENSO.cont,enso.data$ENSO,pch=16,
     xlim=c(minval,2020),bty="n",axes=FALSE,ylab="ENSO Index",xlab="",col=add.alpha(ifelse(enso.data$ENSO>0,"darkred","darkblue"),0.1))
points(enso.yearly$Year,enso.yearly$ENSOmean,pch=16,col="grey20",cex=1.2)
points(enso.yearly$Year,enso.yearly$ENSOmean,pch=16,col="grey60",cex=0.9)

points(ENSO.long$age_AD,ENSO.long$sst.anom,xlim=c(2016,minval),pch=16,col=add.alpha("grey70",0.6))
points(ENSO.long$age_AD,ENSO.long$ENSO_Index_5yr,type="l",lwd=2,col=add.alpha("black",0.3))
points(ENSO.long$age_AD,ENSO.long$ENSO_loess,type="l",col="dodgerblue",lwd=2)

axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
dev.off()



  
## biomarker plot


minval <- min(c(clim.p68$Age..a.AD.CE.,clim.p56$Age..a.AD.CE.))

pdf("figures/climate/clim.biomarkers.pdf",width=12,height=4)
plot(clim.p68$Age..a.AD.CE.,clim.p68$PZIPSO25..PIPSO25.calculated.after.Vorr....,
     pch=16,
     col=coredata$col[coredata$Site=="P68"],
       xlim=c(minval,2020),ylim=c(0,1),
     xlab="year",ylab="BiomarkerIceIndex",bty="n",axes=FALSE)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
points(clim.p56$Age..a.AD.CE.,clim.p56$PZIPSO25..PIPSO25.calculated.after.Vorr....,
       pch=16,
       col=coredata$col[coredata$Site=="P56"])
points(clim.SC$year,clim.SC$PZIPSO25index,
       pch=16,
       col=coredata$col[coredata$Site=="SC"])
dev.off()

pdf("figures/climate/clim.ip25.pdf",width=12,height=4)
plot(clim.p68$Age..a.AD.CE.,clim.p68$C25.HBI.diene.TOC..µg.g...IPSO25.TOC..Ice.Proxy.for.the....,
     pch=16,
     col=coredata$col[coredata$Site=="P68"],
     xlim=c(minval,2020),ylim=c(0,7),
     xlab="year",ylab="IP25 ug/g",bty="n",axes=FALSE)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
points(clim.p56$Age..a.AD.CE.,clim.p56$C25.HBI.diene.TOC..µg.g...IPSO25.TOC..Ice.Proxy.for.the....,
       pch=16,
       col=coredata$col[coredata$Site=="P56"])
points(clim.SC$year,(clim.SC$IPSO25..ng.g./1000),
       pch=16,
       col=coredata$col[coredata$Site=="SC"])
dev.off()

pdf("figures/climate/clim.HBI-Z.pdf",width=12,height=4)
plot(clim.p68$Age..a.AD.CE.,clim.p68$C25.HBI..Z..triene.TOC..µg.g...Gas.chromatography...Mass.spe....,
     pch=16,
     col=coredata$col[coredata$Site=="P68"],
     xlim=c(minval,2020),ylim=c(0,1.7),
     xlab="year",ylab="HBI-Z ug/g",bty="n",axes=FALSE)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
points(clim.p56$Age..a.AD.CE.,clim.p56$C25.HBI..Z..triene.TOC..µg.g...Gas.chromatography...Mass.spe....,
       pch=16,
       col=coredata$col[coredata$Site=="P56"])
points(clim.SC$year,(clim.SC$IPSO25..ng.g./1000),
       pch=16,
       col=coredata$col[coredata$Site=="SC"])

dev.off()



##### Here we are extracting Sea Ice Coverage (SIC) from https://figshare.com/collections/_/5709767 paper - https://www.nature.com/articles/s41558-021-01254-9

# --- Seasonal SIC at sites using geodesic-nearest valid grid cell (per season) ---

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

pdf("mapping/SIC.sites_and_nearest_by_season.pdf", width = 8, height = 12)

basemap(limits = lims, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  
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


dev.off()


### Now get rocking with the sea ice concentration (SIC) 

stopifnot(all(c("Site","Season","Year","SIC") %in% names(sic_sites)))

d <- subset(
  sic_sites,
  is.finite(Year) & is.finite(SIC) & !is.na(Site) & !is.na(Season)
)

d$Site   <- factor(d$Site)
d$Season <- factor(d$Season, levels = c("DJF","MAM","JJA","SON"))

season_cols <- c(DJF="#fee8c9", MAM="#fdbb84", SON="#74a9cf", JJA="#045a8d")
pt_alpha <- 0.50

sites    <- levels(d$Site)
xlim_all <- c(1900, 2020)
ylim_all <- range(d$SIC, na.rm = TRUE)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

pdf("figures/climate/PerSiteSeasons.pdf", width = 3.5, height = 5)
par(mfrow = c(5, 1), mar = c(2, 2, 0.2, 1), xaxs = "i", yaxs = "i", bty = "n")

for (s in sites) {
  ds <- d[d$Site == s, ]
  
  plot(NA, xlim = xlim_all, ylim = ylim_all, xlab = "", ylab = "", main = "", las = 1)
  
  for (sea in levels(d$Season)) {
    dss <- ds[ds$Season == sea, ]
    if (!nrow(dss)) next
    
    
    points(1903,80,pch=16,col=coredata$col[match(s,coredata$Site)],cex=2)
    
    points(
      dss$Year, dss$SIC,
      pch = 16, cex = 0.7,
      col = adjustcolor(season_cols[sea], alpha.f = pt_alpha)
    )
    
    ok <- is.finite(dss$Year) & is.finite(dss$SIC)
    if (sum(ok) >= 10) {
      lo <- loess(
        SIC ~ Year, data = dss[ok, ],
        span = 0.35, degree = 2,
        control = loess.control(surface = "direct")
      )
      xx <- seq(min(dss$Year[ok]), max(dss$Year[ok]), length.out = 200)
      yy <- predict(lo, newdata = data.frame(Year = xx))
      lines(xx, yy, col = season_cols[sea], lwd = 2)
    }
  }
  
  box(col = "grey40")
}

dev.off()


#### What about the average across seasons
out <- aggregate(SIC ~ Site + Year, data = d, FUN = mean, na.rm = TRUE)


pdf("figures/climate/PerSiteSIC.pdf", width = 5, height = 4)
par(mfrow = c(1, 1), mar = c(5, 5, 0.2, 0.2), xaxs = "i", yaxs = "i", bty = "n")

plot(out$Year,out$SIC,col=add.alpha(coredata$col[match(out$Site,coredata$Site)],0.3),pch=16,xlim=c(1900,2026),ylab="Sea Ice Coverage (%)",xlab="Year")

for (site in coredata$Site){
loopdat <- out[out$Site==site,]
looploess <- loess(SIC~Year,data=loopdat,
                   span = 0.35, degree = 2,
                   control = loess.control(surface = "direct"))
xx <- seq(1900, 2020,0.1)
yy <- predict(looploess,newdata = data.frame(Year = xx))  

lines(xx, yy, col = coredata$col[coredata$Site==site], lwd = 3)

}

dev.off()






### Now we take a look at SST netCDF

sites <- read.csv("metadata/sitemetadata.csv")

sst <-nc_open("rawdata/climate/HadISST_sst.nc")

lat <- ncvar_get(sst,"latitude")
lon <- ncvar_get(sst,"longitude")
sst_out <- ncvar_get(sst,"sst")
summary(lat)
summary(lon)
summary(sst_out)
nc_close(sst)



# lat, lon, sst_out, sites already exist

t_index <- 2

# mask of ocean points (non-NA) for this time
sea_mask <- !is.na(sst_out[,, t_index])



idx <- which(sea_mask, arr.ind = TRUE)   # row=lon index, col=lat index

lon_sea <- lon[idx[,"row"]]
lat_sea <- lat[idx[,"col"]]

nearest_sea <- function(phi, lam) {
  d2 <- (lat_sea - phi)^2 + (lon_sea - lam)^2
  k  <- which.min(d2)
  c(newlat = lat_sea[k], newlon = lon_sea[k])
}

res <- t(mapply(nearest_sea, sites$Lat, sites$Lon))

sst_out[sst_out==-1000] <- NA




sites$newlat <- res[,"newlat"]
sites$newlon <- res[,"newlon"]

## ok now lets plot something out 


pdf("figures/climate/clim.sst.pdf",height=9,width = 5)
par(mfrow = c(5,1), bty = "n", mar = c(3,4,1,1))

year <- 1870 + (0:(1871-1))/12  # 1871 monthly points

for (i in 1:5) {
  test <- sst_out[
    match(sites$newlon[i], lon),
    match(sites$newlat[i], lat),
  ]
  
  plot(year, test,
       pch = 16, ylab = "SST", xlab = "year",
       col = add.alpha(coredata$col[i], 0.2))
  
  ok <- is.finite(test) & is.finite(year)
  
  fit <- loess(test[ok] ~ year[ok], span = 0.25, degree = 1,
               family = "symmetric", na.action = na.exclude)
  
  lines(year[ok], predict(fit, year[ok]), lwd = 2, col = coredata$col[i])
}

dev.off()

library(strucchange)

time <- (1:1871)/12 + 1870
mon  <- rep(1:12, length.out=1871)

par(mfrow=c(5,1), bty = "n")

for (i in 1:5) {
  
  sst <- sst_out[
    match(sites$newlon[i], lon),
    match(sites$newlat[i], lat),
  ]
  
  clim <- tapply(sst, mon, mean, na.rm=TRUE)
  anom <- sst - clim[mon]
  
  ok <- is.finite(anom) & is.finite(time)
  
  # fit breakpoints on complete cases only
  bp <- breakpoints(anom[ok] ~ time[ok], h=120)
  
  # plotting
  plot(time, anom, type="l", col=add.alpha("grey", 0.3),
       xlab="Year", ylab="SST anomaly (°C)")
  
  # fitted values only exist for ok points; put them back on full timeline
  fit <- rep(NA_real_, length(anom))
  fit[ok] <- fitted(bp)
  
  lines(time, fit, lwd=3)
}


### Where are the SST sites?
pdf("mapping/SST.sites.pdf",height=7,width =5 )
basemap(limits = lims, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),
                                data = sites, color = "black",
                                cex = 4.2, crs = 4326) +
  
  ggspatial::geom_spatial_point(aes(x = newlon, y = newlat, ),
                                data = sites, cex = 3.2, alpha = 0.9, crs = 4326) +

  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()




#### AI generated basement for wavelet and correlations


## wavelet analysis 

## --- Wavelets for sea-ice proxies vs SAM & ENSO -----------------------------

if (!requireNamespace("WaveletComp", quietly = TRUE)) install.packages("WaveletComp")
library(WaveletComp)

dir.create("figures", showWarnings = FALSE)

z <- function(x) as.numeric(scale(x))  # z-score

# Regularize two series to annual grid over overlap (no extrapolation)
reg_pair <- function(t1, y1, t2, y2){
  ok1 <- is.finite(t1) & is.finite(y1)
  ok2 <- is.finite(t2) & is.finite(y2)
  t1 <- t1[ok1]; y1 <- y1[ok1]
  t2 <- t2[ok2]; y2 <- y2[ok2]
  start <- ceiling(max(min(t1), min(t2)))
  end   <- floor(min(max(t1), max(t2)))
  if (end - start < 10) return(NULL) # too short
  grid  <- seq(start, end, by = 1)
  y1i <- approx(t1, y1, xout = grid)$y
  y2i <- approx(t2, y2, xout = grid)$y
  data.frame(t = grid, a = z(y1i), b = z(y2i))
}

# One-stop wavelet power + coherence and plots
do_wavelets <- function(df, site, driver, upperPeriod = NULL, dj = 1/12, nsim = 500){
  if (is.null(df) || nrow(df) < 32) return(invisible(NULL))
  if (is.null(upperPeriod)) upperPeriod <- floor(nrow(df)/2)
  
  dat <- data.frame(year = df$t, proxy = df$a, driver = df$b)
  
  w.prox <- analyze.wavelet(dat, "proxy",
                            dt = 1, dj = dj,
                            lowerPeriod = 2, upperPeriod = upperPeriod,
                            loess.span = 0, make.pval = TRUE, n.sim = nsim)
  
  pdf(sprintf("figures/wavelet.power_%s.pdf", site), width = 8, height = 5)
  wt.image(w.prox,
           main = sprintf("%s PZIPSO25 – wavelet power", site),
           timelab = "Year",
           legend.params = list(lab = "Power"))
  dev.off()
  
  coh <- analyze.coherency(dat, my.pair = c("proxy","driver"),
                           dt = 1, dj = dj,
                           lowerPeriod = 2, upperPeriod = upperPeriod,
                           loess.span = 0, make.pval = TRUE, n.sim = nsim)
  
  pdf(sprintf("figures/wavelet.coherence_%s_vs_%s.pdf", site, driver), width = 8, height = 5)
  wc.image(coh,
           main = sprintf("%s PZIPSO25 vs %s – coherence", site, driver),
           timelab = "Year",
           legend.params = list(lab = "R²"))
  dev.off()
  invisible(coh)
}
## --- Build series from your objects ----------------------------------------

# Proxies
p68 <- list(t = clim.p68$Age..a.AD.CE.,
            v = clim.p68$PZIPSO25..PIPSO25.calculated.after.Vorr....)
p56 <- list(t = clim.p56$Age..a.AD.CE.,
            v = clim.p56$PZIPSO25..PIPSO25.calculated.after.Vorr....)
sc  <- list(t = clim.SC$year,
            v = clim.SC$PZIPSO25index)

# Drivers (annual)
sam  <- list(t = SAM.long$age_CE,  v = SAM.long$SAM_Index)
enso <- list(t = ENSO.long$age_AD, v = ENSO.long$sst.anom)

sites   <- list(P68 = p68, P56 = p56, SC = sc)
drivers <- list(SAM = sam, ENSO = enso)

## --- Run coherence analyses (writes PDFs) ----------------------------------

for (s in names(sites)){
  for (d in names(drivers)){
    df <- reg_pair(sites[[s]]$t, sites[[s]]$v, drivers[[d]]$t, drivers[[d]]$v)
    do_wavelets(df, site = s, driver = d, nsim = 10000)  # increase nsim for tougher sig
  }
}

## simpler

dir.create("figures", showWarnings = FALSE)

z <- function(x) as.numeric(scale(x))  # standardize

plot_combo <- function(df, site, driver, lag.max = 20){
  if (is.null(df) || nrow(df) < 12) return(invisible(NULL))
  y <- z(df$a); x <- z(df$b); tt <- df$t
  lm1 <- lm(y ~ x); r <- cor(x,y, use="pairwise"); rho <- cor(x,y, method="spearman", use="pairwise")
  r2 <- summary(lm1)$r.squared; b1 <- coef(lm1)[2]
  
  cc <- ccf(y, x, lag.max = lag.max, plot = FALSE, na.action = na.pass)  # y vs x
  k  <- which.max(abs(cc$acf)); best.lag <- cc$lag[k]; best.cc <- cc$acf[k]
  
  pdf(sprintf("figures/driver_links_%s_vs_%s.pdf", site, driver), width=9, height=6)
  op <- par(mfrow=c(2,2), mar=c(3.5,4,2,1), oma=c(0,0,0,0), xaxs="i")
  
  # (1) time series (z-scaled)
  plot(tt, x, type="l", col="grey40", lwd=2, xlab="Year", ylab="z-score",
       main=sprintf("%s vs %s", site, driver))
  lines(tt, y, col="black", lwd=2)
  legend("topleft", bty="n", lty=1, col=c("grey40","black"), lwd=2,
         legend=c(driver,"PZIPSO25"))
  
  # (2) scatter + OLS
  plot(x, y, pch=16, col=rgb(0,0,0,0.5), xlab=paste0(driver," (z)"), ylab="PZIPSO25 (z)")
  abline(lm1, lwd=2)
  legend("topleft", bty="n",
         legend=c(sprintf("Pearson r=%.2f", r),
                  sprintf("Spearman ρ=%.2f", rho),
                  sprintf("β=%.2f,  R²=%.2f", b1, r2)))
  
  # (3) cross-correlation
  plot(cc, main=sprintf("CCF: %s (y) vs %s (x)", "PZIPSO25", driver), ylab="corr")
  abline(v=best.lag, lty=2)
  mtext(sprintf("best lag = %d yr, corr = %.2f", best.lag, best.cc), side=3, line=0.2)
  
  # (4) residuals quick check (optional)
  plot(lm1, which=1)  # residuals vs fitted
  par(op)
  dev.off()
  invisible(list(pearson=r, spearman=rho, beta=b1, r2=r2, best_lag=best.lag, best_cc=best.cc))
}

# run for all site x driver pairs
out <- list()
for (s in names(sites)){
  for (d in names(drivers)){
    df <- reg_pair(sites[[s]]$t, sites[[s]]$v, drivers[[d]]$t, drivers[[d]]$v)
    out[[paste(s,d,sep="_")]] <- plot_combo(df, site=s, driver=d, lag.max=20)
  }
}






