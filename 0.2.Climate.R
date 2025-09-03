########################################################
####=.     Analysis of WP5  Antarctic climate      =####
####==== Luke E. Holman====25.08.2025====################
########################################################

####====0.0 Packages====####

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

minval <- min(c(clim.p68$Age..a.AD.CE.,clim.p56$Age..a.AD.CE.))


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
pdf("figures/clim.SAM.pdf",width=12,height=4)
plot(sam.data$SAM.cont,sam.data$SAM,type="l",col="grey90",
     xlim=c(2016,minval),ylim=c(-4,4),bty="n",axes=FALSE,ylab="SAM Index",xlab="")
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
pdf("figures/clim.ENSO.pdf",width=12,height=4)
plot(enso.data$ENSO.cont,enso.data$ENSO,pch=16,
     xlim=c(2016,minval),bty="n",axes=FALSE,ylab="ENSO Index",xlab="",col=add.alpha(ifelse(enso.data$ENSO>0,"darkred","darkblue"),0.1))
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

pdf("figures/clim.biomarkers.pdf",width=12,height=4)
plot(clim.p68$Age..a.AD.CE.,clim.p68$PZIPSO25..PIPSO25.calculated.after.Vorr....,
     pch=16,
     col=coredata$col[coredata$Site=="P68"],
       xlim=c(2016,minval),ylim=c(0,1),
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






