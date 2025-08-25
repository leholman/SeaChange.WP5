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

## SAM plot
pdf("figures/clim.SAM.pdf",width=12,height=4)
plot(sam.data$SAM.cont,sam.data$SAM,type="l",col="grey90",
     xlim=c(2016,minval),bty="n",axes=FALSE,ylab="SAM Index",xlab="")
points(yearlymean$Year,yearlymean$SAMmean,pch=16,col="grey20",cex=1.2)
points(yearlymean$Year,yearlymean$SAMmean,pch=16,col="grey60",cex=0.9)
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

## ENSO plot
pdf("figures/clim.ENSO.pdf",width=12,height=4)
plot(enso.data$ENSO.cont,enso.data$ENSO,pch=16,
     xlim=c(2016,minval),bty="n",axes=FALSE,ylab="ENSO Index",xlab="",col=add.alpha(ifelse(enso.data$ENSO>0,"darkred","darkblue"),0.1))
points(enso.yearly$Year,enso.yearly$ENSOmean,pch=16,col="grey20",cex=1.2)
points(enso.yearly$Year,enso.yearly$ENSOmean,pch=16,col="grey60",cex=0.9)
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







