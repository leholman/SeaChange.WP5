## age depth
library(rplum)
library(rbacon)

# data in

datain <- read.csv("rawdata/pb210.csv")
ages <- read.csv("metadata/dates.csv")

## data out

ages.out<-cbind(ages[,1:3],
                "upr95"=rep("",length(ages$Depth)),
                "lwr95"=rep("",length(ages$Depth)),
                "median"=rep("",length(ages$Depth)), 
                "mean"=rep("",length(ages$Depth)))


## 68

data68 <- datain[datain$Core=="p68",]
data68.1 <- data68[!is.na(data68$MidDepthOfDatedCore),]

newdat <- data.frame("lab"=data68.1$ID,
                     "depth.cm."=data68.1$DepthTo, 
                     "density.g.cm.3."=data68.1$Density.g.cm., 
                     "X210Pb.Bq.kg."=data68.1$X210Pb.Bq.kg, 
                     "sd.210Pb."=data68.1$sd.210Pb, 
                     "thickness.cm."=data68.1$thickness,
                     "radium"=data68.1$Radium,
                     "radiumSD"=data68.1$Radium.sd.)
dir.create("Plum_runs/indat.68")
write.csv(newdat,file="Plum_runs/indat.68/indat.68.csv",row.names = FALSE)

Plum(core="indat.68",d.min = 0,ra.case=1,d.max = 48,ssize=20000,date.sample = 2016,n.supp=0)
pdf("figures/age.depth/p68.pdf",height = 8,width = 6.5)
agedepth(BCAD=TRUE,d.max = 46,remove.tail=FALSE)
dev.off()

# harvest ages 
for (row in grep("p68",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}

## 56

data56 <- datain[datain$Core=="p56",]
data56.1 <- data56[!is.na(data56$Density.g.cm.)&!is.na(data56$X210Pb.Bq.kg),]

newdat <- data.frame("lab"=data56.1$ID,
                     "depth.cm."=data56.1$DepthTo, 
                     "density.g.cm.3."=data56.1$Density.g.cm., 
                     "X210Pb.Bq.kg."=data56.1$X210Pb.Bq.kg, 
                     "sd.210Pb."=data56.1$sd.210Pb, 
                     "thickness.cm."=data56.1$thickness,
                     "radium"=data56.1$Radium,
                     "radiumSD"=data56.1$Radium.sd.)
dir.create("Plum_runs/indat.56")
write.csv(newdat,file="Plum_runs/indat.56/indat.56.csv",row.names = FALSE)

Plum(core="indat.56",d.min = 0,ra.case=1,d.max = 36,ssize=20000,date.sample = 2016,n.supp=0)
pdf("figures/age.depth/p56.pdf",height = 8,width = 6.5)
agedepth(BCAD=TRUE,d.min=0,d.max = 38,remove.tail=FALSE)
dev.off()

# harvest ages 
for (row in grep("p56",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}



## BB

dataBB <- datain[datain$Core=="BB",]
dataBB <- dataBB[!is.na(dataBB$X210Pb.Bq.kg),]

newdat <- data.frame("lab"=dataBB$ID,
                     "depth.cm."=dataBB$DepthTo, 
                     "density.g.cm.3."=dataBB$Density.g.cm., 
                     "X210Pb.Bq.kg."=dataBB$X210Pb.Bq.kg, 
                     "sd.210Pb."=dataBB$sd.210Pb, 
                     "thickness.cm."=dataBB$thickness,
                     "radium"=dataBB$Radium,
                     "radiumSD"=dataBB$Radium.sd.)
dir.create("Plum_runs/indat.BB")
write.csv(newdat,file="Plum_runs/indat.BB/indat.BB.csv",row.names = FALSE,)

Plum(core="indat.BB",d.min = 0,ra.case=1,d.max = 28,ssize=20000,date.sample = 2019,n.supp=0)
pdf("figures/age.depth/BB.pdf",height = 8,width = 6.5)
agedepth(BCAD=TRUE,d.max = 30,remove.tail=FALSE)
dev.off()

# harvest ages

for (row in grep("BB",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}


## SC

dataSC <- datain[datain$Core=="SC",]
dataSC <- dataSC[!is.na(dataSC$X210Pb.Bq.kg),]

newdat <- data.frame("lab"=dataSC$ID,
                     "depth.cm."=dataSC$DepthTo, 
                     "density.g.cm.3."=dataSC$Density.g.cm., 
                     "X210Pb.Bq.kg."=dataSC$X210Pb.Bq.kg, 
                     "sd.210Pb."=dataSC$sd.210Pb, 
                     "thickness.cm."=dataSC$thickness,
                     "radium"=dataSC$Radium,
                     "radiumSD"=dataSC$Radium.sd.)
dir.create("Plum_runs/indat.SC",)
write.csv(newdat,file="Plum_runs/indat.SC/indat.SC.csv",row.names = FALSE,)

Plum(core="indat.SC",d.min = 0,ra.case=1,d.max = 29,ssize=20000,date.sample = 2019,n.supp=0)
pdf("figures/age.depth/SC.pdf",height = 8,width = 6.5)
agedepth(BCAD=TRUE,d.max = 30,remove.tail=FALSE)
dev.off()

# harvest ages

for (row in grep("SC",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}


## MC

dataMC <- datain[datain$Core=="MC",]
dataMC <- dataMC[!is.na(dataMC$X210Pb.Bq.kg),]

newdat <- data.frame("lab"=dataMC$ID,
                     "depth.cm."=dataMC$DepthTo, 
                     "density.g.cm.3."=dataMC$Density.g.cm., 
                     "X210Pb.Bq.kg."=dataMC$X210Pb.Bq.kg, 
                     "sd.210Pb."=dataMC$sd.210Pb, 
                     "thickness.cm."=dataMC$thickness,
                     "radium"=dataMC$Radium,
                     "radiumSD"=dataMC$Radium.sd.)
dir.create("Plum_runs/indat.MC",)
write.csv(newdat,file="Plum_runs/indat.MC/indat.MC.csv",row.names = FALSE,)

Plum(core="indat.MC",d.min = 0,ra.case=1,d.max = 28,ssize=20000,date.sample = 2019,n.supp=0)
pdf("figures/age.depth/MC.pdf",height = 8,width = 6.5)
agedepth(BCAD=TRUE,d.max = 30,remove.tail=FALSE,pb.lim = c(0,300))
dev.off()

# harvest ages

for (row in grep("MC",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}

Plum(ask=FALSE, suggest=FALSE, ssize=100, coredir=tempdir())

## Age output

write.csv(ages.out,"metadata/dates.out.rPLUM.csv")



############################################################
## 210Pb age–depth models with serac for 5 cores
## Uses existing indat.*.csv files as input
############################################################

## 0. Install and load serac --------------------------------

if (!requireNamespace("serac", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("rosalieb/serac", build_vignettes = FALSE)
}

library(serac)

## 1. Helper: reformat indat.*.csv to serac input ----------
setwd("Serac_runs/")
# core$list element: name, infile, coring_yr
prep_serac_input <- function(core) {
  dat <- read.csv(core$infile, stringsAsFactors = FALSE)
  
  # Columns are assumed in this order in indat.*.csv:
  # 1: labID
  # 2: depth (cm)
  # 3: density (g/cm^3)
  # 4: total 210Pb (Bq/kg)
  # 5: sd(total 210Pb)
  # 6: thickness (cm)
  # 7: 226Ra (Bq/kg)
  # 8: sd(226Ra)
  # 9: settings       – ignored
  
  depth_cm   <- dat[, 2]
  dens_g_cm3 <- dat[, 3]
  pb_tot     <- dat[, 4]
  pb_tot_sd  <- dat[, 5]
  thick_cm   <- dat[, 6]
  ra_Bq_kg   <- dat[, 7]
  ra_sd      <- dat[, 8]
  
  # Excess Pb and propagated sd (independent errors)
  pb_ex      <- pb_tot - ra_Bq_kg
  pb_ex_sd   <- sqrt(pb_tot_sd^2 + ra_sd^2)
  
  # Convert depth/thickness to mm for serac (input_depth_mm = TRUE)
  serac_dat <- data.frame(
    depth_avg = depth_cm * 10,
    thickness = thick_cm * 10,
    density   = dens_g_cm3,
    Pbex      = pb_ex,
    Pbex_er   = pb_ex_sd,
    Ra226     = ra_Bq_kg,
    Ra226_er  = ra_sd
  )
  
  # Create Cores/<name>/<name>.txt for serac
  cores_root <- file.path(getwd(), "Cores")
  dir.create(cores_root, showWarnings = FALSE)
  
  core_dir <- file.path(cores_root, core$name)
  dir.create(core_dir, showWarnings = FALSE)
  
  out_file <- file.path(core_dir, paste0(core$name, ".txt"))
  write.table(
    serac_dat,
    file      = out_file,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
  
  invisible(out_file)
}

## 2. Define the 5 cores ------------------------------------

cores <- list(
  list(name = "indat.68", infile = "indat.68.csv", coring_yr = 2016),
  list(name = "indat.56", infile = "indat.56.csv", coring_yr = 2016),
  list(name = "indat.BB", infile = "indat.BB.csv", coring_yr = 2019),
  list(name = "indat.SC", infile = "indat.SC.csv", coring_yr = 2019),
  list(name = "indat.MC", infile = "indat.MC.csv", coring_yr = 2019)
)

## 3. Run serac for each core -------------------------------

serac_models <- list()  # analogous to storing Plum runs

for (core in cores) {
  
  ## 3.1 Format and write serac input file
  prep_serac_input(core)
  
  ## 3.2 Run serac (simple CFCS + CRS Pb-only model)
  ##     – hemisphere / Cs / Am left out for now
  model_obj <- serac(
    name               = core$name,
    coring_yr          = core$coring_yr,
    model              ="CRS",
    plotphoto          = FALSE,
    plot_Pb            = TRUE,
    plot_Cs            = FALSE,
    plot_Am            = FALSE,
    plot_Pb_inst_deposit = FALSE,
    input_depth_mm     = TRUE,
    plotpdf            = TRUE,   # serac writes its own pdf in Cores/<name>/
    preview            = TRUE,   # set FALSE if you don’t want on-screen plots
    save_code          = FALSE,
    prop_width_fig     = 2.0,     
    prop_height_fig    = 1 
  )
  
  serac_models[[core$name]] <- model_obj
}

## Compare it all 


coredata <- read.csv("metadata/sitemetadata.csv")
coredata$col <- c("#D55E00","#cc79a7","#0072B2","#f0e442","#009e73")

dates.ex <- read.csv("metadata/dates.csv")
dates.plum <- read.csv("metadata/dates.out.rPLUM.csv")

cores <- c("BB","SC","MC","68","56")
dates.serac <- do.call(rbind, lapply(cores, function(core){
  x <- read.table(sprintf("Serac_runs/Cores/indat.%s/indat.%s_CRS_interpolation.txt", core, core), h=1)
  x$core <- ifelse(core=="68"|core=="56",paste0("p",core),core)
  x
}))



pdf("figures/age.depth/comp.pdf",height=4,width = 10)

par(mfrow=c(1,3))
## Expert vs plum 
plot(dates.plum$mean,dates.ex$CalYear,col=coredata$col[match(toupper(dates.plum$Site),coredata$Site)],pch=16,xlab="rPLUM Estimates",ylab="Expert Estimates")
abline(b=1,a=0,lty=2,col="darkred",cex=2)

## Expert vs serac
dates.serac.trim <- dates.serac[match(paste0(dates.ex$Site,dates.ex$Depth),paste0(dates.serac$core,dates.serac$depth_top_mm/10)),]
plot(dates.serac.trim$BestAD_CRS,dates.ex$CalYear,col=coredata$col[match(toupper(dates.serac.trim$core),coredata$Site)],pch=16,xlab="serac Estimates",ylab="Expert Estimates")
abline(b=1,a=0,lty=2,col="darkred",cex=2)

## plum vs serac
plot(dates.serac.trim$BestAD_CRS,dates.plum$mean,col=coredata$col[match(toupper(dates.serac.trim$core),coredata$Site)],pch=16,xlab="serac Estimates",ylab="rPLUM Estimates")
abline(b=1,a=0,lty=2,col="darkred",cex=2)
dev.off()


##Thorbjørn68 
TB <- read.csv("rawdata/68_TJA.csv")

dat.68 <- dates.plum[dates.plum$Site=="p68",]
dates.serac.trim[dates.serac.trim$core=="p68",]


dates.plum$Depth[dates.plum$Site=="p68"]

paste0(dat.68$Site,dat.68$Depth)

paste0(dates.serac.trim$core,(dates.serac.trim$depth_top_mm/10))

dat.68$TB <- TB$Date[match(paste0(dat.68$Site,dat.68$Depth),paste0("p68",TB$Depth))]
dat.68$TB.err <-TB$error.age[match(paste0(dat.68$Site,dat.68$Depth),paste0("p68",TB$Depth))]
dat.68$serac <- dates.serac.trim$BestAD_CRS[match(paste0(dat.68$Site,dat.68$Depth),paste0(dates.serac.trim$core,(dates.serac.trim$depth_top_mm/10)))]
dat.68$exp <- dates$CalYear[match(paste0(dat.68$Site,dat.68$Depth),paste0(dates$Site,dates$Depth))]


par(mfrow=c(1,3))
plot(dat.68$mean,dat.68$TB,xlab = "rplum",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)
plot(dat.68$serac,dat.68$TB,xlab = "serac",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)
plot(dat.68$exp,dat.68$TB,xlab = "expert",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)


##ThorbjørnSC
TB <- read.csv("rawdata/SC_TJA.csv")

dat.SC <- dates.plum[dates.plum$Site=="SC",]
dates.serac.trim[dates.serac.trim$core=="SC",]

SC.interpolated <- spline(TB$Depth,TB$Date,xout = seq(0,13.5,0.5))


dat.SC$TB <- SC.interpolated$y[match(paste0(dat.SC$Site,dat.SC$Depth),paste0("SC",SC.interpolated$x))]
dat.SC$serac <- dates.serac.trim$BestAD_CRS[match(paste0(dat.SC$Site,dat.SC$Depth),paste0(dates.serac.trim$core,(dates.serac.trim$depth_top_mm/10)))]
dat.SC$exp <- dates$CalYear[match(paste0(dat.SC$Site,dat.SC$Depth),paste0(dates$Site,dates$Depth))]


par(mfrow=c(1,3))
plot(dat.SC$mean,dat.SC$TB,xlab = "rplum",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)
plot(dat.SC$serac,dat.SC$TB,xlab = "serac",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)
plot(dat.SC$exp,dat.SC$TB,xlab = "expert",ylab="Thorbjørn",pch=16)
abline(b=1,a=0,lty=2,col="darkred",cex=2)







