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

Plum(core="indat.SC",d.min = 0,ra.case=1,d.max = 28,ssize=20000,date.sample = 2019,n.supp=0)
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
agedepth(BCAD=TRUE,d.max = 30,remove.tail=FALSE)
dev.off()

# harvest ages

for (row in grep("MC",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}








## Age output



write.csv(ages.out,"metadata/dates.out.csv")




