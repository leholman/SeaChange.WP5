## age depth
library(rplum)
library(rbacon)

# data in

datain <- read.csv("rawdata/pb210dat.csv")
ages <- read.csv("metadata/dates.csv")

## data out

ages.out<-cbind(ages[,1:3],
                "upr95"=rep("",length(ages$Depth)),
                "lwr95"=rep("",length(ages$Depth)),
                "median"=rep("",length(ages$Depth)), 
                "mean"=rep("",length(ages$Depth)))




data68 <- datain[datain$Core=="p68",]
data68.1 <- data68[!is.na(data68$MidDepthOfDatedCore),]

newdat <- data.frame("lab"=data68.1$ID,
                     "depth.cm."=data68.1$DepthTo, 
                     "density.g.cm.3."=data68.1$Density..g.cm., 
                     "X210Pb.Bq.kg."=data68.1$X210Pb.Bq.kg, 
                     "sd.210Pb."=data68.1$sd.210Pb, 
                     "thickness.cm."=rep(1,length(data68.1$ID)),
                     "radium"=data68.1$Radium,
                     "radiumSD"=data68.1$Radium.sd.)
write.csv(newdat,file="Plum_runs/indat.68/indat.68.csv",row.names = FALSE)

Plum(core="indat.68",d.min = 0,ra.case=1,d.max = 48,ssize=20000,date.sample = 2016)
agedepth(BCAD=TRUE,d.max = 46,remove.tail=FALSE)

## Age output


####

for (row in grep("p68",ages.out$Site)){
  loopages <- Bacon.Age.d(ages.out$Depth[row])
  ages.out$upr95[row] <- 1950-round(quantile(loopages, 0.95, na.rm=TRUE),digits=0)
  ages.out$lwr95[row] <- 1950-round(quantile(loopages, 0.05, na.rm=TRUE),digits=0)
  ages.out$median[row] <- 1950-round(median(loopages),digits=0)
  ages.out$mean[row] <- 1950-round(mean(loopages),digits=0)
}

write.csv(ages.out,"metadata/dates.out.csv")




