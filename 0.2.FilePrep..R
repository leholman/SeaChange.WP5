
#############################################
####===Data prep before analysis WP5.2===####
####==== Luke E. Holman====01.05.2023====####
#############################################

####====0.0 Packages====####



####====1.0 Make a dataset to play with====####


MAM.data<- read.csv("rawdata/MAM.raw.csv",row.names = 1)

MAM.tax <- read.csv("taxonomy/MAM.parsed.csv",row.names = 1)

colnames(MAM.data)

MAM.all <- cbind(MAM.data,MAM.tax[match(rownames(MAM.data),MAM.tax$OTU),])

sampleIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(MAM.data))

unique(sampleIndex)

MAM.nReps <- data.frame(matrix(0,nrow = length(MAM.data[,1]),ncol=length(unique(sampleIndex))))

colnames(MAM.nReps) <- unique(sampleIndex)
rownames(MAM.nReps) <- rownames(MAM.data)

MAM.binary <-MAM.data
MAM.binary[MAM.binary<11] <- 0
MAM.binary[MAM.binary>10] <- 1
binaryIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(MAM.binary))

for (column in 1:length(MAM.nReps[1,])){
  MAM.nReps[,column] <- rowSums(MAM.binary[,binaryIndex %in% colnames(MAM.nReps)[column]])
  }

MAM.all.nReps <- cbind(MAM.nReps,MAM.tax[match(rownames(MAM.binary),MAM.tax$OTU),])

write.csv(MAM.all,"rawdata/MAM.raw.wTAX.csv")
write.csv(MAM.all.nReps,"rawdata/MAM.raw.nReps.wTAX.csv")



RIZ.data<- read.csv("rawdata/RIZ.raw.csv",row.names = 1)

RIZ.tax <- read.csv("taxonomy/RIZ.parsed.csv",row.names = 1)

colnames(RIZ.data)

RIZ.all <- cbind(RIZ.data,RIZ.tax[match(rownames(RIZ.data),RIZ.tax$OTU),])

sampleIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(RIZ.data))

unique(sampleIndex)

RIZ.nReps <- data.frame(matrix(0,nrow = length(RIZ.data[,1]),ncol=length(unique(sampleIndex))))

colnames(RIZ.nReps) <- unique(sampleIndex)
rownames(RIZ.nReps) <- rownames(RIZ.data)

RIZ.binary <-RIZ.data
RIZ.binary[RIZ.binary<11] <- 0
RIZ.binary[RIZ.binary>10] <- 1
binaryIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(RIZ.binary))

for (column in 1:length(RIZ.nReps[1,])){
  RIZ.nReps[,column] <- rowSums(RIZ.binary[,binaryIndex %in% colnames(RIZ.nReps)[column]])
}

RIZ.all.nReps <- cbind(RIZ.nReps,RIZ.tax[match(rownames(RIZ.binary),RIZ.tax$OTU),])

write.csv(RIZ.all,"rawdata/RIZ.raw.wTAX.csv")
write.csv(RIZ.all.nReps,"rawdata/RIZ.raw.nReps.wTAX.csv")



####====2.0 ====####
####====3.0 ====####


