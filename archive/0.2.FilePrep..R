
#############################################
####===Data prep before analysis WP5.2===####
####==== Luke E. Holman====01.05.2023====####
#############################################

####====0.0 Packages & Parameters====####
library("Biostrings")

#Set some variables 
minreads <- 2
items <- NULL

#Set the seed 
set.seed("123456")

#Read in metadata
metadata<-read.csv("metadata.csv") 
metadata$rep <- gsub(".*-([0-9])$","\\1",metadata$SampleID)



####====1.0 Make a raw dataset to play with/visualise====####


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



####====2.0 Now let's try and do some filtering and make some clean datasets====####

for (dataset in list.files("rawdata",pattern="....raw.csv")){

datasetname <- substr(dataset,1,3)  
  
indata <- read.csv(paste0("rawdata/",datasetname,".raw.csv"),row.names = 1)

expSamples <- indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="Experimental"])),substring(colnames(indata),5)))]

ctlSamples <-indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="ControlN"])),substring(colnames(indata),5)))]

#Filter 1 - minimum number of reads for any ID
expSamples[expSamples< minreads] <- 0
expSamples <- expSamples[rowSums(expSamples) > 0,]

#Filter 2 - within samples OTU must appear in more than one sample (this works because there are lots of reps per site and sample)
filtersam <- expSamples
filtersam[filtersam>0 ] <- 1
filtersam <-filtersam[rowSums(filtersam) > 1,]
expSamples <- expSamples[rownames(expSamples) %in% rownames(filtersam),]

#Filter 3 -Maximum value in neg = 0 value in samples

controlsCONTAM <- ctlSamples[rowSums(ctlSamples) > 0,]
for (contamOTU in 1:length(controlsCONTAM[,1])){
  loopOTU <- row.names(controlsCONTAM[contamOTU,])
  loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
  #loopSum <- sum(as.numeric(controlsCONTAM[contamOTU,]))
  if (any(is.na(expSamples[loopOTU,]))){next}
  expSamples[loopOTU,expSamples[loopOTU,]<loopMax] <- 0
  print(paste("Cleaning contaminants",contamOTU))
}



#Reattach taxonomy and ASVs

rawSeqs <- as.character(readDNAStringSet(paste0("rawdata/ASVs/",datasetname,".DADA2.ASVs.fasta")))
Assignments <- read.csv(paste0("taxonomy/",datasetname,".parsed.csv"),row.names = 1)

CleanedOutput <- cbind(expSamples,
                       unname(rawSeqs)[match(row.names(expSamples),names(rawSeqs))],
                       Assignments[match(row.names(expSamples),Assignments$OTU),])
dir.create("cleanedData")
write.csv(CleanedOutput,paste0("cleanedData/clean.",dataset,".dada2.csv"))

}







####====3.0 ====####



