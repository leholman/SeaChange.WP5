#########################################################
####=An initial analysis of the new Antarctic data  =####
####==== Luke E. Holman====13.06.2023====################
#########################################################

####====0.0 Packages====####

#read in metadata 

metadat <- read.csv("metadata.csv")


####====1.0 Let's look at some whales! ====####

mam <- read.csv("cleanedData/clean.MAM.raw.csv.dada2.csv",row.names = 1)
riz <- read.csv("cleanedData/clean.RIZ.raw.csv.dada2.csv",row.names = 1)

mam.whales <- mam[mam$family=="Balaenopteridae"|
                    mam$family=="Balaenidae"|
                    mam$family=="Delphinidae",]

mam.whales <- mam.whales[!is.na(mam.whales$assignmentQual),]

riz.whales <- riz[riz$family=="Balaenopteridae"|
                    riz$family=="Balaenidae"|
                    riz$family=="Delphinidae",]

riz.whales <- riz.whales[!is.na(riz.whales$assignmentQual),]

write.csv(mam.whales,"cleanedData/clean.MAM.whales.csv")
write.csv(riz.whales,"cleanedData/clean.RIZ.whales.csv")

#Now we check all those whales manually and pull in the checked data  
mam.whales.man <- read.csv("cleanedData/clean.MAM.whales.manual.csv",row.names = 1)
riz.whales.man <- read.csv("cleanedData/clean.RIZ.whales.manual.csv",row.names = 1)


# divide by original read number to get proportion 
riz.whales.man[,1:358] <- sweep(riz.whales.man[,1:358],2,colSums(riz[,1:358]),FUN="/")
riz.whales.man[riz.whales.man=="NaN"] <- 0
riz.whales.man[is.na(riz.whales.man)] <- 0

#same for mam
mam.whales.man[,1:359] <- sweep(mam.whales.man[,1:359],2,colSums(mam[,1:359]),FUN="/")
mam.whales.man[mam.whales.man=="NaN"] <- 0
mam.whales.man[is.na(mam.whales.man)] <- 0

#Aggregate same species only in mam as there are duplicates 
datain <- cbind(mam.whales.man[,1:358],"manualSpecies"=mam.whales.man$manualSpecies)
mam.whales.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(mam.whales.sum) <- mam.whales.sum$manualSpecies
mam.whales.sum <- mam.whales.sum[,-1]


#make a nReps dataset 

#MAM

mam.whales.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum))))),
                                      nrow=length(mam.whales.sum[,1])))
colnames(mam.whales.nReps) <- unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum))))
row.names(mam.whales.nReps) <- row.names(mam.whales.sum)
mam.whales.allReads <- mam.whales.nReps

mam.whales.binary <- mam.whales.sum
mam.whales.binary[mam.whales.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum)))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum)))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum))))[loopcol]
mam.whales.nReps[,loopcol] <- rowSums(mam.whales.binary[index==loopname])
mam.whales.allReads[,loopcol] <- rowSums(mam.whales.sum[index==loopname])
}

#RIZ

riz.whales.sum <- riz.whales.man[,1:357]
rownames(riz.whales.sum) <- riz.whales.man$manualSpecies

riz.whales.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(riz.whales.sum))))),
                                      nrow=length(riz.whales.sum[,1])))
colnames(riz.whales.nReps) <- unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.whales.sum))))
row.names(riz.whales.nReps) <- row.names(riz.whales.sum)
riz.whales.allReads <- riz.whales.nReps

riz.whales.binary <- riz.whales.sum
riz.whales.binary[riz.whales.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.whales.sum)))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.whales.sum)))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.whales.sum))))[loopcol]
  riz.whales.nReps[,loopcol] <- rowSums(riz.whales.binary[index==loopname])
  riz.whales.allReads[,loopcol] <- rowSums(riz.whales.sum[index==loopname])
}



##### NOW WE MELT THE DATA READY FOR PLOTTING 





index <- gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum)))

unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum))))


####====2.0 ====####
####====3.0 ====####




Balaenopteridae OR Otariidae OR Phocidae OR Delphinidae OR Balaenidae