########################################################
####=.     Analysis of WP5  Antarctic data         =####
####==== Luke E. Holman====7.05.2025====################
########################################################

####====0.0 Packages====####
library(tidyverse)

#read in metadata 
metadat <- read.csv("metadata/RepMetaData.csv")
dates <- read.csv("metadata/dates.csv")
coredata <- read.csv("metadata/sitemetadata.csv")

## read in data 
mam <- read.csv("cleandata/MAM.WP5.main.clean.csv",row.names = 1)
mam.tax <- read.csv("cleandata/MAM.WP5.main.clean.tax.csv",row.names = 1)
riz <- read.csv("cleandata/RIZ.WP5.main.clean.csv",row.names = 1)
riz.tax <- read.csv("cleandata/RIZ.WP5.main.clean.tax.csv",row.names = 1)
euk <- read.csv("cleandata/EUK.WP5.main.clean.csv",row.names = 1)
euk.tax <- read.csv("cleandata/EUK.WP5.main.clean.tax.csv",row.names = 1)
euk.pr2 <- read.csv("taxonomy/EUK.pr2.csv",row.names = 1)

## some functions


CountTable <- function(in.taxonomy,in.data,output="Count",some.unassigned=T){
  if(length(in.taxonomy)!=length(in.data[,1])){stop("Dataframe and corresponding taxonomy are not the same length")}
  in.taxonomy[is.na(in.taxonomy)] <- ""
  out.dat <- as.data.frame(matrix(ncol=length(in.data[1,]),nrow=length(unique(in.taxonomy))))
  rownames(out.dat) <- sort(unique(in.taxonomy))
  colnames(out.dat) <- colnames(in.data)    
  out.dat.abundance <- out.dat
  for (sample in 1:length(in.data[1,])){
    out.dat[,sample] <- table(in.taxonomy[in.data[,sample]>0])[match(sort(unique(in.taxonomy)),names(table(in.taxonomy[in.data[,sample]>0])))]
    out.dat.abundance[,sample] <- aggregate(in.data[,sample], by=list(Category=in.taxonomy), FUN=sum)[,2]
  }
  out.dat[is.na(out.dat)] <- 0
  if(some.unassigned==T){rownames(out.dat)[1] <- "Unassigned"}
  if(output=="Count"){return(out.dat)}else if(
    output=="Abundance"){return(out.dat.abundance)}
}

minAbundance <- function(inputtable = NA, minAbun = 0.01) {
  others <- rep(0, ncol(inputtable))
  
  for (col in 1:ncol(inputtable)) {
    threshold = sum(inputtable[, col]) * minAbun
    below_threshold_indices = inputtable[, col] < threshold
    others[col] = sum(inputtable[below_threshold_indices, col])
    inputtable[below_threshold_indices, col] = 0
  }
  
  inputtable <- rbind(inputtable, others)
  rownames(inputtable)[nrow(inputtable)] = "Others"
  
  # Remove rows that sum to zero
  inputtable <- inputtable[rowSums(inputtable) != 0, ]
  
  return(inputtable)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
NrepsMaker <- function(INdataframe,vector){
  ##write these checks
  #check the dataframe is a dataframe
  if(!is.data.frame(INdataframe)){stop("Input dataframe doesn't look like a dataframe")}
  #check the vector is a vector
  if(!is.vector(vector)){stop("Input vector doesn't look like a vector")}
  #check the dataframe contains the vector
  ## TO DO
  #make a new dataframe to captuire the output
  newDataFrame <- data.frame(matrix(0,nrow = length(INdataframe[,1]),ncol=length(unique(vector))))
  #name stuff
  colnames(newDataFrame) <- unique(vector)
  rownames(newDataFrame) <- rownames(INdataframe)
  #make it binary
  INdataframe[INdataframe<1] <- 0
  INdataframe[INdataframe>0] <- 1
  #loop over all the samples with replicates, summing according to the vector
  for (column in 1:length(newDataFrame[1,])){
    # this if statement checks in case there is only one replicate remaining (this sometimes happens in controls)
    if(is.vector(INdataframe[,vector %in% colnames(newDataFrame)[column]])){
      newDataFrame[,column]  <- INdataframe[,vector %in% colnames(newDataFrame)[column]]}else{
        newDataFrame[,column] <- rowSums(INdataframe[,vector %in% colnames(newDataFrame)[column]])}
  }
  return(newDataFrame)

## Taoxnomic overview

getPalette = colorRampPalette(brewer.pal(9, "Set3"))

euk.pr2.2 <- euk.pr2[match(rownames(euk),euk.pr2$X.1),]


#GC1
EUK.tax.a <- minAbundance(CountTable(as.character(euk.pr2.2$Family),euk,output = "Abundance"),minAbun=0.01)
EUK.GC1.tax.c <- minAbundance(CountTable(as.character(EUK.tax.PR2$Family),EUK.GC1.avr,output = "Count"),minAbun=0.01)
row.names(EUK.tax.a)[1] <- "Unknown"
EUK.tax.a.2 <- as.matrix(prop.table(as.matrix(EUK.tax.a[,sort(colnames(EUK.tax.a.2))]),margin = 2))
barplot(EUK.tax.a.2,col=rev(getPalette(dim(EUK.tax.a.2)[1])))

sort(colnames(EUK.tax.a.2))


pdf("figures/fig1/EUK.GC1.tax.Family.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(EUK.GC1.tax.a,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.GC1.tax.a)[1])),ylab="Read Abundance",names.arg=ages$mean[match(colnames(EUK.GC1.tax.a),ages$ID2)])
legend(56,1,rev(rownames(EUK.GC1.tax.a)),fill=getPalette(dim(EUK.GC1.tax.a)[1]),cex=0.5,bty = "n",y.intersp=0.75)
barplot(EUK.GC1.tax.c,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.GC1.tax.c)[1])),ylab="ASV Counts",names.arg=ages$mean[match(colnames(EUK.GC1.tax.c),ages$ID2)])
legend(56,1,rev(rownames(EUK.GC1.tax.c)),fill=getPalette(dim(EUK.GC1.tax.c)[1]),cex=0.5,bty = "n",y.intersp=0.75)
dev.off()





##### don't trust this nREPs plotting - start from scratch here

NrepsMaker()


#### Lets start with some salps and krill

euk.salkril <- euk[euk.tax$family=="Euphausiidae"|
                     euk.tax$family=="Salpidae",]
euk.salkril <- euk.salkril[!is.na(euk.salkril[,1]),]
euk.salkril.tax <- euk.tax[match(rownames(euk.salkril),euk.tax$X.1),]

write.csv(cbind(euk.salkril,euk.salkril.tax),"cleandata/clean.EUK.salkril.csv")


## Pull in manually annotated data 
euk.salkril.man <- read.csv("cleandata/clean.EUK.salkril.byhand.csv",row.names = 1)


# divide by original read number to get proportion 
euk.salkril.man[,1:992] <- sweep(euk.salkril.man[,1:992],2,colSums(euk),FUN="/")
euk.salkril.man[euk.salkril.man=="NaN"] <- 0
euk.salkril.man[is.na(euk.salkril.man)] <- 0


#Aggregate same species
datain <- cbind(euk.salkril.man[,1:992],"manualSpecies"=euk.salkril.man$Assignment)
euk.salkril.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(euk.salkril.sum) <- euk.salkril.sum$manualSpecies
euk.salkril.sum <- euk.salkril.sum[,-1]
euk.salkril.sum <- euk.salkril.sum[!rownames(euk.salkril.sum)=="",]

## change names
colnames(euk.salkril.sum) <- gsub("_", ".", 
                                  gsub("^EUK_", "", 
                                       gsub("^EUK\\.AQ\\.", "", colnames(euk.salkril.sum))
                                  )
)

## Lets make an nREPs dataset 
euk.salkril.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum)))),
                                       nrow=length(euk.salkril.sum[,1])))
colnames(euk.salkril.nReps) <- unique(gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum)))
row.names(euk.salkril.nReps) <- row.names(euk.salkril.sum)
euk.salkril.allReads <- euk.salkril.nReps

euk.salkril.binary <- euk.salkril.sum
euk.salkril.binary[euk.salkril.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum))))[loopcol]
  euk.salkril.nReps[,loopcol] <- rowSums(euk.salkril.binary[index==loopname])
  euk.salkril.allReads[,loopcol] <- rowSums(euk.salkril.sum[index==loopname])
}

#Now make a long dataset 

df <- tibble(spp = rownames(euk.salkril.nReps), euk.salkril.nReps)

euk.salkril.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(euk.salkril.allReads), euk.salkril.allReads)

euk.salkril.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


euk.salkril.long <- cbind(euk.salkril.allReads.long,"nReps"=euk.salkril.nReps.long$value)


#Plots
#change the values a bit so they match the dates
euk.salkril.long$depth[euk.salkril.long$depth==0] <- "0.5"
euk.salkril.long$depth <- as.numeric(euk.salkril.long$depth)
euk.salkril.long$site <- toupper(substr(euk.salkril.long$site,1,2))

#add in cal dates 
paste(euk.salkril.long$site,euk.salkril.long$depth)
paste(substr(dates$Site,1,2),dates$SampleID)

euk.salkril.long$calYear <- dates$CalYear[match(paste(euk.salkril.long$site,euk.salkril.long$depth),
                                                paste(toupper(substr(dates$Site,1,2)),dates$SampleID))]

pdf("figures/euk.salkril.v1.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))
palette( c("#D55E00","#cc79a7","#f0e442","#009e73","#0072B2"))

for (species in unique(euk.salkril.long$spp)){
  
  plot(euk.salkril.long$calYear[euk.salkril.long$spp==species],
       jitter(euk.salkril.long$nReps[euk.salkril.long$spp==species]),
       col=as.factor(euk.salkril.long$site[euk.salkril.long$spp==species]),
       cex=if(species=="Ihlea"){1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*30}else 
         if(species=="Salpa"){1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*300}else 
         {1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*3000},
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=16)
  
  points(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],
         rep(0,length(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0])),
         pch=4,
         col=factor(euk.salkril.long$site[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],levels=c("BB","MC","P5","P6","SC")))
  
  
}

dev.off()





###### OLD




####====1.0 Let's look at some whales! ====####

mam <- read.csv("cleanedData/clean.MAM.raw.csv.dada2.csv",row.names = 1)
riz <- read.csv("cleanedData/clean.RIZ.raw.csv.dada2.csv",row.names = 1)
euk <- read.csv("cleanedData/clean.EUK.raw.csv.dada2.csv",row.names = 1)

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


# Now we manipulate the data to make a long dataset 
#mam



df <- tibble(spp = rownames(mam.whales.nReps), mam.whales.nReps)

mam.whales.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(mam.whales.allReads), mam.whales.allReads)

mam.whales.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


mam.whales.long <- cbind(mam.whales.allReads.long,"nReps"=mam.whales.nReps.long$value)

#Riz

df <- tibble(spp = rownames(riz.whales.nReps), riz.whales.nReps)

riz.whales.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(riz.whales.allReads), riz.whales.allReads)

riz.whales.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


riz.whales.long <- cbind(riz.whales.allReads.long,"nReps"=riz.whales.nReps.long$value)


#Plots
#change the values a bit so they match the dates
mam.whales.long$depth[mam.whales.long$depth==0] <- "0.5"
mam.whales.long$depth <- as.numeric(mam.whales.long$depth)
mam.whales.long$site <- toupper(substr(mam.whales.long$site,1,2))

#add in cal dates 
mam.whales.long$calYear <- dates$CalYear[match(paste(mam.whales.long$site,mam.whales.long$depth),
                paste(dates$Site,dates$Depth))]

pdf("figures/mam.whales.v1.pdf",width = 9,height =6.5)
par(mfrow=c(2,3))

for (species in unique(mam.whales.long$spp)){
  
plot(mam.whales.long$calYear[mam.whales.long$spp==species],
     mam.whales.long$nReps[mam.whales.long$spp==species],
     col=as.factor(mam.whales.long$site[mam.whales.long$spp==species]),
     cex=1+mam.whales.long$value[mam.whales.long$spp==species]/mam.whales.long$nReps[mam.whales.long$spp==species],
     main=species,
     xlab="Date",
     ylab="nPCR reps",
     ylim=c(0,4),
     pch=16)
  abline(v=1984,col="darkred",lty=2)
  
points(mam.whales.long$calYear[mam.whales.long$spp==species&mam.whales.long$nReps==0],
       rep(0,length(mam.whales.long$calYear[mam.whales.long$spp==species&mam.whales.long$nReps==0])),
       pch=4,
       col=as.factor(mam.whales.long$site[mam.whales.long$spp==species&mam.whales.long$nReps==0]))
}

dev.off()


#change the values a bit so they match the dates
riz.whales.long$depth[riz.whales.long$depth==0] <- "0.5"
riz.whales.long$depth <- as.numeric(riz.whales.long$depth)
riz.whales.long$site <- toupper(substr(riz.whales.long$site,1,2))

#add in cal dates 
riz.whales.long$calYear <- dates$CalYear[match(paste(riz.whales.long$site,riz.whales.long$depth),
                                               paste(dates$Site,dates$Depth))]
pdf("figures/riz.whales.v1.pdf",width = 9,height =6.5)

par(mfrow=c(2,3))

for (species in unique(riz.whales.long$spp)){
  
  plot(riz.whales.long$calYear[riz.whales.long$spp==species],
       riz.whales.long$nReps[riz.whales.long$spp==species],
       col=as.factor(riz.whales.long$site[riz.whales.long$spp==species]),
       cex=1+riz.whales.long$value[riz.whales.long$spp==species]/riz.whales.long$nReps[riz.whales.long$spp==species],
       main=species,
       xlab="Date",
       ylab="nPCR",
       ylim=c(0,3),
       pch=16)
  abline(v=1984,col="darkred",lty=2)
  
  points(riz.whales.long$calYear[riz.whales.long$spp==species&riz.whales.long$nReps==0],
         rep(0,length(riz.whales.long$calYear[riz.whales.long$spp==species&riz.whales.long$nReps==0])),
         pch=4,
         col=as.factor(riz.whales.long$site[riz.whales.long$spp==species&riz.whales.long$nReps==0]))
  
  
}

dev.off()


### Now lets combine the datasets 


com.whales.long <- bind_rows(mam.whales.long, riz.whales.long) %>%
  group_by(spp, site, depth, calYear) %>%
  summarise(value = sum(value), nReps = sum(nReps))


pdf("figures/com.mam.whales.v1.pdf",width = 11,height =6.5)

par(mfrow=c(2,3))
palette(c("#009e73","#D55E00","#0072B2"))

for (species in unique(com.whales.long$spp)){
  
  plot(com.whales.long$calYear[com.whales.long$spp==species],
       jitter(com.whales.long$nReps[com.whales.long$spp==species]),
       col=as.factor(com.whales.long$site[com.whales.long$spp==species]),
       cex=1+com.whales.long$value[com.whales.long$spp==species]/com.whales.long$nReps[com.whales.long$spp==species],
       main=species,
       xlab="Date",
       ylab="nPCR",
       ylim=c(0,5),
       pch=16)
  abline(v=1984,col="darkred",lty=2)
  
  points(com.whales.long$calYear[com.whales.long$spp==species&com.whales.long$nReps==0],
         rep(0,length(com.whales.long$calYear[com.whales.long$spp==species&com.whales.long$nReps==0])),
         pch=4,
         col=as.factor(com.whales.long$site[com.whales.long$spp==species&com.whales.long$nReps==0]))
  
  
}

dev.off()


index <- gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum)))

unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.whales.sum))))


####====2.0 What about the seals====####

mam <- read.csv("cleanedData/clean.MAM.raw.csv.dada2.csv",row.names = 1)
riz <- read.csv("cleanedData/clean.RIZ.raw.csv.dada2.csv",row.names = 1)

mam.seals <- mam[mam$family=="Otariidae"|
                    mam$family=="Phocidae",]

mam.seals<- mam.seals[!is.na(mam.seals$assignmentQual),]

riz.seals <- riz[riz$family=="Otariidae"|
                    riz$family=="Phocidae",]

riz.seals <- riz.seals[!is.na(riz.seals$assignmentQual),]

write.csv(mam.seals,"cleanedData/clean.MAM.seals.csv")
write.csv(riz.seals,"cleanedData/clean.RIZ.seals.csv")

#Now we check all those seals manually and pull in the checked data  
mam.seals.man <- read.csv("cleanedData/clean.MAM.seals.manual.csv",row.names = 1)

riz.seals.man <- read.csv("cleanedData/clean.RIZ.seals.manual.csv",row.names = 1)


# divide by original read number to get proportion 
riz.seals.man[,1:358] <- sweep(riz.seals.man[,1:358],2,colSums(riz[,1:358]),FUN="/")
riz.seals.man[riz.seals.man=="NaN"] <- 0
riz.seals.man[is.na(riz.seals.man)] <- 0

#same for mam
mam.seals.man[,1:359] <- sweep(mam.seals.man[,1:359],2,colSums(mam[,1:359]),FUN="/")
mam.seals.man[mam.seals.man=="NaN"] <- 0
mam.seals.man[is.na(mam.seals.man)] <- 0

#Aggregate same species
datain <- cbind(mam.seals.man[,1:358],"manualSpecies"=mam.seals.man$manualSpecies)
mam.seals.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(mam.seals.sum) <- mam.seals.sum$manualSpecies
mam.seals.sum <- mam.seals.sum[,-1]


datain <- cbind(riz.seals.man[,1:358],"manualSpecies"=riz.seals.man$manualSpecies)
riz.seals.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(riz.seals.sum) <- riz.seals.sum$manualSpecies
riz.seals.sum <- riz.seals.sum[,-1]



#make a nReps dataset 

#MAM

mam.seals.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum))))),
                                      nrow=length(mam.seals.sum[,1])))
colnames(mam.seals.nReps) <- unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum))))
row.names(mam.seals.nReps) <- row.names(mam.seals.sum)
mam.seals.allReads <- mam.seals.nReps

mam.seals.binary <- mam.seals.sum
mam.seals.binary[mam.seals.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum)))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum)))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum))))[loopcol]
  mam.seals.nReps[,loopcol] <- rowSums(mam.seals.binary[index==loopname])
  mam.seals.allReads[,loopcol] <- rowSums(mam.seals.sum[index==loopname])
}

#RIZ


riz.seals.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(riz.seals.sum))))),
                                      nrow=length(riz.seals.sum[,1])))
colnames(riz.seals.nReps) <- unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.seals.sum))))
row.names(riz.seals.nReps) <- row.names(riz.seals.sum)
riz.seals.allReads <- riz.seals.nReps

riz.seals.binary <- riz.seals.sum
riz.seals.binary[riz.seals.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.seals.sum)))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.seals.sum)))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("RIZ.AQ.","",colnames(riz.seals.sum))))[loopcol]
  riz.seals.nReps[,loopcol] <- rowSums(riz.seals.binary[index==loopname])
  riz.seals.allReads[,loopcol] <- rowSums(riz.seals.sum[index==loopname])
}


# Now we manipulate the data to make a long dataset 
#mam



df <- tibble(spp = rownames(mam.seals.nReps), mam.seals.nReps)

mam.seals.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(mam.seals.allReads), mam.seals.allReads)

mam.seals.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


mam.seals.long <- cbind(mam.seals.allReads.long,"nReps"=mam.seals.nReps.long$value)

#Riz

df <- tibble(spp = rownames(riz.seals.nReps), riz.seals.nReps)

riz.seals.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(riz.seals.allReads), riz.seals.allReads)

riz.seals.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


riz.seals.long <- cbind(riz.seals.allReads.long,"nReps"=riz.seals.nReps.long$value)


#Plots
#change the values a bit so they match the dates
mam.seals.long$depth[mam.seals.long$depth==0] <- "0.5"
mam.seals.long$depth <- as.numeric(mam.seals.long$depth)
mam.seals.long$site <- toupper(substr(mam.seals.long$site,1,2))

#add in cal dates 
mam.seals.long$calYear <- dates$CalYear[match(paste(mam.seals.long$site,mam.seals.long$depth),
                                               paste(dates$Site,dates$Depth))]

pdf("figures/mam.seals.v1.pdf",width = 9,height =6.5)
par(mfrow=c(2,3))

for (species in unique(mam.seals.long$spp)){
  
  plot(mam.seals.long$calYear[mam.seals.long$spp==species],
       mam.seals.long$nReps[mam.seals.long$spp==species],
       col=as.factor(mam.seals.long$site[mam.seals.long$spp==species]),
       cex=1+mam.seals.long$value[mam.seals.long$spp==species]/mam.seals.long$nReps[mam.seals.long$spp==species],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,4),
       pch=16)
  
  points(mam.seals.long$calYear[mam.seals.long$spp==species&mam.seals.long$nReps==0],
         rep(0,length(mam.seals.long$calYear[mam.seals.long$spp==species&mam.seals.long$nReps==0])),
         pch=4,
         col=as.factor(mam.seals.long$site[mam.seals.long$spp==species&mam.seals.long$nReps==0]))
}

dev.off()


#change the values a bit so they match the dates
riz.seals.long$depth[riz.seals.long$depth==0] <- "0.5"
riz.seals.long$depth <- as.numeric(riz.seals.long$depth)
riz.seals.long$site <- toupper(substr(riz.seals.long$site,1,2))

#add in cal dates 
riz.seals.long$calYear <- dates$CalYear[match(paste(riz.seals.long$site,riz.seals.long$depth),
                                               paste(dates$Site,dates$Depth))]
pdf("figures/riz.seals.v1.pdf",width = 9,height =6.5)

par(mfrow=c(2,3))

for (species in unique(riz.seals.long$spp)){
  
  plot(riz.seals.long$calYear[riz.seals.long$spp==species],
       riz.seals.long$nReps[riz.seals.long$spp==species],
       col=as.factor(riz.seals.long$site[riz.seals.long$spp==species]),
       cex=1+riz.seals.long$value[riz.seals.long$spp==species]/riz.seals.long$nReps[riz.seals.long$spp==species],
       main=species,
       xlab="Date",
       ylab="nPCR",
       ylim=c(0,3),
       pch=16)
  abline(v=1984,col="darkred",lty=2)
  
  points(riz.seals.long$calYear[riz.seals.long$spp==species&riz.seals.long$nReps==0],
         rep(0,length(riz.seals.long$calYear[riz.seals.long$spp==species&riz.seals.long$nReps==0])),
         pch=4,
         col=as.factor(riz.seals.long$site[riz.seals.long$spp==species&riz.seals.long$nReps==0]))
  
  
}

dev.off()


### Now lets combine the datasets 


com.seals.long <- bind_rows(mam.seals.long, riz.seals.long) %>%
  group_by(spp, site, depth, calYear) %>%
  summarise(value = sum(value), nReps = sum(nReps))


pdf("figures/com.mam.seals.v1.pdf",width = 11,height =6.5)

par(mfrow=c(2,3))
palette(c("#009e73","#D55E00","#0072B2"))

for (species in unique(com.seals.long$spp)){
  
  plot(com.seals.long$calYear[com.seals.long$spp==species],
       jitter(com.seals.long$nReps[com.seals.long$spp==species]),
       col=as.factor(com.seals.long$site[com.seals.long$spp==species]),
       cex=1+com.seals.long$value[com.seals.long$spp==species]/com.seals.long$nReps[com.seals.long$spp==species],
       main=species,
       xlab="Date",
       ylab="nPCR",
       ylim=c(0,6),
       pch=16)
  
  points(com.seals.long$calYear[com.seals.long$spp==species&com.seals.long$nReps==0],
         rep(0,length(com.seals.long$calYear[com.seals.long$spp==species&com.seals.long$nReps==0])),
         pch=4,
         col=factor(com.seals.long$site[com.seals.long$spp==species&com.seals.long$nReps==0],levels=c("BB","MC","SC")))
  
  
}

dev.off()


index <- gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum)))

unique(gsub("(.*).[0-9]$","\\1",gsub("MAM.AQ.","",colnames(mam.seals.sum))))


####====3.0 Polarella ===####

euk.polarella <- euk[euk$genus=="Polarella" & euk$assignmentQual=="High",]
euk.polarella <- euk.polarella[!is.na(euk.polarella$EUK.AQ.BB3.0.5.1),]
euk.polarella.clean <- colSums(euk.polarella[,1:358])
euk.polarella.clean.p <- euk.polarella.clean/colSums(euk[,1:358])
euk.dinos.h <- euk[euk$class=="Dinophyceae" & (euk$assignmentQual=="High"|euk$assignmentQual=="High"),]
euk.dinos.h <- euk.dinos.h[!is.na(euk.dinos.h$EUK.AQ.BB3.0.5.1),]
euk.dinos.l <- euk[euk$class=="Dinophyceae",]
euk.dinos.l <- euk.dinos.l[!is.na(euk.dinos.l$EUK.AQ.BB3.0.5.1),]

euk.polarella.dat <- data.frame("raw.reads"=euk.polarella.clean,
                                "prop.reads"=euk.polarella.clean.p,
                                "prop.reads.dino"=euk.polarella.clean/colSums(euk.dinos.h[,1:358]),
                                "prop.reads.dino.l"=euk.polarella.clean/colSums(euk.dinos.l[,1:358]),
                                "site"=toupper(substr(names(euk.polarella.clean),8,9)),
                                "depth"=as.numeric(sapply(strsplit(names(euk.polarella.clean), "\\."), function(x) x[4])))
euk.polarella.dat$depth[euk.polarella.dat$depth=="0"] <- 0.5
euk.polarella.dat$dates <- dates$CalYear[match(paste(euk.polarella.dat$site,euk.polarella.dat$depth),
                    paste(dates$Site,dates$Depth))]

pdf("figures/p.glacialis.pdf",height = 8,width = 7)
par(mfrow=c(3,1))
palette(c("#009e73","#D55E00","#0072B2"))
count <- 1
for (site in unique(euk.polarella.dat$site)){
plot(euk.polarella.dat$dates[euk.polarella.dat$site==site],
     euk.polarella.dat$prop.reads[euk.polarella.dat$site==site],
     col=count,pch=16,ylab="Proportion P.glacialis",xlab="Date (CE)",main=site)
  count <- count+1
}
dev.off()

pdf("figures/p.glacialis.dino.pdf",height = 8,width = 7)
par(mfrow=c(3,1))
palette(c("#009e73","#D55E00","#0072B2"))
count <- 1
for (site in unique(euk.polarella.dat$site)){
  plot(euk.polarella.dat$dates[euk.polarella.dat$site==site],
       euk.polarella.dat$prop.reads.dino[euk.polarella.dat$site==site],
       col=count,pch=16,ylab="Proportion P.glacialis",xlab="Date (CE)",main=site)
  count <- count+1
}
dev.off()

pdf("figures/p.glacialis.dino.l.pdf",height = 8,width = 7)
par(mfrow=c(3,1))
palette(c("#009e73","#D55E00","#0072B2"))
count <- 1
for (site in unique(euk.polarella.dat$site)){
  plot(euk.polarella.dat$dates[euk.polarella.dat$site==site],
       euk.polarella.dat$prop.reads.dino.l[euk.polarella.dat$site==site],
       col=count,pch=16,ylab="Proportion P.glacialis",xlab="Date (CE)",main=site)
  count <- count+1
}
dev.off()

####====3.1 fish====####
euk.fish <- euk["ASV_433",]

euk.fish.dat <- data.frame("raw.reads"=unlist(euk.fish[1,1:358]),
                                "prop.reads"=unlist(euk.fish[1,1:358]/colSums(euk[,1:358])),
                                "site"=toupper(substr(names(euk.fish[,1:358]),8,9)),
                                "depth"=as.numeric(sapply(strsplit(names(euk.fish[,1:358]), "\\."), function(x) x[4])))
euk.fish.dat$binary <- ifelse(euk.fish.dat$raw.reads>0,1,0)

euk.fish.dat$depth[euk.fish.dat$depth=="0"] <- 0.5
euk.fish.dat$dates <- dates$CalYear[match(paste(euk.fish.dat$site,euk.fish.dat$depth),
                                               paste(dates$Site,dates$Depth))]

pdf("figures/fish.pdf",height = 8,width = 7)
par(mfrow=c(3,1))
palette(c("#009e73","#D55E00","#0072B2"))
count <- 1
for (site in unique(euk.fish.dat$site)){
  plot(euk.fish.dat$dates[euk.fish.dat$site==site],
       euk.fish.dat$prop.reads[euk.fish.dat$site==site],
       col=count,pch=16,ylab="Proportion fish",xlab="Date (CE)",main=site)
  count <- count+1
}
dev.off()

summed_binary <-euk.fish.dat %>%
  group_by(site, depth) %>%
  summarise(sum_binary = sum(binary))

plot(dates$CalYear[match(paste(summed_binary$site,summed_binary$depth),
                    paste(dates$Site,dates$Depth))],
     summed_binary$sum_binary,
     col=as.numeric(as.factor(summed_binary$site)),pch=16)


####====3.2 Salps and Krill====####

euk <- read.csv("cleanedData/clean.EUK.raw.csv.dada2.csv",row.names = 1)

euk.salkril <- euk[euk$family=="Euphausiidae"|
                   euk$family=="Salpidae",]

euk.salkril <- euk.salkril[!is.na(euk.salkril$assignmentQual),]
write.csv(euk.salkril,"cleanedData/clean.EUK.salkril.csv")


## Pull in manually annotated data 

euk.salkril.man <- read.csv("cleanedData/clean.EUK.salkril.manual.csv",row.names = 1)


# divide by original read number to get proportion 
euk.salkril.man[,1:358] <- sweep(euk.salkril.man[,1:358],2,colSums(euk[,1:358]),FUN="/")
euk.salkril.man[euk.salkril.man=="NaN"] <- 0
euk.salkril.man[is.na(euk.salkril.man)] <- 0

euk.salkril.man[,1:358] <- sweep(euk.salkril.man[,1:358],2,colSums(euk[,1:358]),FUN="/")
euk.salkril.man[euk.salkril.man=="NaN"] <- 0
euk.salkril.man[is.na(euk.salkril.man)] <- 0




#Aggregate same species
datain <- cbind(euk.salkril.man[,1:358],"manualSpecies"=euk.salkril.man$manualSpecies)
euk.salkril.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(euk.salkril.sum) <- euk.salkril.sum$manualSpecies
euk.salkril.sum <- euk.salkril.sum[,-1]
euk.salkril.sum <- euk.salkril.sum[!rownames(euk.salkril.sum)=="",]

## Lets make an nREPs dataset 
euk.salkril.nReps <- data.frame(matrix(ncol=length(unique(gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum))))),
                                     nrow=length(euk.salkril.sum[,1])))
colnames(euk.salkril.nReps) <- unique(gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum))))
row.names(euk.salkril.nReps) <- row.names(euk.salkril.sum)
euk.salkril.allReads <- euk.salkril.nReps

euk.salkril.binary <- euk.salkril.sum
euk.salkril.binary[euk.salkril.binary>0] <- 1

index <- gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum)))

for (loopcol in 1:length(unique(gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum)))))){
  loopname <- unique(gsub("(.*).[0-9]$","\\1",gsub("EUK.AQ.","",colnames(euk.salkril.sum))))[loopcol]
  euk.salkril.nReps[,loopcol] <- rowSums(euk.salkril.binary[index==loopname])
  euk.salkril.allReads[,loopcol] <- rowSums(euk.salkril.sum[index==loopname])
}

#Now make a long dataset 

df <- tibble(spp = rownames(euk.salkril.nReps), euk.salkril.nReps)

euk.salkril.nReps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


df <- tibble(spp = rownames(euk.salkril.allReads), euk.salkril.allReads)

euk.salkril.allReads.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") %>% 
  separate(column, into = c("site", "depth"), sep = "[.]") 


euk.salkril.long <- cbind(euk.salkril.allReads.long,"nReps"=euk.salkril.nReps.long$value)


#Plots
#change the values a bit so they match the dates
euk.salkril.long$depth[euk.salkril.long$depth==0] <- "0.5"
euk.salkril.long$depth <- as.numeric(euk.salkril.long$depth)
euk.salkril.long$site <- toupper(substr(euk.salkril.long$site,1,2))

#add in cal dates 
euk.salkril.long$calYear <- dates$CalYear[match(paste(euk.salkril.long$site,euk.salkril.long$depth),
                                              paste(dates$Site,dates$Depth))]

pdf("figures/euk.salkril.v1.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))

for (species in unique(euk.salkril.long$spp)){
  
  plot(euk.salkril.long$calYear[euk.salkril.long$spp==species],
       jitter(euk.salkril.long$nReps[euk.salkril.long$spp==species]),
       col=as.factor(euk.salkril.long$site[euk.salkril.long$spp==species]),
       cex=if(species=="Ihlea racovitzai"){1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*30}else 
         if(species=="Salpa thompsoni"){1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*300}else 
           {1+(euk.salkril.long$value[euk.salkril.long$spp==species]/euk.salkril.long$nReps[euk.salkril.long$spp==species])*3000},
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=16)
  
  points(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],
         rep(0,length(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0])),
         pch=4,
         col=factor(euk.salkril.long$site[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],levels=c("BB","MC","SC")))

  
  }

dev.off()

#SALP KRILL RATIO


krill <- euk.salkril.long$value[euk.salkril.long$spp=="Krill"]/euk.salkril.long$nReps[euk.salkril.long$spp=="Krill"]


salps <- euk.salkril.long$value[euk.salkril.long$spp=="Salpa thompsoni"]/euk.salkril.long$nReps[euk.salkril.long$spp=="Salpa thompsoni"]
salps[is.na(salps)] <- 0 
salps2 <- euk.salkril.long$value[euk.salkril.long$spp=="Ihlea racovitzai"]/euk.salkril.long$nReps[euk.salkril.long$spp=="Ihlea racovitzai"]
allsalps <- salps+salps2

pdf("figures/euk.salkril.ratio.v1.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))

plot(euk.salkril.long$calYear[euk.salkril.long$spp==species],
     log10(krill/salps),
     col=as.factor(euk.salkril.long$site[euk.salkril.long$spp==species]),
     main="Salpa thompsoni",
     xlab="Date",
     ylab="Krill-Salp Ratio (log10)",
     pch=16)

points(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],
       rep(0,length(euk.salkril.long$calYear[euk.salkril.long$spp==species&euk.salkril.long$nReps==0])),
       pch=4,
       col=factor(euk.salkril.long$site[euk.salkril.long$spp==species&euk.salkril.long$nReps==0],levels=c("BB","MC","SC")))

plot(euk.salkril.long$calYear[euk.salkril.long$spp==species],
     log10(krill/salps2),
     col=as.factor(euk.salkril.long$site[euk.salkril.long$spp==species]),
     main="Ihlea racovitzai",
     xlab="Date",
     ylab="Krill-Salp Ratio (log10)",
     pch=16)

plot(euk.salkril.long$calYear[euk.salkril.long$spp==species],
     log10(krill/allsalps),
     col=as.factor(euk.salkril.long$site[euk.salkril.long$spp==species]),
     main="Both Salps",
     xlab="Date",
     ylab="Krill-Salp Ratio (log10)",
     pch=16)

dev.off()


##



####====4.0 ====####


