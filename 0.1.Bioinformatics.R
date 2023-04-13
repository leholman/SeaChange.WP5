
## 1.0 Demultiplex Pools & trim primers

##We need to demultiplex the data 

# Some packages we need 
library(metabarTOAD)
packageVersion("metabarTOAD")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(dplyr)
packageVersion("dplyr")


# set the directory 
setwd("1.rawreads/")
# load in the metadata
files <- read.csv("../metadata.csv")

#record th eversion of cutadapt for the analysis
temp <- paste("Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
cat(file="../log.txt", temp, append=T, sep="\n")

#primer information
Euk18S <- c("GTACACACCGCCCGTC","TGATCCTTCTGCAGGTTCACCTAC","18S","EUK")
Riaz12S <-c("TTAGATACCCCACTATGC","TAGAACAGGCTCCTCTAG","12S","RIZ")
Mamm16S <- c("CGGTTGGGGTGACCTCGGA","GCTGTTATCCCTAGGGTAACT","16S","MAM")

primerSets <- as.data.frame(rbind(Euk18S,Riaz12S,Mamm16S))

#a little function to do revcomp 
RevComp <- function(input){
  require(Biostrings)
  dna <- DNAString(input) 
  dna <- reverseComplement(dna)
  return(toString(dna))}

### Generate FASTA files for forward and reverses for demultiplexing  


Forwards <- DNAStringSet(c(paste0(files$TagF_18S[1:40],primerSets$V1[1]),
                           paste0(files$TagF_12S[1:40],primerSets$V1[2]),
                           paste0(files$TagF_16S[1:40],primerSets$V1[3])))

names(Forwards) <- c(paste0("EUK-",formatC(1:40, width = 2, format = "d", flag = "0")),
                     paste0("RIZ-",formatC(1:40, width = 2, format = "d", flag = "0")),
                     paste0("MAM-",formatC(1:40, width = 2, format = "d", flag = "0")))
                 
writeFasta(Forwards,"P.Forwards.fasta")

Reverses <- DNAStringSet(c(paste0(files$TagR_18S[1:40],primerSets$V2[1]),
                           paste0(files$TagR_12S[1:40],primerSets$V2[2]),
                           paste0(files$TagR_16S[1:40],primerSets$V2[3])))

names(Reverses) <- c(paste0("EUK-",formatC(1:40, width = 2, format = "d", flag = "0")),
                     paste0("RIZ-",formatC(1:40, width = 2, format = "d", flag = "0")),
                     paste0("MAM-",formatC(1:40, width = 2, format = "d", flag = "0")))

writeFasta(Reverses,"P.Reverses.fasta")

### Loop across each pool and output demultiplexed data
# Note sometimes you hit file limits so need to run ulimit -S -n 1000

# The entire expression we are after looks like the below 
# cutadapt -e 1 --pair-adapters --trimmed-only -j 7 -g ^file:Forwards.fasta -G ^file:Reverses.fasta -o Pool01-{name}.1.fastq.gz -p Pool01-{name}.2.fastq.gz Pool01_S1_L001_R1_001.fastq.gz Pool01_S1_L001_R2_001.fastq.gz

dir.create("S.demultiplexed",showWarnings = FALSE)
dir.create("A.demultiplexed",showWarnings = FALSE)


for (pool in unique(files$FileName)){

looppoolName <- strsplit(pool, "[_]")[[1]][1]

temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -j 7 -g ^file:P.Forwards.fasta -G ^file:P.Reverses.fasta -o S.demultiplexed/",looppoolName,"-{name}.S.R1.fastq.gz -p S.demultiplexed/",looppoolName,"-{name}.S.R2.fastq.gz ",pool,"1_001.fastq.gz ",pool,"2_001.fastq.gz"), stdout=T, stderr=T)
cat(file="log.txt", temp, append=T, sep="\n")

temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -j 7 -g ^file:P.Reverses.fasta -G ^file:P.Forwards.fasta -o A.demultiplexed/",looppoolName,"-{name}.A.R1.fastq.gz -p A.demultiplexed/",looppoolName,"-{name}.A.R2.fastq.gz ",pool,"1_001.fastq.gz ",pool,"2_001.fastq.gz"), stdout=T, stderr=T)
cat(file="log.txt", temp, append=T, sep="\n")

}

### Rename all files according to metadata 

for (loopIndex in 1:length(files[,1])){

  for (primer in c("EUK","RIZ","MAM")) {
  loopSample <- files[loopIndex,1]
 
  loopTargetfile <- paste0(strsplit(files$FileName[loopIndex], "[_]")[[1]][1],
                           "-",
                           primer,
                           "-",
                           formatC(files$TagID[loopIndex],width = 2, format = "d", flag = "0"),
                           ".")
  
  print(loopSample)
  system2("mv",args=paste0("S.demultiplexed/",loopTargetfile,"S.R1.fastq.gz S.demultiplexed/",primer,".",loopSample,".S.R1.fastq.gz"))
  system2("mv",args=paste0("S.demultiplexed/",loopTargetfile,"S.R2.fastq.gz S.demultiplexed/",primer,".",loopSample,".S.R2.fastq.gz"))
  system2("mv",args=paste0("A.demultiplexed/",loopTargetfile,"A.R1.fastq.gz A.demultiplexed/",primer,".",loopSample,".A.R1.fastq.gz"))
  system2("mv",args=paste0("A.demultiplexed/",loopTargetfile,"A.R2.fastq.gz A.demultiplexed/",primer,".",loopSample,".A.R2.fastq.gz"))
  }
}

### Loop over each marker stripping RevComp primer from F & R
dir.create("S.stripped",showWarnings = FALSE)
dir.create("A.stripped",showWarnings = FALSE)

for (primer in c("EUK","RIZ","MAM")){
  for (file in list.files("S.demultiplexed",pattern=paste0(primer,".*R1.fastq.gz$"))){
  
  filename <- gsub("(.*R)[12].fastq.gz","\\1",file)
  temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -a ",
                                         RevComp(primerSets$V2[primerSets$V4==primer])," -A ",
                                         RevComp(primerSets$V1[primerSets$V4==primer]),
                                         " -j 7 -o S.stripped/",file," -p S.stripped/",paste0(filename,"2.fastq.gz "),
                                         "S.demultiplexed/",file," S.demultiplexed/",paste0(filename,"2.fastq.gz")), stdout=T, stderr=T)
  cat(file="log.txt", temp, append=T, sep="\n")
  
  }
  for (file in list.files("A.demultiplexed",pattern=paste0(primer,".*R1.fastq.gz$"))){
    filename <- gsub("(.*R)[12].fastq.gz","\\1",file)
  temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -a ",
                                         RevComp(primerSets$V1[primerSets$V4==primer])," -A ",
                                         RevComp(primerSets$V2[primerSets$V4==primer]),
                                         " -j 7 -o A.stripped/",file," -p A.stripped/",paste0(filename,"2.fastq.gz "),
                                         "A.demultiplexed/",file," A.demultiplexed/",paste0(filename,"2.fastq.gz")), stdout=T, stderr=T)
  cat(file="log.txt", temp, append=T, sep="\n")
}}



## 2.0 Count Reads
#Before we start the pipeline lets have a look for any big drops in the number of reads retained by trimming. 

#fileList <- list.files("S.demultiplexed",pattern="R1.fastq",full.names = TRUE)
#rawreadcount <- sapply(files,FastqCount)
#names(rawreadcount) <- gsub(".R1.fastq.gz","",sapply(strsplit(basename(names(rawreadcount)),"_"),'[',1))

DeMultiS  <- sapply(list.files("S.demultiplexed",pattern="R1.fastq",full.names = TRUE),FastqCount)
DeMultiA  <- sapply(list.files("A.demultiplexed",pattern="R1.fastq",full.names = TRUE),FastqCount)
StrippedS  <- sapply(list.files("S.stripped",pattern="R1.fastq",full.names = TRUE),FastqCount)
StrippedA  <- sapply(list.files("A.stripped",pattern="R1.fastq",full.names = TRUE),FastqCount)

length(DeMultiS)
length(DeMultiA)
length(StrippedS)
length(StrippedA)

ReadCounts <- data.frame("SenseDemultiplex"=DeMultiS,"SenseStripped"=StrippedS,"AntiSenseDemultiplex"=DeMultiA,"AntSenseStripped"=StrippedA)

rownames(ReadCounts) <- gsub(".R1.fastq.gz","",sapply(strsplit(basename(rownames(ReadCounts)),"_"),'[',1))


write.csv(ReadCounts,file="../Readcount.csv")


### 3.0 DADA2 

## 3.1 EUK.SENSE

path <- "S.stripped"

fnFsSense <- sort(list.files(path, pattern="EUK.*R1.fastq.gz", full.names = TRUE))
fnRsSense <- sort(list.files(path, pattern="EUK.*R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub(".S.R1.fastq.gz","",basename(fnFsSense))

#Lets look at the quality for forwards
pdf("../7.DADA2/S_EUK_forwards.pdf")
plotQualityProfile(fnFsSense[1:4])
dev.off()
#...and reverses
pdf("../7.DADA2/S_EUK_reverse.pdf")
plotQualityProfile(fnRsSense[1:4])
dev.off()
# Place filtered files in filtered/ subdirectory
filtFsSense <- file.path("../7.DADA2/filtered", paste0(sample.names, "_S_F_filt.fastq.gz"))
filtRsSense <- file.path("../7.DADA2/filtered", paste0(sample.names, "_S_R_filt.fastq.gz"))

# Filter expression

out <- filterAndTrim(fnFsSense, filtFsSense, fnRsSense, filtRsSense, minLen = 50,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     matchIDs=TRUE,compress=TRUE, multithread=7)

filtFsSense <- sort(list.files("../7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE))
filtRsSense <- sort(list.files("../7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE))



# Now lets learn the error rates

# Lets try this novel model from https://github.com/benjjneb/dada2/issues/1307
loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}



errF <- learnErrors(filtFsSense, multithread=7,
                    errorEstimationFunction = loessErrfun_mod4,
                    verbose = TRUE)
errF.2 <- learnErrors(filtFsSense, multithread=7)


errR <- learnErrors(filtRsSense, multithread=7,
                    errorEstimationFunction = loessErrfun_mod4,
                    verbose = TRUE)

#lets visualise the error rates

pdf("../7.DADA2/S_EUK_errorF.pdf")
plotErrors(errF)
dev.off()

pdf("../7.DADA2/S_EUK_OLDmodelerrorF.pdf")
plotErrors(errF.2)
dev.off()

pdf("../7.DADA2/S_EUK_errorR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Now lets dereplicate the samples, some of our files have dissapeared so we use a complex expression to restore them

derepFs <- derepFastq(list.files("../7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("../7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(list.files("../7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("../7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

# Now we can pply the DADA2 algorithm to the data
dadaFsSense <- dada(derepFs, err=errF, multithread=TRUE)

dadaRsSense <- dada(derepRs, err=errR, multithread=TRUE)

# We can now merge the seqs
mergers <- mergePairs(dadaFsSense, derepFs, dadaRsSense, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Now we can construct a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=7, verbose=TRUE)
dim(seqtab.nochim)

#Lets count up the reads for the sense section of the pipeline

#How many reads do we lose?
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFsSense, getN), sapply(dadaRsSense, getN), sapply(mergers, getN),rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- paste0("S.",c("denoisedF", "denoisedR", "merged", "lengthtrunc", "nonchim"))
rownames(track) <- names(dadaFsSense)
head(track)


#Let's save the Sense dataset and rerun the pipeline for the Antisense

EUK.Sense <- seqtab.nochim
EUK.SenseTrackedReads <- track




## 3.2 EUK.ANTISENSE

path <- "A.stripped"

fnFsASense <- sort(list.files(path, pattern="EUK.*R1.fastq.gz", full.names = TRUE))
fnRsASense <- sort(list.files(path, pattern="EUK.*R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub(".A.R1.fastq.gz","",basename(fnFsASense))

#Lets look at the quality for forwards
pdf("../7.DADA2/A_EUK_forwards.pdf")
plotQualityProfile(fnFsASense[1:4])
dev.off()
#...and reverses
pdf("../7.DADA2/A_EUK_reverse.pdf")
plotQualityProfile(fnRsASense[1:4])
dev.off()
# Place filtered files in filtered/ subdirectory
filtFsASense <- file.path("../7.DADA2/filtered", paste0(sample.names, "_A_F_filt.fastq.gz"))
filtRsASense <- file.path("../7.DADA2/filtered", paste0(sample.names, "_A_R_filt.fastq.gz"))

# Filter expression

out <- filterAndTrim(fnFsASense, filtFsASense, fnRsASense, filtRsASense, minLen = 50,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     matchIDs=TRUE,compress=TRUE, multithread=7)

#####START FROM HERE 14APRIL


filtFsSense <- sort(list.files("../7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE))
filtRsSense <- sort(list.files("../7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE))



# Now lets learn the error rates

errF <- learnErrors(filtFsASense, multithread=7,
                    errorEstimationFunction = loessErrfun_mod4,
                    verbose = TRUE)

errR <- learnErrors(filtRsASense, multithread=7,
                    errorEstimationFunction = loessErrfun_mod4,
                    verbose = TRUE)

#lets visualise the error rates

pdf("../7.DADA2/A_EUK_errorF.pdf")
plotErrors(errF)
dev.off()

pdf("../7.DADA2/A_EUK_errorR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Now lets dereplicate the samples, some of our files have dissapeared so we use a complex expression to restore them

derepFs <- derepFastq(list.files("../7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("../7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(list.files("../7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("../7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

# Now we can pply the DADA2 algorithm to the data
dadaFsASense <- dada(derepFs, err=errF, multithread=TRUE)

dadaRsASense <- dada(derepRs, err=errR, multithread=TRUE)

# We can now merge the seqs
mergers <- mergePairs(dadaFsASense, derepFs, dadaRsASense, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Now we can construct a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=7, verbose=TRUE)
dim(seqtab.nochim)

#Lets count up the reads for the sense section of the pipeline

#How many reads do we lose?
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFsSense, getN), sapply(dadaRsSense, getN), sapply(mergers, getN),rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- paste0("S.",c("denoisedF", "denoisedR", "merged", "lengthtrunc", "nonchim"))
rownames(track) <- names(dadaFsSense)
head(track)


#Let's save the Sense dataset and rerun the pipeline for the Antisense

EUK.Sense <- seqtab.nochim
EUK.SenseTrackedReads <- track






##Cheeky function from Ben, see https://github.com/benjjneb/dada2/issues/132#issuecomment-255050128

sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- rownames(tables[[1]])
  for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  sams <- unique(sample.names)
  # Make merged table
  rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
  rownames(rval) <- sams
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}

#Merge Tables




dim(Sense)
dim(seqtab.nochim)
AntiSense.RC <- seqtab.nochim 
colnames(AntiSense.RC) <- sapply(colnames(seqtab.nochim),FUN=RevComp)

TotalTable <- sumSequenceTables(Sense,AntiSense.RC)

dim(TotalTable)

## Overall Reads lost 

str(TrackedReads)

TrackedReads <- cbind(TrackedReads,track)

write.csv(TrackedReads,file="TrackedReads.csv")



#Output OTUs for taxonomy 
ASVs <- DNAStringSet(getSequences(TotalTable))
names(ASVs) <- paste0("ASV_",1:length(ASVs))
writeXStringSet(ASVs,"5.OTUs/DADA2.18S.ASVs.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


#Save output with OTU names
TotalTable.named <- TotalTable
colnames(TotalTable.named) <-paste0("ASV_",1:length(ASVs))
write.table(as.data.frame(t(TotalTable.named)),"6.mappings/OTUtabs/DADA2.18S.csv",sep=",")


##Lets also run LULU while we are here
ApplyLulu(seqs="5.OTUs/DADA2.18S.ASVs.fasta",
          table="6.mappings/OTUtabs/DADA2.18S.csv",
          output="8.LULU/DADA2.18S.lulu.csv",
          vsearchdest="/dockeDNAsoftware/vsearch-2.15.1-linux-x86_64")


```

## 6.0 Taxonomic Assignment of ASVs
#Now we need to assign some taxonomy to our big pile of sequences. We're going to do this by splitting our pile of ASVs into a number of smaller files and then sending them all to different nodes for blast assignment, followed by LCa inference using a hacky script. This should give us something to filter and refine further in our analyses.  



#Here I assume you have blast (v2.12.0) & the nt database downloaded, mine is in the directory ~/data/temp_db/nt and I will save the output to ~/data/temp_db/nt/WP5.1 where the input fasta file is.  

#how to update the blast database with the perl script, these are already run so no need to do it again

# the script can be found in the below tar folder
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz 
#tar -xvf ncbi-blast-2.12.0+-x64-linux.tar.gz   

update_blastdb.pl --decompress nt --blastdb_version 5
update_blastdb.pl taxdb
tar xzvf taxdb.tar.gz

#We first make a results folder
cd ~/data/temp_db/nt
mkdir WP5.1/results 

#split fasta into chunks of 100 seqs for fast processing
cd  WP5.1
awk -v size=100 -v pre=split -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' DADA2.18S.ASVs.fasta

##This is the raw expression we will use, the format listed is expected by other scripts/functions, dont mess with it.
#blastn -query WP5.1/split.00001 -db nt -out /WP5.1/results/test.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200


#Rather than submitting each expression we create a bunch of scripts that can then be submitted individually, there is probably a neater way to do this but I like this way as you can add many lines to be evaluated sequentially which I find hard to do with xsbatch

#make little scripts from within nt folder and give them execute permission
cd ~/data/temp_db/nt
for i in WP5.1/split.* ;do bname=$(basename $i); echo -e "cd data/temp_db/nt\n./blastn -query WP5.1/$bname -num_threads 1 -db nt -out WP5.1/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > script.$bname.sh ;done
chmod +x script*

#submit scripts with 160000mb of RAM & 60 min limits
cd ~
for i in ~/data/temp_db/nt/script.* ; do echo "~/data/temp_db/nt/$(basename $i)" | xsbatch -c 1 --time 60 --mem-per-cpu 160000 --;done 

#combine output from all the split files to form one mega taxonomy file
cd ~/data/temp_db/nt/WP5.1/results
cat split* > dada2.raw.taxonomy.txt



```

Now we download the file *dada2.raw.taxonomy.txt* to the local area and run parse taxonomy in R as follows. 

```{r}
rawtax <- read.table("data/taxonomy/dada2.raw.taxonomy.txt",sep="\t")

```







## Links and Extras

Some inspiration for TagSteady and DADA2 from 2bias and Ben [here](https://github.com/benjjneb/dada2/issues/132#issuecomment-255050128)
