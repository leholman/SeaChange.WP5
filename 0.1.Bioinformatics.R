
## 1.0 Demultiplex Pools

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


# set the directory 
setwd("1.rawreads/")
# load in the metadata
files <- read.csv("../metadata.csv")

#record th eversion of cutadapt for the analysis
temp <- paste("Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
cat(file="../log.txt", temp, append=T, sep="\n")

#primer information
Euk18S <- c("GTACACACCGCCCGTC","TGATCCTTCTGCAGGTTCACCTAC","18S","E")
Riaz12S <-c("TTAGATACCCCACTATGC","TAGAACAGGCTCCTCTAG","12S","R")
Mamm16S <- c("CGGTTGGGGTGACCTCGGA","GCTGTTATCCCTAGGGTAACT","16S","M")

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

for (pool in unique(files$FileName)){

looppoolName <- strsplit(pool, "[_]")[[1]][1]

dir.create("demultiplexed")
temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -j 7 -g ^file:P.Forwards.fasta -G ^file:P.Reverses.fasta -o demultiplexed/",looppoolName,"-{name}.R1.fastq.gz -p demultiplexed/",looppoolName,"-{name}.R2.fastq.gz ",pool,"1_001.fastq.gz ",pool,"2_001.fastq.gz"), stdout=T, stderr=T)
cat(file="log.txt", temp, append=T, sep="\n")

}

### Rename all files according to metadata 


################# HERE 11 April 2023 #####


### Loop over each marker stripping RevComp primer from F & R














##### CODE BASEMENT - DRAGONS BELOW #####

#this is a shit way of doing this, it doesnt create an index of the tags so is very slow, can be improved easily but would give all possible tag combinations, see cutadapt documentation https://cutadapt.readthedocs.io/en/v3.1/guide.html#demultiplexing-paired-end-reads-with-combinatorial-dual-indexes 

for (primer in primerSets$V3){
  print(primer)
}

for (loopfile in 1:length(files[,1])){
  
for (loopfile in 46:56){
    
  
loopForward <- paste0(files[loopfile,match(paste0("TagF_",primer),colnames(files))],
                      primerSets$V1[primerSets$V3==primer])  
loopReverse <- paste0(files[loopfile,match(paste0("TagR_",primer),colnames(files))],
                      primerSets$V2[primerSets$V3==primer])   

tempOutF <- paste0(primerSets$V4[primerSets$V3==primer],".",files$SampleID[loopfile],".temp.R1.fastq.gz")
tempOutR <- paste0(primerSets$V4[primerSets$V3==primer],".",files$SampleID[loopfile],".temp.R2.fastq.gz")
OutF <- paste0(primerSets$V4[primerSets$V3==primer],".",files$SampleID[loopfile],".R1.fastq.gz")
OutR <- paste0(primerSets$V4[primerSets$V3==primer],".",files$SampleID[loopfile],".R2.fastq.gz")


temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -g ^",loopForward," -G ^",loopReverse," -j 8 -o ",tempOutF," -p ",tempOutR," ",files$FileName[loopfile],"1_001.fastq.gz ",files$FileName[loopfile],"2_001.fastq.gz"), stdout=T, stderr=T)
cat(file="log.txt", temp, append=T, sep="\n")

temp <- system2("cutadapt",args=paste0("-e 1 --pair-adapters --trimmed-only -a ",RevComp(primerSets$V2[primerSets$V3==primer])," -A ",RevComp(primerSets$V1[primerSets$V3==primer])," -j 7 -o ",OutF," -p ",OutR," ",tempOutF," ",tempOutR), stdout=T, stderr=T)
cat(file="log.txt", temp, append=T, sep="\n")

}
  
  
setwd("../")
```


## 3.0 Strip Primers

Some good ideas for this found [here](https://benjjneb.github.io/dada2/ITS_workflow.html)

This code helps to count the number of times a specific primer is spotted in a sew
```{r,eval=FALSE}
primerF <- "GTACACACCGCCCGTC"
primerR <- "TGATCCTTCTGCAGGTTCACCTAC"
path <- "1.rawreads"  
list.files(path)
fnFs <- sort(list.files(path, pattern = "R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2.fastq.gz", full.names = TRUE))

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(primerF)
REV.orients <- allOrients(primerR)
FWD.orients

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
print(fnFs[1])
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

```

This section strips the primers from the sequence and outputs the sense and antisense strands for each example 

```{r,eval=FALSE}

files <- read.csv("metadata.csv")

loopfile <- 1

files$SampleID[loopfile]


primerF <- "GTACACACCGCCCGTC"
primerR <- "TGATCCTTCTGCAGGTTCACCTAC"
  

forwards <- list.files("1.rawreads/",pattern="R1.fastq.gz")
reverses <- list.files("1.rawreads/",pattern="R2.fastq.gz")


#a little function to do revcomp 
RevComp <- function(input){
require(Biostrings)
dna <- DNAString(input) 
dna <- reverseComplement(dna)
return(toString(dna))}


#little check

if(length(forwards)!=length(reverses)){message("Files look good, crack on!")}else{stop("F & R different lengths")}

loopname <-  gsub("*.R[12].fastq.gz","",forwards)

#make some folders here to take the outputs
dir.create("7.DADA2/trimmed")
dir.create("7.DADA2/trimmed/temp")

##loop over files


for (loopfile in 1:length(files[,1])){


loopname <- files$SampleID[loopfile]

print(loopname)

#Sense
#Trim front anchored primers
  
system2("cutadapt",args=paste0("--pair-adapters --trimmed-only -g ^",primerF," -G ^",primerR," -j 7 -o 7.DADA2/trimmed/temp/",loopname,".T.S.R1.fastq.gz -p 7.DADA2/trimmed/temp/",loopname,".T.S.R2.fastq.gz 1.rawreads/",loopname,".R1.fastq.gz 1.rawreads/",loopname,".R2.fastq.gz"))

#trim primers on second end of amplicon (higher error rate in these so -e 0.2)
system2("cutadapt",args=paste0("-e 0.2 --pair-adapters --trimmed-only -a ",RevComp(primerR)," -A ",RevComp(primerF)," -j 7 -o 7.DADA2/trimmed/",loopname,".TF.S.R1.fastq.gz -p 7.DADA2/trimmed/",loopname,".TF.S.R2.fastq.gz 7.DADA2/trimmed/temp/",loopname,".T.S.R1.fastq.gz 7.DADA2/trimmed/temp/",loopname,".T.S.R2.fastq.gz"))

#AntiSense forward
system2("cutadapt",args=paste0("--pair-adapters --trimmed-only -g ^",primerR," -G ^",primerF," -j 7 -o 7.DADA2/trimmed/temp/",loopname,".T.A.R1.fastq.gz -p 7.DADA2/trimmed/temp/",loopname,".T.A.R2.fastq.gz 1.rawreads/",loopname,".R1.fastq.gz 1.rawreads/",loopname,".R2.fastq.gz"))

#trim primers on second end of antisense amplicon (higher error rate in these so -e 0.2)
system2("cutadapt",args=paste0("-e 0.2 --pair-adapters --trimmed-only -a ",RevComp(primerF)," -A ",RevComp(primerR)," -j 7 -o 7.DADA2/trimmed/",loopname,".TF.A.R1.fastq.gz -p 7.DADA2/trimmed/",loopname,".TF.A.R2.fastq.gz 7.DADA2/trimmed/temp/",loopname,".T.A.R1.fastq.gz 7.DADA2/trimmed/temp/",loopname,".T.A.R2.fastq.gz"))

}
```

## 4.0 Count Reads
Before we start the pipeline lets have a look for any big drops in the number of reads retained by trimming. 

```{r,eval=FALSE}

files <- list.files("1.rawreads",pattern="R1.fastq",full.names = TRUE)
rawreadcount <- sapply(files,FastqCount)
names(rawreadcount) <- gsub(".R1.fastq.gz","",sapply(strsplit(basename(names(rawreadcount)),"_"),'[',1))

Trim1s  <- sapply(list.files("7.DADA2/trimmed/temp",pattern=".T.S.R1.fastq",full.names = TRUE),FastqCount)
Trim2s  <- sapply(list.files("7.DADA2/trimmed",pattern=".TF.S.R1.fastq",full.names = TRUE),FastqCount)
Trim1a  <- sapply(list.files("7.DADA2/trimmed/temp",pattern=".T.A.R1.fastq",full.names = TRUE),FastqCount)
Trim2a  <- sapply(list.files("7.DADA2/trimmed",pattern=".TF.A.R1.fastq",full.names = TRUE),FastqCount)
length(Trim1s)
length(Trim2s)
length(Trim1a)
length(Trim2a)

ReadCounts <- data.frame("RawReads"=rawreadcount,"SenseTrim1"=Trim1s,"SenseTrim2"=Trim2s,"AntiSenseTrim1"=Trim1a,"AntiSenseTrim2"=Trim2a)
row.names(ReadCounts) <- names(rawreadcount)

write.csv(ReadCounts,file="Readcount.csv")

```

## 5.1 DADA2 pipeline 


```{r}
path <- "7.DADA2/trimmed"

fnFsSense <- sort(list.files(path, pattern="TF.S.R1.fastq.gz", full.names = TRUE))
fnRsSense <- sort(list.files(path, pattern="TF.S.R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub(".TF.S.R1.fastq.gz","",basename(fnFsSense))

#Lets look at the quality for forwards
pdf("7.DADA2/S_forwards.pdf")
plotQualityProfile(fnFsSense[1:4])
dev.off()
#...and reverses
pdf("7.DADA2/S_reverses.pdf")
plotQualityProfile(fnRsSense[5:8])
dev.off()
# Place filtered files in filtered/ subdirectory
filtFsSense <- file.path("7.DADA2/filtered", paste0(sample.names, "_S_F_filt.fastq.gz"))
filtRsSense <- file.path("7.DADA2/filtered", paste0(sample.names, "_S_R_filt.fastq.gz"))

# Filter expression

out <- filterAndTrim(fnFsSense, filtFsSense, fnRsSense, filtRsSense, minLen = 50,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     matchIDs=TRUE,compress=TRUE, multithread=7)

filtFsSense <- sort(list.files("7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE))
filtRsSense <- sort(list.files("7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE))



# Now lets learn the error rates

errF <- learnErrors(filtFsSense, multithread=7)

errR <- learnErrors(filtRsSense, multithread=7)

#lets visualise the error rates
pdf("7.DADA2/errorF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("7.DADA2/errorR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Now lets dereplicate the samples, some of our files have dissapeared so we use a complex expression to restore them

derepFs <- derepFastq(list.files("7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="_S_F_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="_S_R_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

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

Sense <- seqtab.nochim
TrackedReads <- track

```


## 5.2 DADA2 AntiSense


```{r, eval=FALSE}

fnFsASense <- sort(list.files(path, pattern="TF.A.R1.fastq.gz", full.names = TRUE))
fnRsASense <- sort(list.files(path, pattern="TF.A.R2.fastq.gz", full.names = TRUE))

sample.names <- gsub(".TF.A.R1.fastq.gz","",basename(fnFsASense))

#Lets look at the quality for forwards
pdf("7.DADA2/A_forwards.pdf")
plotQualityProfile(fnFsASense[1:4])
dev.off()
#...and reverses
pdf("7.DADA2/A_reverses.pdf")
plotQualityProfile(fnRsASense[5:8])
dev.off()
# Place filtered files in filtered/ subdirectory
filtFsASense <- file.path("7.DADA2/filtered", paste0(sample.names, "_A_F_filt.fastq.gz"))
filtRsASense <- file.path("7.DADA2/filtered", paste0(sample.names, "_A_R_filt.fastq.gz"))

# Filter expression


out <- filterAndTrim(fnFsASense, filtFsASense, fnRsASense, filtRsASense, minLen = 50,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     matchIDs=TRUE,compress=TRUE, multithread=7)

filtFsASense <- sort(list.files("7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE))
filtRsASense <- sort(list.files("7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE))



# Now lets learn the error rates

errF <- learnErrors(filtFsASense, multithread=7)

errR <- learnErrors(filtRsASense, multithread=7)

#lets visualise the error rates
pdf("7.DADA2/ASerrorF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("7.DADA2/ASerrorR.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Now lets dereplicate the samples, some of our files have dissapeared so we use a complex expression to restore them

derepFs <- derepFastq(list.files("7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="_A_F_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="_A_R_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

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

#How many reads do we lose?
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFsASense, getN), sapply(dadaRsASense, getN), sapply(mergers, getN),rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- paste0("A.",c("denoisedF", "denoisedR", "merged", "lengthtrunc", "nonchim"))
rownames(track) <- names(dadaFsSense)
head(track)


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
