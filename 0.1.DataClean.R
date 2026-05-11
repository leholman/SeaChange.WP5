### WP5 Antartica data prep ##

# some libaries 
library("seqinr")
library(dada2)
library(Biostrings)

set.seed(1234)

# bring in metadata
dates <-read.csv("metadata/dates.csv") 
metadata <- read.csv("metadata/RepMetaData.csv")

## some cleaning functions
minimumReads <- function(data, minreads) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(data) && !is.matrix(data)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(data, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(data))) {stop("Data contains NA values, which are not allowed.")}
  
  # Store original row names
  original_row_names <- row.names(data)
  
  # Apply thresholding operation based on the structure of the data
  if (is.matrix(data)) {
    # For matrices, use vectorized operation
    data[data <= minreads] <- 0
  } else if (is.data.frame(data)) {
    # For data frames, apply operation column-wise
    data <- as.data.frame(lapply(data, function(x) ifelse(x <= minreads, 0, x)))
  }
  
  
  # Get rid of empty rows
  # For both data frames and matrices, rowSums works and subsetting by rowSums > 0 removes empty rows
  outdata <- data[rowSums(data > 0) > 0, ]
  row.names(outdata) <- original_row_names[rowSums(data > 0) > 0]
  
  
  # Return the cleaned data
  return(outdata)
}
minimumReps <- function(data,min_obs) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(data) && !is.matrix(data)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(data, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(data))) {stop("Data contains NA values, which are not allowed.")}
  
  # Store original row names
  original_row_names <- row.names(data)
  
  # Determine the input type to return the same type
  input_type <- ifelse(is.data.frame(data), "data.frame", "matrix")
  
  # Calculate the number of non-zero observations per row
  total_obs <- apply(data, 1, function(row) sum(row > 0, na.rm = TRUE))
  
  # Filter rows based on the user-specified minimum number of non-zero observations
  filtered_data <- data[total_obs >= min_obs, ]
  
  # Reapply the original row names to the filtered data
  row.names(filtered_data) <- original_row_names[total_obs >= min_obs]
  
  # Ensure the output is of the same type as the input
  if (input_type == "matrix") {
    return(as.matrix(filtered_data))
  } else {
    return(as.data.frame(filtered_data))
  }
}


NS1.clean.1 <- dataCleanBy(NS1.exp,NS1.neg,method = "max",removeRows = FALSE)

inputdata <- NS1.exp
cleaningdata <- NS1.neg
method <- "max"

dataCleanBy <- function(inputdata, cleaningdata, method, removeRows=TRUE) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(inputdata) && !is.matrix(inputdata)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(inputdata, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(inputdata))) {stop("Data contains NA values, which are not allowed.")}
  # Validate method parameter
  if (!method %in% c("max", "min", "avr")) {
    stop("Method must be 'max', 'min', or 'avr'.")
  }
  
  # Preserve original format of inputdata
  inputIsDataFrame <- is.data.frame(inputdata)
  
  # Convert dataframes to matrices for uniform processing
  inputdata <- as.matrix(inputdata)
  cleaningdata <- as.matrix(cleaningdata)
  
  # Ensure both inputs contain only numeric values
  if (!all(is.numeric(inputdata)) || !all(is.numeric(cleaningdata))) {
    stop("Both inputs must contain only numeric values.")
  }
  
  # Ensure both inputs have the same number of rows
  if (nrow(inputdata) != nrow(cleaningdata)) {
    stop("Both inputs must have the same number of rows.")
  }
  
  # Calculate a metric for each row in cleaningdata based on the specified method
  metric_values <- switch(method,
                          "max" = apply(cleaningdata, 1, max),
                          "min" = apply(cleaningdata, 1, min),
                          "avr" = apply(cleaningdata, 1, mean))
  
  # Vectorized update of inputdata based on the calculated metrics
  metric_matrix <- matrix(rep(metric_values, each = ncol(inputdata)), 
                          nrow = nrow(inputdata), 
                          byrow = TRUE)
  
  # Identify positions to be zeroed that are currently >0
  to_zero <- (inputdata > 0) & (inputdata <= metric_matrix)
  num_removed <- sum(to_zero)
  num_total <- sum(inputdata > 0)
  
  #change results
  inputdata[inputdata <= metric_matrix] <- 0
  
  # Report
  message(paste0("Removed ", round(100 * num_removed / num_total, 3), "% of positive ASV/OTU observations < criteria"))
  message(paste0(num_removed, " out of ", num_total, " positive observations set to zero"))
  
  # Remove rows in inputdata with only zero values
  inputdata2 <- inputdata[rowSums(inputdata != 0) > 0, ]
  
  # Return data with empy rows if requested
  if (removeRows==FALSE){
  return(inputdata)  
  } else {
  
  # Return inputdata in its original format
  if (inputIsDataFrame) {
    return(as.data.frame(inputdata2))
  } else {
    return(inputdata2)
  }}
}

### Global settings
min_reads <- 2
min_reps <- 2
comp_method <- "avr"

## Lets start with EUK 
EUK <- read.csv("rawdata/EUK.combined.csv.gz",row.names = 1)
EUK.tax <- read.csv("taxonomy/EUK.parsedtax.csv")
EUK.pr2 <- read.csv("taxonomy/EUK.pr2.csv")

## we have three sequencing lanes to clean using the negatives, combine into 2 datasets - the miseq comparison and the ecological analyses
##then apply minimum reads and reps

## first the MiSeq
neg.m1 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="MiSeq" & metadata$SampleType=="ControlN"]  
exp.m1 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="MiSeq" & metadata$SampleType=="Experimental"]  

## some blanks use hyphen not dot so a weird to lower and gsub
MS1.neg <- EUK[,na.omit(match(tolower(gsub("-", ".", neg.m1)), tolower(colnames(EUK))))]
MS1.exp <- EUK[,na.omit(match(tolower(gsub("-", ".", exp.m1)), tolower(colnames(EUK))))]
  
MS1.clean.1 <- dataCleanBy(MS1.exp,MS1.neg,method = "max",removeRows = FALSE)

#write out the negative results
MS1.neg.2 <- MS1.neg[rowSums(MS1.neg)>0,]
output <- cbind(MS1.neg.2 ,EUK.pr2[match(rownames(MS1.neg.2),EUK.pr2$X.1),],EUK.tax[match(rownames(MS1.neg.2),EUK.tax$OTU),])
write.csv(output,"cleandata/negatives/MS1.negatives.csv")

## now the NovaSeq 1
neg.n1 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="Nova1" & metadata$SampleType=="ControlN"]  
exp.n1 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="Nova1" & metadata$SampleType=="Experimental"]  

NS1.neg <- EUK[,na.omit(match(tolower(gsub("-", ".", neg.n1)), tolower(colnames(EUK))))]
NS1.exp <- EUK[,na.omit(match(tolower(gsub("-", ".", exp.n1)), tolower(colnames(EUK))))]

NS1.clean.1 <- dataCleanBy(NS1.exp,NS1.neg,method = "max",removeRows = FALSE)

#write out the negative results
NS1.neg.2 <- NS1.neg[rowSums(NS1.neg)>0,]
output <- cbind(NS1.neg.2 ,EUK.pr2[match(rownames(NS1.neg.2),EUK.pr2$X.1),],EUK.tax[match(rownames(NS1.neg.2),EUK.tax$OTU),])
write.csv(output,"cleandata/negatives/NS1.negatives.csv")


## now the NovaSeq 2
neg.n2 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="Nova2" & metadata$SampleType=="ControlN"]  
exp.n2 <- metadata$UniqueID[metadata$PrimerID=="EUK" & metadata$SeqExp=="Nova2" & metadata$SampleType=="Experimental"]  

NS2.neg <- EUK[,na.omit(match(tolower(gsub("-", ".", neg.n2)), tolower(colnames(EUK))))]
NS2.exp <- EUK[,na.omit(match(tolower(gsub("-", ".", exp.n2)), tolower(colnames(EUK))))]

NS2.clean.2 <- dataCleanBy(NS2.exp,NS2.neg,method = "max",removeRows = FALSE)

#write out the negative results
NS2.neg.2 <- NS2.neg[rowSums(NS2.neg)>0,]
output <- cbind(NS2.neg.2 ,EUK.pr2[match(rownames(NS2.neg.2),EUK.pr2$X.1),],EUK.tax[match(rownames(NS2.neg.2),EUK.tax$OTU),])
write.csv(output,"cleandata/negatives/NS2.negatives.csv")


### great now we create the combined dataset 
WP5.comp <- cbind(MS1.clean.1,NS1.clean.1)
WP5.comp.2 <- minimumReads(WP5.comp,min_reads)
WP5.comp.3 <- minimumReps(WP5.comp.2,min_reps)
WP5.comp.4 <- WP5.comp.3[nchar(EUK.pr2$X[match(rownames(WP5.comp.3),EUK.pr2$X.1)]) > 90 & nchar(EUK.pr2$X[match(rownames(WP5.comp.3),EUK.pr2$X.1)]) < 200,]
WP5.comp.tax <- cbind(EUK.pr2[match(rownames(WP5.comp.4),EUK.pr2$X.1),],EUK.tax[match(rownames(WP5.comp.4),EUK.tax$OTU),])

write.csv(WP5.comp.4,"cleandata/EUK.WP5.comp.clean.csv")
write.csv(WP5.comp.tax,"cleandata/EUK.WP5.comp.clean.tax.csv")

## now the main dataset 

WP5.main <- cbind(NS1.clean.1,NS2.clean.2)
WP5.main.2 <- minimumReads(WP5.main,min_reads)
WP5.main.3 <- minimumReps(WP5.main.2,min_reps)
WP5.main.4 <- WP5.main.3[nchar(EUK.pr2$X[match(rownames(WP5.main.3),EUK.pr2$X.1)]) > 90 & nchar(EUK.pr2$X[match(rownames(WP5.main.3),EUK.pr2$X.1)]) < 200,]
WP5.main.tax <- cbind(EUK.pr2[match(rownames(WP5.main.4),EUK.pr2$X.1),],EUK.tax[match(rownames(WP5.main.4),EUK.tax$OTU),])

write.csv(WP5.main.4,"cleandata/EUK.WP5.main.clean.csv")
write.csv(WP5.main.tax,"cleandata/EUK.WP5.main.clean.tax.csv")



### Now the RIZ marker

RIZ <- read.csv("rawdata/RIZ.combined.csv.gz", row.names = 1)
RIZ.tax <- read.csv("taxonomy/RIZ.parsedtax.csv")
RIZ.asvs <- read.fasta("rawdata/RIZ.asvs.fasta",as.string = TRUE)

# NovaSeq 1
neg.n1 <- metadata$UniqueID[metadata$PrimerID=="RIZ" & metadata$SeqExp=="Nova1" & metadata$SampleType=="ControlN"]
exp.n1 <- metadata$UniqueID[metadata$PrimerID=="RIZ" & metadata$SeqExp=="Nova1" & metadata$SampleType=="Experimental"]

NS1.neg <- RIZ[, na.omit(match(tolower(gsub("-", ".", neg.n1)), tolower(colnames(RIZ))))]
NS1.exp <- RIZ[, na.omit(match(tolower(gsub("-", ".", exp.n1)), tolower(colnames(RIZ))))]

NS1.clean.1 <- dataCleanBy(NS1.exp, NS1.neg, method = "max", removeRows = FALSE)

# Write negatives
NS1.neg.2 <- NS1.neg[rowSums(NS1.neg) > 0, ]
output <- cbind(NS1.neg.2, RIZ.tax[match(rownames(NS1.neg.2), RIZ.tax$ASV), ])
write.csv(output, "cleandata/negatives/NS1.RIZ.negatives.csv")

# NovaSeq 2
neg.n2 <- metadata$UniqueID[metadata$PrimerID=="RIZ" & metadata$SeqExp=="Nova2" & metadata$SampleType=="ControlN"]
exp.n2 <- metadata$UniqueID[metadata$PrimerID=="RIZ" & metadata$SeqExp=="Nova2" & metadata$SampleType=="Experimental"]

NS2.neg <- RIZ[, na.omit(match(tolower(gsub("-", ".", neg.n2)), tolower(colnames(RIZ))))]
NS2.exp <- RIZ[, na.omit(match(tolower(gsub("-", ".", exp.n2)), tolower(colnames(RIZ))))]

NS2.clean.2 <- dataCleanBy(NS2.exp, NS2.neg, method = "max", removeRows = FALSE)

# Write negatives
NS2.neg.2 <- NS2.neg[rowSums(NS2.neg) > 0, ]
output <- cbind(NS2.neg.2, RIZ.tax[match(rownames(NS2.neg.2), RIZ.tax$ASV), ])
write.csv(output, "cleandata/negatives/NS2.RIZ.negatives.csv")

# Main dataset
WP5.RIZ.main <- cbind(NS1.clean.1, NS2.clean.2)
WP5.RIZ.main.2 <- minimumReads(WP5.RIZ.main, min_reads)
WP5.RIZ.main.3 <- minimumReps(WP5.RIZ.main.2, min_reps)
WP5.RIZ.main.4 <- WP5.RIZ.main.3[nchar(RIZ.tax$ASVseq[match(rownames(WP5.RIZ.main.3),RIZ.tax$ASV)]) > 89 & nchar(RIZ.tax$ASVseq[match(rownames(WP5.RIZ.main.3),RIZ.tax$ASV)]) < 110,]
WP5.RIZ.main.tax <- RIZ.tax[match(rownames(WP5.RIZ.main.4), RIZ.tax$ASV),]

write.csv(WP5.RIZ.main.4, "cleandata/RIZ.WP5.main.clean.csv")
write.csv(WP5.RIZ.main.tax, "cleandata/RIZ.WP5.main.clean.tax.csv")


### Now the MAM marker

MAM <- read.csv("rawdata/MAM.combined.csv.gz", row.names = 1)
MAM.tax <- read.csv("taxonomy/MAM.parsedtax.csv")

# NovaSeq 1
neg.n1 <- metadata$UniqueID[metadata$PrimerID=="MAM" & metadata$SeqExp=="Nova1" & metadata$SampleType=="ControlN"]
exp.n1 <- metadata$UniqueID[metadata$PrimerID=="MAM" & metadata$SeqExp=="Nova1" & metadata$SampleType=="Experimental"]

NS1.neg <- MAM[, na.omit(match(tolower(gsub("-", ".", neg.n1)), tolower(colnames(MAM))))]
NS1.exp <- MAM[, na.omit(match(tolower(gsub("-", ".", exp.n1)), tolower(colnames(MAM))))]

NS1.clean.1 <- dataCleanBy(NS1.exp, NS1.neg, method = "max", removeRows = FALSE)

# Write negatives
NS1.neg.2 <- NS1.neg[rowSums(NS1.neg) > 0, ]
output <- cbind(NS1.neg.2, MAM.tax[match(rownames(NS1.neg.2), MAM.tax$ASV), ])
write.csv(output, "cleandata/negatives/NS1.MAM.negatives.csv")

# NovaSeq 2
neg.n2 <- metadata$UniqueID[metadata$PrimerID=="MAM" & metadata$SeqExp=="Nova2" & metadata$SampleType=="ControlN"]
exp.n2 <- metadata$UniqueID[metadata$PrimerID=="MAM" & metadata$SeqExp=="Nova2" & metadata$SampleType=="Experimental"]

NS2.neg <- MAM[, na.omit(match(tolower(gsub("-", ".", neg.n2)), tolower(colnames(MAM))))]
NS2.exp <- MAM[, na.omit(match(tolower(gsub("-", ".", exp.n2)), tolower(colnames(MAM))))]

NS2.clean.2 <- dataCleanBy(NS2.exp, NS2.neg, method = "max", removeRows = FALSE)

# Write negatives
NS2.neg.2 <- NS2.neg[rowSums(NS2.neg) > 0, ]
output <- cbind(NS2.neg.2, MAM.tax[match(rownames(NS2.neg.2), MAM.tax$ASV), ])
write.csv(output, "cleandata/negatives/NS2.MAM.negatives.csv")

# Main dataset
WP5.MAM.main <- cbind(NS1.clean.1, NS2.clean.2)
WP5.MAM.main.2 <- minimumReads(WP5.MAM.main, min_reads)
WP5.MAM.main.3 <- minimumReps(WP5.MAM.main.2, min_reps)
WP5.MAM.main.4 <- WP5.MAM.main.3[nchar(MAM.tax$ASVseq[match(rownames(WP5.MAM.main.3),MAM.tax$ASV)]) > 89 & nchar(MAM.tax$ASVseq[match(rownames(WP5.MAM.main.3),MAM.tax$ASV)]) < 120,]
WP5.MAM.main.tax <- MAM.tax[match(rownames(WP5.MAM.main.4), MAM.tax$ASV), ]

write.csv(WP5.MAM.main.4, "cleandata/MAM.WP5.main.clean.csv")
write.csv(WP5.MAM.main.tax, "cleandata/MAM.WP5.main.clean.tax.csv")


##### Lets do some read counting 
#let's count up the reads
readCounts <- read.csv("rawdata/readcounts/ReadCounts.csv")

cat_map <- c(EXP = "Experimental", NTC = "ControlN", PTC = "ControlP")

out <- do.call(rbind, lapply(unique(readCounts$Marker), function(m) {
  loopDat <- readCounts[readCounts$Marker == m, ]
  do.call(rbind, lapply(names(cat_map), function(cat) {
    sub <- loopDat[loopDat$Cat == cat, ]
    totals <- sub$SenseDemultiplex + sub$AntiSenseDemultiplex
    data.frame(stage = "demultiplexed", marker = m, sample_type = cat_map[[cat]],
               n = length(totals), mean = mean(totals), sd = sd(totals))
  }))
}))

write.table(out, "summaryfiles/read_counts_summary.csv",
            sep = ",", row.names = FALSE, col.names = FALSE,
            append = TRUE)


### count up the post bioinformatvics pre control Qc'd data  data 

markers <- list(EUK = EUK, RIZ = RIZ, MAM = MAM)
out <- do.call(rbind, lapply(names(markers), function(m) {
  d <- markers[[m]]
  do.call(rbind, lapply(c("Experimental", "ControlN", "ControlP"), function(s) {
    ids <- metadata$UniqueID[metadata$PrimerID == m & metadata$SampleType == s]
    cols <- na.omit(match(tolower(gsub("-", ".", ids)), tolower(colnames(d))))
    sums <- colSums(d[, cols])
    data.frame(stage = "filtered", marker = m, sample_type = s,
               n = length(sums), mean = mean(sums), sd = sd(sums))
  }))
}))
write.csv(out, "summaryfiles//read_counts_summary.csv", row.names = FALSE)

### Now we count up the infal QC data and provide some taxonomic summaries


