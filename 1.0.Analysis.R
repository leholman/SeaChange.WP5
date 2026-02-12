########################################################
####=.     Analysis of WP5  Antarctic data         =####
####==== Luke E. Holman====7.05.2025====################
########################################################

####====0.0 Packages====####
library(tidyverse)
library(RColorBrewer)
library(iNEXT)
library(vegan)  
library(pbmcapply)


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

coredata$col <- c("#D55E00","#cc79a7","#0072B2","#f0e442","#009e73")

## climate data 

clim <- read.csv("cleandata/climate.data.PCA.csv",row.names = 1)
#get rid of bad sample
clim <- clim[!paste(clim$Site,clim$SampleID)=="p56 1",] 
whales <- read.csv("cleandata/whale.model.year.csv",row.names = 1)


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
  ##If a single vector 
  if(is.vector(INdataframe) && !is.list(INdataframe)){
    INdataframe <- as.data.frame(t(INdataframe))
    colnames(INdataframe) <- names(INdataframe)        
  }
  if(is.matrix(INdataframe)){
    INdataframe <- as.data.frame(INdataframe)
  }
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
}  

average_by_reps <- function (df,sampleIndex){
  # Check if the dataframe and cols have compatible sizes
  if (ncol(df) != length(sampleIndex)) {
    stop("The number of columns in the dataframe must match the length of the vector")
  }
  
  # Get unique sample names
  unique_samples <- unique(sampleIndex)
  
  # Initialize the output dataframe
  avg_df <- data.frame(matrix(ncol = length(unique_samples), nrow = nrow(df)))
  colnames(avg_df) <- unique_samples
  rownames(avg_df) <- rownames(df) 
  
  # Compute averages for each sample
  for (sample in unique_samples) {
    indices <- which(sampleIndex == sample)  # Get indices for columns corresponding to this sample
    if (length(indices) == 1) {
      # If only one column, use it directly
      avg_df[, sample] <- df[, indices, drop = FALSE]
    } else {
      # Compute row means for multiple columns
      avg_df[, sample] <- rowMeans(df[, indices, drop = FALSE], na.rm = TRUE)
    }
  }
  
  return(avg_df)
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


## Data Handling

euk.reps <- average_by_reps(euk,gsub("\\.[0-9]+$", "",colnames(euk)))
rownames(euk.reps) <- rownames(euk)

getPalette = colorRampPalette(brewer.pal(9, "Set3"))

euk.pr2.2 <- euk.pr2[match(rownames(euk.reps),euk.pr2$X.1),]

## ASV subsets

ASVs.met <- euk.pr2.2$X.1[euk.pr2.2$Subdivision=="Metazoa" & euk.pr2.2$Subdivision.1 >80]
ASVs.pro <- euk.pr2.2$X.1[
  euk.pr2.2$Subdivision %in% c(
    "Gyrista", "Cercozoa", "Apicomplexa", "Dinoflagellata", "Ciliophora", "Bigyra",
    "Chlorophyta_X", "Haptophyta_X", "Apusomonada_X", "Prasinodermophyta_X",
    "Radiolaria", "Breviatea_X", "Eurhodophytina", "Fornicata", "Cryptophyta_X",
    "Choanoflagellata", "Tubulinea_X", "Picozoa_subdivision_1", "Evosea_X",
    "Foraminifera", "Euglenozoa", "Telonemia_X", "Chrompodellids", "Perkinsea",
    "Picozoa_subdivision_2", "Discoba_X", "Discosea_X", "Hemimastigophora_X",
    "Centroplasthelida_X", "Metamonada_X", "Ancyromonadida_X", "Nebulidia_X",
    "Rotosphaerida", "Kathablepharida", "Cryptophyta_X:nucl", "Rhodophyta_X",
    "Alveolata_X"
  ) & euk.pr2.2$Subdivision.1 > 80
]
ASVs.bac <- euk.pr2.2$X.1[
  euk.pr2.2$Subdivision %in% c(
    "Proteobacteria_X", "Firmicutes_X", "Verrucomicrobia_X", "Planctomycetes_X", 
    "Acidobacteria_X", "Nitrospirae_X", "Dadabacteria_X", "Thaumarchaeota_X", 
    "Kiritimatiellaeota_X", "Nitrospinae_X", "Actinobacteria_X", "Bacteroidetes_X", 
    "Hydrogenedentes_X"
  ) & euk.pr2.2$Subdivision.1 > 80
]

##pull in this little class classifier - note check this if it becomes important 
classes <- read.csv("taxonomy/class_category.csv")

ASVs.diatom <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Diatom"] & euk.pr2.2$Class.1 > 80]
ASVs.dino <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Dinoflagellate"] & euk.pr2.2$Class.1 > 80]
ASVs.smolphy <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Small phytoplankton"] & euk.pr2.2$Class.1 > 80]


## Taxonomic overview

#EUK taxonomy all family 
EUK.tax.a <- minAbundance(CountTable(as.character(euk.pr2.2$Family),euk.reps,output = "Abundance"),minAbun=0.01)
row.names(EUK.tax.a)[1] <- "Unknown"
EUK.tax.a.2 <- as.matrix(prop.table(as.matrix(EUK.tax.a[,sort(colnames(EUK.tax.a))]),margin = 2))
barplot(EUK.tax.a.2,col=rev(getPalette(dim(EUK.tax.a.2)[1])))


pdf("figures/EUK.all.fam.pdf",width = 14,height = 9)
par(mfrow=c(1,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(EUK.tax.a.2,col=rev(getPalette(dim(EUK.tax.a.2)[1])))
legend(150,1,rev(rownames(EUK.tax.a.2)),fill=getPalette(dim(EUK.tax.a.2)[1]),cex=0.5,bty = "n",y.intersp=0.75)
dev.off()


## use this little script to seperate everything 
site.map <- c(BB3 = "BB", scAØ2 = "SC", MC3 = "MC", p68 = "p68", p56 = "p56")

expandedDates <- merge(
  data.frame(SampleName = cn <- colnames(euk.reps),
             Site       = site.map[regmatches(cn, regexpr("BB3|scAØ2|MC3|p68|p56", cn))],
             SampleID   = as.numeric(sub(".*(?:BB3|scAØ2|MC3|p68|p56)[._]", "", cn, perl = TRUE)),
             stringsAsFactors = FALSE),
  dates, by = c("Site", "SampleID"), all.x = TRUE)

##loop over sites 
unique(expandedDates$Site)
for (site in unique(expandedDates$Site)){
  loopdat <- EUK.tax.a.2[,colnames(EUK.tax.a.2)%in%(expandedDates$SampleName[expandedDates$Site==site])]
  pdf(paste0("figures/EUK.",site,".fam.pdf"),width = 7.5,height = 5)
  par(mfrow=c(1,1),mar=c(3.1, 3.1, 1.1, 1.1),xpd=TRUE)
  barplot(loopdat,col=rev(getPalette(dim(EUK.tax.a.2)[1])),las=2,names.arg = expandedDates$CalYear[match(colnames(loopdat),expandedDates$SampleName)])  
  dev.off()
}

## print out a legend

pal <- rev(getPalette(nrow(EUK.tax.a.2))) 
pdf("figures/EUK.family.legend.pdf", width = 4.5, height = 9)
par(mar = c(0, 0, 0, 0))  
plot.new()
legend("center",
       legend = rev(rownames(EUK.tax.a.2)),       
       fill   = rev(pal),                         
       bty    = "n",
       border  = NA,
       cex    = 0.6,
       y.intersp = 0.7)

dev.off()

### attempt of a shared x axis


# Common, ordered time axis across all sites
all_years <- min(unique(na.omit(expandedDates$CalYear))):max(unique(na.omit(expandedDates$CalYear)))
sites     <- unique(expandedDates$Site)

# Helper: make a site x year matrix of *proportions* (averaging if multiple samples per year)
# rename the matrix arg so it won't collide with lapply's X
make_site_year <- function(site, M, map, years){
  out <- sapply(years, function(y){
    smps <- map$SampleName[map$Site == site & map$CalYear == y]
    if (length(smps) == 0) {
      rep(0, nrow(M))
    } else if (length(smps) == 1) {
      M[, smps, drop = FALSE][,1]
    } else {
      rowMeans(M[, smps, drop = FALSE], na.rm = TRUE)
    }
  })
  rownames(out) <- rownames(M)
  colnames(out) <- years
  out
}

### now a helper to plot segments between all the bars
# Draw light "ribbons" connecting non-zero stacked bars across zero gaps

connect_observed <- function(M, bp, cols, alpha = 0.35, minval = 0){
  M[is.na(M)] <- 0
  obs <- which(colSums(M) > 0)
  if (length(obs) < 2) return(invisible(NULL))
  cs <- apply(M, 2, cumsum)
  bt <- rbind(0, cs[-nrow(M), , drop=FALSE])
  tp <- cs
  for (r in seq_len(nrow(M))){
    for (k in seq_len(length(obs) - 1)){
      i <- obs[k]; j <- obs[k + 1]
      if (M[r, i] <= minval && M[r, j] <= minval) next
      polygon(
        x = c(bp[i], bp[i], bp[j], bp[j]),
        y = c(bt[r, i], tp[r, i], tp[r, j], bt[r, j]),
        col = adjustcolor(cols[r], alpha.f = alpha),
        border = NA
      )
    }
  }
  invisible(NULL)
}


# now this works (no name collision)
site_mats <- setNames(
  lapply(sites, make_site_year, M = EUK.tax.a.2, map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 1)
pal <- rev(getPalette(nrow(EUK.tax.a.2)))  # you already defined getPalette

pdf("figures/bigplot.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.5, 4.1, 1.1, 1.1),  # small bottom margins for all but last
    oma = c(6, 0, 0, 0),          # outer bottom margin for axis labels room
    xpd = NA)

for (i in c(2,4,3,1,5)) {
  site <- sites[i]
  mat  <- site_mats[[site]]
  
  # Draw the stacked barplot without x-axis ticks
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.35)
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t",add = TRUE)
  
  # Add a site label on the left
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
 

  
  # Only the bottom panel gets the shared x-axis
  if (i == length(sites)) {
    axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
    mtext("Calibrated year", side = 1, line = 4)
  }
}

# legend
legend(0,3,
       legend = rownames(EUK.tax.a.2),
       fill   = pal,
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()




### metazoan subset 
## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only Metazoan ASVs (use AFTER relative abundance)
euk.met.only <- euk.prop[rownames(euk.prop) %in% ASVs.met, , drop = FALSE]
euk.pr2.2.met <- euk.pr2.2[match(rownames(euk.met.only),euk.pr2.2$X.1),]

EUK.met <- minAbundance(CountTable(as.character(euk.pr2.2.met$Family),euk.met.only,output = "Abundance"),minAbun=0.01)
row.names(EUK.met)[1] <- "Unknown"

barplot(as.matrix(EUK.met),col=rev(getPalette(dim(EUK.met)[1])))
legend(50,0.1,rev(rownames(EUK.met)),fill=getPalette(dim(EUK.met)[1]),cex=0.5,bty = "n",y.intersp=0.75)

## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only Metazoan ASVs (use AFTER relative abundance)
euk.met.only <- euk.prop[rownames(euk.prop) %in% ASVs.met, , drop = FALSE]

## ---- Get Family-level annotations
euk.pr2.2.met <- euk.pr2.2[match(rownames(euk.met.only), euk.pr2.2$X.1),]
met.fam <- euk.pr2.2.met$Family
met.fam[is.na(met.fam) | met.fam == ""] <- "Unknown"

## ---- Collapse to family-level abundance
EUK.met <- minAbundance(CountTable(as.character(euk.pr2.2.met$Family),euk.met.only,output = "Abundance"),minAbun=0.01)
rownames(EUK.met)[1] <- "Unknown"
EUK.met.2 <- as.matrix(EUK.met)

for (site in unique(expandedDates$Site)) {
  # Identify relevant sample names for this site
  site_samples <- expandedDates$SampleName[expandedDates$Site == site]
  
  # Match to calendar years and sort descending
  cal_years <- expandedDates$CalYear[match(site_samples, expandedDates$SampleName)]
  names(cal_years) <- site_samples
  sorted_samples <- names(sort(cal_years, decreasing = TRUE))
  
  # Subset and reorder columns in loopdat
  loopdat <- EUK.met.2[, colnames(EUK.met.2) %in% sorted_samples, drop = FALSE]
  loopdat <- loopdat[, sorted_samples, drop = FALSE]
  
  # Barplot
  pdf(paste0("figures/EUK.metazoa.", site, ".fam.pdf"), width = 7.5, height = 5)
  par(mfrow = c(1,1), mar = c(3.1, 3.1, 1.1, 1.1), xpd = TRUE)
  
  barplot(loopdat,
          col = rev(getPalette(nrow(EUK.met.2))),
          las = 2,
          names.arg = expandedDates$CalYear[match(colnames(loopdat), expandedDates$SampleName)])
  
  dev.off()
}


## ---- Legend for Metazoa only
pdf("figures/EUK.metazoa.family.legend.pdf", width = 4.5, height = 9)
par(mar = c(0, 0, 0, 0))  
plot.new()
legend("center",
       legend = rev(rownames(EUK.met.2)), 
       fill = getPalette(nrow(EUK.met.2)),
       bty = "n",
       border = NA,
       cex = 0.6,
       y.intersp = 0.7)
dev.off()

## what about a big plot 

# now this works (no name collision)
site_mats_met <- setNames(
  lapply(sites, make_site_year, M = EUK.met.2, map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 0.2)
pal <- rev(getPalette(nrow(EUK.met.2)))  # you already defined getPalette

pdf("figures/EUK.large.met.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.5, 4.1, 1.1, 1.1),  # small bottom margins for all but last
    oma = c(6, 0, 0, 0),          # outer bottom margin for axis labels room
    xpd = NA)

for (i in c(2,4,3,1,5)) {
  site <- sites[i]
  mat  <- site_mats_met[[site]]
  
  # Draw the stacked barplot without x-axis ticks
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
   #             ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.35)
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                #             ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t",add=TRUE)
  
  # Add a site label on the left
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
  
  
  
  # Only the bottom panel gets the shared x-axis
  if (i == length(sites)) {
    axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
    mtext("Calibrated year", side = 1, line = 4)
  }
}

# legend
legend(0,0.4,
       legend = rownames(EUK.met.2),
       fill   = pal,
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()





### lets make some other subsets! 


### diatoms
## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only diatoms ASVs (use AFTER relative abundance)
euk.dia.only <- euk.prop[rownames(euk.prop) %in% ASVs.diatom, , drop = FALSE]
euk.pr2.2.dia <- euk.pr2.2[match(rownames(euk.dia.only),euk.pr2.2$X.1),]

EUK.dia <- minAbundance(CountTable(as.character(euk.pr2.2.dia$Order),euk.dia.only,output = "Abundance"),minAbun=0.01)
row.names(EUK.dia)[1] <- "Unknown"

barplot(as.matrix(EUK.dia),col=rev(getPalette(dim(EUK.dia)[1])))
legend(50,0.1,rev(rownames(EUK.dia)),fill=getPalette(dim(EUK.dia)[1]),cex=0.5,bty = "n",y.intersp=0.75)

## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only diaazoan ASVs (use AFTER relative abundance)
euk.dia.only <- euk.prop[rownames(euk.prop) %in% ASVs.diatom, , drop = FALSE]

## ---- Get Family-level annotations
euk.pr2.2.dia <- euk.pr2.2[match(rownames(euk.dia.only), euk.pr2.2$X.1),]
dia.fam <- euk.pr2.2.dia$Family
dia.fam[is.na(dia.fam) | dia.fam == ""] <- "Unknown"

## ---- Collapse to family-level abundance
EUK.dia <- minAbundance(CountTable(dia.fam,euk.dia.only,output = "Abundance"),minAbun=0.01)
EUK.dia.2 <- log10(as.matrix(EUK.dia))
EUK.dia.2[EUK.dia.2 =="-Inf"] <- 0
EUK.dia.3<-as.matrix(EUK.dia)


for (site in unique(expandedDates$Site)) {
  # Identify relevant sample names for this site
  site_samples <- expandedDates$SampleName[expandedDates$Site == site]
  
  # Match to calendar years and sort descending
  cal_years <- expandedDates$CalYear[match(site_samples, expandedDates$SampleName)]
  names(cal_years) <- site_samples
  sorted_samples <- names(sort(cal_years, decreasing = TRUE))
  
  # Subset and reorder columns in loopdat
  loopdat <- EUK.dia.2[, colnames(EUK.dia.2) %in% sorted_samples, drop = FALSE]
  loopdat <- loopdat[, sorted_samples, drop = FALSE]
  
  # Barplot
  pdf(paste0("figures/EUK.dia.", site, ".fam.pdf"), width = 7.5, height = 5)
  par(mfrow = c(1,1), mar = c(3.1, 3.1, 1.1, 1.1), xpd = TRUE)
  
  barplot(loopdat,
          col = rev(getPalette(nrow(EUK.dia.2))),
          las = 2,
          names.arg = expandedDates$CalYear[match(colnames(loopdat), expandedDates$SampleName)])
  
  dev.off()
}

## ---- Legend for diatoms only
pdf("figures/EUK.diatom.family.legend.pdf", width = 4.5, height = 9)
par(mar = c(0, 0, 0, 0))  
plot.new()
legend("center",
       legend = rownames(EUK.dia.2), 
       fill = getPalette(nrow(EUK.dia.2)),
       bty = "n",
       border = NA,
       cex = 0.6,
       y.intersp = 0.7)
dev.off()

## what about a big plot 

# now this works (no name collision)
site_mats_dia <- setNames(
  lapply(sites, make_site_year, M = EUK.dia.3, map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 1)
pal <- rev(getPalette(nrow(EUK.dia.3)))  # you already defined getPalette

pdf("figures/EUK.large.dia.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.5, 4.1, 1.1, 1.1),  # small bottom margins for all but last
    oma = c(6, 0, 0, 0),          # outer bottom margin for axis labels room
    xpd = NA)

for (i in c(2,4,3,1,5)) {
  site <- sites[i]
  mat  <- site_mats_dia[[site]]
  
  # Draw the stacked barplot without x-axis ticks
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.35)
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t",add = TRUE)
  
  
  # Add a site label on the left
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
  
  
  
  # Only the bottom panel gets the shared x-axis
  if (i == length(sites)) {
    axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
    mtext("Calibrated year", side = 1, line = 4)
  }
}

# legend
legend(0,0.4,
       legend = rownames(EUK.dia.3),
       fill   = pal,
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()





## now we do the dinos 


## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only diatoms ASVs (use AFTER relative abundance)
euk.din.only <- euk.prop[rownames(euk.prop) %in% ASVs.dino, , drop = FALSE]
euk.pr2.2.din <- euk.pr2.2[match(rownames(euk.din.only),euk.pr2.2$X.1),]

EUK.din <- minAbundance(CountTable(as.character(euk.pr2.2.din$Order),euk.din.only,output = "Abundance"),minAbun=0.01)
row.names(EUK.din)[1] <- "Unknown"

barplot(as.matrix(EUK.din),col=rev(getPalette(dim(EUK.din)[1])))
legend(50,0.1,rev(rownames(EUK.din)),fill=getPalette(dim(EUK.din)[1]),cex=0.5,bty = "n",y.intersp=0.75)

## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only dinazoan ASVs (use AFTER relative abundance)
euk.din.only <- euk.prop[rownames(euk.prop) %in% ASVs.dino, , drop = FALSE]

## ---- Get Family-level annotations
euk.pr2.2.din <- euk.pr2.2[match(rownames(euk.din.only), euk.pr2.2$X.1),]
din.fam <- euk.pr2.2.din$Family
din.fam[is.na(din.fam) | din.fam == ""] <- "Unknown"

## ---- Collapse to family-level abundance
EUK.din <- minAbundance(CountTable(din.fam,euk.din.only,output = "Abundance"),minAbun=0.01)
EUK.din.2 <- as.matrix(EUK.din)


for (site in unique(expandedDates$Site)) {
  # Identify relevant sample names for this site
  site_samples <- expandedDates$SampleName[expandedDates$Site == site]
  
  # Match to calendar years and sort descending
  cal_years <- expandedDates$CalYear[match(site_samples, expandedDates$SampleName)]
  names(cal_years) <- site_samples
  sorted_samples <- names(sort(cal_years, decreasing = TRUE))
  
  # Subset and reorder columns in loopdat
  loopdat <- EUK.din.2[, colnames(EUK.din.2) %in% sorted_samples, drop = FALSE]
  loopdat <- loopdat[, sorted_samples, drop = FALSE]
  
  # Barplot
  pdf(paste0("figures/EUK.din.", site, ".fam.pdf"), width = 7.5, height = 5)
  par(mfrow = c(1,1), mar = c(3.1, 3.1, 1.1, 1.1), xpd = TRUE)
  
  barplot(loopdat,
          col = rev(getPalette(nrow(EUK.din.2))),
          las = 2,
          names.arg = expandedDates$CalYear[match(colnames(loopdat), expandedDates$SampleName)])
  
  dev.off()
}

## ---- Legend for dintoms only
pdf("figures/EUK.dino.family.legend.pdf", width = 4.5, height = 9)
par(mar = c(0, 0, 0, 0))  
plot.new()
legend("center",
       legend = rev(rownames(EUK.din.2)), 
       fill = getPalette(nrow(EUK.din.2)),
       bty = "n",
       border = NA,
       cex = 0.6,
       y.intersp = 0.7)
dev.off()


## what about a big plot 

# now this works (no name collision)
site_mats_din <- setNames(
  lapply(sites, make_site_year, M = EUK.din.2, map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 1)
pal <- rev(getPalette(nrow(EUK.din.2)))  # you already defined getPalette

pdf("figures/EUK.large.din.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.5, 4.1, 1.1, 1.1),  # small bottom margins for all but last
    oma = c(6, 0, 0, 0),          # outer bottom margin for axis labels room
    xpd = NA)

for (i in c(2,4,3,1,5)) {
  site <- sites[i]
  mat  <- site_mats_din[[site]]
  
  # Draw the stacked barplot without x-axis ticks
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
    #            ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.35)
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t",add = TRUE)
  
  # Add a site label on the left
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
  
  
  
  # Only the bottom panel gets the shared x-axis
  if (i == length(sites)) {
    axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
    mtext("Calibrated year", side = 1, line = 4)
  }
}

# legend
legend(0,0.12,
       legend = rownames(EUK.din.2),
       fill   = pal,
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()




## finally small phyto
## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only diatoms ASVs (use AFTER relative abundance)
euk.smolphy.only <- euk.prop[rownames(euk.prop) %in% ASVs.smolphy, , drop = FALSE]
euk.pr2.2.smolphy <- euk.pr2.2[match(rownames(euk.smolphy.only),euk.pr2.2$X.1),]

EUK.smolphy <- minAbundance(CountTable(as.character(euk.pr2.2.smolphy$Order),euk.smolphy.only,output = "Abundance"),minAbun=0.01)


barplot(as.matrix(EUK.smolphy),col=rev(getPalette(dim(EUK.smolphy)[1])))
legend(50,0.1,rev(rownames(EUK.smolphy)),fill=getPalette(dim(EUK.smolphy)[1]),cex=0.5,bty = "n",y.intersp=0.75)

## ---- Relative abundance first
euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## ---- Subset only smolphyazoan ASVs (use AFTER relative abundance)
euk.smolphy.only <- euk.prop[rownames(euk.prop) %in% ASVs.smolphy, ,drop = FALSE]

## ---- Get Family-level annotations
euk.pr2.2.smolphy <- euk.pr2.2[match(rownames(euk.smolphy.only), euk.pr2.2$X.1),]
smolphy.fam <- euk.pr2.2.smolphy$Family
smolphy.fam[is.na(smolphy.fam) | smolphy.fam == ""] <- "Unknown"

## ---- Collapse to family-level abundance
EUK.smolphy <- minAbundance(CountTable(smolphy.fam,euk.smolphy.only,output = "Abundance"),minAbun=0.01)
EUK.smolphy.2 <- as.matrix(EUK.smolphy)


for (site in unique(expandedDates$Site)) {
  # Identify relevant sample names for this site
  site_samples <- expandedDates$SampleName[expandedDates$Site == site]
  
  # Match to calendar years and sort descensmolphyg
  cal_years <- expandedDates$CalYear[match(site_samples, expandedDates$SampleName)]
  names(cal_years) <- site_samples
  sorted_samples <- names(sort(cal_years, decreasing = TRUE))
  
  # Subset and reorder columns in loopdat
  loopdat <- EUK.smolphy.2[, colnames(EUK.smolphy.2) %in% sorted_samples, drop = FALSE]
  loopdat <- loopdat[, sorted_samples, drop = FALSE]
  
  # Barplot
  pdf(paste0("figures/EUK.smolphy.", site, ".fam.pdf"), width = 7.5, height = 5)
  par(mfrow = c(1,1), mar = c(3.1, 3.1, 1.1, 1.1), xpd = TRUE)
  
  barplot(loopdat,
          col = rev(getPalette(nrow(EUK.smolphy.2))),
          las = 2,
          names.arg = expandedDates$CalYear[match(colnames(loopdat), expandedDates$SampleName)])
  
  dev.off()
}

## ---- Legend for smolphytoms only
pdf("figures/EUK.smolphyo.family.legend.pdf", width = 4.5, height = 9)
par(mar = c(0, 0, 0, 0))  
plot.new()
legend("center",
       legend = rev(rownames(EUK.smolphy.2)), 
       fill = getPalette(nrow(EUK.smolphy.2)),
       bty = "n",
       border = NA,
       cex = 0.6,
       y.intersp = 0.7)
dev.off()

## what about a big plot 

# now this works (no name collision)
site_mats_smol <- setNames(
  lapply(sites, make_site_year, M = EUK.smolphy.2, map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 1)
pal <- rev(getPalette(nrow(EUK.smolphy.2)))  # you already defined getPalette

pdf("figures/EUK.large.smol.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.5, 4.1, 1.1, 1.1),  # small bottom margins for all but last
    oma = c(6, 0, 0, 0),          # outer bottom margin for axis labels room
    xpd = NA)

for (i in c(2,4,3,1,5)) {
  site <- sites[i]
  mat  <- site_mats_smol[[site]]
  
  # Draw the stacked barplot without x-axis ticks
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                #            ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.35)
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "t",add = TRUE)
  
  # Add a site label on the left
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
  
  
  
  # Only the bottom panel gets the shared x-axis
  if (i == length(sites)) {
    axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
    mtext("Calibrated year", side = 1, line = 4)
  }
}

# legend
legend(0,0.04,
       legend = rownames(EUK.smolphy.2),
       fill   = pal,
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()

#### what about a ratio???

Dia.smol.ratio <- colSums(EUK.dia.3)/colSums(EUK.smolphy.2)

pdf("figures/EUK.dia.Smol.ratio.pdf",width = 9,height = 4)
plot(expandedDates$CalYear[match(names(Dia.smol.ratio),expandedDates$SampleName)],
     log10(Dia.smol.ratio),
     col=coredata$col[match(toupper(expandedDates$Site[match(names(Dia.smol.ratio),expandedDates$SampleName)]),coredata$Site)],
     pch=16,
     xlim=c(2016,1780),
     xlab="Date",ylab="log10 Diatom-SmallPhytoplankton Ratio",bty="n",axes=FALSE)
axis(1,lwd=0,lwd.ticks = 1,at = seq(2000,1750,-25),labels = seq(2000,1750,-25) )  
axis(2)
dev.off()



# change point analysis 

# packages
library(rioja)
library(vegan)


unique(expandedDates$Site)

for (dataset in c("met","dia","din","smolphy"))

  
coniss.helper <- function(site,data){
  data.sites <- expandedDates$Site[match(colnames(data),expandedDates$SampleName)]
  data.sub <- data[,grepl(site,data.sites)]
  data.sub.years <- expandedDates$CalYear[match(colnames(data.sub),expandedDates$SampleName)]
  data.sub.sort <- data.sub[,order(data.sub.years)]
  data.dist <- vegdist(t(data.sub.sort), method = "bray")
  data.coniss <- chclust(data.dist, method = "coniss")
  bs <- bstick(data.coniss)                                           # expected SS decrease under null
  return(data.coniss)
}

BB.met.con <- coniss.helper("BB",euk.met.only) # 2
cutree(BB.met.con, k = 2)
MC.met.con <- coniss.helper("MC",euk.met.only) # 5
cutree(MC.met.con, k = 5)
p56.met.con <- coniss.helper("p56",euk.met.only) # 3
cutree(p56.met.con, k = 3)
p68.met.con <- coniss.helper("p68",euk.met.only) # 3
cutree(p68.met.con, k = 3)
SC.met.con <- coniss.helper("SC",euk.met.only) # 6
cutree(SC.met.con, k = 1)


BB.dia.con <- coniss.helper("BB",euk.dia.only) # 3
cutree(BB.dia.con, k = 3)
MC.dia.con <- coniss.helper("MC",euk.dia.only) # 3
cutree(MC.dia.con, k = 3)
p56.dia.con <- coniss.helper("p56",euk.dia.only) # 2
cutree(p56.dia.con, k = 2)
p68.dia.con <- coniss.helper("p68",euk.dia.only) # 3
cutree(p68.dia.con, k = 3)
SC.dia.con <- coniss.helper("SC",euk.dia.only) # 6
cutree(SC.dia.con, k = 1)


BB.din.con <- coniss.helper("BB",euk.din.only) # 2
cutree(BB.din.con, k = 2)
MC.din.con <- coniss.helper("MC",euk.din.only) # 1
cutree(MC.din.con, k = 1)
p56.din.con <- coniss.helper("p56",euk.din.only) # 2
cutree(p56.din.con, k = 2)
p68.din.con <- coniss.helper("p68",euk.din.only) # 4
cutree(p68.din.con, k = 4)
SC.din.con <- coniss.helper("SC",euk.din.only) # 1
cutree(SC.din.con, k = 1)


BB.smolphy.con <- coniss.helper("BB",euk.smolphy.only) # 2
cutree(BB.smolphy.con, k = 2)
MC.smolphy.con <- coniss.helper("MC",euk.smolphy.only) # 1
cutree(MC.smolphy.con, k = 1)
p56.smolphy.con <- coniss.helper("p56",euk.smolphy.only) # 3
cutree(p56.smolphy.con, k = 3)
p68.smolphy.con <- coniss.helper("p68",euk.smolphy.only) # 2
cutree(p68.smolphy.con, k = 2)
SC.smolphy.con <- coniss.helper("SC",euk.smolphy.only) # 2
cutree(SC.smolphy.con, k = 2)


## what about the whole timeseries
BB.con <- coniss.helper("BB", euk.prop) # 2
cutree(BB.con, k = 2)
MC.con <- coniss.helper("MC", euk.prop) # 2
cutree(MC.con, k = 2)
p56.con <- coniss.helper("p56", euk.prop) # 4
cutree(p56.con, k = 4)
p68.con <- coniss.helper("p68", euk.prop) # 3
cutree(p68.con, k = 3)
SC.con <- coniss.helper("SC", euk.prop) # 1
cutree(SC.con, k = 1)



  
  
site <- unique(expandedDates$Site)[5]
data <- euk.met.only
data.sites <- expandedDates$Site[match(colnames(data),expandedDates$SampleName)]
data.sub <- data[,grepl(site,data.sites)]
data.sub.years <- expandedDates$CalYear[match(colnames(data.sub),expandedDates$SampleName)]
data.sub.sort <- data.sub[,order(data.sub.years)]
data.dist <- vegdist(t(data.sub.sort), method = "bray")
data.coniss <- chclust(data.dist, method = "coniss")
bs <- bstick(data.coniss)                                           # expected SS decrease under null
# Quick elbow/BS pick; adjust 'k' as needed after checking plot
k <- 3

# 6) zone labels ----
zones <- cutree(data.coniss, k = k)




### Alpha diversity 

#lets look at uneven read depth 
lib_sizes <- colSums(euk) 
min(lib_sizes)                        
quantile(lib_sizes) 
#ignore samples below 30k
n_reads <- 30000
euk_filt <- as.matrix(euk[ , lib_sizes >= n_reads]) 


#some functions to help rarefy & calc
## Hill numbers for one count vector
hill_3 <- function(v){
  v <- v[v > 0]
  p <- v / sum(v)
  c(q0 = length(v),
    q1 = exp(-sum(p * log(p))),
    q2 = 1 / sum(p^2))
}

## Rarefy the whole matrix once and return a data-frame of Hill nums
one_draw <- function(mat, depth){
  rare <- t(rrarefy(t(mat), sample = depth))
  hills <- t(apply(rare, 2, hill_3))
  data.frame(sample = rownames(hills), hills, row.names = NULL)
}
# number of rarefaction replicates
N     <- 100
set.seed(12345)
cores <- parallel::detectCores() - 1     # leave one core free

out_list <- pbmclapply(
  seq_len(N),
  function(i) one_draw(euk_filt, n_reads),
  mc.cores      = cores,
  mc.preschedule = TRUE                  # preschedule identical tasks
)


alpha_all <- bind_rows(out_list, .id = "draw")  # stack the 100 tables

alpha_avg <- alpha_all %>% 
  group_by(sample) %>% 
  summarise(across(q0:q2,
                   list(mean = mean,  sd = sd), .names = "{.col}_{.fn}"))
alpha_avg_reps <- alpha_avg %>% 
  mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
  group_by(replicate) %>% 
  summarise(across(ends_with("_mean"), mean),
            across(ends_with("_sd"),   \(x) sqrt(mean(x^2))))   # pooled SD

##Plot richness
##Hill NUmber Q=0 (richness)
pdf("figures/EUK.Q0.pdf",height =7,width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
     alpha_avg$q0_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
     pch=16,
     xlab="Year",
     ylab="ASV Richness")
points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
       alpha_avg_reps$q0_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16,cex=1.2)
dev.off()

##Hill NUmber Q=1 (inverse shannon)
pdf("figures/EUK.Q1.pdf",height =7,width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
     alpha_avg$q1_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
     pch=16,
     xlab="Year",
     ylab="Exp Shannon Index")
points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
       alpha_avg_reps$q1_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16,cex=1.2)
dev.off()

##Hill NUuber Q=2 (simpson)
pdf("figures/EUK.Q2.pdf",height =7,width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
     alpha_avg$q2_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
     pch=16,
     xlab="Year",
     ylab="Simpson Index")
points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
       alpha_avg_reps$q2_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16,cex=1.2)
dev.off()

#### Now we loop over subsets 

asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

# Rarefaction + Hill function definitions stay unchanged
hill_3 <- function(v){
  v <- v[v > 0]
  p <- v / sum(v)
  c(q0 = length(v),
    q1 = exp(-sum(p * log(p))),
    q2 = 1 / sum(p^2))
}

for (group_name in names(asv_groups)) {
  
  ASVs.keep <- asv_groups[[group_name]]
  
  one_draw <- function(mat, depth){
    rare <- t(rrarefy(t(mat), sample = depth))
    rare_sub <- rare[rownames(rare) %in% ASVs.keep, , drop = FALSE]
    hills <- t(apply(rare_sub, 2, hill_3))
    data.frame(sample = rownames(hills), hills, row.names = NULL)
  }
  
  N <- 100
  set.seed(12345)
  cores <- parallel::detectCores() - 1
  
  out_list <- pbmclapply(
    seq_len(N),
    function(i) one_draw(euk_filt, n_reads),
    mc.cores = cores,
    mc.preschedule = TRUE
  )
  
  alpha_all <- bind_rows(out_list, .id = "draw")
  
  alpha_avg <- alpha_all %>% 
    group_by(sample) %>% 
    summarise(across(q0:q2,
                     list(mean = mean,  sd = sd), .names = "{.col}_{.fn}"))
  
  alpha_avg_reps <- alpha_avg %>% 
    mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
    group_by(replicate) %>% 
    summarise(across(ends_with("_mean"), mean),
              across(ends_with("_sd"),   \(x) sqrt(mean(x^2))))   # pooled SD
  
  ## Plot Q0
  pdf(paste0("figures/", group_name, ".Q0.pdf"), height = 5, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
       alpha_avg$q0_mean,
       col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
       pch=16,
       xlab="Year",
       ylab="ASV Richness")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
         alpha_avg_reps$q0_mean,
         col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
         pch=16,cex=1.2)
  dev.off()
  
  ## Plot Q1
  pdf(paste0("figures/", group_name, ".Q1.pdf"), height = 5, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
       alpha_avg$q1_mean,
       col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
       pch=16,
       xlab="Year",
       ylab="Exp Shannon Index")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
         alpha_avg_reps$q1_mean,
         col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
         pch=16,cex=1.2)
  dev.off()
  
  ## Plot Q2
  pdf(paste0("figures/", group_name, ".Q2.pdf"), height = 5, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)],
       alpha_avg$q2_mean,
       col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_avg$sample),expandedDates$SampleName)]),toupper(coredata$Site))],0.2),
       pch=16,
       xlab="Year",
       ylab="Simpson Index")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate,expandedDates$SampleName)],
         alpha_avg_reps$q2_mean,
         col=coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate,expandedDates$SampleName)]),toupper(coredata$Site))],
         pch=16,cex=1.2)
  dev.off()
}


#### check these stats by genera

valid_genus_ASVs <- euk.pr2.2$X.1[euk.pr2.2$Genus.1 > 60 & !is.na(euk.pr2.2$Genus)]
genus_map <- setNames(euk.pr2.2$Genus, euk.pr2.2$X.1)  # Map ASV ID → Genus
# Filter to valid ASVs
euk_genus <- euk[rownames(euk) %in% valid_genus_ASVs, ]
# Set grouping vector explicitly
genus_names <- genus_map[rownames(euk_genus)]  # This may contain duplicates
# Sum ASV counts by genus
euk_genus_agg <- rowsum(euk_genus, group = genus_names)

n_reads <- 20000
colSums(euk_genus_agg)[colSums(euk_genus_agg)<n_reads]
lib_sizes_genus <- colSums(euk_genus_agg)
euk_genus_filt <- as.matrix(euk_genus_agg[, lib_sizes_genus >= n_reads])

one_draw_genus <- function(mat, depth){
  rare <- t(rrarefy(t(mat), sample = depth))
  hills <- t(apply(rare, 2, hill_3))
  data.frame(sample = rownames(hills), hills, row.names = NULL)
}

set.seed(12345)
out_list_genus <- pbmclapply(
  seq_len(N),
  function(i) one_draw_genus(euk_genus_filt, n_reads),
  mc.cores = cores,
  mc.preschedule = TRUE
)

alpha_genus_all <- bind_rows(out_list_genus, .id = "draw")

alpha_genus_avg <- alpha_genus_all %>% 
  group_by(sample) %>% 
  summarise(across(q0:q2, list(mean = mean,  sd = sd), .names = "{.col}_{.fn}"))

alpha_genus_avg_reps <- alpha_genus_avg %>% 
  mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
  group_by(replicate) %>% 
  summarise(across(ends_with("_mean"), mean),
            across(ends_with("_sd"), \(x) sqrt(mean(x^2))))

pdf("figures/EUK_genus.Q0.pdf", height = 7, width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)],
     alpha_genus_avg$q0_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
     pch=16,
     xlab="Year",
     ylab="Genus Richness")
points(expandedDates$CalYear[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)],
       alpha_genus_avg_reps$q0_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16, cex=1.2)
dev.off()
pdf("figures/EUK_genus.Q1.pdf", height = 7, width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)],
     alpha_genus_avg$q1_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
     pch=16,
     xlab="Year",
     ylab="Genus Richness")
points(expandedDates$CalYear[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)],
       alpha_genus_avg_reps$q1_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16, cex=1.2)
dev.off()

pdf("figures/EUK_genus.Q2.pdf", height = 7, width = 9)
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)],
     alpha_genus_avg$q2_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
     pch=16,
     xlab="Year",
     ylab="Genus Richness")
points(expandedDates$CalYear[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)],
       alpha_genus_avg_reps$q2_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16, cex=1.2)
dev.off()

asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

hill_3 <- function(v){
  v <- v[v > 0]
  p <- v / sum(v)
  c(q0 = length(v),
    q1 = exp(-sum(p * log(p))),
    q2 = 1 / sum(p^2))
}

# Filter genus-assigned ASVs with high bootstrap
valid_genus_ASVs <- euk.pr2.2$X.1[euk.pr2.2$Genus.1 > 60 & !is.na(euk.pr2.2$Genus)]
genus_map <- setNames(euk.pr2.2$Genus, euk.pr2.2$X.1)

for (group_name in names(asv_groups)) {
  
  ASVs.keep <- intersect(asv_groups[[group_name]], valid_genus_ASVs)
  
  one_draw <- function(mat, depth){
    rare <- t(rrarefy(t(mat), sample = depth))
    rare_sub <- rare[rownames(rare) %in% ASVs.keep, , drop = FALSE]
    
    # Map ASVs to genera
    genus_names <- genus_map[rownames(rare_sub)]
    
    # Aggregate to genus
    genus_mat <- rowsum(rare_sub, group = genus_names)
    
    # Hill numbers on genus counts
    hills <- t(apply(genus_mat, 2, hill_3))
    data.frame(sample = rownames(hills), hills, row.names = NULL)
  }
  
  N <- 100
  set.seed(12345)
  cores <- parallel::detectCores() - 1
  
  out_list <- pbmclapply(
    seq_len(N),
    function(i) one_draw(euk_filt, n_reads),
    mc.cores = cores,
    mc.preschedule = TRUE
  )
  
  alpha_all <- bind_rows(out_list, .id = "draw")
  
  alpha_avg <- alpha_all %>%
    group_by(sample) %>%
    summarise(across(q0:q2,
                     list(mean = mean,  sd = sd), .names = "{.col}_{.fn}"))
  
  alpha_avg_reps <- alpha_avg %>%
    mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>%
    group_by(replicate) %>%
    summarise(across(ends_with("_mean"), mean),
              across(ends_with("_sd"),   \(x) sqrt(mean(x^2))))
  
  ## Plot Q0
  pdf(paste0("figures/", group_name, ".genus.Q0.pdf"), height = 7, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)],
       alpha_avg$q0_mean,
       col = add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
       pch = 16,
       xlab = "Year",
       ylab = "Genus Richness")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate, expandedDates$SampleName)],
         alpha_avg_reps$q0_mean,
         col = coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
         pch = 16, cex = 1.2)
  dev.off()
  
  ## Plot Q1
  pdf(paste0("figures/", group_name, ".genus.Q1.pdf"), height = 7, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)],
       alpha_avg$q1_mean,
       col = add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
       pch = 16,
       xlab = "Year",
       ylab = "Genus Exp Shannon Index")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate, expandedDates$SampleName)],
         alpha_avg_reps$q1_mean,
         col = coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
         pch = 16, cex = 1.2)
  dev.off()
  
  ## Plot Q2
  pdf(paste0("figures/", group_name, ".genus.Q2.pdf"), height = 7, width = 9)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)],
       alpha_avg$q2_mean,
       col = add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
       pch = 16,
       xlab = "Year",
       ylab = "Genus Simpson Index")
  points(expandedDates$CalYear[match(alpha_avg_reps$replicate, expandedDates$SampleName)],
         alpha_avg_reps$q2_mean,
         col = coredata$col[match(toupper(expandedDates$Site[match(alpha_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
         pch = 16, cex = 1.2)
  dev.off()
}

####### now let's run a mixed model 

## first we make the dataset 


#we make a new var in clim we can use to match
clim$longname <- expandedDates$SampleName[match(paste(clim$Site,clim$SampleID,sep="_"),paste(expandedDates$Site,expandedDates$SampleID,sep="_"))]

alpha <- data.frame("sample" = paste(clim$Site,clim$SampleID),
                    "site" = as.factor(clim$Site),
                    "year" =clim$CalYear,
                    "ASV_Q0" = alpha_avg_reps$q0_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q1" = alpha_avg_reps$q1_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q2" = alpha_avg_reps$q2_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "genus_Q0" = alpha_genus_avg_reps$q0_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q1" = alpha_genus_avg_reps$q1_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q2" = alpha_genus_avg_reps$q2_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "whale" = whales$biomass[match(round(clim$CalYear),whales$year)],
                    "whale_s" = scale(whales$biomass[match(round(clim$CalYear),whales$year)]),
                    "PC1" = clim$PC1,
                    "PC2" = clim$PC2,
                    "PC1_s" = scale(clim$PC1),
                    "PC2_s" = scale(clim$PC2)
                    
)

library(mgcv)

# ASV models (with year)
a1_A0 <- gam(ASV_Q0 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a1_A1 <- gam(ASV_Q1 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a1_A2 <- gam(ASV_Q2 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")

# ASV models (no year)
a2_A0 <- gam(ASV_Q0 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a2_A1 <- gam(ASV_Q1 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a2_A2 <- gam(ASV_Q2 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")

# Genus models (with year)
a1_G0 <- gam(genus_Q0 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a1_G1 <- gam(genus_Q1 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a1_G2 <- gam(genus_Q2 ~ s(year, k = 10) + s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")

# Genus models (no year)
a2_G0 <- gam(genus_Q0 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a2_G1 <- gam(genus_Q1 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")
a2_G2 <- gam(genus_Q2 ~ s(site, bs = "re") + PC1_s + PC2_s + whale_s, data = alpha, method = "REML")

summary(a1_A0)
summary(a1_A1)
summary(a1_A2)

summary(a2_A0)
summary(a2_A1)
summary(a2_A2)

summary(a1_G0)
summary(a1_G1)
summary(a1_G2)

summary(a2_G0)
summary(a2_G1)
summary(a2_G2)


## beta diversity 





#### now lets run a HMSC model 

library(Hmsc)

# binary indicidence
comm.b <- (as.matrix(comm) > 0) * 1

prev <- colMeans(comm.b)
comm.b.smol <- comm.b[, prev >= 0.20, drop = FALSE]


X <- data.frame(
  PC1   = as.numeric(scale(beta$PC1_s)),
  whale = as.numeric(scale(beta$whale_s))
)

studyDesign <- data.frame(site = beta$site)
rL <- HmscRandomLevel(units = studyDesign$site)


m <- Hmsc(
  Y = comm.b.smol,
  XData = X,
  XFormula = ~ PC1 + whale,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = list(site = rL)
)


m <- sampleMcmc(
  m,
  thin = 5,
  samples = 600,
  transient = 3000,
  nChains = 2,
  nParallel = 2,
  verbose=1)


mpost <- convertToCodaObject(m)
effectiveSize(mpost)
gelman.diag(mpost)


summary(mpost$Beta)

preds <- computePredictedValues(m)

MF <- evaluateModelFit(hM=m,predy = preds)



group <- c(0, 1, 2)   # 0=Intercept, 1=Climate, 2=Whaling
VP <- computeVariancePartitioning(m)

## variance parititonwarn 
group <- c(1, 2)  # 1 = Climate, 2 = Whaling
VP <- computeVariancePartitioning(m, group = group)
colMeans(VP$vals)

# first a little test with UMAP
test <- umap(t(prop.table(as.matrix(euk.reps),2)))
plot(test$layout,pch=16, col=coredata$col[match(toupper(expandedDates$Site[match(row.names(test$layout), expandedDates$SampleName)]),toupper(coredata$Site))])
text(test$layout, labels = expandedDates$CalYear[match(row.names(test$layout), expandedDates$SampleName)],cex=0.3)

# now the NMDS
EUK.nMDS.j <- metaMDS(vegdist(t(prop.table(as.matrix(euk.reps),2)),binary = TRUE,method = "jaccard"),trymax = 200,parallel=7)
EUK.nMDS.b <- metaMDS(vegdist(t(prop.table(as.matrix(euk.reps),2)),method = "bray"),trymax = 200,parallel=7)

pdf(paste0("figures/NMDS_ALL_jaccard.pdf"), height = 7, width = 9)
plot(EUK.nMDS.j$points[,1], 
EUK.nMDS.j$points[,2],
col=coredata$col[match(toupper(expandedDates$Site[match(row.names(EUK.nMDS.j$points), expandedDates$SampleName)]),toupper(coredata$Site))],
pch=16)

text(EUK.nMDS.j$points[,1],
     EUK.nMDS.j$points[,2],
     labels = expandedDates$CalYear[match(rownames(EUK.nMDS.j$points), expandedDates$SampleName)],
     pos = 3, cex = 0.6)

dev.off()


pdf(paste0("figures/NMDS_ALL_bray.pdf"), height = 7, width = 9)
plot(EUK.nMDS.b$points[,1], 
     EUK.nMDS.b$points[,2],
     col=coredata$col[match(toupper(expandedDates$Site[match(row.names(EUK.nMDS.b$points), expandedDates$SampleName)]),toupper(coredata$Site))],
     pch=16)
text(EUK.nMDS.b$points[,1],
     EUK.nMDS.b$points[,2],
     labels = expandedDates$CalYear[match(rownames(EUK.nMDS.j$points), expandedDates$SampleName)],
     pos = 3, cex = 0.6)
dev.off()

### beta by group 

# Define ASV groups
asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

# Loop through each group
for (group_name in names(asv_groups)) {
  
  # Subset count matrix
  group_asvs <- asv_groups[[group_name]]
  mat <- euk.reps[rownames(euk.reps) %in% group_asvs, , drop = FALSE]
  
  # Transpose + convert to relative abundance
  mat_rel <- prop.table(as.matrix(mat), 2)
  
  # NMDS - Jaccard (presence/absence)
  set.seed(42)
  nmds_j <- metaMDS(vegdist(t(mat_rel), binary = TRUE, method = "jaccard"), trymax = 200, parallel = 7)
  
  # NMDS - Bray-Curtis
  set.seed(42)
  nmds_b <- metaMDS(vegdist(t(mat_rel), method = "bray"), trymax = 200, parallel = 7)
  
  # Get matching metadata
  sample_names <- rownames(nmds_j$points)
  match_idx <- match(sample_names, expandedDates$SampleName)
  years <- expandedDates$CalYear[match_idx]
  site_col <- coredata$col[match(toupper(expandedDates$Site[match_idx]), toupper(coredata$Site))]
  
  # Plot Jaccard
  pdf(paste0("figures/NMDS_", group_name, "_jaccard.pdf"), height = 7, width = 9)
  plot(nmds_j$points[,1], nmds_j$points[,2],
       col = site_col,
       pch = 16,
       main = paste(group_name, "– Jaccard"))
  text(nmds_j$points[,1], nmds_j$points[,2],
       labels = years, pos = 3, cex = 0.6)
  dev.off()
  
  # Plot Bray
  pdf(paste0("figures/NMDS_", group_name, "_bray.pdf"), height = 7, width = 9)
  plot(nmds_b$points[,1], nmds_b$points[,2],
       col = site_col,
       pch = 16,
       main = paste(group_name, "– Bray-Curtis"))
  text(nmds_b$points[,1], nmds_b$points[,2],
       labels = years, pos = 3, cex = 0.6)
  dev.off()
}

#### now let's run some tests 


beta <- data.frame("sample" = paste(clim$Site,clim$SampleID),
                   "site" = as.factor(clim$Site),
                   "year" =clim$CalYear,
                   "whale" = whales$biomass[match(round(clim$CalYear),whales$year)],
                   "whale_s" = scale(whales$biomass[match(round(clim$CalYear),whales$year)]),
                   "PC1" = clim$PC1,
                   "PC2" = clim$PC2,
                   "PC1_s" = scale(clim$PC1),
                   "PC2_s" = scale(clim$PC2))



## distance based rDA

library(vegan)

comm.raw <- prop.table(as.matrix(euk.reps),2)

comm <- t(comm.raw[,match(clim$longname,colnames(comm.raw))])


## distance basd redundnacy analysis 
mod_bray <- dbrda(
  comm ~ whale_s + PC1_s + Condition(site),
  data = beta,
  distance = "bray"
)

mod_jacc <- dbrda(
  comm ~ whale_s +PC1_s  + Condition(site),
  data = beta,
  distance = "jacc"
)

##  anova with perms within cores
perm <- how(blocks = alpha$site,nperm = 1000)

anova(mod_bray, permutations = perm)
anova(mod_bray, by = "terms", permutations = perm)

anova(mod_jacc, permutations = perm)
anova(mod_jacc, by = "terms", permutations = perm)

## variance partititoning 

var_bray <- varpart(vegdist(comm, "bray"), ~ PC1_s, ~ whale_s, data = beta)
var_jacc <- varpart(vegdist(comm, "jacc",binary = TRUE), ~ PC1_s, ~ whale_s, data = beta)

var_bray
plot(var_bray)
var_jacc
plot(var_jacc)


## PERMANOVA
perm_bray <- adonis2(vegdist(comm,method="bray")~site+PC1_s+whale_s,
                     data=beta, permutations = 10000,
                     by="margin")
perm_bray


perm_jacc <- adonis2(vegdist(comm,method="jacc",binary = TRUE)~site+PC1_s+whale_s,
                     data=beta, permutations = 10000,
                     by="margin")
perm_jacc



## write results to file
sink("results.txt")

cat("=== Distance-based RDA and PERMANOVA results by taxonomic group ===\n\n")

library(vegan)

# define groups (must exist in workspace)
asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

for (group_name in names(asv_groups)) {
  
  cat("############################################################\n")
  cat("Group:", group_name, "\n")
  cat("############################################################\n\n")
  
  # subset ASVs
  group_asvs <- asv_groups[[group_name]]
  mat <- euk.reps[rownames(euk.reps) %in% group_asvs, , drop = FALSE]
  
  # relative abundance + match ordering
  comm.raw <- prop.table(as.matrix(mat), 2)
  comm <- t(comm.raw[, match(clim$longname, colnames(comm.raw))])
  
  ## dbRDA
  cat("---- dbRDA: Bray ----\n")
  mod_bray <- dbrda(
    comm ~ whale_s + PC1_s + Condition(site),
    data = beta,
    distance = "bray"
  )
  
  perm <- how(blocks = alpha$site, nperm = 1000)
  
  print(anova(mod_bray, permutations = perm))
  print(anova(mod_bray, by = "terms", permutations = perm))
  
  cat("\n---- dbRDA: Jaccard ----\n")
  mod_jacc <- dbrda(
    comm ~ whale_s + PC1_s + Condition(site),
    data = beta,
    distance = "jacc"
  )
  
  print(anova(mod_jacc, permutations = perm))
  print(anova(mod_jacc, by = "terms", permutations = perm))
  
  ## variance partitioning
  cat("\n---- Variance partitioning: Bray ----\n")
  var_bray <- varpart(
    vegdist(comm, "bray"),
    ~ PC1_s,
    ~ whale_s,
    data = beta
  )
  print(var_bray)
  
  cat("\n---- Variance partitioning: Jaccard ----\n")
  var_jacc <- varpart(
    vegdist(comm, "jacc", binary = TRUE),
    ~ PC1_s,
    ~ whale_s,
    data = beta
  )
  print(var_jacc)
  
  ## PERMANOVA
  cat("\n---- PERMANOVA: Bray ----\n")
  perm_bray <- adonis2(
    vegdist(comm, method = "bray") ~ site + PC1_s + whale_s,
    data = beta,
    permutations = 10000,
    by = "margin"
  )
  print(perm_bray)
  
  cat("\n---- PERMANOVA: Jaccard ----\n")
  perm_jacc <- adonis2(
    vegdist(comm, method = "jacc", binary = TRUE) ~ site + PC1_s + whale_s,
    data = beta,
    permutations = 10000,
    by = "margin"
  )
  print(perm_jacc)
  
  cat("\n\n")
}

sink()











### Now do individual SPP

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

## Lets make an nREPs dataset 

euk.salkril.nReps <- NrepsMaker(euk.salkril.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum)))
euk.salkril.avr <- average_by_reps(euk.salkril.sum,gsub("(.*).[0-9]$","\\1",colnames(euk.salkril.sum)))



nr_long <- tibble(spp = rownames(euk.salkril.nReps), euk.salkril.nReps) |>
  pivot_longer(-spp, names_to = "column", values_to = "value")

# long averages
avr_long <- tibble(spp = rownames(euk.salkril.avr), euk.salkril.avr) |>
  pivot_longer(-spp, names_to = "column", values_to = "avg")

# add avg as a new column
euk.salkril.nReps.long <- nr_long |>
  left_join(avr_long, by = c("spp","column"))


euk.salkril.nReps.long$site <- expandedDates$Site[match(euk.salkril.nReps.long$column,expandedDates$SampleName)]
euk.salkril.nReps.long$dates <- expandedDates$CalYear[match(euk.salkril.nReps.long$column,expandedDates$SampleName)]


#Plots
#change the values a bit so they match the dates

pdf("figures/EUK.salkril.v1.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))

for (species in unique(euk.salkril.nReps.long$spp)){
  
  loopdat <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))

  
}

dev.off()
pdf("figures/EUK.salkril.v1.trends.pdf", width = 9, height = 3.5)
par(mfrow = c(1, 3))

for (species in unique(euk.salkril.nReps.long$spp)) {
  
  loopdat <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)],0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
       ylim = c(0, 8),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    # sort by date to ensure smooth lines
    subdat <- subdat[order(subdat$dates), ]
    # only fit if enough points
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit), 
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()

####lets plot the salps and krill 


## now lets combine them into one

sp_list <- unique(euk.salkril.nReps.long$spp)

for (sp in sp_list) {
  pdf(paste0("figures/EUK.big.",sp,".pdf"), width=9, height=6)
  par(mfrow=c(2,1), mar=c(3.5,4,2.5,1), oma=c(0,0,0,0), xaxs="i")
  
  # --- data per species
  d1 <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp==sp, ]
  d1 <- d1[!is.na(d1$dates), ]
  
  # site-color map
  sites <- sort(unique(toupper(d1$site)))
  cols  <- coredata$col[match(sites, coredata$Site)]
  colv1 <- cols[match(toupper(d1$site), sites)]
  
  # --- TOP: nREps
  set.seed(1234)
  plot(d1$dates, jitter(d1$value, 0.8),
       xlim=c(2024,1760),
       xlab="", ylab="nPCR reps", main=sp,
       pch=ifelse(d1$value>0,16,4), col=colv1, bty="n")
  for (si in sites) {
    d1s <- d1[toupper(d1$site)==si, ]
    if (nrow(d1s) > 2) {
      fit <- loess(value ~ dates, data=d1s, span=0.6, degree=1)
      xx  <- seq(min(d1s$dates), max(d1s$dates), length.out=300)
      yy  <- predict(fit, data.frame(dates=xx))
      lines(xx, yy, col=cols[sites==si], lwd=2)
    }
  }
  legend("topright", legend=sites, col=cols, pch=16, cex=0.8, bty="n")
  
  # --- MID: average proportion of reads
  plot(d1$dates, d1$avg,
       xlim=c(2024,1760), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d1$avg>0,16,4), col=colv1, bty="n")
  dev.off()
  
}








## Lets focus on some ice associated diatoms 

iceSPP <-c("Polarella_glacialis","Islandinium_minutum","Islandinium_tricingulatum")

euk.ice <- euk[euk.tax$Species%in%iceSPP & 
                 euk.tax$assignmentQual%in%c("High","High-MH","High_MH-S"),]
euk.ice.2 <- euk.ice

# divide by original read number to get proportion 
euk.ice[,1:992] <- sweep(euk.ice[,1:992],2,colSums(euk),FUN="/")
euk.ice[euk.ice=="NaN"] <- 0
euk.ice[is.na(euk.ice)] <- 0


euk.ice.2[,1:992] <- sweep(euk.ice.2[,1:992],2,colSums(euk[rownames(euk)%in% ASVs.dino, , drop = FALSE]),FUN="/")
euk.ice.2[euk.ice.2=="NaN"] <- 0
euk.ice.2[is.na(euk.ice.2)] <- 0


#Aggregate same species
datain <- cbind(euk.ice[,1:992],"Species"=euk.tax$species[match(rownames(euk.ice),euk.tax$OTU)])
euk.ice.sum <- aggregate(.~Species,data=datain,FUN=sum)
rownames(euk.ice.sum) <- euk.ice.sum$Species
euk.ice.sum<- euk.ice.sum[,-1]

datain <- cbind(euk.ice.2[,1:992],"Species"=euk.tax$species[match(rownames(euk.ice.2),euk.tax$OTU)])
euk.ice.2.sum <- aggregate(.~Species,data=datain,FUN=sum)
rownames(euk.ice.2.sum) <- euk.ice.2.sum$Species
euk.ice.2.sum<- euk.ice.2.sum[,-1]





## Lets make an nREPs dataset 

euk.ice.sum.nreps <- NrepsMaker(euk.ice.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(euk.ice.sum)))

df <- tibble(spp = rownames(euk.ice.sum), euk.ice.sum.nreps)

euk.ice.sum.nreps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

euk.ice.sum.nreps.long$site <- expandedDates$Site[match(euk.ice.sum.nreps.long$column,expandedDates$SampleName)]
euk.ice.sum.nreps.long$dates <- expandedDates$CalYear[match(euk.ice.sum.nreps.long$column,expandedDates$SampleName)]

## lets make a avr dataset

euk.ice.sum.avr <- average_by_reps(euk.ice.sum,gsub("(.*).[0-9]$","\\1",colnames(euk.ice.sum)))


df <- tibble(spp = rownames(euk.ice.sum), euk.ice.sum.avr)

euk.ice.sum.avr.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

euk.ice.sum.avr.long$site <- expandedDates$Site[match(euk.ice.sum.avr.long$column,expandedDates$SampleName)]
euk.ice.sum.avr.long$dates <- expandedDates$CalYear[match(euk.ice.sum.avr.long$column,expandedDates$SampleName)]

#finally a dino normalised version 
euk.ice.2.sum.avr <- average_by_reps(euk.ice.2.sum,gsub("(.*).[0-9]$","\\1",colnames(euk.ice.2.sum)))

df <- tibble(spp = rownames(euk.ice.2.sum), euk.ice.2.sum.avr)

euk.ice.2.sum.avr.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

euk.ice.2.sum.avr.long$site <- expandedDates$Site[match(euk.ice.2.sum.avr.long$column,expandedDates$SampleName)]
euk.ice.2.sum.avr.long$dates <- expandedDates$CalYear[match(euk.ice.2.sum.avr.long$column,expandedDates$SampleName)]





#change the values a bit so they match the dates

pdf("figures/EUK.ice.v1.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))

for (species in unique(euk.ice.sum.nreps.long$spp)){
  
  loopdat <- euk.ice.sum.nreps.long[euk.ice.sum.nreps.long$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))
  
  
}

dev.off()
pdf("figures/EUK.ice.v1.trends.pdf", width = 9, height = 3.5)
par(mfrow = c(1, 3))

for (species in unique(euk.ice.sum.nreps.long$spp)) {
  
  loopdat <- euk.ice.sum.nreps.long[euk.ice.sum.nreps.long$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)],0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
       ylim = c(0, 9),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    # sort by date to ensure smooth lines
    subdat <- subdat[order(subdat$dates), ]
    # only fit if enough points
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit), 
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()



##### an idea here of a new plot

par(mfrow = c(3, 1))

for (species in unique(euk.ice.sum.avr.long$spp)) {
  
  loopdat <- euk.ice.sum.avr.long[euk.ice.sum.avr.long$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       loopdat$value,
       col =coredata$col[match(toupper(loopdat$site), coredata$Site)],
       main = species,
       xlab = "Date",
       xlim=c(2020,1780),
       pch = ifelse(loopdat$value > 0, 16, 4),bty="n")

}

## now lets combine them into one

sp_list <- intersect(unique(euk.ice.sum.nreps.long$spp),
                     unique(euk.ice.sum.avr.long$spp))

for (sp in sp_list) {
  pdf(paste0("figures/EUK.ice.big.",sp,".pdf"), width=9, height=6)
  par(mfrow=c(3,1), mar=c(3.5,4,2.5,1), oma=c(0,0,0,0), xaxs="i")
  
  # --- data per species
  d1 <- euk.ice.sum.nreps.long[euk.ice.sum.nreps.long$spp==sp, ]
  d2 <- euk.ice.sum.avr.long  [euk.ice.sum.avr.long$spp==sp, ]
  d3 <- euk.ice.2.sum.avr.long  [euk.ice.2.sum.avr.long$spp==sp, ]
  d1 <- d1[!is.na(d1$dates), ]
  d2 <- d2[!is.na(d2$dates), ]
  d3 <- d3[!is.na(d3$dates), ]
  
  # site-color map
  sites <- sort(unique(toupper(c(d1$site, d2$site))))
  cols  <- coredata$col[match(sites, coredata$Site)]
  colv1 <- cols[match(toupper(d1$site), sites)]
  colv2 <- cols[match(toupper(d2$site), sites)]
  colv3 <- cols[match(toupper(d3$site), sites)]
  
  # --- TOP: nREps
  set.seed(1234)
  plot(d1$dates, jitter(d1$value, 0.8),
       xlim=c(2024,1760),
       xlab="", ylab="nPCR reps", main=sp,
       pch=ifelse(d1$value>0,16,4), col=colv1, bty="n")
  for (si in sites) {
    d1s <- d1[toupper(d1$site)==si, ]
    if (nrow(d1s) > 2) {
      fit <- loess(value ~ dates, data=d1s, span=0.6, degree=1)
      xx  <- seq(min(d1s$dates), max(d1s$dates), length.out=300)
      yy  <- predict(fit, data.frame(dates=xx))
      lines(xx, yy, col=cols[sites==si], lwd=2)
    }
  }
  legend("topright", legend=sites, col=cols, pch=16, cex=0.8, bty="n")
  
  # --- MID: average proportion of reads
  plot(d2$dates, d2$value,
       xlim=c(2024,1760), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d2$value>0,16,4), col=colv2, bty="n")


plot(d3$dates, d3$value,
     xlim=c(2024,1760), 
     xlab="Calibrated year", ylab="Avg proportion in Dinos",
     pch=ifelse(d3$value>0,16,4), col=colv3, bty="n")
dev.off()
}

## now the phaeocystis

new_asvs  <- c("ASV_95","ASV_1748")

euk.pha <- euk[rownames(euk)%in%new_asvs,]
# divide by original read number to get proportion 
euk.pha[,1:992] <- sweep(euk.pha[,1:992],2,colSums(euk),FUN="/")
euk.pha[euk.pha=="NaN"] <- 0
euk.pha[is.na(euk.pha)] <- 0


#Aggregate same species
euk.pha.sum<- colSums(euk.pha)

## Lets make an nREPs dataset 

euk.pha.sum.nreps <- NrepsMaker(euk.pha.sum*100000000000000,gsub("(.*).[0-9]$","\\1",names(euk.pha.sum)))

df <- tibble(spp = rownames(euk.pha.sum), euk.pha.sum.nreps)

euk.pha.sum.nreps.long <- data.frame("Site1"=names(euk.pha.sum.nreps),
                                     "nREPs"=unlist(unname(euk.pha.sum.nreps)),
                                     "site"=expandedDates$Site[match(names(euk.pha.sum.nreps),expandedDates$SampleName)],
                                     "dates"=expandedDates$CalYear[match(names(euk.pha.sum.nreps),expandedDates$SampleName)])
## lets make a avr dataset
euk.pha.sum.avr <- tapply(euk.pha.sum, gsub("(.*).[0-9]$","\\1",names(euk.pha.sum)), mean, na.rm=TRUE)

euk.pha.sum.avr.long <- data.frame("Site1"=names(euk.pha.sum.avr),
                                     "value"=unlist(unname(euk.pha.sum.avr)),
                                     "site"=expandedDates$Site[match(names(euk.pha.sum.avr),expandedDates$SampleName)],
                                     "dates"=expandedDates$CalYear[match(names(euk.pha.sum.avr),expandedDates$SampleName)])


  # --- data per species
  d1 <- euk.pha.sum.nreps.long
  d2 <- euk.pha.sum.avr.long  

 
  
  # site-color map
  sites <- sort(unique(toupper(c(d1$site, d2$site))))
  cols  <- coredata$col[match(sites, coredata$Site)]
  colv1 <- cols[match(toupper(d1$site), sites)]
  colv2 <- cols[match(toupper(d2$site), sites)]
  
  # --- TOP: nREps
  pdf("figures/EUK.p.antarctica.pdf",, width=9, height=6)
  par(mfrow=c(3,1), mar=c(3.5,4,2.5,1), oma=c(0,0,0,0), xaxs="i")
  set.seed(1234)
  par(mfrow=c(2,1))
  plot(d1$dates, jitter(d1$nREPs, 0.8),
       xlim=c(2024,1760),
       xlab="", ylab="nPCR reps", main="P antarctica",
       pch=ifelse(d1$nREPs>0,16,4), col=colv1, bty="n")
  for (si in sites) {
    d1s <- d1[toupper(d1$site)==si, ]
    if (nrow(d1s) > 2) {
      fit <- loess(nREPs ~ dates, data=d1s, span=0.6, degree=1)
      xx  <- seq(min(d1s$dates), max(d1s$dates), length.out=300)
      yy  <- predict(fit, data.frame(dates=xx))
      lines(xx, yy, col=cols[sites==si], lwd=2)
    }
  }
  # --- MID: average proportion of reads
  plot(d2$dates, d2$value,
       xlim=c(2024,1760), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d2$value>0,16,4), col=colv2, bty="n")
  legend(1780,0.007, legend=sites, col=cols, pch=16, cex=0.8, bty="n")
  dev.off()







### How aboiut this HAB species 
spp <- "Karlodinium veneficum" 


euk.kar <- euk[euk.tax$species%in%spp & 
                 euk.tax$assignmentQual%in%c("High","High-MH","High_MH-S"),]

# divide by original read number to get proportion 
euk.kar[,1:992] <- sweep(euk.kar[,1:992],2,colSums(euk),FUN="/")
euk.kar[euk.kar=="NaN"] <- 0
euk.kar[is.na(euk.kar)] <- 0
euk.kar.sum <- colSums(euk.kar)

## Lets make an nREPs dataset 

euk.kar.sum.nreps <- NrepsMaker(euk.kar.sum*100000000000000,gsub("(.*).[0-9]$","\\1",names(euk.kar.sum)))

df <- tibble(spp = rownames(euk.kar.sum.nreps), euk.kar.sum.nreps)

euk.kar.sum.nreps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

euk.kar.sum.nreps.long$site <- expandedDates$Site[match(euk.kar.sum.nreps.long$column,expandedDates$SampleName)]
euk.kar.sum.nreps.long$dates <- expandedDates$CalYear[match(euk.kar.sum.nreps.long$column,expandedDates$SampleName)]
loopdat <- euk.kar.sum.nreps.long

par(mfrow = c(1, 1))
set.seed(1234)
plot(loopdat$dates,
     jitter(loopdat$value, 0.8),
     col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)],0.5),
     main = spp,
     xlab = "Date",
     ylab = "nPCR reps",
     ylim = c(0, 9),
     pch = ifelse(loopdat$value > 0, 16, 4))

# Add LOESS smooth per site
for (s in unique(loopdat$site)) {
  subdat <- loopdat[loopdat$site == s, ]
  # sort by date to ensure smooth lines
  subdat <- subdat[order(subdat$dates), ]
  # only fit if enough points
  if (nrow(subdat) > 4) {
    fit <- loess(value ~ dates, data = subdat, span = 0.8)
    lines(subdat$dates, predict(fit), 
          col = coredata$col[match(toupper(s), coredata$Site)],
          lwd = 2)
  }
}

## how about these genera
genera <- c("Sanguina","Ascoseira","Deepstaria")


euk.gen <- euk[euk.tax$genus%in%genera & 
                 euk.tax$assignmentQual%in%c("High","High-MH","High_MH-S"),]

# divide by original read number to get proportion 
euk.gen[,1:992] <- sweep(euk.gen[,1:992],2,colSums(euk),FUN="/")
euk.gen[euk.gen=="NaN"] <- 0
euk.gen[is.na(euk.gen)] <- 0



#Aggregate same species
datain <- cbind(euk.gen[,1:992],"Genus"=euk.tax$genus[match(rownames(euk.gen),euk.tax$OTU)])
euk.gen.sum <- aggregate(.~Genus,data=datain,FUN=sum)
rownames(euk.gen.sum) <- euk.gen.sum$Genus
euk.gen.sum<- euk.gen.sum[,-1]

## Lets make an nREPs dataset 

euk.gen.sum.nreps <- NrepsMaker(euk.gen.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(euk.gen.sum)))

df <- tibble(spp = rownames(euk.gen.sum.nreps), euk.gen.sum.nreps)

euk.gen.sum.nreps.long <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

euk.gen.sum.nreps.long$site <- expandedDates$Site[match(euk.gen.sum.nreps.long$column,expandedDates$SampleName)]
euk.gen.sum.nreps.long$dates <- expandedDates$CalYear[match(euk.gen.sum.nreps.long$column,expandedDates$SampleName)]


#change the values a bit so they match the dates

pdf("figures/EUK.gen.pdf",width = 9,height =3.5)
par(mfrow=c(1,3))

for (species in unique(euk.gen.sum.nreps.long$spp)){
  
  loopdat <- euk.gen.sum.nreps.long[euk.gen.sum.nreps.long$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))
  
  
}

dev.off()
pdf("figures/EUK.gen.trends.pdf", width = 9, height = 3.5)
par(mfrow = c(1, 3))

for (species in unique(euk.gen.sum.nreps.long$spp)) {
  
  loopdat <- euk.gen.sum.nreps.long[euk.gen.sum.nreps.long$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)],0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
       ylim = c(0, 9),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    # sort by date to ensure smooth lines
    subdat <- subdat[order(subdat$dates), ]
    # only fit if enough points
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit), 
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()




### Now let's examine the whales
mam <- read.csv("cleandata/MAM.WP5.main.clean.csv",row.names = 1)
mam.tax <- read.csv("cleandata/MAM.WP5.main.clean.tax.csv",row.names = 1)
mam.tax.whales.asvs <- na.omit(mam.tax$ASV[mam.tax$family=="Balaenopteridae"|
                    mam.tax$family=="Balaenidae"|
                    mam.tax$family=="Delphinidae"])

mam.whales <- cbind(mam[match(mam.tax.whales.asvs,rownames(mam)),],mam.tax[match(mam.tax.whales.asvs,mam.tax$ASV),])

riz <- read.csv("cleandata/riz.WP5.main.clean.csv",row.names = 1)
riz.tax <- read.csv("cleandata/riz.WP5.main.clean.tax.csv",row.names = 1)
riz.tax.whales.asvs <- na.omit(riz.tax$ASV[riz.tax$family=="Balaenopteridae"|
                                             riz.tax$family=="Balaenidae"|
                                             riz.tax$family=="Delphinidae"])

riz.whales <- cbind(riz[match(riz.tax.whales.asvs,rownames(riz)),],riz.tax[match(riz.tax.whales.asvs,riz.tax$ASV),])


write.csv(mam.whales,"cleandata/MAM.whales.clean.csv")
write.csv(riz.whales,"cleandata/RIZ.whales.clean.csv")

#Now we check all those whales manually and pull in the checked data  
mam.whales.man <- read.csv("cleandata/MAM.whales.clean.manual.csv",row.names = 1)
riz.whales.man <- read.csv("cleandata/RIZ.whales.clean.manual.csv",row.names = 1)


# divide by original read number to get proportion 
mam.whales.man[1:983] <- sweep(mam.whales.man[1:983],2,colSums(mam),FUN="/")
mam.whales.man[mam.whales.man=="NaN"] <- 0
mam.whales.man[is.na(mam.whales.man)] <- 0

riz.whales.man[1:813] <- sweep(riz.whales.man[1:813],2,colSums(riz),FUN="/")
riz.whales.man[riz.whales.man=="NaN"] <- 0
riz.whales.man[is.na(riz.whales.man)] <- 0

#Aggregate same species
datain <- cbind(mam.whales.man[,1:983],"manualSpecies"=mam.whales.man$ManualSpecies)
mam.whales.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(mam.whales.sum) <- mam.whales.sum$manualSpecies
mam.whales.sum <- mam.whales.sum[,-1]
mam.whales.sum <- mam.whales.sum[!rownames(mam.whales.sum)=="",]
colnames(mam.whales.sum) <- gsub("^MAM[._]", "", colnames(mam.whales.sum))

datain <- cbind(riz.whales.man[,1:813],"manualSpecies"=riz.whales.man$ManualSpecies)
riz.whales.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(riz.whales.sum) <- riz.whales.sum$manualSpecies
riz.whales.sum <- riz.whales.sum[,-1]
riz.whales.sum <- riz.whales.sum[!rownames(riz.whales.sum)=="",]
colnames(riz.whales.sum) <- gsub("^RIZ[._]", "", colnames(riz.whales.sum))


### summing observations

# Ensure both are matrices
#riz <- as.matrix(riz.whales.sum)
#mam <- as.matrix(mam.whales.sum)

riz <- as.matrix(NrepsMaker(riz.whales.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(riz.whales.sum))))
mam <- as.matrix(NrepsMaker(mam.whales.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(mam.whales.sum))))


# Align dimensions via padding with 0s
all_species <- union(rownames(riz), rownames(mam))
all_samples <- union(colnames(riz), colnames(mam))

pad_matrix <- function(mat, all_rows, all_cols) {
  out <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
                dimnames = list(all_rows, all_cols))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

riz_pad <- pad_matrix(riz, all_species, all_samples)
mam_pad <- pad_matrix(mam, all_species, all_samples)

# Final merged result
whales.sum.nreps <-as.data.frame(riz_pad + mam_pad)


df <- tibble(spp = rownames(whales.sum.nreps), whales.sum.nreps)

whales.sum.nreps <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

whales.sum.nreps$site <- expandedDates$Site[match(whales.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]
whales.sum.nreps$dates <- expandedDates$CalYear[match(whales.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]


#Plots

pdf("figures/whales.pdf",width = 9,height =6.5)
par(mfrow=c(2,3))

for (species in unique(whales.sum.nreps$spp)){
  
  loopdat <- whales.sum.nreps[whales.sum.nreps$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))
  
  
}

dev.off()
pdf("figures/whales.loess.pdf", width = 9, height = 6.5)
par(mfrow = c(2, 3))

for (species in unique(whales.sum.nreps$spp)) {
  
  loopdat <- whales.sum.nreps[whales.sum.nreps$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
  ylim = c(0, 12),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    subdat <- subdat[order(subdat$dates), ]
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit),
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()



### Now let's examine the seals
mam <- read.csv("cleandata/MAM.WP5.main.clean.csv",row.names = 1)
mam.tax <- read.csv("cleandata/MAM.WP5.main.clean.tax.csv",row.names = 1)
mam.tax.seals.asvs <- na.omit(mam.tax$ASV[mam.tax$family=="Otariidae"|
                                             mam.tax$family=="Phocidae"])
mam.seals <- cbind(mam[match(mam.tax.seals.asvs,rownames(mam)),],mam.tax[match(mam.tax.seals.asvs,mam.tax$ASV),])

riz <- read.csv("cleandata/riz.WP5.main.clean.csv",row.names = 1)
riz.tax <- read.csv("cleandata/riz.WP5.main.clean.tax.csv",row.names = 1)
riz.tax.seals.asvs <- na.omit(riz.tax$ASV[riz.tax$family=="Otariidae"|
                                             riz.tax$family=="Phocidae"])

riz.seals <- cbind(riz[match(riz.tax.seals.asvs,rownames(riz)),],riz.tax[match(riz.tax.seals.asvs,riz.tax$ASV),])


write.csv(mam.seals,"cleandata/MAM.seals.clean.csv")
write.csv(riz.seals,"cleandata/RIZ.seals.clean.csv")

#Now we check all those seals manually and pull in the checked data  
mam.seals.man <- read.csv("cleandata/MAM.seals.clean.manual.csv",row.names = 1)
riz.seals.man <- read.csv("cleandata/RIZ.seals.clean.manual.csv",row.names = 1)


# divide by original read number to get proportion 
mam.seals.man[1:983] <- sweep(mam.seals.man[1:983],2,colSums(mam),FUN="/")
mam.seals.man[mam.seals.man=="NaN"] <- 0
mam.seals.man[is.na(mam.seals.man)] <- 0

riz.seals.man[1:813] <- sweep(riz.seals.man[1:813],2,colSums(riz),FUN="/")
riz.seals.man[riz.seals.man=="NaN"] <- 0
riz.seals.man[is.na(riz.seals.man)] <- 0

#Aggregate same species
datain <- cbind(mam.seals.man[,1:983],"manualSpecies"=mam.seals.man$ManualSpecies)
mam.seals.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(mam.seals.sum) <- mam.seals.sum$manualSpecies
mam.seals.sum <- mam.seals.sum[,-1]
mam.seals.sum <- mam.seals.sum[!rownames(mam.seals.sum)=="",]
colnames(mam.seals.sum) <- gsub("^MAM[._]", "", colnames(mam.seals.sum))

datain <- cbind(riz.seals.man[,1:813],"manualSpecies"=riz.seals.man$ManualSpecies)
riz.seals.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(riz.seals.sum) <- riz.seals.sum$manualSpecies
riz.seals.sum <- riz.seals.sum[,-1]
riz.seals.sum <- riz.seals.sum[!rownames(riz.seals.sum)=="",]
colnames(riz.seals.sum) <- gsub("^RIZ[._]", "", colnames(riz.seals.sum))


### summing observations

# Ensure both are nreps matrices
riz <- as.matrix(NrepsMaker(riz.seals.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(riz.seals.sum))))
mam <- as.matrix(NrepsMaker(mam.seals.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(mam.seals.sum))))


# Align dimensions via padding with 0s
all_species <- union(rownames(riz), rownames(mam))
all_samples <- union(colnames(riz), colnames(mam))

pad_matrix <- function(mat, all_rows, all_cols) {
  out <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
                dimnames = list(all_rows, all_cols))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

riz_pad <- pad_matrix(riz, all_species, all_samples)
mam_pad <- pad_matrix(mam, all_species, all_samples)

# Final merged result
seals.sum.nreps <-as.data.frame(riz_pad + mam_pad)


df <- tibble(spp = rownames(seals.sum.nreps), seals.sum.nreps)

seals.sum.nreps <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

seals.sum.nreps$site <- expandedDates$Site[match(seals.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]
seals.sum.nreps$dates <- expandedDates$CalYear[match(seals.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]




#Plots
pdf("figures/seals.pdf",width = 9,height =6.5)
par(mfrow=c(2,3))

for (species in unique(seals.sum.nreps$spp)){
  
  loopdat <- seals.sum.nreps[seals.sum.nreps$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))
  
  
}

dev.off()
pdf("figures/seals.loess.pdf", width = 9, height = 6.5)
par(mfrow = c(2, 3))

for (species in unique(seals.sum.nreps$spp)) {
  
  loopdat <- seals.sum.nreps[seals.sum.nreps$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
       ylim = c(0, 8),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    subdat <- subdat[order(subdat$dates), ]
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit),
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()

#### Now fish 

mam <- read.csv("cleandata/MAM.WP5.main.clean.csv",row.names = 1)
mam.tax <- read.csv("cleandata/MAM.WP5.main.clean.tax.csv",row.names = 1)
mam.tax.fish.asvs <- na.omit(mam.tax$ASV[mam.tax$class=="Actinopteri"])

mam.fish <- cbind(mam[match(mam.tax.fish.asvs,rownames(mam)),],mam.tax[match(mam.tax.fish.asvs,mam.tax$ASV),])

riz <- read.csv("cleandata/riz.WP5.main.clean.csv",row.names = 1)
riz.tax <- read.csv("cleandata/riz.WP5.main.clean.tax.csv",row.names = 1)
riz.tax.fish.asvs <- na.omit(riz.tax$ASV[riz.tax$class=="Actinopteri"])

riz.fish <- cbind(riz[match(riz.tax.fish.asvs,rownames(riz)),],riz.tax[match(riz.tax.fish.asvs,riz.tax$ASV),])


write.csv(mam.fish,"cleandata/MAM.fish.clean.csv")
write.csv(riz.fish,"cleandata/RIZ.fish.clean.csv")

#Now we check all those fish manually and pull in the checked data  (none for MAM)
riz.fish.man <- read.csv("cleandata/RIZ.fish.clean.manual.csv",row.names = 1)


# divide by original read number to get proportion 
riz.fish.man[1:813] <- sweep(riz.fish.man[1:813],2,colSums(riz),FUN="/")
riz.fish.man[riz.fish.man=="NaN"] <- 0
riz.fish.man[is.na(riz.fish.man)] <- 0

#Aggregate same species
datain <- cbind(riz.fish.man[,1:813],"manualSpecies"=riz.fish.man$ManualSpecies)
riz.fish.sum <- aggregate(.~manualSpecies,data=datain,FUN=sum)
rownames(riz.fish.sum) <- riz.fish.sum$manualSpecies
riz.fish.sum <- riz.fish.sum[,-1]
riz.fish.sum <- riz.fish.sum[!rownames(riz.fish.sum)=="",]
colnames(riz.fish.sum) <- gsub("^RIZ[._]", "", colnames(riz.fish.sum))




fish.sum.nreps <- NrepsMaker(riz.fish.sum*100000000000000,gsub("(.*).[0-9]$","\\1",colnames(riz.fish.sum)))

df <- tibble(spp = rownames(fish.sum.nreps), fish.sum.nreps)

fish.sum.nreps <- df %>% 
  pivot_longer(cols = -spp, 
               names_to = "column", 
               values_to = "value") 

fish.sum.nreps$site <- expandedDates$Site[match(fish.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]
fish.sum.nreps$dates <- expandedDates$CalYear[match(fish.sum.nreps$column,gsub("^EUK[._]", "",expandedDates$SampleName))]


#Plots
pdf("figures/fish.pdf",width = 12,height =8)
par(mfrow=c(3,5))

for (species in unique(fish.sum.nreps$spp)){
  
  loopdat <- fish.sum.nreps[fish.sum.nreps$spp==species,]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value,0.8),
       col=coredata$col[match(toupper(loopdat$site),coredata$Site)],
       main=species,
       xlab="Date",
       ylab="nPCR reps",
       ylim=c(0,8),
       pch=ifelse(loopdat$value>0,16,4))
  
  
}

dev.off()
pdf("figures/fish.loess.pdf", width = 12, height = 8)
par(mfrow = c(3, 5))

for (species in unique(fish.sum.nreps$spp)) {
  
  loopdat <- fish.sum.nreps[fish.sum.nreps$spp == species, ]
  set.seed(1234)
  plot(loopdat$dates,
       jitter(loopdat$value, 0.8),
       col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.5),
       main = species,
       xlab = "Date",
       ylab = "nPCR reps",
       ylim = c(0, 8),
       pch = ifelse(loopdat$value > 0, 16, 4))
  
  # Add LOESS smooth per site
  for (s in unique(loopdat$site)) {
    subdat <- loopdat[loopdat$site == s, ]
    subdat <- subdat[order(subdat$dates), ]
    if (nrow(subdat) > 4) {
      fit <- loess(value ~ dates, data = subdat, span = 0.8)
      lines(subdat$dates, predict(fit),
            col = coredata$col[match(toupper(s), coredata$Site)],
            lwd = 2)
    }
  }
}

dev.off()





####====3.0 Polarella ===####

euk.polarella <- euk[euk.tax$genus=="Polarella" & euk.tax$assignmentQual=="High",]
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


