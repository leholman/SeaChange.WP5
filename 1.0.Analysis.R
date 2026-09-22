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
library(rioja)
library(mgcv)


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
site_cols <- setNames(coredata$col, tolower(coredata$Site))
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
## Prokaryotic (bacterial + archaeal) subdivisions recovered as 18S off-target
## amplification. Defined once here and reused wherever prokaryotes are grouped.
prok.subdiv <- c(
  "Proteobacteria_X", "Firmicutes_X", "Verrucomicrobia_X", "Planctomycetes_X",
  "Acidobacteria_X", "Nitrospirae_X", "Dadabacteria_X", "Thaumarchaeota_X",
  "Kiritimatiellaeota_X", "Nitrospinae_X", "Actinobacteria_X", "Bacteroidetes_X",
  "Hydrogenedentes_X", "Euryarchaeota_X", "Epsilonbacteraeota_X"
)

ASVs.bac <- euk.pr2.2$X.1[
  euk.pr2.2$Subdivision %in% prok.subdiv & euk.pr2.2$Subdivision.1 > 80
]

## Coverage check: which subdivisions fall into none of the three subsets?
local({
  present  <- unique(na.omit(euk.pr2.2$Subdivision[
    match(rownames(euk.reps), euk.pr2.2$X.1)]))
  assigned <- unique(na.omit(euk.pr2.2$Subdivision[
    euk.pr2.2$X.1 %in% c(ASVs.met, ASVs.pro, ASVs.bac)]))
  miss <- setdiff(present, assigned)
  if (length(miss)) cat("Subdivisions in no subset:", paste(miss, collapse = ", "), "\n")
  ghost <- setdiff(c(prok.subdiv), present)
  if (length(ghost)) cat("Named but absent from data:", paste(ghost, collapse = ", "), "\n")
})

##pull in this little class classifier - note check this if it becomes important 
classes <- read.csv("taxonomy/class_category.csv")

ASVs.diatom <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Diatom"] & euk.pr2.2$Class.1 > 80]
ASVs.dino <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Dinoflagellate"] & euk.pr2.2$Class.1 > 80]
ASVs.smolphy <- euk.pr2.2$X.1[euk.pr2.2$Class %in%  classes$Class[classes$Category=="Pico/nanophytoplankton"] & euk.pr2.2$Class.1 > 80]


## Taxonomic overview

#EUK taxonomy all family 
EUK.tax.a <- minAbundance(CountTable(as.character(euk.pr2.2$Subdivision),euk.reps,output = "Abundance"),minAbun=0.01)
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

expandedDates$originalCalYear <- expandedDates$CalYear

expandedDates$CalYear <- round(expandedDates$CalYear)

write.csv(expandedDates,"metadata/expandeddates.csv")

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
  lapply(sites, make_site_year, M = EUK.tax.a.2[dim(EUK.tax.a.2)[1]:1,], map = expandedDates, years = all_years),
  sites
)

# Consistent spacing and limits
sp <- 0.2
ylim <- c(0, 1)
pal <- getPalette(nrow(EUK.tax.a.2)) # you already defined getPalette

pdf("figures/EUK.bigplot.pdf",width = 14,height = 10)
par(mfrow = c(length(sites), 1),
    mar = c(0.2, 4.1, 0.5, 1.1),  # small bottom margins for all but last
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
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.70)
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
       fill   = rev(pal),
       border = NA,
       cex    = 0.6,
       bty    = "n",
       inset  = c(-0.01, 0), xpd = NA)


dev.off()


##### How about this 
# Common, ordered time axis across all sites
all_years <- min(unique(na.omit(expandedDates$CalYear))):max(unique(na.omit(expandedDates$CalYear)))
sites     <- unique(expandedDates$Site)

# Helper: make a site x year matrix of *proportions* (averaging if multiple samples per year)
make_site_year <- function(site, M, map, years){
  out <- sapply(years, function(y){
    smps <- map$SampleName[map$Site == site & map$CalYear == y]
    if (length(smps) == 0) {
      rep(0, nrow(M))
    } else if (length(smps) == 1) {
      M[, smps, drop = FALSE][, 1]
    } else {
      rowMeans(M[, smps, drop = FALSE], na.rm = TRUE)
    }
  })
  rownames(out) <- rownames(M)
  colnames(out) <- years
  out
}

# Draw light "ribbons" connecting non-zero stacked bars across zero gaps
connect_observed <- function(M, bp, cols, alpha = 0.35, minval = 0){
  M[is.na(M)] <- 0
  obs <- which(colSums(M) > 0)
  if (length(obs) < 2) return(invisible(NULL))
  cs <- apply(M, 2, cumsum)
  bt <- rbind(0, cs[-nrow(M), , drop = FALSE])
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

site_mats <- setNames(
  lapply(sites, make_site_year,
         M = EUK.tax.a.2[dim(EUK.tax.a.2)[1]:1, ],
         map = expandedDates,
         years = all_years),
  sites
)

# --- Pastel, colourblind-friendly palette for up to ~15 categories ---
pastelize <- function(hex, mix = 0.40){
  rgb_mat <- col2rgb(hex) / 255
  white   <- matrix(1, nrow = 3, ncol = ncol(rgb_mat))
  out     <- rgb_mat * (1 - mix) + white * mix
  rgb(out[1, ], out[2, ], out[3, ])
}


# --- Soft palette with more contrast between taxa ---
soft_pal <- c(
  "Unknown"          = "#BDB7AE",  # darker warm grey — recedes
  "Apicomplexa"      = "#9B7EBD",  # deeper lilac/purple
  "Bigyra"           = "#5FA8A3",  # teal — saturated but soft
  "Cercozoa"         = "#E89F71",  # warm peach-orange
  "Chlorophyta_X"    = "#B8D87A",  # bright sage — pops
  "Ciliophora"       = "#6FA8DC",  # clear powder blue
  "Dinoflagellata"   = "#E87A8A",  # rosy pink — saturated
  "Fungi"            = "#8B6F47",  # warm brown — earthy, distinct
  "Gyrista"          = "#4A9D6C",  # deeper diatom green
  "Ichthyosporea"    = "#F2C9A0",  # pale apricot
  "Metazoa"          = "#C75B93",  # magenta — animal pop
  "Proteobacteria_X" = "#3D5A80",  # deep slate blue
  "Radiolaria"       = "#F0D060",  # golden yellow
  "Others"           = "#8A8A7A"   # darker olive-grey — recedes
)

pal <- rev(soft_pal[rownames(EUK.tax.a.2)])

# Consistent spacing and limits
sp   <- 0.2
ylim <- c(0, 1)

# --- X-axis tick positions: every ~25 years, on round numbers ---
tick_step  <- 25
tick_years <- seq(1775,
                  2025,
                  by = tick_step)


pdf("figures/EUK.bigplot.2.pdf", width = 10, height = 6.5)
par(mfrow = c(length(sites), 1),
    mar = c(0.2, 4.1, 0.5, 1.1),
    oma = c(6, 0, 0, 0),
    xpd = NA,
    las = 1)

site_order <- c(2, 4, 3, 1, 5)

for (idx in seq_along(site_order)) {
  i    <- site_order[idx]
  site <- sites[i]
  mat  <- site_mats[[site]]
  mat[, colSums(mat) == 0] <- NA
  
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width  = 30,
                ylim   = ylim,
                border = NA,
                xaxt   = "n",
                yaxt   = "s")
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.70)
  barplot(mat,
          col    = pal,
          space  = 0,
          width  = 30,
          ylim   = ylim,
          border = NA,
          xaxt   = "n",
          yaxt   = "n",
          add    = TRUE)
  
  mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
  
  # Only the bottom panel gets the shared x-axis
  if (idx == length(site_order)) {
    # Map calendar years to bar-plot positions (bp is aligned with all_years)
    tick_pos <- bp[match(tick_years, all_years)]
    keep     <- !is.na(tick_pos)
    axis(1,
         at      = tick_pos[keep],
         labels  = tick_years[keep],
         las     = 1,
         cex.axis = 1.5)
    mtext("", side = 1, line = 4)
  }
}

# legend
legend(0,2.4,
       legend = gsub("_X","",rownames(EUK.tax.a.2)),
       fill   = rev(pal),            # match plotting order (no rev needed)
       border = NA,
       cex    = 1.2,
       bty    = "n",
       inset  = c(-0.01, 0),
       xpd    = NA)

dev.off()

#### NEW TEST 
pdf("figures/DraftFigures/Fig.2/EUK.bigplot.2v2.pdf", width = 10, height = 6.5)

par(mfrow = c(length(sites), 1),
    mar   = c(0.2, 4.1, 0.5, 1.1),
    oma   = c(6, 0, 0, 0),
    xpd   = NA,
    las   = 1)

site_order <- c(2, 4, 3, 1, 5)

## shared y-axis ticks: drop the top one so adjacent panels don't collide
yticks <- pretty(ylim, n = 4)
yticks <- yticks[yticks >= ylim[1] & yticks < ylim[2]]

## right-hand padding, in bar-width units (bars are 30 wide, space = 0)
xpad <- 60

for (idx in seq_along(site_order)) {
  i    <- site_order[idx]
  site <- sites[i]
  mat  <- site_mats[[site]]
  mat[, colSums(mat) == 0] <- NA
  
  xr <- c(0, ncol(mat) * 30 + xpad)
  
  bp <- barplot(mat,
                col    = pal,
                space  = 0,
                width  = 30,
                ylim   = ylim,
                xlim   = xr,
                border = NA,
                xaxt   = "n",
                yaxt   = "n")
  
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.70)
  
  barplot(mat,
          col    = pal,
          space  = 0,
          width  = 30,
          ylim   = ylim,
          border = NA,
          xaxt   = "n",
          yaxt   = "n",
          add    = TRUE)
  
  ## y-axis on the right, sitting in the padded space
  axis(4,
       at       = c(yticks,1),
       pos      = ncol(mat) * 30 + 6,
       las      = 1,
       tcl      = -0.25,
       mgp      = c(3, 0.4, 0),
       cex.axis = 0.8)
  
  ## site name stays in the left margin
  mtext(site, side = 4, line = 0, las = 0, cex = 0.9)
  
  ## bottom panel only: shared x-axis
 # if (idx == length(site_order)) {
#    tick_pos <- bp[match(tick_years, all_years)]
#    keep     <- !is.na(tick_pos)
#    axis(1,
#         at       = tick_pos[keep],
#         labels   = tick_years[keep],
#         las      = 1,
#         cex.axis = 1.5)
#  }
}

## legend — fill order now matches the label order
legend(0, 2.4,
       legend = gsub("_X", "", rownames(EUK.tax.a.2)),
       fill   = pal,
       border = NA,
       cex    = 1.2,
       bty    = "n",
       xpd    = NA)

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
  connect_observed(M = mat, bp = bp, cols = pal, alpha = 0.70)
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
cutree(SC.met.con, k = 6)


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


####====Taxonomic overview + CONISS zoned big plots====####
## Assumes everything up to and including '## Taxonomic overview' has been run

####----0. Build the EUK Subdivision-level matrix used by the first big plot----####
EUK.tax.a <- minAbundance(
  CountTable(as.character(euk.pr2.2$Subdivision), euk.reps, output = "Abundance"),
  minAbun = 0.01
)
row.names(EUK.tax.a)[1] <- "Unknown"
EUK.tax.a.2 <- as.matrix(prop.table(
  as.matrix(EUK.tax.a[, sort(colnames(EUK.tax.a))]), margin = 2
))


####----1. Helpers shared by all big plots----####

make_site_year <- function(site, M, map, years){
  out <- sapply(years, function(y){
    smps <- map$SampleName[map$Site == site & map$CalYear == y]
    if (length(smps) == 0) {
      rep(0, nrow(M))
    } else if (length(smps) == 1) {
      M[, smps, drop = FALSE][, 1]
    } else {
      rowMeans(M[, smps, drop = FALSE], na.rm = TRUE)
    }
  })
  rownames(out) <- rownames(M)
  colnames(out) <- years
  out
}

connect_observed <- function(M, bp, cols, alpha = 0.35, minval = 0){
  M[is.na(M)] <- 0
  obs <- which(colSums(M) > 0)
  if (length(obs) < 2) return(invisible(NULL))
  cs <- apply(M, 2, cumsum)
  bt <- rbind(0, cs[-nrow(M), , drop = FALSE])
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

coniss.helper <- function(site, data){
  data.sites <- expandedDates$Site[match(colnames(data), expandedDates$SampleName)]
  data.sub   <- data[, grepl(site, data.sites), drop = FALSE]
  yrs        <- expandedDates$CalYear[match(colnames(data.sub), expandedDates$SampleName)]
  data.sort  <- data.sub[, order(yrs), drop = FALSE]
  d   <- vegdist(t(data.sort), method = "bray")
  con <- chclust(d, method = "coniss")
  con
}

coniss_boundaries <- function(con, k, map = expandedDates){
  if (is.null(con) || k < 2) return(numeric(0))
  cl  <- cutree(con, k = k)
  yrs <- map$CalYear[match(names(cl), map$SampleName)]
  br  <- which(diff(cl) != 0)
  if (length(br) == 0) return(numeric(0))
  (yrs[br] + yrs[br + 1]) / 2
}

draw_zone_lines <- function(boundary_years, bp, all_years,
                            col = "darkred", lwd = 2, lty = 5){
  if (length(boundary_years) == 0) return(invisible(NULL))
  x <- approx(x = all_years, y = bp, xout = boundary_years, rule = 2)$y
  abline(v = x, col = col, lwd = lwd, lty = lty, xpd = FALSE)
}


####----2. Build per-taxon community matrices----####

euk.prop <- prop.table(as.matrix(euk.reps), margin = 2)

## Metazoa
euk.met.only  <- euk.prop[rownames(euk.prop) %in% ASVs.met, , drop = FALSE]
euk.pr2.2.met <- euk.pr2.2[match(rownames(euk.met.only), euk.pr2.2$X.1), ]
met.fam <- euk.pr2.2.met$Family
met.fam[is.na(met.fam) | met.fam == ""] <- "Unknown"
EUK.met   <- minAbundance(CountTable(met.fam, euk.met.only, output = "Abundance"), minAbun = 0.01)
EUK.met.2 <- as.matrix(EUK.met)

## Diatoms
euk.dia.only  <- euk.prop[rownames(euk.prop) %in% ASVs.diatom, , drop = FALSE]
euk.pr2.2.dia <- euk.pr2.2[match(rownames(euk.dia.only), euk.pr2.2$X.1), ]
dia.fam <- euk.pr2.2.dia$Family
dia.fam[is.na(dia.fam) | dia.fam == ""] <- "Unknown"
EUK.dia   <- minAbundance(CountTable(dia.fam, euk.dia.only, output = "Abundance"), minAbun = 0.01)
EUK.dia.3 <- as.matrix(EUK.dia)

## Dinoflagellates
euk.din.only  <- euk.prop[rownames(euk.prop) %in% ASVs.dino, , drop = FALSE]
euk.pr2.2.din <- euk.pr2.2[match(rownames(euk.din.only), euk.pr2.2$X.1), ]
din.fam <- euk.pr2.2.din$Family
din.fam[is.na(din.fam) | din.fam == ""] <- "Unknown"
EUK.din   <- minAbundance(CountTable(din.fam, euk.din.only, output = "Abundance"), minAbun = 0.01)
EUK.din.2 <- as.matrix(EUK.din)

## Small phytoplankton
euk.smolphy.only  <- euk.prop[rownames(euk.prop) %in% ASVs.smolphy, , drop = FALSE]
euk.pr2.2.smolphy <- euk.pr2.2[match(rownames(euk.smolphy.only), euk.pr2.2$X.1), ]
smolphy.fam <- euk.pr2.2.smolphy$Family
smolphy.fam[is.na(smolphy.fam) | smolphy.fam == ""] <- "Unknown"
EUK.smolphy   <- minAbundance(CountTable(smolphy.fam, euk.smolphy.only, output = "Abundance"), minAbun = 0.01)
EUK.smolphy.2 <- as.matrix(EUK.smolphy)


####----3. Site/year scaffolding----####

site.map <- c(BB3 = "BB", scAØ2 = "SC", MC3 = "MC", p68 = "p68", p56 = "p56")

if (!exists("expandedDates")) {
  cn <- colnames(euk.reps)
  expandedDates <- merge(
    data.frame(SampleName = cn,
               Site       = site.map[regmatches(cn, regexpr("BB3|scAØ2|MC3|p68|p56", cn))],
               SampleID   = as.numeric(sub(".*(?:BB3|scAØ2|MC3|p68|p56)[._]", "", cn, perl = TRUE)),
               stringsAsFactors = FALSE),
    dates, by = c("Site", "SampleID"), all.x = TRUE
  )
  expandedDates$originalCalYear <- expandedDates$CalYear
  expandedDates$CalYear         <- round(expandedDates$CalYear)
}

all_years  <- min(unique(na.omit(expandedDates$CalYear))):max(unique(na.omit(expandedDates$CalYear)))
sites      <- unique(expandedDates$Site)
site_order <- c(2, 4, 3, 1, 5)


####----4. CONISS for every site x dataset combination----####

k_all  <- c(BB = 2, MC = 2, p56 = 4, p68 = 3, SC = 1)
k_met  <- c(BB = 2, MC = 5, p56 = 3, p68 = 3, SC = 6)
k_dia  <- c(BB = 3, MC = 3, p56 = 2, p68 = 3, SC = 1)
k_din  <- c(BB = 2, MC = 1, p56 = 2, p68 = 4, SC = 1)
k_smol <- c(BB = 2, MC = 1, p56 = 3, p68 = 2, SC = 2)

run_coniss_all_sites <- function(M){
  setNames(lapply(sites, function(s) coniss.helper(s, M)), sites)
}

con_all  <- run_coniss_all_sites(euk.prop)
con_met  <- run_coniss_all_sites(euk.met.only)
con_dia  <- run_coniss_all_sites(euk.dia.only)
con_din  <- run_coniss_all_sites(euk.din.only)
con_smol <- run_coniss_all_sites(euk.smolphy.only)

zone_years <- function(con_list, k_vec){
  setNames(lapply(names(con_list),
                  function(s) coniss_boundaries(con_list[[s]], k_vec[s])),
           names(con_list))
}

zones_all  <- zone_years(con_all,  k_all)
zones_met  <- zone_years(con_met,  k_met)
zones_dia  <- zone_years(con_dia,  k_dia)
zones_din  <- zone_years(con_din,  k_din)
zones_smol <- zone_years(con_smol, k_smol)


####----5. Pre-compute site x year matrices and shared ylim per dataset----####

build_site_mats <- function(M_full){
  setNames(
    lapply(sites, make_site_year, M = M_full,
           map = expandedDates, years = all_years),
    sites
  )
}

# largest stacked-bar total observed across all sites for this dataset
shared_ylim_max <- function(site_mats){
  vals <- unlist(lapply(site_mats, function(mat) {
    mat[is.na(mat)] <- 0
    colSums(mat)
  }))
  max(vals, na.rm = TRUE)
}

site_mats_all  <- build_site_mats(EUK.tax.a.2[nrow(EUK.tax.a.2):1, ])
site_mats_met  <- build_site_mats(EUK.met.2)
site_mats_dia  <- build_site_mats(EUK.dia.3)
site_mats_din  <- build_site_mats(EUK.din.2)
site_mats_smol <- build_site_mats(EUK.smolphy.2)

ymax_all  <- shared_ylim_max(site_mats_all)    # ~1 (proportions sum to 1)
ymax_met  <- shared_ylim_max(site_mats_met)
ymax_dia  <- shared_ylim_max(site_mats_dia)
ymax_din  <- shared_ylim_max(site_mats_din)
ymax_smol <- shared_ylim_max(site_mats_smol)


####----6. Core plotting function----####
# Same y-axis (0, ymax) across all 5 site panels of a dataset.
# Legend placed on the LAST (bottom) panel using user-supplied coords,
# matching the carefully placed legends from the original plots.

big_plot_with_zones <- function(site_mats,
                                zones_by_site,
                                pal,
                                ymax,
                                pdf_path,
                                legend_x, legend_y,
                                legend_labels,
                                legend_fill,
                                legend_cex   = 0.6,
                                ribbon_alpha = 0.70,
                                width = 14, height = 10){
  pdf(pdf_path, width = width, height = height)
  par(mfrow = c(length(sites), 1),
      mar = c(0.5, 4.1, 1.1, 1.1),
      oma = c(6, 0, 0, 0),
      xpd = NA)
  
  ylim <- c(0, ymax)
  
  for (idx in seq_along(site_order)) {
    i    <- site_order[idx]
    site <- sites[i]
    mat  <- site_mats[[site]]
    
    bp <- barplot(mat,
                  col    = pal,
                  space  = 0,
                  width  = 30,
                  ylim   = ylim,
                  border = NA,
                  xaxt   = "n",
                  yaxt   = "t")
    connect_observed(M = mat, bp = bp, cols = pal, alpha = ribbon_alpha)
    barplot(mat,
            col    = pal,
            space  = 0,
            width  = 30,
            ylim   = ylim,
            border = NA,
            xaxt   = "n",
            yaxt   = "t",
            add    = TRUE)
    
    draw_zone_lines(zones_by_site[[site]], bp, all_years)
    
    mtext(site, side = 2, line = 3, las = 0, cex = 0.9)
    
    if (idx == length(site_order)) {
      axis(1, at = bp, labels = all_years, las = 2, cex.axis = 0.8)
      mtext("Calibrated year", side = 1, line = 4)
      
      # Legend on the bottom panel; xpd = NA lets it draw outside the panel
      legend(legend_x, legend_y,
             legend = legend_labels,
             fill   = legend_fill,
             border = NA,
             cex    = legend_cex,
             bty    = "n",
             inset  = c(-0.01, 0),
             xpd    = NA)
    }
  }
  dev.off()
}


####----7. Produce all five zoned big plots, with original legend placement----####

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# 7a. Whole-community Subdivision plot
#     Original (EUK.bigplot.pdf): legend(0, 3, ..., fill = rev(pal), cex = 0.6)
pal_all <- getPalette(nrow(EUK.tax.a.2))
big_plot_with_zones(
  site_mats     = site_mats_all,
  zones_by_site = zones_all,
  pal           = pal_all,
  ymax          = ymax_all,
  pdf_path      = "figures/EUK.bigplot.zoned.pdf",
  legend_x      = 0, legend_y = 3,
  legend_labels = rownames(EUK.tax.a.2),
  legend_fill   = rev(pal_all),
  legend_cex    = 0.6,
  ribbon_alpha  = 0.70
)

# 7b. Metazoa
#     Original (EUK.large.met.pdf): legend(0, 0.4, ..., fill = pal, cex = 0.6)
pal_met <- rev(getPalette(nrow(EUK.met.2)))
big_plot_with_zones(
  site_mats     = site_mats_met,
  zones_by_site = zones_met,
  pal           = pal_met,
  ymax          = ymax_met,
  pdf_path      = "figures/EUK.large.met.zoned.pdf",
  legend_x      = 0, legend_y = 0.4,
  legend_labels = rownames(EUK.met.2),
  legend_fill   = pal_met,
  legend_cex    = 0.6,
  ribbon_alpha  = 0.70
)

# 7c. Diatoms
#     Original (EUK.large.dia.pdf): legend(0, 0.4, ..., fill = pal, cex = 0.6)
pal_dia <- rev(getPalette(nrow(EUK.dia.3)))
big_plot_with_zones(
  site_mats     = site_mats_dia,
  zones_by_site = zones_dia,
  pal           = pal_dia,
  ymax          = ymax_dia,
  pdf_path      = "figures/EUK.large.dia.zoned.pdf",
  legend_x      = 0, legend_y = 0.4,
  legend_labels = rownames(EUK.dia.3),
  legend_fill   = pal_dia,
  legend_cex    = 0.6,
  ribbon_alpha  = 0.35
)

# 7d. Dinoflagellates
#     Original (EUK.large.din.pdf): legend(0, 0.12, ..., fill = pal, cex = 0.6)
pal_din <- rev(getPalette(nrow(EUK.din.2)))
big_plot_with_zones(
  site_mats     = site_mats_din,
  zones_by_site = zones_din,
  pal           = pal_din,
  ymax          = ymax_din,
  pdf_path      = "figures/EUK.large.din.zoned.pdf",
  legend_x      = 0, legend_y = 0.12,
  legend_labels = rownames(EUK.din.2),
  legend_fill   = pal_din,
  legend_cex    = 0.6,
  ribbon_alpha  = 0.35
)

# 7e. Small phytoplankton
#     Original (EUK.large.smol.pdf): legend(0, 0.04, ..., fill = pal, cex = 0.6)
pal_smol <- rev(getPalette(nrow(EUK.smolphy.2)))
big_plot_with_zones(
  site_mats     = site_mats_smol,
  zones_by_site = zones_smol,
  pal           = pal_smol,
  ymax          = ymax_smol,
  pdf_path      = "figures/EUK.large.smol.zoned.pdf",
  legend_x      = 0, legend_y = 0.2,
  legend_labels = rownames(EUK.smolphy.2),
  legend_fill   = pal_smol,
  legend_cex    = 0.6,
  ribbon_alpha  = 0.35
)

####----8. PDF → PNG via sips (one line each)----####

library(parallel)

pdfs <- c("figures/EUK.bigplot.zoned.pdf",
          "figures/EUK.large.met.zoned.pdf",
          "figures/EUK.large.dia.zoned.pdf",
          "figures/EUK.large.din.zoned.pdf",
          "figures/EUK.large.smol.zoned.pdf")

invisible(mclapply(pdfs, function(f) {
  system2("sips", c("-s","format","png","-z",6000,8400, f, "--out", sub("\\.pdf$",".png",f)))
}, mc.cores = 5))




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

alpha_all <- bind_rows(out_list, .id = "draw")

alpha_avg <- alpha_all %>% 
  group_by(sample) %>% 
  summarise(across(q0:q2,
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

alpha_avg_reps <- alpha_avg %>% 
  mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
  group_by(replicate) %>% 
  summarise(across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
            across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))))

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

ASV.alpha_avg_reps <- alpha_avg_reps

#### Now we loop over subsets 

asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

# Rarefaction + Hill function definitions stay unchanged
hill_3 <- function(v){
  v <- v[v > 0]
  if(length(v) == 0 || sum(v) == 0) return(c(q0 = NA, q1 = NA, q2 = NA))
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
                     list(mean = \(x) mean(x, na.rm = TRUE),
                          sd   = \(x) sd(x,   na.rm = TRUE)),
                     .names = "{.col}_{.fn}"))
  
  alpha_avg_reps <- alpha_avg %>% 
    mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
    group_by(replicate) %>% 
    summarise(across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
              across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))))
  
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

n_reads <- 30000
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
  summarise(across(q0:q2,
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

alpha_genus_avg_reps <- alpha_genus_avg %>% 
  mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>% 
  group_by(replicate) %>% 
  summarise(across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
            across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))))

pdf("figures/EUK_genus.Q0.pdf", height = 4.5, width = 12)
par(mar = c(5.1, 2.1, 4.1, 4.1))
plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)],
     alpha_genus_avg$q0_mean,
     col=add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "",alpha_genus_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
     pch=16,
     bty = "n",
 xlab="",
     ylab="",
     yaxt="n",
     xaxt="n")
axis(1,at=seq(1775,2000,by=25),labels=seq(1775,2000,by=25),cex.axis=1.5)
axis(4,at=c(50,100,150,200),labels=c(50,100,150,200))
mtext("Genus Richness", side = 4, line = 3)
points(expandedDates$CalYear[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)],
       alpha_genus_avg_reps$q0_mean,
       col=coredata$col[match(toupper(expandedDates$Site[match(alpha_genus_avg_reps$replicate, expandedDates$SampleName)]),toupper(coredata$Site))],
       pch=16, cex=1.2)
dev.off()

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
                     list(mean = \(x) mean(x, na.rm = TRUE),
                          sd   = \(x) sd(x,   na.rm = TRUE)),
                     .names = "{.col}_{.fn}"))
  
  alpha_avg_reps <- alpha_avg %>%
    mutate(replicate = sub("(.*)\\.[0-9]+$", "\\1", sample)) %>%
    group_by(replicate) %>%
    summarise(across(ends_with("_mean"), \(x) mean(x, na.rm = TRUE)),
              across(ends_with("_sd"),   \(x) sqrt(mean(x^2, na.rm = TRUE))))
  
  
  
  ## Plot Q0
  pdf(paste0("figures/", group_name, ".genus.Q0.pdf"),  height = 4.5, width = 12)
  plot(expandedDates$CalYear[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)],
       alpha_avg$q0_mean,
       col = add.alpha(coredata$col[match(toupper(expandedDates$Site[match(gsub("\\.[0-9]+$", "", alpha_avg$sample), expandedDates$SampleName)]),toupper(coredata$Site))], 0.2),
       pch = 16,
       bty="n",
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
                    "ASV_Q0" = ASV.alpha_avg_reps$q0_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q1" = ASV.alpha_avg_reps$q1_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q2" = ASV.alpha_avg_reps$q2_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "genus_Q0" = alpha_genus_avg_reps$q0_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q1" = alpha_genus_avg_reps$q1_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q2" = alpha_genus_avg_reps$q2_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "whale" = whales$biomass[match(round(clim$CalYear),whales$year)],
                    "whale_s" = scale(whales$biomass[match(round(clim$CalYear),whales$year)]),
                    "SAM" = clim$SAM_Index_3yr,
                    "ENSO" = clim$ENSO_3yr,
                    "SeaIce" = clim$SIC_mean,
                    "SAM_s" = scale(clim$SAM_Index_3yr),
                    "ENSO_s" = scale(clim$ENSO_3yr),
                    "SeaIce_s" = scale(clim$SIC_mean)
                    
)

## Lets make a little plot of analysed values 

site_cols <- coredata$col[match(toupper(alpha$site), toupper(coredata$Site))]
ylim <- range(c(alpha$whale_s, alpha$SAM_s), na.rm = TRUE)


pdf("figures/climate/ScaledWahles&Clim.pdf",width = 11,height = 5.5)
plot(alpha$year, alpha$whale_s,
     bg = site_cols, col = "black", pch = 24, cex = 1.2,
     ylim = ylim,
     xlab = "", ylab = "Scaled value",
     main = "",
     bty = "n", xaxt = "n")
points(alpha$year, alpha$SAM_s,
       bg = site_cols, col = "black", pch = 25, cex = 1.2)
axis(3, at = axTicks(1), labels = axTicks(1))
dev.off()


library(lme4)
library(lmerTest)

# Fit models
a2_A0 <- lmer(ASV_Q0   ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)
a2_A1 <- lmer(ASV_Q1   ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)
a2_A2 <- lmer(ASV_Q2   ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)

a2_G0 <- lmer(genus_Q0 ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)
a2_G1 <- lmer(genus_Q1 ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)
a2_G2 <- lmer(genus_Q2 ~ SAM_s + ENSO_s + SeaIce_s + whale_s + (1|site), data = alpha, REML = TRUE)

mods <- list(
  "ASV Q0"   = a2_A0, "ASV Q1"   = a2_A1, "ASV Q2"   = a2_A2,
  "Genus Q0" = a2_G0, "Genus Q1" = a2_G1, "Genus Q2" = a2_G2
)

predictors <- c("SAM_s", "ENSO_s", "SeaIce_s", "whale_s")

# Check singular fits
singular_status <- sapply(mods, isSingular)
print(singular_status)

# Variance components
lapply(mods, VarCorr)

# Summaries
lapply(mods, summary)

# ---- Diagnostics ----
for (m in names(mods)) {
  mod <- mods[[m]]
  par(mfrow = c(2, 2))
  plot(mod, main = m)
  qqnorm(resid(mod), main = paste("QQ:", m))
  qqline(resid(mod))
  hist(resid(mod), breaks = 20, main = paste("Residuals:", m))
  plot(fitted(mod), resid(mod), main = paste("Fitted vs Resid:", m))
  abline(h = 0, lty = 2)
}
par(mfrow = c(1, 1))

# ---- Extract coefficients safely ----
coef_df <- do.call(rbind, lapply(names(mods), function(m) {
  mod <- mods[[m]]
  sm  <- coef(summary(mod))
  ci  <- confint(mod, method = "Wald")
  
  out <- lapply(predictors, function(p) {
    if (!p %in% rownames(sm)) return(NULL)
    data.frame(
      model = m,
      predictor = p,
      estimate = sm[p, "Estimate"],
      se = sm[p, "Std. Error"],
      lo = ci[p, 1],
      hi = ci[p, 2],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}))

# Response SD from actual model frame
resp_sd <- sapply(names(mods), function(m) {
  mf <- model.frame(mods[[m]])
  sd(model.response(mf), na.rm = TRUE)
})

coef_df$est_std <- coef_df$estimate / resp_sd[coef_df$model]
coef_df$lo_std  <- coef_df$lo / resp_sd[coef_df$model]
coef_df$hi_std  <- coef_df$hi / resp_sd[coef_df$model]

# ---- Coefficient plot ----
xpos <- c("ASV Q0" = 1, "ASV Q1" = 2, "ASV Q2" = 3,
          "Genus Q0" = 5, "Genus Q1" = 6, "Genus Q2" = 7)

pred_cols <- c("SAM_s" = "#141d23", "ENSO_s" = "#31424e", "SeaIce_s" = "#5a6c76", "whale_s" = "#7d9bb0")

offsets   <- c("SAM_s" = -0.18, "ENSO_s" = -0.06, "SeaIce_s" = 0.06, "whale_s" = 0.18)

ylim <- range(c(coef_df$lo_std, coef_df$hi_std), na.rm = TRUE)


pdf("figures/stat.alpha.coeff.pdf",width = 5.5,height = 9)
par(mar = c(5, 5, 3, 1))

plot(NULL, xlim = c(0.4, 7.6), ylim = ylim,
     xaxt = "n", bty = "n",
     xlab = "", ylab = "Response-standardised coefficient (± 95% CI)",
     main = "LMM fixed-effect coefficients")
abline(h = 0, lty = 2, col = "grey40")
abline(v = 4, lty = 3, col = "grey70")

for (p in predictors) {
  sub <- coef_df[coef_df$predictor == p, ]
  x <- xpos[sub$model] + offsets[p]
  segments(x, sub$lo_std, x, sub$hi_std, col = pred_cols[p], lwd = 2)
  points(x, sub$est_std, col = pred_cols[p], pch = 16, cex = 1.4)
}

axis(1, at = c(1, 2, 3, 5, 6, 7),
     labels = c("Q0", "Q1", "Q2", "Q0", "Q1", "Q2"),
     tick = FALSE, line = 0)

mtext("ASV", side = 1, at = 2, line = 2.5, cex = 1.1, font = 2)
mtext("Genus", side = 1, at = 6, line = 2.5, cex = 1.1, font = 2)

legend("topleft",
       legend = c("SAM", "ENSO", "Sea ice", "Whale biomass"),
       col = pred_cols, pch = 16, lwd = 2,
       cex = 0.8, bty = "n")


dev.off()

####====Export LMM results to text file====####

mods <- list(
  "ASV Q0"   = a2_A0, "ASV Q1"   = a2_A1, "ASV Q2"   = a2_A2,
  "Genus Q0" = a2_G0, "Genus Q1" = a2_G1, "Genus Q2" = a2_G2
)

predictors <- c("SAM_s", "ENSO_s", "SeaIce_s", "whale_s")

sink("summaryfiles/LMM_alpha_results.txt")

cat("======================================================================\n")
cat("Linear Mixed Model results: alpha diversity ~ SAM_s + ENSO_s + SeaIce_s + whale_s\n")
cat("Random effect: (1 | site)\n")
cat("Fit by REML using lme4 / lmerTest\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("======================================================================\n\n")

# Sample sizes
cat("---- Sample sizes ----\n")
for (m in names(mods)) {
  n_obs <- nobs(mods[[m]])
  n_grp <- ngrps(mods[[m]])
  cat(sprintf("  %-10s  n obs = %d   n sites = %d\n", m, n_obs, n_grp))
}
cat("\n")

# Singular fits
cat("---- Singular fit check ----\n")
sing <- sapply(mods, isSingular)
for (m in names(sing)) cat(sprintf("  %-10s  isSingular = %s\n", m, sing[m]))
cat("\n")

# Per-model detail
for (m in names(mods)) {
  mod <- mods[[m]]
  cat("======================================================================\n")
  cat("Model:", m, "\n")
  cat("Formula:", deparse(formula(mod)), "\n")
  cat("----------------------------------------------------------------------\n")
  
  # Variance components
  cat("Random effects (variance components):\n")
  vc <- as.data.frame(VarCorr(mod))
  print(vc[, c("grp", "var1", "vcov", "sdcor")], row.names = FALSE)
  cat("\n")
  
  # Fixed effects with df + p (lmerTest summary)
  cat("Fixed effects (Satterthwaite df, lmerTest):\n")
  sm <- coef(summary(mod))
  print(round(sm, 4))
  cat("\n")
  
  # Wald 95% CIs for fixed effects
  ci <- confint(mod, method = "Wald")
  ci_fx <- ci[rownames(ci) %in% rownames(sm), , drop = FALSE]
  cat("95% Wald confidence intervals (fixed effects):\n")
  print(round(ci_fx, 4))
  cat("\n")
  
  # Marginal / conditional R^2 if MuMIn is around
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2 <- tryCatch(MuMIn::r.squaredGLMM(mod), error = function(e) NULL)
    if (!is.null(r2)) {
      cat("R^2 (MuMIn):\n")
      cat(sprintf("  Marginal R2  = %.4f\n", r2[1, "R2m"]))
      cat(sprintf("  Conditional R2 = %.4f\n", r2[1, "R2c"]))
      cat("\n")
    }
  }
  
  # Response SD (for context)
  mf <- model.frame(mod)
  cat(sprintf("Response SD (for scaling coefficients): %.4f\n",
              sd(model.response(mf), na.rm = TRUE)))
  cat(sprintf("Residual SD: %.4f\n",
              sigma(mod)))
  cat("\n")
}

# Compact summary table across all models
cat("======================================================================\n")
cat("Compact summary table: estimate (SE), 95% CI, p-value\n")
cat("======================================================================\n\n")

compact <- do.call(rbind, lapply(names(mods), function(m) {
  mod <- mods[[m]]
  sm  <- coef(summary(mod))
  ci  <- confint(mod, method = "Wald")
  do.call(rbind, lapply(predictors, function(p) {
    if (!p %in% rownames(sm)) return(NULL)
    data.frame(
      model     = m,
      predictor = p,
      estimate  = unname(sm[p, "Estimate"]),
      SE        = unname(sm[p, "Std. Error"]),
      df        = unname(sm[p, "df"]),
      t         = unname(sm[p, "t value"]),
      CI_lo     = unname(ci[p, 1]),
      CI_hi     = unname(ci[p, 2]),
      p_value   = unname(sm[p, "Pr(>|t|)"]),
      sig       = ifelse(sm[p, "Pr(>|t|)"] < 0.001, "***",
                         ifelse(sm[p, "Pr(>|t|)"] < 0.01,  "**",
                                ifelse(sm[p, "Pr(>|t|)"] < 0.05,  "*",
                                       ifelse(sm[p, "Pr(>|t|)"] < 0.1,   ".", "")))),
      stringsAsFactors = FALSE
    )
  }))
}))

compact[, c("estimate","SE","df","t","CI_lo","CI_hi","p_value")] <-
  round(compact[, c("estimate","SE","df","t","CI_lo","CI_hi","p_value")], 4)

print(compact, row.names = FALSE)
cat("\nSignificance codes:  *** < 0.001   ** < 0.01   * < 0.05   . < 0.1\n")

sink()

# Also write the compact table as CSV for convenience
write.csv(compact, "summaryfiles/LMM_alpha_results_compact.csv", row.names = FALSE)

cat("Wrote: figures/LMM_alpha_results.txt\n")
cat("Wrote: figures/LMM_alpha_results_compact.csv\n")




## beta diversity 



# first a little test with UMAP
#test <- umap(t(prop.table(as.matrix(euk.reps),2)))
#plot(test$layout,pch=16, col=coredata$col[match(toupper(expandedDates$Site[match(row.names(test$layout), expandedDates$SampleName)]),toupper(coredata$Site))])
#text(test$layout, labels = expandedDates$CalYear[match(row.names(test$layout), expandedDates$SampleName)],cex=0.7)

# now the NMDS
EUK.nMDS.j <- metaMDS(vegdist(t(prop.table(as.matrix(euk.reps),2)),binary = TRUE,method = "jaccard"),trymax = 200,parallel=7)
EUK.nMDS.b <- metaMDS(vegdist(t(prop.table(as.matrix(euk.reps),2)),method = "bray"),trymax = 200,parallel=7)

pdf(paste0("figures/NMDS_ALL_jaccard.pdf"), height = 7, width = 7.5)
plot(EUK.nMDS.j$points[,1], 
     EUK.nMDS.j$points[,2],
     col=coredata$col[match(toupper(expandedDates$Site[match(row.names(EUK.nMDS.j$points), expandedDates$SampleName)]),toupper(coredata$Site))],
     pch=16,ylab="",xlab="")

text(EUK.nMDS.j$points[,1],
     EUK.nMDS.j$points[,2],
     labels = expandedDates$CalYear[match(rownames(EUK.nMDS.j$points), expandedDates$SampleName)],
     pos = 3, cex = 0.6)

dev.off()


pdf(paste0("figures/NMDS_ALL_bray.pdf"), height = 7, width = 7.5)
plot(EUK.nMDS.b$points[,1], 
     EUK.nMDS.b$points[,2],
     col=coredata$col[match(toupper(expandedDates$Site[match(row.names(EUK.nMDS.b$points), expandedDates$SampleName)]),toupper(coredata$Site))],
     pch=16,ylab="",xlab="")
text(EUK.nMDS.b$points[,1],
     EUK.nMDS.b$points[,2],
     labels = expandedDates$CalYear[match(rownames(EUK.nMDS.j$points), expandedDates$SampleName)],
     pos = 3, cex = 0.6)
dev.off()


# A new version with connecting lines
# --- helper: mix a color toward white or black ---
shade_col <- function(col, amount) {
  # amount in [-1, 1]: positive = toward white, negative = toward black
  rgb_vals <- col2rgb(col) / 255
  if (amount >= 0) {
    mixed <- rgb_vals + (1 - rgb_vals) * amount
  } else {
    mixed <- rgb_vals * (1 + amount)
  }
  rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

# --- helper: draw trajectories on an existing plot ---
draw_trajectories <- function(points, sites, years, site_col,
                              smooth = FALSE, spar = 0.7,
                              year_step = 25,
                              base_lighten = 0.5,   # how much to lighten site col -> trajectory midtone
                              shade_offset = 0.18) { # +/- around midtone for alternating bands
  
  for (s in unique(sites)) {
    idx <- which(sites == s)
    ord <- idx[order(years[idx])]
    
    # build the two alternating shades from the site color
    midtone     <- shade_col(site_col[idx[1]], base_lighten)
    shade_light <- shade_col(midtone,  shade_offset)
    shade_dark  <- shade_col(midtone, -shade_offset)
    
    if (smooth && length(ord) >= 4) {
      x <- points[ord, 1]; y <- points[ord, 2]; yr <- years[ord]
      t <- c(0, cumsum(sqrt(diff(x)^2 + diff(y)^2)))
      tt <- seq(min(t), max(t), length.out = 400)
      sx <- predict(smooth.spline(t, x, spar = spar), tt)$y
      sy <- predict(smooth.spline(t, y, spar = spar), tt)$y
      
      # year at every point along tt: invert (yr -> t) mapping, then back to year via tt
      tt_to_yr <- approx(t, yr, xout = tt, rule = 2)$y
      
      # alternation index: floor((year - oldest) / year_step) %% 2
      band_idx <- floor((tt_to_yr - min(yr)) / year_step) %% 2
      seg_cols <- ifelse(band_idx[-length(band_idx)] == 0, shade_dark, shade_light)
      
      segments(sx[-length(sx)], sy[-length(sy)],
               sx[-1],           sy[-1],
               col = seg_cols, lwd = 2)
      
      # arrowhead in whichever shade the final segment uses
      final_col <- if (tail(band_idx, 1) == 0) shade_dark else shade_light
      arrows(sx[length(sx) - 1], sy[length(sy) - 1],
             sx[length(sx)],     sy[length(sy)],
             col = final_col, length = 0.08, lwd = 2)
      
    } else if (length(ord) >= 2) {
      x <- points[ord, 1]; y <- points[ord, 2]; yr <- years[ord]
      # alternation by year of the starting point of each segment
      band_idx <- floor((yr[-length(yr)] - min(yr)) / year_step) %% 2
      seg_cols <- ifelse(band_idx == 0, shade_dark, shade_light)
      for (i in seq_len(length(ord) - 1)) {
        arrows(x[i], y[i], x[i+1], y[i+1],
               col = seg_cols[i], length = 0.08, lwd = 2)
      }
    }
  }
  
  # redraw points on top in their original (unlightened) site color
  points(points[, 1], points[, 2], col = site_col, pch = 16, cex = 0.8)
  
  # label first and last year per site
  for (s in unique(sites)) {
    idx <- which(sites == s)
    ord <- idx[order(years[idx])]
    endpts <- ord[c(1, length(ord))]
    text(points[endpts, 1], points[endpts, 2],
         labels = years[endpts], pos = 3, cex = 0.6,col=adjustcolor(coredata$col[match(toupper(s), toupper(coredata$Site))],red.f = 0.5, green.f = 0.5, blue.f = 0.5))
  }
}

# --- helper: full trajectory plot to PDF ---
trajectory_pdf <- function(file, points, sites, years, site_col,
                           main, smooth = FALSE, year_step = 20,
                           base_lighten = 0.5, shade_offset = 0.3) {
  pdf(file, height = 5.5, width = 5)
  plot(points[,1], points[,2],
       col = site_col, pch = 16, cex=0.8,
       main = main, xlab = "", ylab = "")
  draw_trajectories(points, sites, years, site_col,
                    smooth = smooth, year_step = year_step,
                    base_lighten = base_lighten, shade_offset = shade_offset)
  dev.off()
}

# ======================================================
# ALL eukaryotes
# ======================================================
nmds_names_j <- rownames(EUK.nMDS.j$points)
nmds_sites_j <- expandedDates$Site[match(nmds_names_j, expandedDates$SampleName)]
nmds_years_j <- expandedDates$CalYear[match(nmds_names_j, expandedDates$SampleName)]
nmds_cols_j  <- coredata$col[match(toupper(nmds_sites_j), toupper(coredata$Site))]

nmds_names_b <- rownames(EUK.nMDS.b$points)
nmds_sites_b <- expandedDates$Site[match(nmds_names_b, expandedDates$SampleName)]
nmds_years_b <- expandedDates$CalYear[match(nmds_names_b, expandedDates$SampleName)]
nmds_cols_b  <- coredata$col[match(toupper(nmds_sites_b), toupper(coredata$Site))]

# Jaccard – straight + wiggly
trajectory_pdf("figures/NMDS_ALL_jaccard_trajectories.pdf",
               EUK.nMDS.j$points, nmds_sites_j, nmds_years_j, nmds_cols_j,
               main = "",
               smooth = FALSE)

trajectory_pdf("figures/NMDS_ALL_jaccard_wiggly_trajectories.pdf",
               EUK.nMDS.j$points, nmds_sites_j, nmds_years_j, nmds_cols_j,
               main = "",
               smooth = TRUE)

# Bray – straight + wiggly
trajectory_pdf("figures/NMDS_ALL_bray_trajectories.pdf",
               EUK.nMDS.b$points, nmds_sites_b, nmds_years_b, nmds_cols_b,
               main = "",
               smooth = FALSE)

trajectory_pdf("figures/NMDS_ALL_bray_wiggly_trajectories.pdf",
               EUK.nMDS.b$points, nmds_sites_b, nmds_years_b, nmds_cols_b,
               main = "",
               smooth = TRUE)


# ======================================================
# By group
# ======================================================
asv_groups <- list(
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

for (group_name in names(asv_groups)) {
  
  group_asvs <- asv_groups[[group_name]]
  mat <- euk.reps[rownames(euk.reps) %in% group_asvs, , drop = FALSE]
  mat_rel <- prop.table(as.matrix(mat), 2)
  
  set.seed(42)
  nmds_j <- metaMDS(vegdist(t(mat_rel), binary = TRUE, method = "jaccard"),
                    trymax = 200, parallel = 7)
  set.seed(42)
  nmds_b <- metaMDS(vegdist(t(mat_rel), method = "bray"),
                    trymax = 200, parallel = 7)
  
  sample_names <- rownames(nmds_j$points)
  match_idx <- match(sample_names, expandedDates$SampleName)
  sites <- expandedDates$Site[match_idx]
  years <- expandedDates$CalYear[match_idx]
  site_col <- coredata$col[match(toupper(sites), toupper(coredata$Site))]
  
  # Plain plots (no trajectories)
  pdf(paste0("figures/NMDS_", group_name, "_jaccard.pdf"), height = 7, width = 7.5)
  plot(nmds_j$points[,1], nmds_j$points[,2],
       col = site_col, pch = 16,
       main = paste(group_name, "– Jaccard"))
  text(nmds_j$points[,1], nmds_j$points[,2],
       labels = years, pos = 3, cex = 0.6)
  dev.off()
  
  pdf(paste0("figures/NMDS_", group_name, "_bray.pdf"), height = 7, width = 7.5)
  plot(nmds_b$points[,1], nmds_b$points[,2],
       col = site_col, pch = 16,
       main = paste(group_name, "– Bray-Curtis"))
  text(nmds_b$points[,1], nmds_b$points[,2],
       labels = years, pos = 3, cex = 0.6)
  dev.off()
  
  # Jaccard trajectories – straight + wiggly
  trajectory_pdf(paste0("figures/NMDS_", group_name, "_jaccard_trajectories.pdf"),
                 nmds_j$points, sites, years, site_col,
                 main = paste(group_name, "– Jaccard (trajectories)"),
                 smooth = FALSE)
  
  trajectory_pdf(paste0("figures/NMDS_", group_name, "_jaccard_wiggly_trajectories.pdf"),
                 nmds_j$points, sites, years, site_col,
                 main = paste(group_name, "– Jaccard (smoothed trajectories)"),
                 smooth = TRUE)
  
  # Bray trajectories – straight + wiggly
  trajectory_pdf(paste0("figures/NMDS_", group_name, "_bray_trajectories.pdf"),
                 nmds_b$points, sites, years, site_col,
                 main = paste(group_name, "– Bray-Curtis (trajectories)"),
                 smooth = FALSE)
  
  trajectory_pdf(paste0("figures/NMDS_", group_name, "_bray_wiggly_trajectories.pdf"),
                 nmds_b$points, sites, years, site_col,
                 main = paste(group_name, "– Bray-Curtis (smoothed trajectories)"),
                 smooth = TRUE)
}



#### now let's run some tests 
beta <- data.frame("sample" = paste(clim$Site,clim$SampleID),
                   "site" = as.factor(clim$Site),
                   "year" =clim$CalYear,
                   "whale" = whales$biomass[match(round(clim$CalYear),whales$year)],
                   "whale_s" = scale(whales$biomass[match(round(clim$CalYear),whales$year)]),
                   "SAM" = clim$SAM_Index_3yr,
                   "ENSO" = clim$ENSO_3yr,
                   "SeaIce" = clim$SIC_mean,
                   "SAM_s" = scale(clim$SAM_Index_3yr),
                   "ENSO_s" = scale(clim$ENSO_3yr),
                   "SeaIce_s" = scale(clim$SIC_mean))



## distance based rDA

library(vegan)

comm.raw <- prop.table(as.matrix(euk.reps),2)

comm <- t(comm.raw[,match(clim$longname,colnames(comm.raw))])


## distance basd redundnacy analysis 
mod_bray <- dbrda(
  vegdist(comm, "bray") ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
  data = beta
)

mod_jacc <- dbrda(
  vegdist(comm, "jacc",binary = TRUE) ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
  data = beta
)




plot_dbrda <- function(mod, title, trajectories = FALSE) {
  
  sc_sites <- scores(mod, display = "sites", choices = 1:2)
  sc_bp <- scores(mod, display = "bp", choices = 1:2)
  
  # beta rows are in the same order as sc_sites, so use position directly
  sample_sites <- beta$site
  site_cols <- coredata$col[match(toupper(sample_sites), toupper(coredata$Site))]
  
  eig <- eigenvals(mod)
  eig_pct <- 100 * eig / sum(eig)
  xlab <- paste0("dbRDA1 (", round(eig_pct[1], 1), "%)")
  ylab <- paste0("dbRDA2 (", round(eig_pct[2], 1), "%)")
  
  arrow_scale <- 0.5 * max(abs(sc_sites)) / max(abs(sc_bp))
  sc_bp_scaled <- sc_bp * arrow_scale
  
  plot(sc_sites[,1], sc_sites[,2],
       col = site_cols, pch = 16, cex = 1.2,
       xlab = xlab, ylab = ylab,
       main = title)
  abline(h = 0, v = 0, lty = 2, col = "grey60")
  
  if (trajectories) {
    years <- beta$year
    for (s in unique(sample_sites)) {
      idx <- which(sample_sites == s)
      ord <- idx[order(years[idx])]
      if (length(ord) >= 3) {
        x <- sc_sites[ord, 1]
        y <- sc_sites[ord, 2]
        t <- c(0, cumsum(sqrt(diff(x)^2 + diff(y)^2)))
        tt <- seq(min(t), max(t), length.out = 100)
        sx <- predict(smooth.spline(t, x, spar = 0.5), tt)$y
        sy <- predict(smooth.spline(t, y, spar = 0.5), tt)$y
        lines(sx, sy, col = add.alpha(site_cols[idx[1]],alpha = 0.4), lwd = 1.5)
        arrows(sx[length(sx) - 1], sy[length(sy) - 1],
               sx[length(sx)],     sy[length(sy)],
               col = add.alpha(site_cols[idx[1]],alpha = 0.4), length = 0.08, lwd = 1.5)
      }
    }
    
    points(sc_sites[,1], sc_sites[,2],
           col = site_cols, pch = 16, cex = 1.2)
    
    # year labels (unchanged)
    for (s in unique(sample_sites)) {
      idx <- which(sample_sites == s)
      ord <- idx[order(years[idx])]
      thin <- ord[seq(1, length(ord), by = 5)]
      thin <- union(thin, ord[c(1, length(ord))])
      text(sc_sites[thin, 1], sc_sites[thin, 2],
           labels = years[thin], pos = 3, cex = 0.6)
    }
  }
  
  # match arrow colours to predictor names via the biplot rownames
  arrow_cols <- pred_cols[rownames(sc_bp_scaled)]
  
  arrows(0, 0, sc_bp_scaled[,1], sc_bp_scaled[,2],
         col = arrow_cols, lwd = 2, length = 0.1)
  
  #legend("topleft",
   #      legend = gsub("whale","WhaleBiomass",gsub("_s", "", names(arrow_cols))),
    #     col    = arrow_cols,
     #    lwd    = 2, cex = 0.85, bty = "n")
  
  u_sites <- unique(sample_sites)
  u_cols <- coredata$col[match(toupper(u_sites), toupper(coredata$Site))]
  # legend("topright", legend = u_sites, col = u_cols,
  #      pch = 16, cex = 0.75, bty = "n")
}
pdf("figures/stat.beta.dbRDA.b.pdf",height = 8,width = 8)
plot_dbrda(mod_bray, "",trajectories =FALSE )
dev.off()
pdf("figures/stat.beta.dbRDA.j.pdf",height = 8,width = 8)
plot_dbrda(mod_jacc, "",trajectories =FALSE )
dev.off()




##  anova with perms within cores
perm <- how(blocks = alpha$site,nperm = 100000)

anova(mod_bray, permutations = perm)
anova(mod_bray, by = "margin", permutations = perm)

anova(mod_jacc, permutations = perm)
anova(mod_jacc, by = "margin", permutations = perm)

## variance partititoning 

var_bray <- varpart(vegdist(comm, "bray"),
                    ~ SAM_s + ENSO_s + SeaIce_s,
                    ~ whale_s,
                    ~ site,
                    data = beta)
var_jacc <- varpart(vegdist(comm, "jacc",binary = TRUE),
                    ~ SAM_s + ENSO_s + SeaIce_s,
                    ~ whale_s,
                    ~ site,
                    data = beta)

var_bray
plot(var_bray)
var_jacc
plot(var_jacc)



## PERMANOVA
perm_bray <- adonis2(vegdist(comm, method = "bray") ~ SAM_s + ENSO_s + SeaIce_s + whale_s,
                     data = beta,
                     permutations = how(blocks = beta$site, nperm = 100000),
                     by = "margin")
perm_bray

perm_jacc <- adonis2(vegdist(comm, method = "jacc",binary = TRUE) ~ SAM_s + ENSO_s + SeaIce_s + whale_s,
                     data = beta,
                     permutations = how(blocks = beta$site, nperm = 100000),
                     by = "margin")

perm_jacc




## write results to file
sink("summaryfiles/Betaresults.txt")

####---- Global analysis (all eukaryotes) ----####

cat("############################################################\n")
cat("Group: All Taxa (global)\n")
cat("############################################################\n\n")

comm.raw <- prop.table(as.matrix(euk.reps), 2)
comm     <- t(comm.raw[, match(clim$longname, colnames(comm.raw))])

perm <- how(blocks = beta$site, nperm = 100000)

cat("---- dbRDA: Bray ----\n")
mod_bray_glob <- dbrda(vegdist(comm, "bray") ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
                       data = beta)
print(anova(mod_bray_glob, permutations = perm))
print(anova(mod_bray_glob, by = "margin", permutations = perm))

cat("\n---- dbRDA: Jaccard ----\n")
mod_jacc_glob <- dbrda(vegdist(comm, "jacc",binary = TRUE) ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
                       data = beta)
print(anova(mod_jacc_glob, permutations = perm))
print(anova(mod_jacc_glob, by = "margin", permutations = perm))

cat("\n---- Variance partitioning: Bray ----\n")
var_bray_glob <- varpart(vegdist(comm, "bray"),
                         ~ SAM_s + ENSO_s + SeaIce_s,
                         ~ whale_s,
                         ~ site,
                         data = beta)
print(var_bray_glob)

cat("\n---- Variance partitioning: Jaccard ----\n")
var_jacc_glob <- varpart(vegdist(comm, "jacc",binary = TRUE),
                         ~ SAM_s + ENSO_s + SeaIce_s,
                         ~ whale_s,
                         ~ site,
                         data = beta)
print(var_jacc_glob)

cat("\n\n")





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
    vegdist(comm, "bray") ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
    data = beta
    
  )
  
  perm <- how(blocks = alpha$site, nperm = 100000)
  
  print(anova(mod_bray, permutations = perm))
  print(anova(mod_bray, by = "margin", permutations = perm))
  
  cat("\n---- dbRDA: Jaccard ----\n")
  mod_jacc <- dbrda(
    vegdist(comm, "jacc",binary = TRUE) ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
    data = beta
  )
  
  print(anova(mod_jacc, permutations = perm))
  print(anova(mod_jacc, by = "margin", permutations = perm))
  
  ## variance partitioning
  cat("\n---- Variance partitioning: Bray ----\n")
  var_bray <- varpart(
    vegdist(comm, "bray"),
    ~ SAM_s + ENSO_s + SeaIce_s,
    ~ whale_s,
    ~ site,
    data = beta
  )
  print(var_bray)
  
  cat("\n---- Variance partitioning: Jaccard ----\n")
  var_jacc <- varpart(
    vegdist(comm, "jacc", binary = TRUE),
    ~ SAM_s + ENSO_s + SeaIce_s,
    ~ whale_s,
    ~ site,
    data = beta
  )
  print(var_jacc)
  
  ## PERMANOVA
  cat("\n---- PERMANOVA: Bray ----\n")
  perm_bray <- adonis2(
    vegdist(comm, method = "bray") ~  SAM_s + ENSO_s + SeaIce_s + whale_s,
    data = beta,
    permutations = how(blocks = beta$site, nperm = 100000),
    by = "margin"
  )
  print(perm_bray)
  
  cat("\n---- PERMANOVA: Jaccard ----\n")
  perm_jacc <- adonis2(
    vegdist(comm, method = "jacc", binary = TRUE) ~  SAM_s + ENSO_s + SeaIce_s + whale_s,
    data = beta,
    permutations = how(blocks = beta$site, nperm = 100000),
    by = "margin"
  )
  print(perm_jacc)
  
  cat("\n\n")
}

sink()

####---- Build tidy summary tables ----####

groups <- list(
  Global   = NULL,                 # NULL means "use everything"
  Metazoa  = ASVs.met,
  Protists = ASVs.pro,
  Bacteria = ASVs.bac
)

# Tidy a single dbRDA marginal anova into a data frame
tidy_dbrda <- function(mod, group, distance){
  ov  <- anova(mod, permutations = perm)
  mar <- anova(mod, by = "margin", permutations = perm)
  eig <- eigenvals(mod)
  R2_constr <- sum(eig[grep("^dbRDA", names(eig))]) / sum(eig)
  data.frame(
    Group     = group,
    Distance  = distance,
    Term      = c("Overall model", rownames(mar)[rownames(mar) != "Residual"]),
    Df        = c(ov$Df[1], mar$Df[rownames(mar) != "Residual"]),
    SumOfSqs  = c(ov$SumOfSqs[1], mar$SumOfSqs[rownames(mar) != "Residual"]),
    F_value   = c(ov$F[1], mar$F[rownames(mar) != "Residual"]),
    p_value   = c(ov$`Pr(>F)`[1], mar$`Pr(>F)`[rownames(mar) != "Residual"]),
    R2_constrained = c(R2_constr, rep(NA, sum(rownames(mar) != "Residual"))),
    stringsAsFactors = FALSE
  )
}

# Tidy a single varpart object into a data frame
# X1 = Climate (SAM+ENSO+SeaIce), X2 = Whale, X3 = Site
# vegan returns 8 fractions for 3 explanatory tables:
#   [a] X1|X2+X3  [b] X2|X1+X3  [c] X3|X1+X2
#   [d] X1nX2     [e] X2nX3     [f] X1nX3     [g] X1nX2nX3   [h] Residual
tidy_varpart <- function(vp, group, distance){
  parts <- vp$part$indfract
  stopifnot(nrow(parts) == 8)   # guard: fails loudly if a call reverts to 2 tables
  data.frame(
    Group    = group,
    Distance = distance,
    Fraction = c("Climate (SAM+ENSO+SeaIce) | Whale + Site",
                 "Whale | Climate + Site",
                 "Site | Climate + Whale",
                 "Shared: Climate & Whale",
                 "Shared: Whale & Site",
                 "Shared: Climate & Site",
                 "Shared: Climate & Whale & Site",
                 "Residual"),
    Df       = c(parts$Df[1:3], rep(NA, 5)),
    Adj_R2 = parts[["Adj.R.square"]],
    stringsAsFactors = FALSE
  )
}

dbrda_rows  <- list()
varpart_rows <- list()

for (g in names(groups)) {
  asvs <- groups[[g]]
  
  if (is.null(asvs)) {
    mat <- euk.reps                                    # global = everything
  } else {
    mat <- euk.reps[rownames(euk.reps) %in% asvs, , drop = FALSE]
  }
  
  comm.raw <- prop.table(as.matrix(mat), 2)
  comm     <- t(comm.raw[, match(clim$longname, colnames(comm.raw))])
  
  perm <- how(blocks = beta$site, nperm = 100000)
  bray_dist <- vegan::vegdist(comm, method = "bray")
  
  jacc_dist <- vegan::vegdist(
    comm,
    method = "jaccard",
    binary = TRUE
  )
  
  mb <- vegan::dbrda(
    bray_dist ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
    data = beta
  )
  
  mj <- vegan::dbrda(
    jacc_dist ~ whale_s + SAM_s + ENSO_s + SeaIce_s + Condition(site),
    data = beta
  )
  
  vb <- vegan::varpart(
    bray_dist,
    ~ SAM_s + ENSO_s + SeaIce_s,
    ~ whale_s,
    ~ site,
    data = beta
  )
  
  vj <- vegan::varpart(
    jacc_dist,
    ~ SAM_s + ENSO_s + SeaIce_s,
    ~ whale_s,
    ~ site,
    data = beta
  )
  dbrda_rows[[paste(g, "bray")]]    <- tidy_dbrda(mb, g, "Bray-Curtis")
  dbrda_rows[[paste(g, "jaccard")]] <- tidy_dbrda(mj, g, "Jaccard")
  
  
  varpart_rows[[paste(g, "bray")]]    <- tidy_varpart(vb, g, "Bray-Curtis")
  varpart_rows[[paste(g, "jaccard")]] <- tidy_varpart(vj, g, "Jaccard")
}

dbrda_summary  <- do.call(rbind, dbrda_rows)
varpart_summary <- do.call(rbind, varpart_rows)

# Round for display
# was: round(...4 columns including p_value...)
dbrda_summary[c("SumOfSqs","F_value","R2_constrained")] <-
  round(dbrda_summary[c("SumOfSqs","F_value","R2_constrained")], 4)
dbrda_summary$p_value <- signif(dbrda_summary$p_value, 3)
varpart_summary$Adj_R2 <- round(varpart_summary$Adj_R2, 4)

write.csv(dbrda_summary,   "summaryfiles/dbRDA_summary.csv",   row.names = FALSE)
write.csv(varpart_summary, "summaryfiles/varpart_summary.csv", row.names = FALSE)

print(dbrda_summary)
print(varpart_summary)

#### now lets run a HMSC model with all taxa
library(Hmsc)

## Beta coefficient gradient: red = negative, green = positive
beta.grad <- colorRampPalette(c("#c13d2a","#db964d","#FFFFFF","#6ab79f","#244e51"))(100)


# -- Data prep ----------------------------------------------------------------
comm.raw.t <- t(as.matrix(euk.reps))
comm.raw.t <- comm.raw.t[match(clim$longname, rownames(comm.raw.t)), ]

# Normalise on full community FIRST
comm.prop <- comm.raw.t / rowSums(comm.raw.t)

# High-confidence Subdivision assignments only (>=60%)
tax.filt <- euk.tax[euk.tax$Subdivision.1 >= 60 & !is.na(euk.tax$Subdivision), ]

# Intersect ASVs, subset, collapse to Subdivision
shared.asvs <- intersect(colnames(comm.prop), tax.filt$ASV)
comm.filt   <- comm.prop[, shared.asvs, drop = FALSE]
tax.filt    <- tax.filt[match(shared.asvs, tax.filt$ASV), ]
comm.sub    <- t(rowsum(t(comm.filt), group = tax.filt$Subdivision, na.rm = TRUE))

# Zero replacement and CLR transformation
zero_replace <- function(x, delta = NULL) {
  if (is.null(delta)) delta <- min(x[x > 0]) * 0.65
  x[x == 0] <- delta
  x <- x / rowSums(x)
  x
}

clr <- function(x) {
  lx <- log(x)
  lx - rowMeans(lx)
}

zero_pct <- colMeans(comm.sub == 0)
#comm.sub <- comm.sub[, zero_pct < 0.70]
comm.sub <- zero_replace(comm.sub)
Y <- clr(comm.sub)
cat("Subdivisions retained:", ncol(Y), "\n")
for (p in c("SAM", "ENSO", "SeaIce")) {
  kk <- grep(p, m.fam$covNames, ignore.case = TRUE)
  ps <- sapply(post.g, function(s) {
    b <- s$Beta[kk, ]; mean(outer(b[wc.grp], b[!wc.grp], ">")) })
  cat(sprintf("%s: P(sup) = %.2f (%.2f-%.2f)\n",
      p, mean(ps), quantile(ps, .025), quantile(ps, .975)))
}

comm.sub <- zero_replace(comm.sub)
Y <- clr(comm.sub)

cat("Samples:", nrow(Y), "| Subdivisions:", ncol(Y), "\n")

# Predictors
X <- data.frame(
  SAM    = as.numeric(scale(beta$SAM_s)),
  ENSO   = as.numeric(scale(beta$ENSO_s)),
  SeaIce = as.numeric(scale(beta$SeaIce_s)),
  whale = as.numeric(scale(beta$whale_s))
)

# -- Model --------------------------------------------------------------------

studyDesign <- data.frame(site = factor(beta$site))
rL.site     <- HmscRandomLevel(units = levels(studyDesign$site))

m <- Hmsc(
  Y           = Y,
  XData       = X,
  XFormula    = ~ SAM + ENSO + SeaIce + whale,
  distr       = "normal",
  studyDesign = studyDesign,
  ranLevels   = list(site = rL.site)
)

#m <- sampleMcmc(  m,
#                  thin      = 100,
#                  samples   = 4000,
#                  transient = 100000,
#                  nChains   = 6,
#                  nParallel = 6,
#                  verbose   = 5000
#)
#saveRDS(m, file = "cleandata/MCMCchains/SubDiv.rds")
m <- readRDS(file = "cleandata/MCMCchains/SubDiv.rds")

# -- Convergence --------------------------------------------------------------

mpost    <- convertToCodaObject(m)
beta_ess <- effectiveSize(mpost$Beta)
beta_gel <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

cat("ESS — Min:", round(min(beta_ess)),
    "| Median:", round(median(beta_ess)),
    "| % > 200:", round(mean(beta_ess > 200) * 100, 1), "%\n")
cat("G-R — Max:", round(max(beta_gel[, 1]), 3),
    "| % > 1.05:", round(mean(beta_gel[, 1] > 1.05) * 100, 1), "%\n")

# -- Results ------------------------------------------------------------------

# Variance partitioning — design matrix cols: 1=Intercept, 2:4=Climate, 5=whale
VP <- computeVariancePartitioning(m,
                                  group = c(1, 2, 2, 2, 3), 
                                  groupnames = c("Intercept", "Climate", "Whaling"))
print(round(rowMeans(VP$vals), 3))
plotVariancePartitioning(m, VP = VP,col=c("darkred","darkgreen","dodgerblue"))

# Posterior support for fixed effects
postBeta   <- getPostEstimate(m, parName = "Beta")
support_df <- data.frame(
  Subdivision = rep(colnames(Y), each = 4),
  predictor   = rep(c("SAM", "ENSO", "SeaIce", "Whaling"), times = ncol(Y)),
  mean_beta   = c(rbind(postBeta$mean[2, ],       postBeta$mean[3, ],       postBeta$mean[4, ], postBeta$mean[5, ])),
  P_pos       = c(rbind(postBeta$support[2, ],    postBeta$support[3, ],    postBeta$support[4, ], postBeta$support[5, ])),
  P_neg       = c(rbind(postBeta$supportNeg[2, ], postBeta$supportNeg[3, ], postBeta$supportNeg[4, ], postBeta$supportNeg[5, ]))
)
support_df$direction <- ifelse(support_df$mean_beta > 0, "positive", "negative")
support_df$strong    <- support_df$P_pos > 0.90 | support_df$P_neg > 0.90
print(support_df[support_df$strong, c("Subdivision", "predictor", "mean_beta",
                                      "P_pos", "P_neg", "direction")])


## -- Prokaryotic subdivisions: retained in the fitted model (they form part of
## -- the compositional denominator) but dropped from all reporting figures,
## -- since 18S gives no quantitative measure of prokaryotic abundance.
## -- prok.subdiv is defined once near the top, alongside ASVs.bac.
euk.keep <- setdiff(colnames(Y), prok.subdiv)
cat("Subdivisions plotted:", length(euk.keep), "of", ncol(Y),
    "(", ncol(Y) - length(euk.keep), "prokaryotic dropped from figures )\n")


## Nice plot

# Extract and reshape VP values
vp_mat <- VP$vals  # rows = components, cols = subdivisions
vp_mat_rev <- vp_mat[,dim(vp_mat)[2]:1] 
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
vp_mat_rev <- vp_mat_rev[, colnames(vp_mat_rev) %in% euk.keep, drop = FALSE]
site_var <- vp_mat_rev["Random: site", ]
vp_mat_rev <- vp_mat_rev[, order(site_var, decreasing = TRUE)]

# Colours matching your original call
cols <- c("#31424e","#7d9bb0","#d6d2c9") 


# Margins: bottom, left, top, right — left margin wide for taxon names

pdf("figures/HMSC/subdivision.var.pdf",height=9,width = 5)
par(mar = c(4, 10, 2, 8))

bp <- barplot(
  vp_mat_rev,
  horiz    = TRUE,
  col      = cols,
  border   = NA,
  las      = 1,
  xlim     = c(0, 1),
  xlab     = "Variance explained",
  names.arg = gsub("_X","",colnames(vp_mat_rev)),
  cex.names = 0.75,
  axes = FALSE)

# x-axis as percentages
axis(1, at = seq(0, 1, 0.25),
     labels = paste0(seq(0, 100, 25), "%"),
     cex.axis = 0.85)

# Legend outside right margin
legend(1.01,48,
       inset   = c(-0.35, 0),
       legend  = rownames(vp_mat_rev),
       fill    = cols,
       border  = NA,
       bty     = "n",
       cex     = 0.85,
       xpd     = TRUE)

dev.off()




#Now we plot what subdivisions respond to what 

# Identify subdivisions with at least one significant response
sup_pos <- postBeta$support[2:5, ]
sup_neg <- postBeta$supportNeg[2:5, ]
keep <- apply(sup_pos > 0.95 | sup_neg > 0.95, 2, any)
keep <- keep & (colnames(Y) %in% euk.keep)

beta_mat <- postBeta$mean[2:5, keep]
sup_pos  <- sup_pos[, keep]
sup_neg  <- sup_neg[, keep]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")

# Mark significance at 0.95 threshold
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

# Order subdivisions by whaling response for readability

ord <- order(beta_mat["WhaleBiomass", ])
#ord <- rev(order(colnames(beta_mat)))
beta_plot <- t(beta_mat[, ord])
sig_plot  <- t(sig_mat[, ord])

# Plot
pdf("figures/HMSC/subdivision.var.sig.pdf",height=9,width=6)
par(mar = c(5, 12, 3, 2))
n_sub <- nrow(beta_plot)
n_pred <- ncol(beta_plot)

image(1:n_pred, 1:n_sub, t(beta_plot),
      col = beta.grad,
      zlim = c(-max(abs(beta_plot)), max(abs(beta_plot))),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:n_pred, labels = colnames(beta_plot), las = 1)
axis(2, at = 1:n_sub, labels = gsub("_X","",rownames(beta_plot)), las = 2, cex.axis = 0.6)

# Add significance markers
for (i in 1:n_pred) {
  for (j in 1:n_sub) {
    if (sig_plot[j, i] != "") {
      text(i, j, sig_plot[j, i], cex = 0.7)
    }
  }
}
#title("Subdivision responses to environmental drivers",cex=0.6)
dev.off()

###Everything together plot 

beta_to_col <- function(x, zlim, cols) {
  if (is.na(x)) return("grey90")
  idx <- round((x - zlim[1]) / diff(zlim) * (length(cols) - 1)) + 1
  idx <- max(1, min(length(cols), idx))
  cols[idx]
}
heat_cols <- beta.grad


# -- Combined: Subdivision -----------------------------------------------------
vp_mat <- VP$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]

# Order by whaling beta 
whale_beta <- postBeta$mean[5, ]  # row 5 = whale
whale_order <- names(sort(whale_beta))  # ascending = most positive at top in horiz barplot
whale_order <- whale_order[whale_order %in% euk.keep]
vp_mat_rev <- vp_mat_rev[, whale_order]
taxon_order <- whale_order

# Beta and significance in same order
beta_mat <- postBeta$mean[2:5, , drop = FALSE]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")
sup_pos <- postBeta$support[2:5, ]
sup_neg <- postBeta$supportNeg[2:5, ]

sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

beta_full <- beta_mat[, taxon_order]
sig_full  <- sig_mat[, match(taxon_order, colnames(beta_mat))]

zlim <- c(-max(abs(beta_full)), max(abs(beta_full)))
cols_vp <- c("#31424e","#7d9bb0","#d6d2c9") 
pred_names <- rownames(beta_full)


pdf("figures/HMSC/subdivision.combined.pdf", height = 9, width = 7)
par(mar = c(4, 10, 3, 12),xpd = TRUE)

bp <- barplot(
  vp_mat_rev,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "", 
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = gsub("_X", "", colnames(vp_mat_rev)),
  cex.names = 0.7,
  axes      = FALSE
)
axis(1, at = seq(0, 1, 0.25),
     labels = paste0(seq(0, 100, 25), "%"), cex.axis = 0.85)
mtext("Variance explained", side = 1, line = 2.5, at = 0.5, cex = 0.85)
legend(-1,-2,
       legend = rownames(vp_mat_rev),
       fill   = cols_vp,
       border = NA, bty = "n", cex = 0.8)

# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.30
gap     <- 0.04
x_start <- 1.15

par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order)) {
    cell_col <- beta_to_col(beta_full[k, j], zlim, heat_cols)
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_full[k, j] != "") {
      text((xl + xr) / 2, bp[j], sig_full[k, j], cex = 0.6)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}
# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5

for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}

# Numbers above the bar
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))

# Beta label below the bar
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)

dev.off()
## now just the sig taxa

# Keep only taxa with at least one significant predictor
keep <- apply(sig_full, 2, function(x) any(x != ""))

vp_mat_rev   <- vp_mat_rev[, keep, drop = FALSE]
beta_full    <- beta_full[, keep, drop = FALSE]
sig_full     <- sig_full[, keep, drop = FALSE]
taxon_order  <- taxon_order[keep]

## ---- R2 scaling for VP bars -------------------------------------------------
## Bars show share of TOTAL variance (VP proportions x per-taxon R2), normalised
## so the largest single bar spans the full bar panel; axis relabelled to match.
## NA R2 (Tjur undefined for taxa present at all sites) -> zero-length bar.
scale_vp <- function(vp, r2) {
  r2[is.na(r2)] <- 0
  sc <- max(r2[colnames(vp)])
  list(mat = sweep(vp, 2, r2[colnames(vp)] / sc, `*`), sc = sc)
}
vp_axis <- function(sc) {
  tk <- pretty(c(0, sc)); tk <- tk[tk <= sc + 1e-9]   # drop this filter to allow a tick past the bar end
  axis(1, at = tk / sc, labels = paste0(round(100 * tk), "%"), cex.axis = 0.85)
}

fit.sub <- evaluateModelFit(m, predY = computePredictedValues(m))
r2.sub  <- setNames(fit.sub$R2, m$spNames)
sv <- scale_vp(vp_mat_rev, r2.sub)

pdf("figures/HMSC/subdivision.combined.sig.pdf", height = 9, width = 7)
par(mar = c(4, 10, 3, 12),xpd = TRUE)

bp <- barplot(
  sv$mat,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "", 
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = gsub("_X", "", colnames(vp_mat_rev)),
  cex.names = 0.7,
  axes      = FALSE
)
vp_axis(sv$sc)
mtext("Variance explained (% of total)", side = 1, line = 2.5, at = 0.5, cex = 0.85)
legend(-1,-2,
       legend = rownames(vp_mat_rev),
       fill   = cols_vp,
       border = NA, bty = "n", cex = 0.8)

# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.30
gap     <- 0.04
x_start <- 1.15

par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order)) {
    cell_col <- if (sig_full[k, j] != "") {
      beta_to_col(beta_full[k, j], zlim, heat_cols)
    } else "#f2f0ec"
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_full[k, j] == "**") {
      rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
           border = "black", lwd = 0.9)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}
# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5

for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}

# Numbers above the bar
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))

# Beta label below the bar
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)

dev.off()


# -- Metazoa family-level HMSC ------------------------------------------------
# Filter to Metazoa with perfect family assignments only
tax.meta <- euk.tax[euk.tax$Subdivision == "Metazoa" & 
                      !is.na(euk.tax$family) &
                      euk.tax$Family.1 == 100, ]
tax.meta <- tax.meta[tax.meta$family != "", ]
cat("ASVs retained:", nrow(tax.meta), "\n")
cat("Unique families:", length(unique(tax.meta$family)), "\n")
print(sort(table(tax.meta$family), decreasing = TRUE))

# Intersect with proportion matrix (normalised on full community)
shared.meta <- intersect(colnames(comm.prop), tax.meta$ASV)
comm.meta   <- comm.prop[, shared.meta, drop = FALSE]
tax.meta    <- tax.meta[match(shared.meta, tax.meta$ASV), ]

# Collapse to family
comm.fam <- t(rowsum(t(comm.meta), group = tax.meta$family, na.rm = TRUE))

# Prevalence filter
prev <- colMeans(comm.fam > 0)
cat("\nPrevalence distribution:\n")
print(round(sort(prev, decreasing = TRUE), 2))
comm.fam <- comm.fam[, prev > 0.10]
cat("Families after prevalence filter:", ncol(comm.fam), "\n")

# Zero replacement and CLR
comm.fam <- zero_replace(comm.fam)
Y.fam    <- clr(comm.fam)

cat("Samples:", nrow(Y.fam), "| Families:", ncol(Y.fam), "\n")
cat("Any NA:", anyNA(Y.fam), "\n")

# -- Model --------------------------------------------------------------------
m.fam <- Hmsc(
  Y           = Y.fam,
  XData       = X,
  XFormula    = ~ SAM + ENSO + SeaIce + whale,
  distr       = "normal",
  studyDesign = studyDesign,
  ranLevels   = list(site = rL.site)
)
#m.fam <- sampleMcmc(m.fam,
#                    thin      = 100,
#                    samples   = 4000,
#                    transient = 100000,
#                    nChains   = 6,
#                    nParallel = 6,
#                    verbose   = 5000
#)
#saveRDS(m.fam, file = "cleandata/MCMCchains/Meta.rds")
m.fam <- readRDS(file = "cleandata/MCMCchains/Meta.rds")

# -- Convergence --------------------------------------------------------------
mpost.fam    <- convertToCodaObject(m.fam)
beta_ess.fam <- effectiveSize(mpost.fam$Beta)
beta_gel.fam <- gelman.diag(mpost.fam$Beta, multivariate = FALSE)$psrf
cat("ESS — Min:", round(min(beta_ess.fam)),
    "| Median:", round(median(beta_ess.fam)),
    "| % > 200:", round(mean(beta_ess.fam > 200) * 100, 1), "%\n")
cat("G-R — Max:", round(max(beta_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(beta_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

omega_ess.fam <- effectiveSize(mpost.fam$Omega[[1]])
omega_gel.fam <- gelman.diag(mpost.fam$Omega[[1]], multivariate = FALSE)$psrf
cat("Omega ESS — Min:", round(min(omega_ess.fam)),
    "| Median:", round(median(omega_ess.fam)),
    "| % > 200:", round(mean(omega_ess.fam > 200) * 100, 1), "%\n")
cat("Omega G-R — Max:", round(max(omega_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(omega_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

sigma_ess.fam <- effectiveSize(mpost.fam$Sigma)
sigma_gel.fam <- gelman.diag(mpost.fam$Sigma, multivariate = FALSE)$psrf
cat("Sigma ESS — Min:", round(min(sigma_ess.fam)),
    "| Median:", round(median(sigma_ess.fam)),
    "| % > 200:", round(mean(sigma_ess.fam > 200) * 100, 1), "%\n")
cat("Sigma G-R — Max:", round(max(sigma_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(sigma_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

# -- Results ------------------------------------------------------------------
VP.fam <- computeVariancePartitioning(m.fam,
                                      group = c(1, 2, 2, 2, 3),
                                      groupnames = c("Intercept", "Climate", "Whaling"))
print(round(rowMeans(VP.fam$vals), 3))
plotVariancePartitioning(m.fam, VP = VP.fam)

postBeta.fam <- getPostEstimate(m.fam, parName = "Beta")
support.fam <- data.frame(
  family    = rep(colnames(Y.fam), each = 4),
  predictor = rep(c("SAM", "ENSO", "SeaIce", "Whaling"), times = ncol(Y.fam)),
  mean_beta = c(rbind(postBeta.fam$mean[2, ],       postBeta.fam$mean[3, ],       postBeta.fam$mean[4, ], postBeta.fam$mean[5, ])),
  P_pos     = c(rbind(postBeta.fam$support[2, ],    postBeta.fam$support[3, ],    postBeta.fam$support[4, ], postBeta.fam$support[5, ])),
  P_neg     = c(rbind(postBeta.fam$supportNeg[2, ], postBeta.fam$supportNeg[3, ], postBeta.fam$supportNeg[4, ], postBeta.fam$supportNeg[5, ]))
)
support.fam$direction <- ifelse(support.fam$mean_beta > 0, "positive", "negative")
support.fam$strong    <- support.fam$P_pos > 0.95 | support.fam$P_neg > 0.95
print(support.fam[support.fam$strong, c("family", "predictor", "mean_beta",
                                        "P_pos", "P_neg", "direction")])

## Nice plot

# Extract and reshape VP values
vp_mat <- VP.fam$vals  # rows = families, cols = components
vp_mat_rev <- vp_mat[,dim(vp_mat)[2]:1] 
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
site_var <- vp_mat_rev["Random: site", ]
vp_mat_rev <- vp_mat_rev[, order(site_var, decreasing = TRUE)]

# Colours matching your original call
cols <- c("#31424e","#7d9bb0","#d6d2c9") 

# Margins: bottom, left, top, right — left margin wide for taxon names

pdf("figures/HMSC/metazoa.var.pdf",height=9,width = 5)
par(mar = c(4, 10, 2, 8))

bp <- barplot(
  vp_mat_rev,
  horiz    = TRUE,
  col      = cols,
  border   = NA,
  las      = 1,
  xlim     = c(0, 1),
  xlab     = "Variance explained",
  names.arg = gsub("_X","",colnames(vp_mat_rev)),
  cex.names = 0.75,
  axes = FALSE)

# x-axis as percentages
axis(1, at = seq(0, 1, 0.25),
     labels = paste0(seq(0, 100, 25), "%"),
     cex.axis = 0.85)

# Legend outside right margin
legend(1.01,36,
       inset   = c(-0.35, 0),
       legend  = rownames(vp_mat_rev),
       fill    = cols,
       border  = NA,
       bty     = "n",
       cex     = 0.85,
       xpd     = TRUE)

dev.off()


# -- heatmap ------------------------------------------------------------------


# Extract posterior means and support
postBeta.fam <- getPostEstimate(m.fam, parName = "Beta")

beta_mat <- postBeta.fam$mean[2:5, ]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")

sup_pos <- postBeta.fam$support[2:5, ]
sup_neg <- postBeta.fam$supportNeg[2:5, ]

# Filter to families with at least one significant response
keep <- apply(sup_pos > 0.95 | sup_neg > 0.95, 2, any)
beta_mat <- beta_mat[, keep]
sup_pos  <- sup_pos[, keep]
sup_neg  <- sup_neg[, keep]

cat("Retained:", sum(keep), "of", length(keep), "families\n")

# Significance markers
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

# Alphabetical order
ord <- rev(order(colnames(beta_mat)))
beta_plot <- t(beta_mat[, ord])
sig_plot  <- t(sig_mat[, ord])

# Plot
pdf("figures/HMSC/metazoa.var.sig.pdf",height=9,width=6)
par(mar = c(5, 12, 3, 2))
n_fam  <- nrow(beta_plot)
n_pred <- ncol(beta_plot)

image(1:n_pred, 1:n_fam, t(beta_plot),
      col = beta.grad,
      zlim = c(-max(abs(beta_plot)), max(abs(beta_plot))),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:n_pred, labels = colnames(beta_plot), las = 1)
axis(2, at = 1:n_fam, labels = rownames(beta_plot), las = 2, cex.axis = 0.7)

for (i in 1:n_pred) {
  for (j in 1:n_fam) {
    if (sig_plot[j, i] != "") {
      text(i, j, sig_plot[j, i], cex = 0.7)
    }
  }
}
title("Metazoan family responses to environmental drivers")
dev.off()

# -- Combined: Metazoa CLR ----------------------------------------------------
vp_mat <- VP.fam$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
# Order by whaling beta
whale_beta <- postBeta.fam$mean[5, ]
whale_order <- names(sort(whale_beta))
vp_mat_rev <- vp_mat_rev[, whale_order]
taxon_order <- whale_order
# Beta and significance in same order
beta_mat <- postBeta.fam$mean[2:5, , drop = FALSE]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")
sup_pos <- postBeta.fam$support[2:5, ]
sup_neg <- postBeta.fam$supportNeg[2:5, ]
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"
beta_full <- beta_mat[, taxon_order]
sig_full  <- sig_mat[, match(taxon_order, colnames(beta_mat))]
zlim <- c(-max(abs(beta_full)), max(abs(beta_full)))
cols_vp <- c("#31424e","#7d9bb0","#d6d2c9")
pred_names <- rownames(beta_full)
pdf("figures/HMSC/metazoa.combined.pdf", height = 9, width = 7)
par(mar = c(4, 10, 3, 12), xpd = TRUE)
bp <- barplot(
  vp_mat_rev,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "",
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = gsub("_X", "", colnames(vp_mat_rev)),
  cex.names = 0.7,
  axes      = FALSE
)
axis(1, at = seq(0, 1, 0.25),
     labels = paste0(seq(0, 100, 25), "%"), cex.axis = 0.85)
mtext("Variance explained", side = 1, line = 2.5, at = 0.5, cex = 0.85)
legend(-1, -2,
       legend = rownames(vp_mat_rev),
       fill   = cols_vp,
       border = NA, bty = "n", cex = 0.8)
# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.30
gap     <- 0.04
x_start <- 1.15
par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order)) {
    cell_col <- beta_to_col(beta_full[k, j], zlim, heat_cols)
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_full[k, j] != "") {
      text((xl + xr) / 2, bp[j], sig_full[k, j], cex = 0.6)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}
# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5
for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)
dev.off()

# -- Combined: Metazoa CLR ----------------------------------------------------
vp_mat <- VP.fam$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
# Order by whaling beta
whale_beta <- postBeta.fam$mean[5, ]
whale_order <- names(sort(whale_beta))
vp_mat_rev <- vp_mat_rev[, whale_order]
taxon_order <- whale_order
# Beta and significance in same order
beta_mat <- postBeta.fam$mean[2:5, , drop = FALSE]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")
sup_pos <- postBeta.fam$support[2:5, ]
sup_neg <- postBeta.fam$supportNeg[2:5, ]
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"
beta_full <- beta_mat[, taxon_order]
sig_full  <- sig_mat[, match(taxon_order, colnames(beta_mat))]
zlim <- c(-max(abs(beta_full)), max(abs(beta_full)))
cols_vp <- c("#31424e","#7d9bb0","#d6d2c9")
pred_names <- rownames(beta_full)

# Keep only taxa with ≥1 significant predictor

keep <- colSums(sig_full != "") > 0

vp_mat_sig   <- vp_mat_rev[, keep, drop = FALSE]
beta_sig     <- beta_full[, keep, drop = FALSE]
sig_sig      <- sig_full[, keep, drop = FALSE]
taxon_order_sig <- taxon_order[keep]

# Functional groupings + RAMS-screened names (rawdata/metazoa_functional_groups_RAMS.csv)
met.groups <- read.csv("rawdata/metazoa_functional_groups_RAMS.csv")
guild_pal <- c("Pelagic water-column feeder"       = "#24507F",
               "Benthic suspension feeder"          = "#8FB8DE",
               "Soft-sediment deposit feeder"       = "#7A5C3D",
               "Epibenthic grazer"                  = "#B09265",
               "Mobile benthic predator/scavenger"  = "#8F7F6B",
               "Meiofauna"                          = "#D4C2A0")
guild_cols <- setNames(guild_pal[met.groups$group], met.groups$family)
disp     <- setNames(gsub(" \\(fam\\. indet\\.\\)", "", met.groups$display_name),
                     met.groups$family)
disp_dag <- setNames(grepl("fam\\. indet", met.groups$display_name),
                     met.groups$family)
add_daggers <- function(fams, x, y, dy, cex = 0.45) {
  i <- which(disp_dag[fams])
  if (length(i)) text(x, y[i] + dy, "\u2020", adj = 0, cex = cex, xpd = TRUE)
}

## Scale VP stacks by per-family R2 (share of TOTAL variance), normalised to
## the largest bar. fit.fam$R2 is unnamed - stamp names first!
fit.fam <- evaluateModelFit(m.fam, predY = computePredictedValues(m.fam))
r2.fam  <- setNames(fit.fam$R2, m.fam$spNames)
sv <- scale_vp(vp_mat_sig, r2.fam)

pdf("figures/HMSC/metazoa.combined.sig.pdf", height = 9, width = 7)
par(mar = c(7, 10, 3, 12), xpd = TRUE)
bp <- barplot(
  sv$mat,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "",
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = rep("", ncol(vp_mat_sig)),
  cex.names = 0.7,
  axes      = FALSE
)
vp_axis(sv$sc)
mtext("Variance explained (% of total)", side = 1, line = 2.5, at = 0.5, cex = 0.85)
legend(-1.2, -6.5,
       legend = rownames(vp_mat_sig),
       fill   = cols_vp,
       border = NA, bty = "n", cex = 0.5)
# Taxon labels (italic = resolved above family rank; see legend) + guild dots
text(-0.1, bp, labels = disp[taxon_order_sig], adj = 1, cex = 0.7)
add_daggers(taxon_order_sig, -0.1, bp, dy = 0.28)
points(rep(-0.03, length(bp)), bp, pch = 21,
       bg = guild_cols[taxon_order_sig], col = "white", lwd = 0.4, cex = 1.8)
legend(-1.2, -1.5, legend = names(guild_pal), pch = 21, pt.bg = guild_pal,
       col = "white", pt.cex = 1.3, bty = "n", cex = 0.5)

# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.30
gap     <- 0.04
x_start <- 1.15
par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order_sig)) {
    cell_col <- if (sig_sig[k, j] != "") {
      beta_to_col(beta_sig[k, j], zlim, heat_cols)
    } else "#f2f0ec"
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_sig[k, j] == "**") {
      rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
           border = "black", lwd = 0.9)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}

# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5
for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)
dev.off()

#### DIFFERENT ORDER

# -- Combined: Metazoa CLR ----------------------------------------------------
vp_mat <- VP.fam$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
fit.fam <- evaluateModelFit(m.fam, predY = computePredictedValues(m.fam))
r2.fam  <- setNames(fit.fam$R2, m.fam$spNames)
sv <- scale_vp(vp_mat_sig, r2.fam)

# Ordering
order <- names(sort(colSums(sv$mat[1:2,])))
vp_mat_rev <- vp_mat_rev[, order]
taxon_order <- order
# Beta and significance in same order
beta_mat <- postBeta.fam$mean[2:5, , drop = FALSE]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")
sup_pos <- postBeta.fam$support[2:5, ]
sup_neg <- postBeta.fam$supportNeg[2:5, ]
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"
beta_full <- beta_mat[, taxon_order]
sig_full  <- sig_mat[, match(taxon_order, colnames(beta_mat))]
pred_ord  <- c(4, 1, 2, 3)            # WhaleBiomass, SAM, ENSO, SeaIce
beta_full <- beta_full[pred_ord, , drop = FALSE]
sig_full  <- sig_full[pred_ord, , drop = FALSE]
zlim <- c(-max(abs(beta_full)), max(abs(beta_full)))
cols_vp <- c("#31424e","#7d9bb0","#d6d2c9")
pred_names <- rownames(beta_full)

# Keep only taxa with ≥1 significant predictor

keep <- colSums(sig_full != "") > 0

vp_mat_sig   <- vp_mat_rev[, keep, drop = FALSE]
beta_sig     <- beta_full[, keep, drop = FALSE]
sig_sig      <- sig_full[, keep, drop = FALSE]
taxon_order_sig <- taxon_order[keep]

disp     <- setNames(gsub(" \\(fam\\. indet\\.\\)", "", met.groups$display_name),
                     met.groups$family)
disp_dag <- setNames(grepl("fam\\. indet", met.groups$display_name),
                     met.groups$family)
add_daggers <- function(fams, x, y, dy, cex = 0.45) {
  i <- which(disp_dag[fams])
  if (length(i)) text(x, y[i] + dy, "\u2020", adj = 0, cex = cex, xpd = TRUE)
}

## Scale VP stacks by per-family R2 (share of TOTAL variance), normalised to
## the largest bar. fit.fam$R2 is unnamed - stamp names first!
fit.fam <- evaluateModelFit(m.fam, predY = computePredictedValues(m.fam))
r2.fam  <- setNames(fit.fam$R2, m.fam$spNames)
sv <- scale_vp(vp_mat_sig, r2.fam)

pdf("figures/HMSC/metazoa.combined.sig.o2.pdf", height = 9, width = 7)
par(mar = c(7, 10, 3, 12), xpd = TRUE)
bp <- barplot(
  sv$mat,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "",
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = rep("", ncol(vp_mat_sig)),
  cex.names = 0.7,
  axes      = FALSE
)
vp_axis(sv$sc)
mtext("Variance explained (%)", side = 1, line = 2.5, at = 0.5, cex = 0.85)
#legend(-1.2, -6.5,
 #      legend = rownames(vp_mat_sig),
  #     fill   = cols_vp,
   #    border = NA, bty = "n", cex = 0.5)
# Taxon labels (italic = resolved above family rank; see legend) + guild dots
text(-0.1, bp, labels = disp[taxon_order_sig], adj = 1, cex = 0.7)
add_daggers(taxon_order_sig, -0.1, bp, dy = 0.28)
#points(rep(-0.03, length(bp)), bp, pch = 21,
   #    bg = guild_cols[taxon_order_sig], col = "white", lwd = 0.4, cex = 1.8)
#legend(-1.2, -1.5, legend = names(guild_pal), pch = 21, pt.bg = guild_pal,
   #    col = "white", pt.cex = 1.3, bty = "n", cex = 0.5)

# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.15
gap     <- 0.04
x_start <- 1.15
par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order_sig)) {
    cell_col <- if (sig_sig[k, j] != "") {
      beta_to_col(beta_sig[k, j], zlim, heat_cols)
    } else "#f2f0ec"
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_sig[k, j] == "**") {
      rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
           border = "black", lwd = 0.9)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}

# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5
for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)
dev.off()





## ===== Guild contrast: whaling beta ~ feeding ecology (Metazoa CLR) ==========
## Water-column feeders (pelagic + benthic suspension) vs all benthic taxa.
## All modelled families included (no significance filtering). Families treated
## as independent units  descriptive association.

wc.grp <- met.groups$group[match(colnames(m.fam$Y), met.groups$family)] %in%
  c("Pelagic water-column feeder", "Benthic suspension feeder")
k.wh   <- grep("whale", m.fam$covNames, ignore.case = TRUE)
post.g <- Hmsc::poolMcmcChains(m.fam$postList)

gdiff <- sapply(post.g, function(s)                       # effect size (in-sample)
  mean(s$Beta[k.wh, wc.grp]) - mean(s$Beta[k.wh, !wc.grp]))
psup  <- sapply(post.g, function(s) {                     # generalising statistic
  b <- s$Beta[k.wh, ]; mean(outer(b[wc.grp], b[!wc.grp], ">")) })
bhat  <- Hmsc::getPostEstimate(m.fam, "Beta")$mean[k.wh, ]
wt    <- wilcox.test(bhat[wc.grp], bhat[!wc.grp])         # frequentist companion

sink("summaryfiles/HMSC_guild_tests.txt")
cat("==== Guild contrasts: whaling beta ~ feeding ecology ====\n\n")
cat("---- Metazoa CLR ----\n")
cat(sprintf("Water-column n = %d; benthic n = %d (all modelled families)\n",
            sum(wc.grp), sum(!wc.grp)))
cat(sprintf("Group-mean difference: %.3f (95%% CI %.3f to %.3f); P(diff > 0) %s\n",
            mean(gdiff), quantile(gdiff, .025), quantile(gdiff, .975),
            if (all(gdiff > 0)) sprintf("> %.4f (all %d draws positive)",
                                        1 - 1/length(gdiff), length(gdiff))
            else sprintf("= %.4f", mean(gdiff > 0))))
cat(sprintf("P(superiority): %.3f (95%% CI %.3f to %.3f)\n",
            mean(psup), quantile(psup, .025), quantile(psup, .975)))
cat(sprintf("Wilcoxon on posterior means: W = %.0f, p = %.4g\n\n",
            wt$statistic, wt$p.value))
sink()

## Figure: portrait - pooled compartments on top, families by descending beta
fam.draws <- sapply(post.g, function(s) s$Beta[k.wh, ])
rownames(fam.draws) <- m.fam$spNames
ord <- rev(taxon_order_sig)                     # largest beta at top
fam.draws <- fam.draws[ord, ]
gcol <- guild_pal[met.groups$group[match(ord, met.groups$family)]]

grp.draws <- rbind(                             # full 55-family pools
  benthic = sapply(post.g, function(s) mean(s$Beta[k.wh, !wc.grp])),
  water   = sapply(post.g, function(s) mean(s$Beta[k.wh,  wc.grp])))
grp.cols <- c(benthic = "#B09265", water = "#3F6F9F")
grp.labs <- c("Benthic taxa", "Water-column feeders")


pdf("figures/HMSC/metazoa.whaling.guild.posteriors.pdf", width = 6, height = 10)
par(mar = c(4, 9, 1, 1))
n  <- nrow(fam.draws)
yf <- n:1                                        # family rows, top -> bottom
yg <- n + c(1.8, 3.2)                            # pooled rows above (benthic, water)
xr <- range(fam.draws, grp.draws)
plot(NA, xlim = xr, ylim = c(0.4, max(yg) + 0.9), axes = FALSE,
     xlab = expression("Whaling "*beta*" (posterior)"), ylab = "")
abline(v = 0, col = "grey70", lwd = 0.8)
for (i in 1:n) {
  d <- density(fam.draws[i, ]); s <- d$y / max(d$y) * 0.42
  polygon(c(d$x, rev(d$x)), c(yf[i] - s, rev(yf[i] + s)),
          col = adjustcolor(gcol[i], 0.75), border = NA)
  segments(mean(fam.draws[i, ]), yf[i] - 0.1,
           mean(fam.draws[i, ]), yf[i] + 0.1, lwd = 1.2)
}
for (j in 1:2) {
  d <- density(grp.draws[j, ]); s <- d$y / max(d$y) * 0.65
  polygon(c(d$x, rev(d$x)), c(yg[j] - s, rev(yg[j] + s)),
          col = adjustcolor(grp.cols[j], 0.9), border = grp.cols[j], lwd = 0.8)
  ci <- quantile(grp.draws[j, ], c(0.025, 0.5, 0.975))
  segments(ci[1], yg[j], ci[3], yg[j], lwd = 1.4)
  segments(ci[2], yg[j] - 0.15, ci[2], yg[j] + 0.15, lwd = 2)
}
axis(1, cex.axis = 0.8)
text(xr[1] - 0.03 * diff(xr), yf, labels = disp[ord], adj = 1,
     cex = 0.55, xpd = TRUE)
text(xr[1] - 0.03 * diff(xr), yg, labels = grp.labs, adj = 1,
     cex = 0.65, font = 2, xpd = TRUE)
text(xr[2], mean(yg), sprintf("P(sup) = %.2f", mean(psup)),
     adj = 1, cex = 0.65)
legend("bottomright", legend = names(guild_pal), pt.bg = guild_pal, pch = 21,
       col = "white", pt.cex = 1.2, bty = "n", cex = 0.55)

stopifnot(identical(rownames(fam.draws), ord),
          !any(is.na(gcol)), !any(is.na(disp[ord])))
dev.off()

## ===== Figure: whaling beta by family, pooled guild violins on top ==========

## Families shown = those with >=1 predictor at >95% support.
## Defined here directly so it doesn't depend on taxon_order_sig from earlier blocks.
pb.fam  <- Hmsc::getPostEstimate(m.fam, "Beta")
#sig.any <- colSums(pb.fam$support[2:5, , drop = FALSE] > 0.95 |
#                     pb.fam$supportNeg[2:5, , drop = FALSE] > 0.95) > 0
#fam.show <- m.fam$spNames[sig.any]

fam.draws <- sapply(post.g, function(s) s$Beta[k.wh, ])
rownames(fam.draws) <- m.fam$spNames
fam.draws <- fam.draws[fam.show, , drop = FALSE]


## Families shown = every family in the model, i.e. exactly those pooled in the stats
fam.draws <- sapply(post.g, function(s) s$Beta[k.wh, ])
rownames(fam.draws) <- m.fam$spNames

## (3) order by posterior mean whaling beta, largest at top
ord       <- names(sort(rowMeans(fam.draws), decreasing = TRUE))
fam.draws <- fam.draws[ord, , drop = FALSE]
gcol      <- guild_pal[met.groups$group[match(ord, met.groups$family)]]

stopifnot(identical(rownames(fam.draws), ord),
          !any(is.na(gcol)), !any(is.na(disp[ord])), !any(is.na(disp_dag[ord])))

## (4) labels as in the combined plot: RAMS display names, dagger = fam. indet.
fam.labs <- as.expression(lapply(ord, function(f)
  if (disp_dag[[f]]) bquote(.(disp[[f]])^"\u2020") else disp[[f]]))

## pooled guild draws: per MCMC draw, mean whaling beta over ALL modelled families
grp.draws <- rbind(
  benthic = sapply(post.g, function(s) mean(s$Beta[k.wh, !wc.grp])),
  water   = sapply(post.g, function(s) mean(s$Beta[k.wh,  wc.grp])))
grp.cols <- c(benthic = "#B09265", water = "#3F6F9F")
grp.labs <- c(sprintf("Benthic taxa (n = %d)",        sum(!wc.grp)),
              sprintf("Water-column feeders (n = %d)", sum(wc.grp)))

## layout
n  <- nrow(fam.draws)
yf <- n:1                          # family rows, top -> bottom
yg <- n + c(2.2, 4.4)              # (1) pooled rows, more space for bigger violins
xr <- range(fam.draws, grp.draws)

pdf("figures/HMSC/metazoa.whaling.guild.posteriors.pdf", width = 6, height = 10)
par(mar = c(7, 10, 1, 1))
plot(NA, xlim = rev(xr),           # (3) reversed: left = higher beta
     ylim = c(0.4, max(yg) + 2.4), axes = FALSE,
     xlab = expression("Whaling " * beta * " (posterior)"), ylab = "")
abline(v = 0, col = "grey70", lwd = 0.8)
abline(h = n + 0.9, col = "grey85", lwd = 0.8)   # separates summary from families

## family violins (posterior mean tick kept)
for (i in 1:n) {
  d <- density(fam.draws[i, ]); s <- d$y / max(d$y) * 0.42
  polygon(c(d$x, rev(d$x)), c(yf[i] - s, rev(yf[i] + s)),
          col = adjustcolor(gcol[i], 0.75), border = NA)
  m <- mean(fam.draws[i, ])
  segments(m, yf[i] - 0.1, m, yf[i] + 0.1, lwd = 1.2)
}

## (1)+(2) larger pooled violins, no CI bar or median tick
for (j in 1:2) {
  d <- density(grp.draws[j, ]); s <- d$y / max(d$y) * 1.0
  polygon(c(d$x, rev(d$x)), c(yg[j] - s, rev(yg[j] + s)),
          col = adjustcolor(grp.cols[j], 0.9), border = grp.cols[j], lwd = 0.8)
}



## axes and labels (axis() copes with the reversed x-range, so no manual offsets)
axis(1, cex.axis = 0.8)
axis(2, at = yf, labels = fam.labs, las = 1, tick = FALSE,
     cex.axis = 0.55, mgp = c(3, 0.3, 0))
axis(2, at = yg, labels = grp.labs, las = 1, tick = FALSE,
     cex.axis = 0.65, font.axis = 2, mgp = c(3, 0.3, 0))

## guild legend below the axis so it can't collide with violins
legend("bottom", inset = c(0, -0.13), xpd = NA, ncol = 2,
       legend = names(guild_pal), pt.bg = guild_pal, pch = 21,
       col = "white", pt.cex = 1.2, bty = "n", cex = 0.55)
dev.off()





# -- Metazoa family-level HMSC (presence/absence) -----------------------------
tax.meta <- euk.tax[euk.tax$Subdivision == "Metazoa" & 
                      !is.na(euk.tax$family) &
                      euk.tax$Family.1 == 100, ]
tax.meta <- tax.meta[tax.meta$family != "", ]
cat("ASVs retained:", nrow(tax.meta), "\n")
cat("Unique families:", length(unique(tax.meta$family)), "\n")

# Intersect with raw count matrix (not proportions — we just need presence)
shared.meta <- intersect(colnames(comm.raw.t), tax.meta$ASV)
comm.meta.pa <- (comm.raw.t[, shared.meta, drop = FALSE] > 0) * 1
tax.meta     <- tax.meta[match(shared.meta, tax.meta$ASV), ]

# Collapse to family — present if any ASV detected
comm.fam.pa <- t(rowsum(t(comm.meta.pa), group = tax.meta$family, na.rm = TRUE))
comm.fam.pa <- (comm.fam.pa > 0) * 1

# Prevalence filter
prev <- colMeans(comm.fam.pa)
cat("\nPrevalence distribution:\n")
print(round(sort(prev, decreasing = TRUE), 2))

# Remove families at 100% or 0% prevalence (no variation to model)
comm.fam.pa <- comm.fam.pa[, prev > 0.1 & prev < 1.00]
cat("Families after filter:", ncol(comm.fam.pa), "\n")

Y.fam.pa <- comm.fam.pa

cat("Samples:", nrow(Y.fam.pa), "| Families:", ncol(Y.fam.pa), "\n")

# -- Model --------------------------------------------------------------------
m.fam.pa <- Hmsc(
  Y           = Y.fam.pa,
  XData       = X,
  XFormula    = ~ SAM + ENSO + SeaIce + whale,
  distr       = "probit",
  studyDesign = studyDesign,
  ranLevels   = list(site = rL.site)
)
#m.fam.pa <- sampleMcmc(m.fam.pa,
#                       thin      = 100,
#                       samples   = 4000,
#                       transient = 100000,
#                       nChains   = 6,
#                       nParallel = 6,
#                       verbose   = 5000
#)
#saveRDS(m.fam.pa, file = "cleandata/MCMCchains/Meta.PA.rds")
m.fam.pa <- readRDS(file = "cleandata/MCMCchains/Meta.PA.rds")


# -- Convergence --------------------------------------------------------------
mpost.fam.pa    <- convertToCodaObject(m.fam.pa)
beta_ess.fam.pa <- effectiveSize(mpost.fam.pa$Beta)
beta_gel.fam.pa <- gelman.diag(mpost.fam.pa$Beta, multivariate = FALSE)$psrf
cat("ESS — Min:", round(min(beta_ess.fam.pa)),
    "| Median:", round(median(beta_ess.fam.pa)),
    "| % > 200:", round(mean(beta_ess.fam.pa > 200) * 100, 1), "%\n")
cat("G-R — Max:", round(max(beta_gel.fam.pa[, 1]), 3),
    "| % > 1.05:", round(mean(beta_gel.fam.pa[, 1] > 1.05) * 100, 1), "%\n")

omega_ess.fam.pa <- effectiveSize(mpost.fam.pa$Omega[[1]])
omega_gel.fam.pa <- gelman.diag(mpost.fam.pa$Omega[[1]], multivariate = FALSE)$psrf
cat("Omega ESS — Min:", round(min(omega_ess.fam.pa)),
    "| Median:", round(median(omega_ess.fam.pa)),
    "| % > 200:", round(mean(omega_ess.fam.pa > 200) * 100, 1), "%\n")
cat("Omega G-R — Max:", round(max(omega_gel.fam.pa[, 1]), 3),
    "| % > 1.05:", round(mean(omega_gel.fam.pa[, 1] > 1.05) * 100, 1), "%\n")

sigma_ess.fam.pa <- effectiveSize(mpost.fam.pa$Sigma)
sigma_gel.fam.pa <- gelman.diag(mpost.fam.pa$Sigma, multivariate = FALSE)$psrf
cat("Sigma ESS — Min:", round(min(sigma_ess.fam.pa)),
    "| Median:", round(median(sigma_ess.fam.pa)),
    "| % > 200:", round(mean(sigma_ess.fam.pa > 200) * 100, 1), "%\n")
cat("Sigma G-R — Max:", round(max(sigma_gel.fam.pa[, 1]), 3),
    "| % > 1.05:", round(mean(sigma_gel.fam.pa[, 1] > 1.05) * 100, 1), "%\n")

# -- Results ------------------------------------------------------------------
VP.fam.pa <- computeVariancePartitioning(m.fam.pa,
                                         group = c(1, 2, 2, 2, 3),
                                         groupnames = c("Intercept", "Climate", "Whaling"))
print(round(rowMeans(VP.fam.pa$vals), 3))
plotVariancePartitioning(m.fam.pa, VP = VP.fam.pa)

postBeta.fam.pa <- getPostEstimate(m.fam.pa, parName = "Beta")
support.fam.pa <- data.frame(
  family    = rep(colnames(Y.fam.pa), each = 4),
  predictor = rep(c("SAM", "ENSO", "SeaIce", "Whaling"), times = ncol(Y.fam.pa)),
  mean_beta = c(rbind(postBeta.fam.pa$mean[2, ],       postBeta.fam.pa$mean[3, ],       postBeta.fam.pa$mean[4, ], postBeta.fam.pa$mean[5, ])),
  P_pos     = c(rbind(postBeta.fam.pa$support[2, ],    postBeta.fam.pa$support[3, ],    postBeta.fam.pa$support[4, ], postBeta.fam.pa$support[5, ])),
  P_neg     = c(rbind(postBeta.fam.pa$supportNeg[2, ], postBeta.fam.pa$supportNeg[3, ], postBeta.fam.pa$supportNeg[4, ], postBeta.fam.pa$supportNeg[5, ]))
)
support.fam.pa$direction <- ifelse(support.fam.pa$mean_beta > 0, "positive", "negative")
support.fam.pa$strong    <- support.fam.pa$P_pos > 0.95 | support.fam.pa$P_neg > 0.95
print(support.fam.pa[support.fam.pa$strong, c("family", "predictor", "mean_beta",
                                              "P_pos", "P_neg", "direction")])

## Nice VP plot
vp_mat <- VP.fam.pa$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]

cols <- c("#31424e","#7d9bb0","#d6d2c9")

pdf("figures/HMSC/metazoa.pa.var.pdf", height = 9, width = 5)
par(mar = c(4, 10, 2, 8))

bp <- barplot(
  vp_mat_rev,
  horiz     = TRUE,
  col       = cols,
  border    = NA,
  las       = 1,
  xlim      = c(0, 1),
  xlab      = "Variance explained",
  names.arg = gsub("_X", "", colnames(vp_mat_rev)),
  cex.names = 0.75,
  axes      = FALSE)

axis(1, at = seq(0, 1, 0.25),
     labels = paste0(seq(0, 100, 25), "%"),
     cex.axis = 0.85)

legend(1.01, 36,
       inset  = c(-0.35, 0),
       legend = rownames(vp_mat_rev),
       fill   = cols,
       border = NA,
       bty    = "n",
       cex    = 0.85,
       xpd    = TRUE)

dev.off()

# -- Heatmap ------------------------------------------------------------------
beta_mat <- postBeta.fam.pa$mean[2:5, ]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")

sup_pos <- postBeta.fam.pa$support[2:5, ]
sup_neg <- postBeta.fam.pa$supportNeg[2:5, ]

keep <- apply(sup_pos > 0.95 | sup_neg > 0.95, 2, any)
beta_mat <- beta_mat[, keep]
sup_pos  <- sup_pos[, keep]
sup_neg  <- sup_neg[, keep]

cat("Retained:", sum(keep), "of", length(keep), "families\n")

sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

ord <- rev(order(colnames(beta_mat)))
beta_plot <- t(beta_mat[, ord])
sig_plot  <- t(sig_mat[, ord])

pdf("figures/HMSC/metazoa.pa.sig.pdf", height = 9, width = 6)
par(mar = c(5, 12, 3, 2))
n_fam  <- nrow(beta_plot)
n_pred <- ncol(beta_plot)

image(1:n_pred, 1:n_fam, t(beta_plot),
      col = beta.grad,
      zlim = c(-max(abs(beta_plot)), max(abs(beta_plot))),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:n_pred, labels = colnames(beta_plot), las = 1)
axis(2, at = 1:n_fam, labels = rownames(beta_plot), las = 2, cex.axis = 0.7)

for (i in 1:n_pred) {
  for (j in 1:n_fam) {
    if (sig_plot[j, i] != "") {
      text(i, j, sig_plot[j, i], cex = 0.7)
    }
  }
}
title("Metazoan family responses (presence/absence)")
dev.off()

# -- Combined: Metazoa PA -----------------------------------------------------
vp_mat <- VP.fam.pa$vals
vp_mat_rev <- vp_mat[, dim(vp_mat)[2]:1]
vp_mat_rev <- vp_mat_rev[rownames(vp_mat_rev) != "Intercept", ]
# Order by whaling beta
whale_beta <- postBeta.fam.pa$mean[5, ]
whale_order <- names(sort(whale_beta))
vp_mat_rev <- vp_mat_rev[, whale_order]
taxon_order <- whale_order
# Beta and significance in same order
beta_mat <- postBeta.fam.pa$mean[2:5, , drop = FALSE]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")
sup_pos <- postBeta.fam.pa$support[2:5, ]
sup_neg <- postBeta.fam.pa$supportNeg[2:5, ]
sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"
beta_full <- beta_mat[, taxon_order]
sig_full  <- sig_mat[, match(taxon_order, colnames(beta_mat))]
zlim <- c(-max(abs(beta_full)), max(abs(beta_full)))
cols_vp <- c("#31424e","#7d9bb0","#d6d2c9")
pred_names <- rownames(beta_full)

# Keep only taxa with at least one significant predictor
keep <- apply(sig_full, 2, function(x) any(x != ""))

vp_mat_rev   <- vp_mat_rev[, keep, drop = FALSE]
beta_full    <- beta_full[, keep, drop = FALSE]
sig_full     <- sig_full[, keep, drop = FALSE]
taxon_order  <- taxon_order[keep]


fit.pa <- evaluateModelFit(m.fam.pa, predY = computePredictedValues(m.fam.pa))
r2.pa  <- setNames(fit.pa$TjurR2, m.fam.pa$spNames)
sv <- scale_vp(vp_mat_rev, r2.pa)

pdf("figures/HMSC/metazoa.pa.combined.sig.pdf", height = 10, width = 7)
par(mar = c(4, 10, 3, 12), xpd = TRUE)
bp <- barplot(
  sv$mat,
  horiz     = TRUE,
  col       = cols_vp,
  border    = NA,
  xlab      = "",
  las       = 1,
  xlim      = c(0, 1.8),
  names.arg = rep("", ncol(vp_mat_rev)),
  cex.names = 0.7,
  axes      = FALSE
)
# Taxon labels (italic = resolved above family rank; see legend); no guild dots
text(-0.02, bp, labels = disp[colnames(vp_mat_rev)], adj = 1, cex = 0.7)
vp_axis(sv$sc)
mtext(expression("Variance explained (% of total, Tjur "*R^2*")"),
      side = 1, line = 2.5, at = 0.5, cex = 0.85)
legend(-1, -1,
       legend = rownames(vp_mat_rev),
       fill   = cols_vp,
       border = NA, bty = "n", cex = 0.8)
# Heatmap strips
bar_h   <- diff(bp[1:2]) * 0.85
strip_w <- 0.30
gap     <- 0.04
x_start <- 1.15
par(xpd = TRUE)
for (k in 1:4) {
  xl <- x_start + (k - 1) * (strip_w + gap)
  xr <- xl + strip_w
  for (j in seq_along(taxon_order)) {
    cell_col <- if (sig_full[k, j] != "") {
      beta_to_col(beta_full[k, j], zlim, heat_cols)
    } else "#f2f0ec"
    rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
         col = cell_col, border = NA)
    if (sig_full[k, j] == "**") {
      rect(xl, bp[j] - bar_h / 2, xr, bp[j] + bar_h / 2,
           border = "black", lwd = 0.9)
    }
  }
  text((xl + xr) / 2, max(bp) + bar_h, pred_names[k],
       cex = 0.6, srt = 45, adj = c(0, 0.5))
}
# Colour bar legend for beta
n_leg <- 100
leg_xl <- x_start
leg_xr <- x_start + 4 * strip_w + 3 * gap
leg_w  <- (leg_xr - leg_xl) / n_leg
leg_yb <- min(bp) - bar_h * 2.7
leg_yt <- leg_yb + bar_h * 0.5
for (i in 1:n_leg) {
  rect(leg_xl + (i - 1) * leg_w, leg_yb,
       leg_xl + i * leg_w, leg_yt,
       col = heat_cols[i], border = NA)
}
text(leg_xl, leg_yt + bar_h * 0.2,
     round(zlim[1], 2), cex = 0.6, adj = c(0, 0))
text(leg_xr, leg_yt + bar_h * 0.2,
     round(zlim[2], 2), cex = 0.6, adj = c(1, 0))
text((leg_xl + leg_xr) / 2, leg_yt + bar_h * 0.2,
     "0", cex = 0.6, adj = c(0.5, 0))
text((leg_xl + leg_xr) / 2, leg_yb - bar_h * 0.4,
     expression(beta), cex = 0.7, adj = c(0.5, 1))
par(xpd = FALSE)
dev.off()

## ===== Guild contrast: whaling beta ~ feeding ecology (Metazoa PA) ===========
wc.pa <- met.groups$group[match(colnames(m.fam.pa$Y), met.groups$family)] %in%
  c("Pelagic water-column feeder", "Benthic suspension feeder")
k.pa    <- grep("whale", m.fam.pa$covNames, ignore.case = TRUE)
post.pa <- Hmsc::poolMcmcChains(m.fam.pa$postList)

gdiff.pa <- sapply(post.pa, function(s)
  mean(s$Beta[k.pa, wc.pa]) - mean(s$Beta[k.pa, !wc.pa]))
psup.pa  <- sapply(post.pa, function(s) {
  b <- s$Beta[k.pa, ]; mean(outer(b[wc.pa], b[!wc.pa], ">")) })
bhat.pa  <- Hmsc::getPostEstimate(m.fam.pa, "Beta")$mean[k.pa, ]
wt.pa    <- wilcox.test(bhat.pa[wc.pa], bhat.pa[!wc.pa])

sink("summaryfiles/HMSC_guild_tests.txt", append = TRUE)
cat("---- Metazoa PA (probit scale) ----\n")
cat(sprintf("Water-column n = %d; benthic n = %d (all modelled families)\n",
            sum(wc.pa), sum(!wc.pa)))
cat(sprintf("Group-mean difference: %.3f (95%% CI %.3f to %.3f); P(diff > 0) %s\n",
            mean(gdiff.pa), quantile(gdiff.pa, .025), quantile(gdiff.pa, .975),
            if (all(gdiff.pa > 0)) sprintf("> %.4f (all %d draws positive)",
                                           1 - 1/length(gdiff.pa), length(gdiff.pa))
            else sprintf("= %.4f", mean(gdiff.pa > 0))))
cat(sprintf("P(superiority): %.3f (95%% CI %.3f to %.3f)\n",
            mean(psup.pa), quantile(psup.pa, .025), quantile(psup.pa, .975)))
cat(sprintf("Wilcoxon on posterior means: W = %.0f, p = %.4g\n",
            wt.pa$statistic, wt.pa$p.value))
sink()

## Figure: portrait, PA - binary colours, largest beta at top
fam.draws <- sapply(post.pa, function(s) s$Beta[k.pa, ])
rownames(fam.draws) <- m.fam.pa$spNames
ord <- rev(taxon_order)                          # PA figure order, flipped
fam.draws <- fam.draws[ord, ]
bcol <- ifelse(wc.pa[match(ord, colnames(m.fam.pa$Y))], "#46309B", "#B09265")

grp.draws <- rbind(
  benthic = sapply(post.pa, function(s) mean(s$Beta[k.pa, !wc.pa])),
  water   = sapply(post.pa, function(s) mean(s$Beta[k.pa,  wc.pa])))
grp.cols <- c(benthic = "#B09265", water = "#46309B")
grp.labs <- c("Benthic taxa", "Water-column feeders")

pdf("figures/HMSC/metazoa.pa.whaling.guild.posteriors.pdf", width = 6, height = 10)
par(mar = c(4, 9, 1, 1))
n  <- nrow(fam.draws)
yf <- n:1
yg <- n + c(1.8, 3.2)
xr <- range(fam.draws, grp.draws)
plot(NA, xlim = xr, ylim = c(0.4, max(yg) + 0.9), axes = FALSE,
     xlab = expression("Whaling "*beta*" (posterior, probit scale)"), ylab = "")
abline(v = 0, col = "grey70", lwd = 0.8)
for (i in 1:n) {
  d <- density(fam.draws[i, ]); s <- d$y / max(d$y) * 0.42
  polygon(c(d$x, rev(d$x)), c(yf[i] - s, rev(yf[i] + s)),
          col = adjustcolor(bcol[i], 0.75), border = NA)
  segments(mean(fam.draws[i, ]), yf[i] - 0.1,
           mean(fam.draws[i, ]), yf[i] + 0.1, lwd = 1.2)
}
for (j in 1:2) {
  d <- density(grp.draws[j, ]); s <- d$y / max(d$y) * 0.65
  polygon(c(d$x, rev(d$x)), c(yg[j] - s, rev(yg[j] + s)),
          col = adjustcolor(grp.cols[j], 0.9), border = grp.cols[j], lwd = 0.8)
  ci <- quantile(grp.draws[j, ], c(0.025, 0.5, 0.975))
  segments(ci[1], yg[j], ci[3], yg[j], lwd = 1.4)
  segments(ci[2], yg[j] - 0.15, ci[2], yg[j] + 0.15, lwd = 2)
}
axis(1, cex.axis = 0.8)
text(xr[1] - 0.03 * diff(xr), yf, labels = disp[ord], adj = 1,
     cex = 0.55, xpd = TRUE)
text(xr[1] - 0.03 * diff(xr), yg, labels = grp.labs, adj = 1,
     cex = 0.65, font = 2, xpd = TRUE)
text(xr[2], mean(yg), sprintf("P(sup) = %.2f", mean(psup.pa)),
     adj = 1, cex = 0.65)
legend("bottomright", legend = c("Water-column feeder", "Benthic"),
       pt.bg = c("#46309B", "#B09265"), pch = 21, col = "white",
       pt.cex = 1.2, bty = "n", cex = 0.55)

stopifnot(identical(rownames(fam.draws), ord),
          !any(is.na(bcol)), !any(is.na(disp[ord])))
dev.off()

####====Export HMSC results to text file====####
## Writes a single summary file for all three HMSC models:
##   1. Subdivision-level CLR model        (m,        VP,        postBeta,        Y)
##   2. Metazoa family-level CLR model     (m.fam,    VP.fam,    postBeta.fam,    Y.fam)
##   3. Metazoa family P/A probit model    (m.fam.pa, VP.fam.pa, postBeta.fam.pa, Y.fam.pa)
##
## Assumes the HMSC section of the main script has been run so all of the
## objects above exist in the workspace, along with: beta, clim, studyDesign, X,
## mpost*, *_ess, *_gel for Beta/Omega/Sigma of each model.

if (!dir.exists("summaryfiles")) dir.create("summaryfiles", recursive = TRUE)

## ---- helpers ---------------------------------------------------------------

# Per-model summariser. mod = Hmsc object; Y = response matrix; VP = varpart;
# postBeta = posterior estimates; mpost = coda object; ess/gel lists for the
# three parameter blocks.
hmsc_summary_block <- function(name,
                               mod, Y, X, VP, postBeta, mpost,
                               beta_ess,  beta_gel,
                               omega_ess, omega_gel,
                               sigma_ess, sigma_gel,
                               predictor_labels = c("SAM", "ENSO", "SeaIce", "Whaling"),
                               support_threshold = 0.95) {
  
  cat("######################################################################\n")
  cat("Model:", name, "\n")
  cat("######################################################################\n\n")
  
  ## ---- Input data summary --------------------------------------------------
  cat("---- Input data ----\n")
  cat(sprintf("  Samples (rows of Y):        %d\n", nrow(Y)))
  cat(sprintf("  Response groups (cols of Y): %d\n", ncol(Y)))
  cat(sprintf("  Distribution:               %s\n",
              paste(mod$distr[, "family"], collapse = ", ")))
  cat(sprintf("  Fixed-effect formula:       %s\n",
              deparse(mod$XFormula)))
  cat(sprintf("  Random levels:              %s\n",
              paste(names(mod$ranLevels), collapse = ", ")))
  cat(sprintf("  Predictors in XData:        %s\n",
              paste(colnames(X), collapse = ", ")))
  cat("\n")
  
  cat("  Response groups included:\n")
  resp_names <- colnames(Y)
  print(resp_names)
  cat("\n")
  
  ## Response moments — useful sanity check for CLR vs PA
  cat("  Y column summary (mean, sd, min, max):\n")
  y_summary <- data.frame(
    mean = round(colMeans(Y, na.rm = TRUE), 4),
    sd   = round(apply(Y, 2, sd, na.rm = TRUE), 4),
    min  = round(apply(Y, 2, min, na.rm = TRUE), 4),
    max  = round(apply(Y, 2, max, na.rm = TRUE), 4)
  )
  print(y_summary)
  cat("\n")
  
  ## ---- MCMC settings -------------------------------------------------------
  cat("---- MCMC settings ----\n")
  cat(sprintf("  Chains:    %d\n", length(mod$postList)))
  cat(sprintf("  Samples per chain: %d\n", mod$samples))
  cat(sprintf("  Thinning:  %d\n", mod$thin))
  cat(sprintf("  Transient: %d\n", mod$transient))
  cat(sprintf("  Total post-burn-in iterations per chain: %d\n",
              mod$samples * mod$thin))
  cat("\n")
  
  ## ---- Convergence ---------------------------------------------------------
  cat("---- Convergence diagnostics ----\n")
  
  diag_block <- function(label, ess, gel) {
    cat(sprintf("  %s:\n", label))
    cat(sprintf("    ESS  — min = %s, median = %s, mean = %s, %% > 200 = %.1f%%\n",
                round(min(ess, na.rm = TRUE)),
                round(median(ess, na.rm = TRUE)),
                round(mean(ess, na.rm = TRUE)),
                mean(ess > 200, na.rm = TRUE) * 100))
    cat(sprintf("    PSRF — max = %.3f, median = %.3f, %% > 1.05 = %.1f%%\n",
                max(gel[, 1], na.rm = TRUE),
                median(gel[, 1], na.rm = TRUE),
                mean(gel[, 1] > 1.05, na.rm = TRUE) * 100))
  }
  
  diag_block("Beta (fixed effects)", beta_ess,  beta_gel)
  if (!is.null(omega_ess) && !is.null(omega_gel))
    diag_block("Omega (residual associations)", omega_ess, omega_gel)
  if (!is.null(sigma_ess) && !is.null(sigma_gel))
    diag_block("Sigma (residual variance)",     sigma_ess, sigma_gel)
  cat("\n")
  
  ## ---- Variance partitioning ----------------------------------------------
  cat("---- Variance partitioning ----\n")
  cat("Group means across response groups (proportion of variance):\n")
  print(round(rowMeans(VP$vals), 4))
  cat("\n")
  
  cat("Variance partitioning per response group:\n")
  vp_tab <- t(VP$vals)
  vp_tab <- round(vp_tab, 4)
  print(vp_tab)
  cat("\n")
  
  ## ---- Fixed-effect posterior summaries -----------------------------------
  cat(sprintf("---- Fixed-effect posterior summaries (support threshold = %.2f) ----\n",
              support_threshold))
  
  # Rows 2:4 = the three environmental predictors (row 1 = intercept).
  pred_rows  <- 2:(1 + length(predictor_labels))
  beta_mean  <- postBeta$mean[pred_rows, , drop = FALSE]
  P_pos      <- postBeta$support[pred_rows, , drop = FALSE]
  P_neg      <- postBeta$supportNeg[pred_rows, , drop = FALSE]
  rownames(beta_mean) <- rownames(P_pos) <- rownames(P_neg) <- predictor_labels
  
  support_df <- data.frame(
    response  = rep(colnames(Y), each = length(predictor_labels)),
    predictor = rep(predictor_labels, times = ncol(Y)),
    mean_beta = as.vector(beta_mean),
    P_pos     = as.vector(P_pos),
    P_neg     = as.vector(P_neg),
    stringsAsFactors = FALSE
  )
  support_df$direction <- ifelse(support_df$mean_beta > 0, "positive", "negative")
  support_df$strong    <- support_df$P_pos > support_threshold |
    support_df$P_neg > support_threshold
  support_df[, c("mean_beta", "P_pos", "P_neg")] <-
    round(support_df[, c("mean_beta", "P_pos", "P_neg")], 4)
  
  cat("Full table of posterior support for fixed effects:\n")
  print(support_df[, c("response", "predictor", "mean_beta",
                       "P_pos", "P_neg", "direction", "strong")],
        row.names = FALSE)
  cat("\n")
  
  strong_df <- support_df[support_df$strong, ]
  if (nrow(strong_df) > 0) {
    cat(sprintf("Strongly supported responses (P_pos or P_neg > %.2f):\n",
                support_threshold))
    print(strong_df[, c("response", "predictor", "mean_beta",
                        "P_pos", "P_neg", "direction")],
          row.names = FALSE)
  } else {
    cat(sprintf("No responses with posterior support > %.2f.\n",
                support_threshold))
  }
  cat("\n")
  
  invisible(support_df)
}


## ---- Write the file --------------------------------------------------------

out_file <- "summaryfiles/HMSC_results.txt"
sink(out_file)

cat("======================================================================\n")
cat("Hmsc model results: community ~ SAM + ENSO + SeaIce + whale + (1 | site)\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("======================================================================\n\n")

cat("Sites used as random effect:\n")
print(table(beta$site))
cat("\n")

cat("Predictor summary (XData):\n")
print(round(sapply(X, function(z)
  c(mean = mean(z, na.rm = TRUE),
    sd   = sd(z,   na.rm = TRUE),
    min  = min(z,  na.rm = TRUE),
    max  = max(z,  na.rm = TRUE))), 4))
cat("\n")


## One support threshold for all models
support_threshold <- 0.95

## 1. Subdivision-level CLR model
support_sub <- hmsc_summary_block(
  name      = "Subdivision-level CLR (m)",
  mod       = m,
  Y         = Y,
  X         = X,
  VP        = VP,
  postBeta  = postBeta,
  mpost     = mpost,
  beta_ess  = beta_ess,                   beta_gel  = beta_gel,
  omega_ess = effectiveSize(mpost$Omega[[1]]),
  omega_gel = gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf,
  sigma_ess = effectiveSize(mpost$Sigma),
  sigma_gel = gelman.diag(mpost$Sigma, multivariate = FALSE)$psrf,
  support_threshold = support_threshold
)

## 2. Metazoa family-level CLR model
support_met <- hmsc_summary_block(
  name      = "Metazoa family CLR (m.fam)",
  mod       = m.fam,
  Y         = Y.fam,
  X         = X,
  VP        = VP.fam,
  postBeta  = postBeta.fam,
  mpost     = mpost.fam,
  beta_ess  = beta_ess.fam,  beta_gel  = beta_gel.fam,
  omega_ess = omega_ess.fam, omega_gel = omega_gel.fam,
  sigma_ess = sigma_ess.fam, sigma_gel = sigma_gel.fam,
  support_threshold = support_threshold
)

## 3. Metazoa family-level presence/absence probit model
support_met_pa <- hmsc_summary_block(
  name      = "Metazoa family P/A probit (m.fam.pa)",
  mod       = m.fam.pa,
  Y         = Y.fam.pa,
  X         = X,
  VP        = VP.fam.pa,
  postBeta  = postBeta.fam.pa,
  mpost     = mpost.fam.pa,
  beta_ess  = beta_ess.fam.pa,  beta_gel  = beta_gel.fam.pa,
  omega_ess = omega_ess.fam.pa, omega_gel = omega_gel.fam.pa,
  sigma_ess = sigma_ess.fam.pa, sigma_gel = sigma_gel.fam.pa,
  support_threshold = support_threshold
)

## 4. Explanatory power per taxon (residual = 1 - R2)
cat("\n======================================================================\n")
cat("Explanatory power per taxon (residual = 1 - R2)\n")
cat("======================================================================\n\n")
for (blk in list(list("Subdivision CLR",        fit.sub$R2,     m$spNames),
                 list("Metazoa family CLR",     fit.fam$R2,     m.fam$spNames),
                 list("Metazoa family PA (Tjur R2)", fit.pa$TjurR2, m.fam.pa$spNames))) {
  r <- setNames(round(blk[[2]], 3), blk[[3]])
  cat(sprintf("---- %s: mean R2 = %.3f (range %.3f-%.3f); mean residual = %.3f ----\n",
              blk[[1]], mean(r, na.rm = TRUE), min(r, na.rm = TRUE),
              max(r, na.rm = TRUE), 1 - mean(r, na.rm = TRUE)))
  print(sort(r, decreasing = TRUE))
  cat("\n")
}
cat(sprintf("Metazoa family PA: mean AUC = %.3f\n",
            mean(fit.pa$AUC, na.rm = TRUE)))
if (any(is.na(fit.pa$TjurR2)))
  cat("NOTE: Tjur R2 is NA for taxa present in all samples; shown as NA above\n",
      "and plotted with zero-length bars in the PA figure.\n")

sink()


## ---- Also dump compact CSVs of fixed-effect support ------------------------

support_sub$model    <- "Subdivision_CLR"
support_met$model    <- "Metazoa_family_CLR"
support_met_pa$model <- "Metazoa_family_PA"

compact <- rbind(support_sub, support_met, support_met_pa)
compact <- compact[, c("model", "response", "predictor",
                       "mean_beta", "P_pos", "P_neg",
                       "direction", "strong")]
write.csv(compact, "summaryfiles/HMSC_fixedeffects_compact.csv",
          row.names = FALSE)

cat("Wrote: ", out_file, "\n",
    "Wrote:  summaryfiles/HMSC_fixedeffects_compact.csv\n", sep = "")

### build tables for output 

## ===== Supplementary tables: per-taxon beta, posterior support, and =========
## ===== variance shares (of total) for each HMSC model =======================
## One CSV per model: beta and posterior support (probability of the estimated
## sign) for each predictor, plus climate / whaling / site shares of TOTAL
## variance (VP proportion x taxon R2) and residual (1 - R2). Rows ordered by
## whaling beta (descending, matching the figures); mean row appended.
## PA model uses Tjur R2 for the conversion (approximation, as in Fig. S$);
## taxa present in all samples have undefined Tjur R2 and are reported NA.

hmsc_supp_table <- function(model, r2, outfile, display = NULL) {
  pb   <- getPostEstimate(model, "Beta")
  beta <- pb$mean;    supp <- pmax(pb$support, 1 - pb$support)
  rownames(beta) <- rownames(supp) <- model$covNames
  colnames(beta) <- colnames(supp) <- model$spNames
  vp <- computeVariancePartitioning(model, group = c(1, 2, 2, 2, 3),
                                    groupnames = c("Intercept", "Climate", "Whaling"))$vals
  colnames(vp) <- model$spNames
  r2 <- r2[model$spNames]
  
  preds <- model$covNames[grep("Intercept", model$covNames, invert = TRUE)]
  tab <- data.frame(taxon = if (is.null(display)) model$spNames
                    else display[model$spNames])
  for (p in preds) {
    tab[[paste0("beta_", p)]]    <- round(beta[p, ], 3)
    tab[[paste0("support_", p)]] <- round(supp[p, ], 3)
  }
  tab$var_climate_pct  <- round(100 * vp["Climate", ]      * r2, 1)
  tab$var_whaling_pct  <- round(100 * vp["Whaling", ]      * r2, 1)
  tab$var_site_pct     <- round(100 * vp["Random: site", ] * r2, 1)
  tab$var_residual_pct <- round(100 * (1 - r2), 1)
  
  k <- grep("whale", preds, ignore.case = TRUE, value = TRUE)
  tab <- tab[order(-tab[[paste0("beta_", k)]]), ]
  means <- c(taxon = "Mean",
             lapply(tab[-1], function(x) round(mean(x, na.rm = TRUE), 3)))
  tab <- rbind(tab, means)
  
  ## tripwire: variance columns must sum to ~100% for every taxon
  sums <- rowSums(tab[seq_len(nrow(tab) - 1),
                      c("var_climate_pct", "var_whaling_pct",
                        "var_site_pct", "var_residual_pct")])
  stopifnot(all(is.na(sums) | abs(sums - 100) < 0.5))
  
  write.csv(tab, outfile, row.names = FALSE)
  cat("wrote", outfile, "-", nrow(tab) - 1, "taxa\n")
  invisible(tab)
}

disp_tab <- setNames(paste0(disp, ifelse(disp_dag, "\u2020", "")),
                     names(disp))                       # dagger as plain char

hmsc_supp_table(m.fam,    setNames(fit.fam$R2, m.fam$spNames),
                "summaryfiles/TableS_HMSC_metazoa_CLR.csv",  disp_tab)
hmsc_supp_table(m,        setNames(fit.sub$R2, m$spNames),
                "summaryfiles/TableS_HMSC_subdivision_CLR.csv",
                setNames(gsub("_X", "", m$spNames), m$spNames))
hmsc_supp_table(m.fam.pa, setNames(fit.pa$TjurR2, m.fam.pa$spNames),
                "summaryfiles/TableS_HMSC_metazoa_PA.csv",   disp_tab)


## ===== Supplementary table: guild contrasts for every predictor =============
## Water-column feeders vs all other taxa, per predictor, per model. All
## modelled families included (no significance filtering). Difference and
## superiority computed within each posterior sample; Wilcoxon on posterior
## means as frequentist companion. Backs the results sentence and Fig. S$.

guild_contrast_table <- function(model, label) {
  wc   <- met.groups$group[match(model$spNames, met.groups$family)] %in%
    c("Pelagic water-column feeder", "Benthic suspension feeder")
  post <- Hmsc::poolMcmcChains(model$postList)
  bhat <- Hmsc::getPostEstimate(model, "Beta")$mean
  rownames(bhat) <- model$covNames
  preds <- grep("Intercept", model$covNames, invert = TRUE, value = TRUE)
  
  out <- do.call(rbind, lapply(preds, function(p) {
    k  <- match(p, model$covNames)
    dd <- sapply(post, function(s) mean(s$Beta[k, wc]) - mean(s$Beta[k, !wc]))
    ps <- sapply(post, function(s) {
      b <- s$Beta[k, ]; mean(outer(b[wc], b[!wc], ">")) })
    wt <- wilcox.test(bhat[p, wc], bhat[p, !wc])
    data.frame(model = label, predictor = p,
               n_water = sum(wc), n_other = sum(!wc),
               diff_mean   = round(mean(dd), 3),
               diff_lo95   = round(quantile(dd, .025), 3),
               diff_hi95   = round(quantile(dd, .975), 3),
               P_diff_gt0  = round(mean(dd > 0), 4),
               superiority = round(mean(ps), 3),
               sup_lo95    = round(quantile(ps, .025), 3),
               sup_hi95    = round(quantile(ps, .975), 3),
               wilcoxon_W  = unname(wt$statistic),
               wilcoxon_p  = signif(wt$p.value, 3))
  }))
  rownames(out) <- NULL
  out
}

tab <- rbind(guild_contrast_table(m.fam,    "Metazoa CLR"),
             guild_contrast_table(m.fam.pa, "Metazoa PA"))
write.csv(tab, "summaryfiles/TableS_guild_contrasts.csv", row.names = FALSE)

sink("summaryfiles/HMSC_guild_tests.txt", append = TRUE)
cat("\n==== Guild contrasts, all predictors (Table S$) ====\n")
print(tab, row.names = FALSE)
sink()


##### Now we make a composite figure for specific taxa

## First we pull in whales
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

###### SALPS KRILL
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



#####BUILD THE PLOT

margins <- c(2.1, 4.1, 1.1, 2.1)

loopdat <- whales.sum.nreps[whales.sum.nreps$spp == "Balaenoptera musculus", ]

pdf("figures/DraftFigures/Fig.4/whale.pdf", width = 12, height = 3)
par(xpd = NA,mar=margins)  
set.seed(1234)
plot(loopdat$dates, jitter(loopdat$value, 0.8)/16,
     xlim = c(1750, 2026), ylim = c(0, 1),
     col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.6),
     pch = ifelse(loopdat$value > 0, 16, 4),
     main = "", xlab = "", ylab = "",axes = FALSE)
axis(2,las=2)

for (s in unique(loopdat$site)) {
  sub <- loopdat[loopdat$site == s, ]
  sub <- sub[order(sub$dates), ]
  if (nrow(sub) > 4) {
    fit <- loess(value/16 ~ dates, data = sub, span = 0.8)
    lines(sub$dates, predict(fit),
          col = coredata$col[match(toupper(s), coredata$Site)], lwd = 2)
  }
}
dev.off()

loopdat <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp == "Euphausiidae", ]

pdf("figures/DraftFigures/Fig.4/krill.pdf", width = 12, height = 3)
par(xpd = NA,mar=margins)  
set.seed(1234)
plot(loopdat$dates, jitter(loopdat$value, 0.8)/8,
     xlim = c(1750, 2026), ylim = c(0, 1),
     col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.6),
     pch = ifelse(loopdat$value > 0, 16, 4),
     main = "", xlab = "", ylab = "", bty = "n",axes = FALSE)
axis(4,las=2)


for (s in unique(loopdat$site)) {
  sub <- loopdat[loopdat$site == s, ]
  sub <- sub[order(sub$dates), ]
  if (nrow(sub) > 4) {
    fit <- loess(value/8 ~ dates, data = sub, span = 0.8)
    lines(sub$dates, predict(fit),
          col = coredata$col[match(toupper(s), coredata$Site)], lwd = 2)
  }
}
dev.off()

loopdat <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp == "Salpa", ]

pdf("figures/DraftFigures/Fig.4/salpa.pdf", width = 12, height = 3)
par(xpd = NA,mar=margins)  
set.seed(1234)
plot(loopdat$dates, jitter(loopdat$value, 0.8)/8,
     xlim = c(1750, 2026), ylim = c(0, 1),
     col = add.alpha(coredata$col[match(toupper(loopdat$site), coredata$Site)], 0.6),
     pch = ifelse(loopdat$value > 0, 16, 4),
     main = "", xlab = "", ylab = "", bty = "n",axes = FALSE)
axis(2,las=2)

for (s in unique(loopdat$site)) {
  sub <- loopdat[loopdat$site == s, ]
  sub <- sub[order(sub$dates), ]
  if (nrow(sub) > 4) {
    fit <- loess(value/8 ~ dates, data = sub, span = 0.8)
    lines(sub$dates, predict(fit),
          col = coredata$col[match(toupper(s), coredata$Site)], lwd = 2)
  }
}
dev.off()



Dia.smol.ratio <- colSums(EUK.dia.3)/colSums(EUK.smolphy.2)

pdf("figures/DraftFigures/Fig.4/EUK.dia.Smol.ratio.pdf", width = 12, height = 3)
par(xpd = NA, mar = margins)

dates <- expandedDates$CalYear[match(names(Dia.smol.ratio), expandedDates$SampleName)]
sites <- expandedDates$Site[match(names(Dia.smol.ratio), expandedDates$SampleName)]
vals  <- log10(Dia.smol.ratio)

plot(dates, vals,
     col = coredata$col[match(toupper(sites), coredata$Site)],
     pch = 16,
     xlim = c(1750, 2026),
     xlab = "Date", ylab = "", bty = "n")

for (s in unique(sites)) {
  keep <- sites == s & is.finite(vals) & !is.na(dates)
  if (sum(keep) > 4) {
    ord <- order(dates[keep])
    d <- dates[keep][ord]
    v <- vals[keep][ord]
    fit <- loess(v ~ d, span = 0.8)
    lines(d, predict(fit),
          col = coredata$col[match(toupper(s), coredata$Site)], lwd = 2)
  }
}
dev.off()


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
       xlim=c(1760,2024),
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
       xlim=c(1760,2024), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d1$avg>0,16,4), col=colv1, bty="n")
  dev.off()
  
}

### Let's make some GAMS

for (sp in unique(euk.salkril.nReps.long$spp)){
  
  
  # --- data per species
  d1 <- euk.salkril.nReps.long[euk.salkril.nReps.long$spp==sp, ]
  d1 <- d1[!is.na(d1$dates), ]
  
  d2 <- cbind(d1,beta[match(paste(d1$site,d1$dates,sep=" "),paste(beta$site,round(beta$year),sep=" ")),])
  
  pdf(paste0("figures/test/EUK.big.",sp,".pdf"), width=9, height=6)
  par(mfrow=c(2,1), mar=c(3.5,4,2.5,1), oma=c(0,0,0,0), xaxs="i",xpd = TRUE)
  
  
  # site-color map
  sites <- sort(unique(toupper(d2$site)))
  cols  <- coredata$col[match(sites, coredata$Site)]
  colv1 <- cols[match(toupper(d2$site), sites)]
  
  # --- TOP: nREps
  set.seed(1234)
  plot(d2$dates, jitter(d2$value, 0.8),
       xlim=c(1760,2024),
       xlab="", ylab="nPCR reps", main=sp,
       pch=ifelse(d1$value>0,16,4), col=colv1, bty="n")
  #for (si in sites) {
  #  d2s <- d2[toupper(d2$site)==si, ]
  #  if (nrow(d2s) > 2) {
  #    fit <- loess(value ~ dates, data=d2s, span=0.6, degree=1)
  #    xx  <- seq(min(d2s$dates), max(d2s$dates), length.out=300)
  #    yy  <- predict(fit, data.frame(dates=xx))
  #    lines(xx, yy, col=cols[sites==si], lwd=2)
  # }
  #}
  legend(1750,-0.5, legend=sites, col=cols, pch=16, cex=0.6, bty="n")
  
  # --- MID: average proportion of reads
  plot(d2$dates, d2$avg,
       xlim=c(1760,2024), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d1$avg>0,16,4), col=colv1, bty="n")
  dev.off()
  
  
}



fit <- glm(cbind(value, 8 - value) ~ SAM_s + ENSO_s + SeaIce_s + whale_s + site,
           family = binomial, data = d2)
plot(fit)

summary(fit)



sp_list <- unique(euk.salkril.nReps.long$spp)

for (sp in sp_list) {
  pdf(paste0("figures/test/EUK.big.",sp,".pdf"), width=9, height=6)
  par(mfrow=c(2,1), mar=c(3.5,4,2.5,1), oma=c(0,0,0,0), xaxs="i")
  
  # --- data per species
  d1 <- d2[d2$spp==sp, ]
  
  # site-color map
  sites <- sort(unique(toupper(d1$site)))
  cols  <- coredata$col[match(sites, coredata$Site)]
  colv1 <- cols[match(toupper(d1$site), sites)]
  
  # --- TOP: nREps
  set.seed(1234)
  plot(d1$dates, jitter(d1$value, 0.8),
       xlim=c(1760,2024),
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
       xlim=c(1760,2024), 
       xlab="", ylab="Avg proportion",
       pch=ifelse(d1$avg>0,16,4), col=colv1, bty="n")
  dev.off()
  
}



plotall <- function(tax_filter,
                    euk_data = euk,
                    euk_taxonomy = euk.tax,
                    manual_csv = NULL,
                    group_name = "group",
                    expanded_dates = expandedDates,
                    core_data = coredata,
                    beta_data = beta,
                    out_prefix = NULL,
                    n_cols_data = 992) {
  
  
  plotall <- function(tax_filter,
                      euk_data = euk,
                      euk_taxonomy = euk.tax,
                      manual_csv = NULL,
                      group_name = "group",
                      expanded_dates = expandedDates,
                      core_data = coredata,
                      beta_data = beta,
                      out_prefix = NULL,
                      n_cols_data = 992) {
    
    tax_filter='family=="Euphausiidae"'
    euk_data = euk
    euk_taxonomy = euk.tax
    expanded_dates = expandedDates
    core_data = coredata
    beta_data = beta
    
    # --- 1. Subset by the taxonomy filter ---------------------------------
    # tax_filter should be a logical vector matching rows of euk_taxonomy
    sub <- euk_taxonomy[family=="Euphausiidae", ]
    euk_taxonomy[euk_taxonomy$family=="Euphausiidae",]
    
    
    # --- 1. Evaluate the filter expression inside euk_taxonomy ------------
    keep <- eval(str2lang(tax_filter), euk_taxonomy)
    
    sub_tax <- euk_taxonomy[keep, ]
    sub_tax <- sub_tax[!is.na(sub_tax[, 1]), ]
    sub     <- euk_data[match(sub_tax$X.1, rownames(euk_data)), ]
    sub     <- sub[!is.na(sub[, 1]), ]
    sub_tax <- euk_taxonomy[match(rownames(sub), euk_taxonomy$X.1), ]
    
    if (nrow(sub) == 0) stop("No rows match the taxonomy filter.")
    sub <- sub[!is.na(sub[, 1]), ]
    sub_tax <- euk_taxonomy[match(rownames(sub), euk_taxonomy$X.1), ]
    
    if (nrow(sub) == 0) {
      stop("No rows match the taxonomy filter.")
    }
    
    # --- 2. Either load manual annotation, or use raw taxonomy ------------
    # normalize by original library size
    
    ###################### WE ARE HERE!!!!!
    
    
    dat[, 1:n_cols_data] <- sweep(dat[, 1:n_cols_data], 2,
                                  colSums(euk_data), FUN = "/")
    dat[is.na(dat) | dat == "NaN"] <- 0
    species_col <- dat$Assignment
  } else {
    # fall back to using the taxonomy column that was filtered on
    dat <- sub
    dat[, 1:n_cols_data] <- sweep(dat[, 1:n_cols_data], 2,
                                  colSums(euk_data), FUN = "/")
    dat[is.na(dat)] <- 0
    # use the finest taxonomy level that is non-empty; default: species if available
    species_col <- if (!is.null(sub_tax$species)) sub_tax$species else sub_tax$genus
  }
  
  # --- 3. Aggregate by species ------------------------------------------
  datain <- cbind(dat[, 1:n_cols_data], manualSpecies = species_col)
  agg <- aggregate(. ~ manualSpecies, data = datain, FUN = sum)
  rownames(agg) <- agg$manualSpecies
  agg <- agg[, -1]
  agg <- agg[rownames(agg) != "", , drop = FALSE]
  
  if (nrow(agg) == 0) stop("Aggregation produced no rows.")
  
  # --- 4. Build nReps + averages ----------------------------------------
  sample_groups <- gsub("(.*).[0-9]$", "\\1", colnames(agg))
  nReps <- NrepsMaker(agg * 1e14, sample_groups)
  avr   <- average_by_reps(agg, sample_groups)
  
  # --- 5. Long-format (base R, no tidyverse) ----------------------------
  to_long <- function(mat, valname) {
    m <- as.matrix(mat)
    data.frame(
      spp    = rep(rownames(m), ncol(m)),
      column = rep(colnames(m), each = nrow(m)),
      v      = as.vector(m),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) |> setNames(c("spp", "column", valname))
  }
  
  nr_long  <- to_long(nReps, "value")
  avr_long <- to_long(avr,   "avg")
  
  long <- merge(nr_long, avr_long, by = c("spp", "column"), all.x = TRUE)
  
  long$site  <- expanded_dates$Site[match(long$column, expanded_dates$SampleName)]
  long$dates <- expanded_dates$CalYear[match(long$column, expanded_dates$SampleName)]
  
  # --- 6. Plot per species ----------------------------------------------
  tag <- if (!is.null(out_prefix)) out_prefix else group_name
  
  for (sp in unique(long$spp)) {
    d1 <- long[long$spp == sp & !is.na(long$dates), ]
    if (nrow(d1) == 0) next
    
    # optional join with beta (if provided)
    d2 <- if (!is.null(beta_data)) {
      cbind(d1, beta_data[match(paste(d1$site, d1$dates),
                                paste(beta_data$site, round(beta_data$year))), ])
    } else d1
    
    sites <- sort(unique(toupper(d2$site)))
    cols  <- core_data$col[match(sites, core_data$Site)]
    colv1 <- cols[match(toupper(d2$site), sites)]
    
    pdf(paste0("figures/EUK.", tag, ".", sp, ".pdf"), width = 9, height = 6)
    par(mfrow = c(2, 1), mar = c(3.5, 4, 2.5, 1), oma = c(0, 0, 0, 0),
        xaxs = "i", xpd = TRUE)
    
    # Top: nREps
    set.seed(1234)
    plot(d2$dates, jitter(d2$value, 0.8),
         xlim = c(1760, 2024),
         xlab = "", ylab = "nPCR reps", main = sp,
         pch = ifelse(d2$value > 0, 16, 4), col = colv1, bty = "n")
    legend(1750, -0.5, legend = sites, col = cols,
           pch = 16, cex = 0.6, bty = "n")
    
    # Bottom: avg proportion
    plot(d2$dates, d2$avg,
         xlim = c(1760, 2024),
         xlab = "", ylab = "Avg proportion",
         pch = ifelse(d2$avg > 0, 16, 4), col = colv1, bty = "n")
    
    dev.off()
  }
  
  invisible(long)
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

####### why we dont use a year smooth for the alpha diversity metrics?



#we make a new var in clim we can use to match
clim$longname <- expandedDates$SampleName[match(paste(clim$Site,clim$SampleID,sep="_"),paste(expandedDates$Site,expandedDates$SampleID,sep="_"))]

alpha <- data.frame("sample" = paste(clim$Site,clim$SampleID),
                    "site" = as.factor(clim$Site),
                    "year" =clim$CalYear,
                    "ASV_Q0" = ASV.alpha_avg_reps$q0_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q1" = ASV.alpha_avg_reps$q1_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "ASV_Q2" = ASV.alpha_avg_reps$q2_mean[match(clim$longname,alpha_avg_reps$replicate)],
                    "genus_Q0" = alpha_genus_avg_reps$q0_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q1" = alpha_genus_avg_reps$q1_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "genus_Q2" = alpha_genus_avg_reps$q2_mean[match(clim$longname,alpha_genus_avg_reps$replicate)],
                    "whale" = whales$biomass[match(round(clim$CalYear),whales$year)],
                    "whale_s" = scale(whales$biomass[match(round(clim$CalYear),whales$year)]),
                    "SAM" = clim$SAM_Index_3yr,
                    "ENSO" = clim$ENSO_3yr,
                    "SeaIce" = clim$SIC_mean,
                    "SAM_s" = scale(clim$SAM_Index_3yr),
                    "ENSO_s" = scale(clim$ENSO_3yr),
                    "SeaIce_s" = scale(clim$SIC_mean)
                    
)

# ASV models (with year)
a1_A0 <- gam(ASV_Q0 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a1_A1 <- gam(ASV_Q1 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a1_A2 <- gam(ASV_Q2 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")

# ASV models (no year)
a2_A0 <- gam(ASV_Q0 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a2_A1 <- gam(ASV_Q1 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a2_A2 <- gam(ASV_Q2 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")

# Genus models (with year)
a1_G0 <- gam(genus_Q0 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a1_G1 <- gam(genus_Q1 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a1_G2 <- gam(genus_Q2 ~ s(year, k = 10) + s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")

# Genus models (no year)
a2_G0 <- gam(genus_Q0 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a2_G1 <- gam(genus_Q1 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")
a2_G2 <- gam(genus_Q2 ~ s(site, bs = "re") + SAM_s + ENSO_s + SeaIce_s + whale_s, data = alpha, method = "REML")

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





#### HMSC models — Hellinger transformation ####

library(Hmsc)

## Beta coefficient gradient: red = negative, green = positive
beta.grad <- colorRampPalette(c("#c13d2a","#db964d","#e1dcb3","#6ab79f","#244e51"))(100)

# -- Data prep ----------------------------------------------------------------
comm.raw.t <- t(as.matrix(euk.reps))
comm.raw.t <- comm.raw.t[match(clim$longname, rownames(comm.raw.t)), ]

# Normalise on full community FIRST
comm.prop <- comm.raw.t / rowSums(comm.raw.t)

# High-confidence Subdivision assignments only (>=60%)
tax.filt <- euk.tax[euk.tax$Subdivision.1 >= 60 & !is.na(euk.tax$Subdivision), ]

# Intersect ASVs, subset, collapse to Subdivision
shared.asvs <- intersect(colnames(comm.prop), tax.filt$ASV)
comm.filt   <- comm.prop[, shared.asvs, drop = FALSE]
tax.filt    <- tax.filt[match(shared.asvs, tax.filt$ASV), ]
comm.sub    <- t(rowsum(t(comm.filt), group = tax.filt$Subdivision, na.rm = TRUE))

# Hellinger transformation
Y <- sqrt(comm.sub)

cat("Samples:", nrow(Y), "| Subdivisions:", ncol(Y), "\n")

# Predictors
X <- data.frame(
  SAM    = as.numeric(scale(beta$SAM_s)),
  ENSO   = as.numeric(scale(beta$ENSO_s)),
  SeaIce = as.numeric(scale(beta$SeaIce_s)),
  whale = as.numeric(scale(beta$whale_s))
)

# -- Subdivision model --------------------------------------------------------
studyDesign <- data.frame(site = factor(beta$site))
rL.site     <- HmscRandomLevel(units = levels(studyDesign$site))

m <- Hmsc(
  Y           = Y,
  XData       = X,
  XFormula    = ~ SAM + ENSO + SeaIce + whale,
  distr       = "normal",
  studyDesign = studyDesign,
  ranLevels   = list(site = rL.site)
)

m <- sampleMcmc(m,
                thin      = 10,
                samples   = 1000,
                transient = 5000,
                nChains   = 4,
                nParallel = 4,
                verbose   = 500
)

# -- Convergence --------------------------------------------------------------
mpost    <- convertToCodaObject(m)
beta_ess <- effectiveSize(mpost$Beta)
beta_gel <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

cat("ESS — Min:", round(min(beta_ess)),
    "| Median:", round(median(beta_ess)),
    "| % > 200:", round(mean(beta_ess > 200) * 100, 1), "%\n")
cat("G-R — Max:", round(max(beta_gel[, 1]), 3),
    "| % > 1.05:", round(mean(beta_gel[, 1] > 1.05) * 100, 1), "%\n")

omega_ess <- effectiveSize(mpost$Omega[[1]])
omega_gel <- gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf
cat("Omega ESS — Min:", round(min(omega_ess)),
    "| Median:", round(median(omega_ess)),
    "| % > 200:", round(mean(omega_ess > 200) * 100, 1), "%\n")
cat("Omega G-R — Max:", round(max(omega_gel[, 1]), 3),
    "| % > 1.05:", round(mean(omega_gel[, 1] > 1.05) * 100, 1), "%\n")

sigma_ess <- effectiveSize(mpost$Sigma)
sigma_gel <- gelman.diag(mpost$Sigma, multivariate = FALSE)$psrf
cat("Sigma ESS — Min:", round(min(sigma_ess)),
    "| Median:", round(median(sigma_ess)),
    "| % > 200:", round(mean(sigma_ess > 200) * 100, 1), "%\n")
cat("Sigma G-R — Max:", round(max(sigma_gel[, 1]), 3),
    "| % > 1.05:", round(mean(sigma_gel[, 1] > 1.05) * 100, 1), "%\n")

# -- Subdivision results ------------------------------------------------------
VP <- computeVariancePartitioning(m,
                                  group = c(1, 2, 2, 2, 3),
                                  groupnames = c("Intercept", "Climate", "Whaling"))
print(round(rowMeans(VP$vals), 3))
plotVariancePartitioning(m, VP = VP)

postBeta <- getPostEstimate(m, parName = "Beta")
support_df <- data.frame(
  Subdivision = rep(colnames(Y), each = 4),
  predictor   = rep(c("SAM", "ENSO", "SeaIce", "Whaling"), times = ncol(Y)),
  mean_beta   = c(rbind(postBeta$mean[2, ],       postBeta$mean[3, ],       postBeta$mean[4, ], postBeta$mean[5, ])),
  P_pos       = c(rbind(postBeta$support[2, ],    postBeta$support[3, ],    postBeta$support[4, ], postBeta$support[5, ])),
  P_neg       = c(rbind(postBeta$supportNeg[2, ], postBeta$supportNeg[3, ], postBeta$supportNeg[4, ], postBeta$supportNeg[5, ]))
)
support_df$direction <- ifelse(support_df$mean_beta > 0, "positive", "negative")
support_df$strong    <- support_df$P_pos > 0.95 | support_df$P_neg > 0.95
print(support_df[support_df$strong, c("Subdivision", "predictor", "mean_beta",
                                      "P_pos", "P_neg", "direction")])

# -- Subdivision heatmap ------------------------------------------------------
sup_pos <- postBeta$support[2:5, ]
sup_neg <- postBeta$supportNeg[2:5, ]
keep <- apply(sup_pos > 0.95 | sup_neg > 0.95, 2, any)

beta_mat <- postBeta$mean[2:5, keep]
sup_pos  <- sup_pos[, keep]
sup_neg  <- sup_neg[, keep]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")

sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

ord <- rev(order(colnames(beta_mat)))
beta_plot <- t(beta_mat[, ord])
sig_plot  <- t(sig_mat[, ord])

par(mar = c(5, 12, 3, 2))
n_sub  <- nrow(beta_plot)
n_pred <- ncol(beta_plot)

image(1:n_pred, 1:n_sub, t(beta_plot),
      col = beta.grad,
      zlim = c(-max(abs(beta_plot)), max(abs(beta_plot))),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:n_pred, labels = colnames(beta_plot), las = 1)
axis(2, at = 1:n_sub, labels = gsub("_X", "", rownames(beta_plot)), las = 2, cex.axis = 0.6)

for (i in 1:n_pred) {
  for (j in 1:n_sub) {
    if (sig_plot[j, i] != "") {
      text(i, j, sig_plot[j, i], cex = 0.7)
    }
  }
}
title("Subdivision responses to environmental drivers")

# =============================================================================
# -- Metazoa family-level HMSC ------------------------------------------------
# =============================================================================

tax.meta <- euk.tax[euk.tax$Subdivision == "Metazoa" &
                      !is.na(euk.tax$family) &
                      euk.tax$Family.1 == 100, ]
tax.meta <- tax.meta[tax.meta$family != "", ]
cat("ASVs retained:", nrow(tax.meta), "\n")
cat("Unique families:", length(unique(tax.meta$family)), "\n")
print(sort(table(tax.meta$family), decreasing = TRUE))

# Intersect with proportion matrix (normalised on full community)
shared.meta <- intersect(colnames(comm.prop), tax.meta$ASV)
comm.meta   <- comm.prop[, shared.meta, drop = FALSE]
tax.meta    <- tax.meta[match(shared.meta, tax.meta$ASV), ]

# Collapse to family
comm.fam <- t(rowsum(t(comm.meta), group = tax.meta$family, na.rm = TRUE))

# Prevalence filter
prev <- colMeans(comm.fam > 0)
cat("\nPrevalence distribution:\n")
print(round(sort(prev, decreasing = TRUE), 2))
comm.fam <- comm.fam[, prev > 0.20]
cat("Families after prevalence filter:", ncol(comm.fam), "\n")

# Hellinger transformation
Y.fam <- sqrt(comm.fam)

cat("Samples:", nrow(Y.fam), "| Families:", ncol(Y.fam), "\n")
cat("Any NA:", anyNA(Y.fam), "\n")

# -- Family model -------------------------------------------------------------
m.fam <- Hmsc(
  Y           = Y.fam,
  XData       = X,
  XFormula    = ~ SAM + ENSO + SeaIce + whale,
  distr       = "normal",
  studyDesign = studyDesign,
  ranLevels   = list(site = rL.site)
)
m.fam <- sampleMcmc(m.fam,
                    thin      = 100,
                    samples   = 1000,
                    transient = 10000,
                    nChains   = 4,
                    nParallel = 4,
                    verbose   = 500
)

# -- Convergence --------------------------------------------------------------
mpost.fam    <- convertToCodaObject(m.fam)
beta_ess.fam <- effectiveSize(mpost.fam$Beta)
beta_gel.fam <- gelman.diag(mpost.fam$Beta, multivariate = FALSE)$psrf
cat("ESS — Min:", round(min(beta_ess.fam)),
    "| Median:", round(median(beta_ess.fam)),
    "| % > 200:", round(mean(beta_ess.fam > 200) * 100, 1), "%\n")
cat("G-R — Max:", round(max(beta_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(beta_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

omega_ess.fam <- effectiveSize(mpost.fam$Omega[[1]])
omega_gel.fam <- gelman.diag(mpost.fam$Omega[[1]], multivariate = FALSE)$psrf
cat("Omega ESS — Min:", round(min(omega_ess.fam)),
    "| Median:", round(median(omega_ess.fam)),
    "| % > 200:", round(mean(omega_ess.fam > 200) * 100, 1), "%\n")
cat("Omega G-R — Max:", round(max(omega_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(omega_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

sigma_ess.fam <- effectiveSize(mpost.fam$Sigma)
sigma_gel.fam <- gelman.diag(mpost.fam$Sigma, multivariate = FALSE)$psrf
cat("Sigma ESS — Min:", round(min(sigma_ess.fam)),
    "| Median:", round(median(sigma_ess.fam)),
    "| % > 200:", round(mean(sigma_ess.fam > 200) * 100, 1), "%\n")
cat("Sigma G-R — Max:", round(max(sigma_gel.fam[, 1]), 3),
    "| % > 1.05:", round(mean(sigma_gel.fam[, 1] > 1.05) * 100, 1), "%\n")

# -- Family results -----------------------------------------------------------
VP.fam <- computeVariancePartitioning(m.fam,
                                      group = c(1, 2, 2, 2, 3),
                                      groupnames = c("Intercept", "Climate", "Whaling"))
print(round(rowMeans(VP.fam$vals), 3))
plotVariancePartitioning(m.fam, VP = VP.fam)

postBeta.fam <- getPostEstimate(m.fam, parName = "Beta")
support.fam <- data.frame(
  family    = rep(colnames(Y.fam), each = 4),
  predictor = rep(c("SAM", "ENSO", "SeaIce", "Whaling"), times = ncol(Y.fam)),
  mean_beta = c(rbind(postBeta.fam$mean[2, ],       postBeta.fam$mean[3, ],       postBeta.fam$mean[4, ], postBeta.fam$mean[5, ])),
  P_pos     = c(rbind(postBeta.fam$support[2, ],    postBeta.fam$support[3, ],    postBeta.fam$support[4, ], postBeta.fam$support[5, ])),
  P_neg     = c(rbind(postBeta.fam$supportNeg[2, ], postBeta.fam$supportNeg[3, ], postBeta.fam$supportNeg[4, ], postBeta.fam$supportNeg[5, ]))
)
support.fam$direction <- ifelse(support.fam$mean_beta > 0, "positive", "negative")
support.fam$strong    <- support.fam$P_pos > 0.95 | support.fam$P_neg > 0.95
print(support.fam[support.fam$strong, c("family", "predictor", "mean_beta",
                                        "P_pos", "P_neg", "direction")])

plotBeta(m.fam, post = postBeta.fam, param = "Support", supportLevel = 0.95,
         split = 0.4, spNamesNumbers = c(TRUE, FALSE))

# -- Family heatmap -----------------------------------------------------------
beta_mat <- postBeta.fam$mean[2:5, ]
rownames(beta_mat) <- c("SAM", "ENSO", "SeaIce", "WhaleBiomass")

sup_pos <- postBeta.fam$support[2:5, ]
sup_neg <- postBeta.fam$supportNeg[2:5, ]

keep <- apply(sup_pos > 0.95 | sup_neg > 0.95, 2, any)
beta_mat <- beta_mat[, keep]
sup_pos  <- sup_pos[, keep]
sup_neg  <- sup_neg[, keep]

cat("Retained:", sum(keep), "of", length(keep), "families\n")

sig_mat <- matrix("", nrow = 4, ncol = ncol(beta_mat))
sig_mat[sup_pos > 0.95 | sup_neg > 0.95] <- "*"
sig_mat[sup_pos > 0.99 | sup_neg > 0.99] <- "**"

ord <- rev(order(colnames(beta_mat)))
beta_plot <- t(beta_mat[, ord])
sig_plot  <- t(sig_mat[, ord])

par(mar = c(5, 12, 3, 2))
n_fam  <- nrow(beta_plot)
n_pred <- ncol(beta_plot)

image(1:n_pred, 1:n_fam, t(beta_plot),
      col = beta.grad,
      zlim = c(-max(abs(beta_plot)), max(abs(beta_plot))),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = 1:n_pred, labels = colnames(beta_plot), las = 1)
axis(2, at = 1:n_fam, labels = rownames(beta_plot), las = 2, cex.axis = 0.7)

for (i in 1:n_pred) {
  for (j in 1:n_fam) {
    if (sig_plot[j, i] != "") {
      text(i, j, sig_plot[j, i], cex = 0.7)
    }
  }
}
title("Metazoan family responses to environmental drivers")


nmds_names <- rownames(EUK.nMDS.j$points)
nmds_sites <- expandedDates$Site[match(nmds_names, expandedDates$SampleName)]
nmds_years <- expandedDates$CalYear[match(nmds_names, expandedDates$SampleName)]
nmds_cols <- coredata$col[match(toupper(nmds_sites), toupper(coredata$Site))]

# --- Jaccard NMDS with trajectories ---
pdf("figures/NMDS_ALL_jaccard_trajectories.pdf", height = 7, width = 7.5)

plot(EUK.nMDS.j$points[,1], 
     EUK.nMDS.j$points[,2],
     col = nmds_cols, pch = 16, xlab = "", ylab = "")

# Add trajectory lines per site
for (s in unique(nmds_sites)) {
  idx <- which(nmds_sites == s)
  ord <- idx[order(nmds_years[idx])]
  for (i in 1:(length(ord) - 1)) {
    arrows(EUK.nMDS.j$points[ord[i], 1], EUK.nMDS.j$points[ord[i], 2],
           EUK.nMDS.j$points[ord[i+1], 1], EUK.nMDS.j$points[ord[i+1], 2],
           col = nmds_cols[idx[1]], length = 0.08, lwd = 2)
  }
}
# Year labels
for (s in unique(nmds_sites)) {
  idx <- which(nmds_sites == s)
  ord <- idx[order(nmds_years[idx])]
  thin <- ord[seq(1, length(ord), by = 5)]
  thin <- union(thin, ord[c(1, length(ord))])  # always include endpoints
  text(EUK.nMDS.j$points[thin, 1], EUK.nMDS.j$points[thin, 2],
       labels = nmds_years[thin], pos = 3, cex = 0.6,)
}

dev.off()


# --- Bray-Curtis NMDS with trajectories ---
nmds_names_b <- rownames(EUK.nMDS.b$points)
nmds_sites_b <- expandedDates$Site[match(nmds_names_b, expandedDates$SampleName)]
nmds_years_b <- expandedDates$CalYear[match(nmds_names_b, expandedDates$SampleName)]
nmds_cols_b <- coredata$col[match(toupper(nmds_sites_b), toupper(coredata$Site))]

pdf("figures/NMDS_ALL_bray_trajectories.pdf", height = 7, width = 7.5)
plot(EUK.nMDS.b$points[,1], 
     EUK.nMDS.b$points[,2],
     col = nmds_cols_b, pch = 16, xlab = "", ylab = "")
# Add trajectory lines per site
for (s in unique(nmds_sites_b)) {
  idx <- which(nmds_sites_b == s)
  ord <- idx[order(nmds_years_b[idx])]
  for (i in 1:(length(ord) - 1)) {
    arrows(EUK.nMDS.b$points[ord[i], 1], EUK.nMDS.b$points[ord[i], 2],
           EUK.nMDS.b$points[ord[i+1], 1], EUK.nMDS.b$points[ord[i+1], 2],
           col = nmds_cols_b[idx[1]], length = 0.08, lwd = 2)
  }
}
# Year labels
for (s in unique(nmds_sites_b)) {
  idx <- which(nmds_sites_b == s)
  ord <- idx[order(nmds_years_b[idx])]
  thin <- ord[seq(1, length(ord), by = 5)]
  thin <- union(thin, ord[c(1, length(ord))])
  text(EUK.nMDS.b$points[thin, 1], EUK.nMDS.b$points[thin, 2],
       labels = nmds_years_b[thin], pos = 3, cex = 0.6)
}
dev.off()