### Whale analysis 

## helper: load, reshape, add pristine=1850
prep <- function(path) {
  m <- t(read.csv(path))
  df <- as.data.frame(m[-1,], stringsAsFactors = FALSE)
  colnames(df) <- m[1,]
  df$year <- as.numeric(gsub("X","", rownames(df)))
  df$`0.025` <- as.numeric(df$`0.025`)
  df$`0.5`   <- as.numeric(df$`0.5`)
  df$`0.975` <- as.numeric(df$`0.975`)
  df
}

df <- prep("rawdata/biomassdatafromMilo/Blue_whale__SH.csv")
add_pristine <- function(df, y0 = 1750) {
  stopifnot(all(c("year","0.025","0.5","0.975") %in% names(df)))
  df <- df[order(df$year), ]
  
  base_cols <- c("0.025","0.5","0.975","year")
  # enforce numeric types
  df$`0.025` <- as.numeric(df$`0.025`)
  df$`0.5`   <- as.numeric(df$`0.5`)
  df$`0.975` <- as.numeric(df$`0.975`)
  df$year    <- as.integer(df$year)
  
  first_year <- df$year[1]
  if (y0 < first_year) {
    n <- first_year - y0
    first_vals <- df[1, c("0.025","0.5","0.975")]
    pad <- data.frame(
      `0.025` = rep(first_vals$`0.025`, n),
      `0.5`   = rep(first_vals$`0.5`,   n),
      `0.975` = rep(first_vals$`0.975`, n),
      year    = seq.int(y0, first_year - 1L),
      check.names = FALSE
    )
    # build from an empty template to lock names/order
    out <- df[0, base_cols]
    out <- rbind(out, pad[, base_cols], df[, base_cols])
  } else {
    out <- df[, base_cols]
  }
  
  rownames(out) <- NULL
  out
}
add_pristine(df)


## data
blue <- add_pristine(prep("rawdata/biomassdatafromMilo/Blue_whale__SH.csv"))
orca <- add_pristine(prep("rawdata/biomassdatafromMilo/Killer_whale__SH.csv"))
hump <- add_pristine(prep("rawdata/biomassdatafromMilo/Humpback_whale__SH.csv"))
fin <- add_pristine(prep("rawdata/biomassdatafromMilo/Fin_whale__SH.csv"))
mink <- add_pristine(prep("rawdata/biomassdatafromMilo/Antarctic_minke_whale__SH.csv"))
sei <- add_pristine(prep("rawdata/biomassdatafromMilo/Sei_whale__SH.csv"))

## ribbon (only for real data, i.e., from first non-1850 year)
ribbon <- function(df, col,yr=1850) {
  i0 <- which(df$year != yr)[1]
  x  <- df$year[i0:nrow(df)]
  yL <- df$`0.025`[i0:nrow(df)]
  yU <- df$`0.975`[i0:nrow(df)]
  polygon(c(x, rev(x)), c(yL, rev(yU)), border = NA, col = adjustcolor(col, alpha.f = 0.18))
}

## plotting
pdf("figures/whale.biomass.pdf",height = 6.5,width = 9)
yl <- c(0, max(blue$`0.975`, orca$`0.975`, hump$`0.975`, fin$`0.975`, na.rm = TRUE))
plot(blue$year, blue$`0.5`, type = "n", ylim = yl, xlab = "Year", ylab = "Biomass")

## ribbons first 
ribbon(blue, "dodgerblue")
ribbon(orca, "red")
ribbon(hump, "grey")
ribbon(fin,  "purple")
ribbon(mink,  "darkgreen")
ribbon(sei,"orange3")

## lines as current (include 1850 point)
lines(blue$year, blue$`0.5`,lwd=2, col = "dodgerblue")
lines(blue$year, blue$`0.975`, lty = 2, col = "dodgerblue")
lines(blue$year, blue$`0.025`, lty = 2, col = "dodgerblue")

lines(orca$year, orca$`0.5`,lwd=2,   col = "red")
lines(orca$year, orca$`0.975`, lty = 2, col = "red")
lines(orca$year, orca$`0.025`, lty = 2, col = "red")

lines(hump$year, hump$`0.5`, lwd=2,  col = "grey")
lines(hump$year, hump$`0.975`, lty = 2, col = "grey")
lines(hump$year, hump$`0.025`, lty = 2, col = "grey")

lines(fin$year, fin$`0.5`, lwd=2,  col = "purple")
lines(fin$year, fin$`0.975`, lty = 2, col = "purple")
lines(fin$year, fin$`0.025`, lty = 2, col = "purple")

lines(mink$year, mink$`0.5`, lwd=2,  col = "darkgreen")
lines(mink$year, mink$`0.975`, lty = 2, col = "darkgreen")
lines(mink$year, mink$`0.025`, lty = 2, col = "darkgreen")

lines(sei$year, sei$`0.5`, lwd=2,  col = "orange3")
lines(sei$year, sei$`0.975`, lty = 2, col = "orange3")
lines(sei$year, sei$`0.025`, lty = 2, col = "orange3")

legend("topright",
       legend = c("Blue", "Orca", "Humpback", "Fin","Minke","Sei"),
       lty = 1, col = c("dodgerblue","red","grey","purple","darkgreen","orange3"), bty = "n")

dev.off()

#### Now lets pull the data out for statistics 
## here we use only blue, hump, fin & sei 

whale.out <- data.frame("year"=blue$year,
                        "biomass"=blue$`0.5`+fin$`0.5`+hump$`0.5`+sei$`0.5`)
plot(whale.out$year,whale.out$biomass,type="l",lwd=2)
write.csv(whale.out,"cleandata/whale.model.year.csv")


### what about a comparison to metabarcoding data?


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

whales.sum.nreps$spp_short <- NA
whales.sum.nreps$spp_short[whales.sum.nreps$spp == "Balaenoptera bonaerensis"]  <- "mink"
whales.sum.nreps$spp_short[whales.sum.nreps$spp == "Balaenoptera musculus"]     <- "blue"
whales.sum.nreps$spp_short[whales.sum.nreps$spp == "Balaenoptera physalus"]     <- "fin"
whales.sum.nreps$spp_short[whales.sum.nreps$spp == "Megaptera novaeangliae"]    <- "hump"
whales.sum.nreps$spp_short[whales.sum.nreps$spp == "Orcinus orca"]              <- "orca"


## link biomass to metabarcoding (median trajectory) ----------------------
mk <- function(df, key) {
  out <- df[, c("year","0.025","0.5","0.975")]
  out$spp_short <- key
  names(out) <- c("year","b025","b50","b975","spp_short")
  out
}

series <- rbind(
  mk(blue, "blue"),
  mk(orca, "orca"),
  mk(hump, "hump"),
  mk(fin,  "fin"),
  mk(mink, "mink")
)

whales.sum.nreps$biomass <- NA_real_
for (line in seq_len(nrow(whales.sum.nreps))) {
  loopdat <- whales.sum.nreps$dates[line]
  loopspp <- whales.sum.nreps$spp_short[line]
  idx <- which(series$year == loopdat & series$spp_short == loopspp)
  if (length(idx) > 0) whales.sum.nreps$biomass[line] <- series$b50[idx[1]]
}

## model covariate scaling ------------------------------------------------
whales.sum.nreps <- whales.sum.nreps %>%
  filter(!is.na(biomass),
         !is.na(value),
         !is.na(dates),
         !is.na(site),
         !is.na(spp_short))

whales.sum.nreps$biomass_z <- as.numeric(
  scale(log1p(whales.sum.nreps$biomass))
)
whales.sum.nreps$dates_z <- as.numeric(
  scale(whales.sum.nreps$dates)
)

whales.sum.nreps$spp_short <- factor(whales.sum.nreps$spp_short)
whales.sum.nreps$site      <- factor(whales.sum.nreps$site)

## fit negative binomial mixed model --------------------------------------
m_nb <- glmmTMB(
  value ~ biomass_z * dates_z + (1 | spp_short) + (1 | site),
  family = nbinom2,
  data   = whales.sum.nreps
)

summary(m_nb)

## prediction grid per species --------------------------------------------
d <- whales.sum.nreps   # shorthand

bmin <- tapply(d$biomass, d$spp_short, min, na.rm = TRUE)
bmax <- tapply(d$biomass, d$spp_short, max, na.rm = TRUE)

grid_list <- lapply(levels(d$spp_short), function(sp) {
  data.frame(
    spp_short = sp,
    biomass   = seq(bmin[[sp]], bmax[[sp]], length.out = 120)
  )
})
newdat <- do.call(rbind, grid_list)

newdat$biomass_z <- (log1p(newdat$biomass) -
                       mean(log1p(d$biomass), na.rm = TRUE)) /
  sd(log1p(d$biomass), na.rm = TRUE)

med_year <- median(d$dates, na.rm = TRUE)
newdat$dates_z <- (med_year - mean(d$dates, na.rm = TRUE)) /
  sd(d$dates, na.rm = TRUE)

newdat$site <- NA

p <- predict(m_nb, newdata = newdat, type = "response",
             se.fit = TRUE, re.form = NA)
newdat$fit   <- as.numeric(p$fit)
newdat$lower <- pmax(0, newdat$fit - 1.96 * p$se.fit)
newdat$upper <- newdat$fit + 1.96 * p$se.fit

## log10(y+1) plot with raw data + fitted curves --------------------------
d$val_p1        <- d$value + 1
newdat$fit_p1   <- newdat$fit   + 1
newdat$low_p1   <- pmax(1e-9, newdat$lower + 1)
newdat$upp_p1   <- newdat$upper + 1

ggplot() +
  geom_point(data = d,
             aes(x = biomass, y = val_p1),
             alpha = 0.4, size = 1.6, color = "grey30") +
  geom_ribbon(data = newdat,
              aes(x = biomass, ymin = low_p1, ymax = upp_p1),
              alpha = 0.2, fill = "steelblue", color = NA) +
  geom_line(data = newdat,
            aes(x = biomass, y = fit_p1),
            color = "steelblue", linewidth = 1.1) +
  facet_wrap(~ spp_short, scales = "free_x") +
  scale_x_continuous(labels = comma, trans = "log10") +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(0:6),
                     labels = function(z) comma(z - 1)) +
  labs(x = "Modeled biomass (log10 scale)",
       y = "Metabarcoding counts (log10(y+1))") +
  theme_classic()



