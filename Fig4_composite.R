
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

### now we pull in the diatom vs flagellete comp

dia.key  <- read.csv("taxonomy/diatom_flagellate_key.csv", stringsAsFactors = FALSE)
dia.boot <- c(Class = 80, Order = 60, Family = 60, Genus = 60)
dia.cat  <- list()
for (col in c("category", "suppl_category")) {
  dia.cat[[col]] <- rep(NA, nrow(euk.pr2.2))
  for (rk in names(dia.boot)) {
    kk <- dia.key[dia.key$rank == rk, ]
    ok <- suppressWarnings(as.numeric(euk.pr2.2[[paste0(rk, ".1")]])) > dia.boot[rk] & euk.pr2.2[[rk]] %in% kk$taxon
    ok[is.na(ok)] <- FALSE
    if (rk != "Class") ok <- ok & !is.na(dia.cat[[col]])   # lower ranks only refine key classes
    dia.cat[[col]][ok] <- kk[[col]][match(euk.pr2.2[[rk]][ok], kk$taxon)]
  }
}

## everything comes from the key: category -> Fig 4c, suppl_category -> Fig S
EUK.dia.key <- euk.prop[rownames(euk.prop) %in% euk.pr2.2$X.1[dia.cat$category %in% "Diatom"], ]
EUK.fla     <- euk.prop[rownames(euk.prop) %in% euk.pr2.2$X.1[dia.cat$category %in% "Flagellate"], ]
EUK.pho     <- euk.prop[rownames(euk.prop) %in% euk.pr2.2$X.1[dia.cat$suppl_category %in% "Phototroph"], ]

## sample years and sites (same samples for every ratio)
rat.yr   <- expandedDates$CalYear[match(colnames(euk.prop), expandedDates$SampleName)]
rat.site <- expandedDates$Site[match(colnames(euk.prop), expandedDates$SampleName)]


#####BUILD THE PLOT

margins <- c(2.1, 4.1, 1.1, 2.1)

## Trend lines: one GAM per site, smoothness chosen by REML from that site's own data,
## with a 95% CI band that widens where samples are sparse.
## nPCR panels use a binomial GAM (positive reps out of 16 or 8), so fit and band stay within 0-1.
library(mgcv)

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
  sub <- sub[!is.na(sub$dates), ]
  if (nrow(sub) > 4 && sum(sub$value) > 0) {
    fit <- gam(cbind(value, 16 - value) ~ s(dates, k = min(10, length(unique(sub$dates)) - 2)),
               family = binomial, data = sub, method = "REML")
    nd  <- data.frame(dates = seq(min(sub$dates), max(sub$dates), length.out = 200))
    pr  <- predict(fit, nd, se.fit = TRUE)
    polygon(c(nd$dates, rev(nd$dates)), plogis(c(pr$fit + 1.96*pr$se.fit, rev(pr$fit - 1.96*pr$se.fit))),
            col = add.alpha(coredata$col[match(toupper(s), coredata$Site)], 0.2), border = NA)
    lines(nd$dates, plogis(pr$fit),
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
  sub <- sub[!is.na(sub$dates), ]
  if (nrow(sub) > 4 && sum(sub$value) > 0) {
    fit <- gam(cbind(value, 8 - value) ~ s(dates, k = min(10, length(unique(sub$dates)) - 2)),
               family = binomial, data = sub, method = "REML")
    nd  <- data.frame(dates = seq(min(sub$dates), max(sub$dates), length.out = 200))
    pr  <- predict(fit, nd, se.fit = TRUE)
    polygon(c(nd$dates, rev(nd$dates)), plogis(c(pr$fit + 1.96*pr$se.fit, rev(pr$fit - 1.96*pr$se.fit))),
            col = add.alpha(coredata$col[match(toupper(s), coredata$Site)], 0.2), border = NA)
    lines(nd$dates, plogis(pr$fit),
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
  sub <- sub[!is.na(sub$dates), ]
  if (nrow(sub) > 4 && sum(sub$value) > 0) {
    fit <- gam(cbind(value, 8 - value) ~ s(dates, k = min(10, length(unique(sub$dates)) - 2)),
               family = binomial, data = sub, method = "REML")
    nd  <- data.frame(dates = seq(min(sub$dates), max(sub$dates), length.out = 200))
    pr  <- predict(fit, nd, se.fit = TRUE)
    polygon(c(nd$dates, rev(nd$dates)), plogis(c(pr$fit + 1.96*pr$se.fit, rev(pr$fit - 1.96*pr$se.fit))),
            col = add.alpha(coredata$col[match(toupper(s), coredata$Site)], 0.2), border = NA)
    lines(nd$dates, plogis(pr$fit),
          col = coredata$col[match(toupper(s), coredata$Site)], lwd = 2)
  }
}
dev.off()


####----Fig 4c: diatom : flagellate----####
Dia.fla.ratio <- colSums(EUK.dia.key)/colSums(EUK.fla)
pdf("figures/DraftFigures/Fig.4/EUK.dia.fla.ratio.pdf",width = 12,height = 3)
par(xpd = NA,mar=margins)
plot(rat.yr, log10(Dia.fla.ratio),
     col=coredata$col[match(toupper(rat.site),toupper(coredata$Site))],
     pch=16,
     xlim = c(1750, 2026),
     xlab="Date",ylab="",bty="n")
for (s in unique(rat.site)) {
  keep <- rat.site == s & is.finite(log10(Dia.fla.ratio)) & !is.na(rat.yr)
  if (sum(keep) > 4) {
    d <- rat.yr[keep][order(rat.yr[keep])]
    v <- log10(Dia.fla.ratio)[keep][order(rat.yr[keep])]
    fit <- gam(v ~ s(d, k = min(10, length(unique(d)) - 2)), method = "REML")
    nd  <- data.frame(d = seq(min(d), max(d), length.out = 200))
    pr  <- predict(fit, nd, se.fit = TRUE)
    polygon(c(nd$d, rev(nd$d)), c(pr$fit + 1.96*pr$se.fit, rev(pr$fit - 1.96*pr$se.fit)),
            col = add.alpha(coredata$col[match(toupper(s), toupper(coredata$Site))], 0.2), border = NA)
    lines(nd$d, pr$fit, col = coredata$col[match(toupper(s), toupper(coredata$Site))], lwd = 2)
  }
}
dev.off()

####----Fig S: diatom : all non-diatom phototrophs (mixed groups excluded)----####
Dia.pho.ratio <- colSums(EUK.dia.key)/colSums(EUK.pho)
pdf("figures/DraftFigures/FigS_EUK.dia.pho.ratio.pdf",width = 12,height = 3)
par(xpd = NA,mar=margins)
plot(rat.yr, log10(Dia.pho.ratio),
     col=coredata$col[match(toupper(rat.site),toupper(coredata$Site))],
     pch=16,
     xlim = c(1750, 2026),
     xlab="Date",ylab="",bty="n")
for (s in unique(rat.site)) {
  keep <- rat.site == s & is.finite(log10(Dia.pho.ratio)) & !is.na(rat.yr)
  if (sum(keep) > 4) {
    d <- rat.yr[keep][order(rat.yr[keep])]
    v <- log10(Dia.pho.ratio)[keep][order(rat.yr[keep])]
    fit <- gam(v ~ s(d, k = min(10, length(unique(d)) - 2)), method = "REML")
    nd  <- data.frame(d = seq(min(d), max(d), length.out = 200))
    pr  <- predict(fit, nd, se.fit = TRUE)
    polygon(c(nd$d, rev(nd$d)), c(pr$fit + 1.96*pr$se.fit, rev(pr$fit - 1.96*pr$se.fit)),
            col = add.alpha(coredata$col[match(toupper(s), toupper(coredata$Site))], 0.2), border = NA)
    lines(nd$d, pr$fit, col = coredata$col[match(toupper(s), toupper(coredata$Site))], lwd = 2)
  }
}
dev.off()