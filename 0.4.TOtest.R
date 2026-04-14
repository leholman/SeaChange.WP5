###############################################################################
#  Whale biomass analysis
#  – reads modelled biomass trajectories (from Milo)
#  – compares with metabarcoding eDNA data (MAM + RIZ markers)
#  – fits a negative-binomial GLMM (biomass ~ eDNA counts)
###############################################################################

library(glmmTMB)
coredata <- read.csv("metadata/sitemetadata.csv")

coredata$col <- c("#D55E00","#cc79a7","#0072B2","#f0e442","#009e73")

# ── 1. Helper functions ─────────────────────────────────────────────────────

#' Read one of Milo's biomass CSVs (transposed layout) and return a tidy
#' data.frame with columns: 0.025, 0.5, 0.975, year
prep <- function(path) {
  m  <- t(read.csv(path))
  df <- as.data.frame(m[-1, ], stringsAsFactors = FALSE)
  colnames(df) <- m[1, ]
  df$year    <- as.numeric(gsub("X", "", rownames(df)))
  df$`0.025` <- as.numeric(df$`0.025`)
  df$`0.5`   <- as.numeric(df$`0.5`)
  df$`0.975` <- as.numeric(df$`0.975`)
  df
}

#' Pad a biomass data.frame backwards to a "pristine" baseline year (y0).
#' Values before the first real year are held constant at that first value.
add_pristine <- function(df, y0 = 1750) {
  stopifnot(all(c("year", "0.025", "0.5", "0.975") %in% names(df)))
  df <- df[order(df$year), ]
  
  base_cols <- c("0.025", "0.5", "0.975", "year")
  df$`0.025` <- as.numeric(df$`0.025`)
  df$`0.5`   <- as.numeric(df$`0.5`)
  df$`0.975` <- as.numeric(df$`0.975`)
  df$year    <- as.integer(df$year)
  
  first_year <- df$year[1]
  
  if (y0 < first_year) {
    n   <- first_year - y0
    val <- df[1, c("0.025", "0.5", "0.975")]
    pad <- data.frame(
      `0.025` = rep(val$`0.025`, n),
      `0.5`   = rep(val$`0.5`,   n),
      `0.975` = rep(val$`0.975`, n),
      year    = seq.int(y0, first_year - 1L),
      check.names = FALSE
    )
    out <- rbind(pad[, base_cols], df[, base_cols])
  } else {
    out <- df[, base_cols]
  }
  
  rownames(out) <- NULL
  out
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

# ── 1.1 Data in ────────────────────────────────────────────
expandedDates <- read.csv("metadata/expandeddates.csv",row.names = 1)

# ── 2. Load biomass trajectories ────────────────────────────────────────────

datadir <- "rawdata/biomassdatafromMilo"

blue <- add_pristine(prep(file.path(datadir, "Blue_whale__SH.csv")))
orca <- add_pristine(prep(file.path(datadir, "Killer_whale__SH.csv")))
hump <- add_pristine(prep(file.path(datadir, "Humpback_whale__SH.csv")))
fin  <- add_pristine(prep(file.path(datadir, "Fin_whale__SH.csv")))
mink <- add_pristine(prep(file.path(datadir, "Antarctic_minke_whale__SH.csv")))
sei  <- add_pristine(prep(file.path(datadir, "Sei_whale__SH.csv")))


# ── 3. Plot biomass trajectories ────────────────────────────────────────────

# Draw a 95 % credible-interval ribbon, starting from the first "real"
# (non-padded) year.
ribbon <- function(df, col, pad_year = 1850) {
  i0 <- which(df$year != pad_year)[1]
  x  <- df$year[i0:nrow(df)]
  yL <- df$`0.025`[i0:nrow(df)]
  yU <- df$`0.975`[i0:nrow(df)]
  polygon(c(x, rev(x)), c(yL, rev(yU)),
          border = NA, col = adjustcolor(col, alpha.f = 0.18))
}

# Convenience: plot median + 95 % dashed lines for one species
add_lines <- function(df, col) {
  lines(df$year, df$`0.5`,   lwd = 2, col = col)
  lines(df$year, df$`0.975`, lty = 2, col = col)
  lines(df$year, df$`0.025`, lty = 2, col = col)
}

# Colour palette (one per species)
spp_cols <- c(blue = "dodgerblue", orca = "red", hump = "grey",
              fin = "purple", mink = "darkgreen", sei = "orange3")

pdf("figures/whale.biomass.pdf", height = 6.5, width = 9)

yl <- c(0, max(blue$`0.975`, orca$`0.975`, hump$`0.975`,
               fin$`0.975`, na.rm = TRUE))

plot(blue$year, blue$`0.5`, type = "n", ylim = yl,
     xlab = "Year", ylab = "Biomass")

# Ribbons first (so lines sit on top)
ribbon(blue, spp_cols["blue"])
ribbon(orca, spp_cols["orca"])
ribbon(hump, spp_cols["hump"])
ribbon(fin,  spp_cols["fin"])
ribbon(mink, spp_cols["mink"])
ribbon(sei,  spp_cols["sei"])

# Median + CI lines
add_lines(blue, spp_cols["blue"])
add_lines(orca, spp_cols["orca"])
add_lines(hump, spp_cols["hump"])
add_lines(fin,  spp_cols["fin"])
add_lines(mink, spp_cols["mink"])
add_lines(sei,  spp_cols["sei"])

legend("topright",
       legend = c("Blue", "Orca", "Humpback", "Fin", "Minke", "Sei"),
       col = spp_cols, lty = 1, lwd = 2, bty = "n")

dev.off()


# ── 4. Combined biomass for modelling ───────────────────────────────────────
# Sum median biomass of the four baleen species used in the model

whale_out <- data.frame(
  year    = blue$year,
  biomass = blue$`0.5` + fin$`0.5` + hump$`0.5` + sei$`0.5`
)

plot(whale_out$year, whale_out$biomass, type = "l", lwd = 2,
     xlab = "Year", ylab = "Summed biomass (blue + fin + hump + sei)")

write.csv(whale_out, "cleandata/whale.model.year.csv")


# ── 5. Extract whale ASVs from metabarcoding data ──────────────────────────
# Two markers: MAM and RIZ.  Keep families Balaenopteridae, Balaenidae,
# Delphinidae, then manually curate species assignments.

whale_families <- c("Balaenopteridae", "Balaenidae", "Delphinidae")

# --- MAM marker ---
mam_otu <- read.csv("cleandata/MAM.WP5.main.clean.csv",     row.names = 1)
mam_tax <- read.csv("cleandata/MAM.WP5.main.clean.tax.csv", row.names = 1)

mam_whale_asvs <- na.omit(
  mam_tax$ASV[mam_tax$family %in% whale_families]
)
mam_whales <- cbind(
  mam_otu[match(mam_whale_asvs, rownames(mam_otu)), ],
  mam_tax[match(mam_whale_asvs, mam_tax$ASV), ]
)

# --- RIZ marker ---
riz_otu <- read.csv("cleandata/riz.WP5.main.clean.csv",     row.names = 1)
riz_tax <- read.csv("cleandata/riz.WP5.main.clean.tax.csv", row.names = 1)

riz_whale_asvs <- na.omit(
  riz_tax$ASV[riz_tax$family %in% whale_families]
)
riz_whales <- cbind(
  riz_otu[match(riz_whale_asvs, rownames(riz_otu)), ],
  riz_tax[match(riz_whale_asvs, riz_tax$ASV), ]
)

write.csv(mam_whales, "cleandata/MAM.whales.clean.csv")
write.csv(riz_whales, "cleandata/RIZ.whales.clean.csv")


# ── 6. Read back manually curated whale tables ─────────────────────────────

mam_whales_man <- read.csv("cleandata/MAM.whales.clean.manual.csv", row.names = 1)
riz_whales_man <- read.csv("cleandata/RIZ.whales.clean.manual.csv", row.names = 1)

# Convert counts to proportions of total reads per sample
n_mam <- 983   # number of sample columns in MAM
n_riz <- 813   # number of sample columns in RIZ

mam_whales_man[1:n_mam] <- sweep(mam_whales_man[1:n_mam], 2,
                                 colSums(mam_otu), FUN = "/")
mam_whales_man[is.na(mam_whales_man) | mam_whales_man == "NaN"] <- 0

riz_whales_man[1:n_riz] <- sweep(riz_whales_man[1:n_riz], 2,
                                 colSums(riz_otu), FUN = "/")
riz_whales_man[is.na(riz_whales_man) | riz_whales_man == "NaN"] <- 0


# ── 7. Aggregate to species level ──────────────────────────────────────────

aggregate_species <- function(dat, n_samp, prefix_re) {
  spp_col <- data.frame(manualSpecies = dat$ManualSpecies)
  tmp     <- cbind(dat[, 1:n_samp], spp_col)
  agg     <- aggregate(. ~ manualSpecies, data = tmp, FUN = sum)
  rownames(agg) <- agg$manualSpecies
  agg <- agg[, -1]
  agg <- agg[rownames(agg) != "", ]
  colnames(agg) <- gsub(prefix_re, "", colnames(agg))
  agg
}

mam_whales_sum <- aggregate_species(mam_whales_man, n_mam, "^MAM[._]")
riz_whales_sum <- aggregate_species(riz_whales_man, n_riz, "^RIZ[._]")


# ── 8. Merge the two markers ───────────────────────────────────────────────
# NrepsMaker collapses technical replicates (assumed available in workspace).
# Scale proportions up to pseudo-counts before merging.

scale_factor <- 1e14   # large multiplier so NrepsMaker sees integers

strip_rep <- function(x) gsub("(.*).[0-9]$", "\\1", x)

riz_mat <- as.matrix(
  NrepsMaker(riz_whales_sum * scale_factor, strip_rep(colnames(riz_whales_sum)))
)
mam_mat <- as.matrix(
  NrepsMaker(mam_whales_sum * scale_factor, strip_rep(colnames(mam_whales_sum)))
)

# Pad both to the same dimensions
all_species <- union(rownames(riz_mat), rownames(mam_mat))
all_samples <- union(colnames(riz_mat), colnames(mam_mat))

pad_matrix <- function(mat, rn, cn) {
  out <- matrix(0, nrow = length(rn), ncol = length(cn),
                dimnames = list(rn, cn))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

whales_merged <- as.data.frame(
  pad_matrix(riz_mat, all_species, all_samples) +
    pad_matrix(mam_mat, all_species, all_samples)
)


# ── 9. Reshape to long format (base R) ─────────────────────────────────────

whales_long <- data.frame(
  spp    = rep(rownames(whales_merged), each = ncol(whales_merged)),
  column = rep(colnames(whales_merged), times = nrow(whales_merged)),
  value  = as.vector(t(as.matrix(whales_merged))),
  stringsAsFactors = FALSE
)

# Match sample metadata (site, calibrated year)
sample_key <- gsub("^EUK[._]", "", expandedDates$SampleName)

whales_long$site  <- expandedDates$Site[match(whales_long$column, sample_key)]
whales_long$dates <- expandedDates$CalYear[match(whales_long$column, sample_key)]

# Short species labels
spp_lookup <- c(
  "Balaenoptera bonaerensis" = "mink",
  "Balaenoptera musculus"    = "blue",
  "Balaenoptera physalus"    = "fin",
  "Megaptera novaeangliae"   = "hump",
  "Orcinus orca"             = "orca"
)
whales_long$spp_short <- spp_lookup[whales_long$spp]


# ── 10. Link modelled biomass to eDNA rows ──────────────────────────────────

# Build a single long table of biomass trajectories
mk <- function(df, key) {
  out <- df[, c("year", "0.025", "0.5", "0.975")]
  out$spp_short <- key
  names(out) <- c("year", "b025", "b50", "b975", "spp_short")
  out
}

series <- rbind(
  mk(blue, "blue"), mk(orca, "orca"), mk(hump, "hump"),
  mk(fin,  "fin"),  mk(mink, "mink")
)

# Match on year × species
whales_long$biomass <- NA_real_
for (i in seq_len(nrow(whales_long))) {
  idx <- which(series$year      == whales_long$dates[i] &
                 series$spp_short == whales_long$spp_short[i])
  if (length(idx) > 0) whales_long$biomass[i] <- series$b50[idx[1]]
}


# ── 11. Fit negative-binomial GLMM ─────────────────────────────────────────

# Drop incomplete cases
whales_long <- whales_long[
  complete.cases(whales_long[, c("biomass", "value", "dates",
                                 "site", "spp_short")]), ]

# Centre & scale covariates
whales_long$biomass_z <- as.numeric(scale(log1p(whales_long$biomass)))
whales_long$dates_z   <- as.numeric(scale(whales_long$dates))

whales_long$spp_short <- factor(whales_long$spp_short)
whales_long$site      <- factor(whales_long$site)

whales_long_noBlue <- whales_long[whales_long$spp_short!="blue",]

m_nb <- glmmTMB(
  value ~ biomass_z + dates_z + spp_short + (1 | site),
  family = poisson,
  data   = whales_long
)

summary(m_nb)


## Output

# ── 11b. Export model summary to text file ──────────────────────────────────

sink("figures/whalecomp/model_summary.txt")

cat("================================================================\n")
cat("  Poisson GLMM: eDNA ~ modelled whale biomass\n")
cat("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("================================================================\n\n")

# ── Data summary ──
cat("── Data summary ──\n\n")
cat("Total observations:", nrow(whales_long), "\n")
cat("Sites:", nlevels(whales_long$site),
    paste0("(", paste(levels(whales_long$site), collapse = ", "), ")"), "\n")
cat("Species:", nlevels(whales_long$spp_short),
    paste0("(", paste(levels(whales_long$spp_short), collapse = ", "), ")"), "\n")
cat("Year range:", paste(range(whales_long$dates), collapse = " – "), "\n")
cat("Response (value) range:", paste(range(whales_long$value), collapse = " – "), "\n\n")

# ── Scaling parameters (needed to back-transform) ──
cat("── Scaling parameters ──\n\n")
cat("biomass_z = (log1p(biomass) - mu) / sd\n")
cat("  mu =", round(mean(log1p(whales_long$biomass)), 4), "\n")
cat("  sd =", round(sd(log1p(whales_long$biomass)), 4), "\n\n")
cat("dates_z = (year - mu) / sd\n")
cat("  mu =", round(mean(whales_long$dates), 4), "\n")
cat("  sd =", round(sd(whales_long$dates), 4), "\n\n")

# ── Full model summary ──
cat("── Model summary ──\n\n")
print(summary(m_nb))
cat("\n")

# ── Fixed effects table ──
cat("── Fixed effects (exponentiated = rate ratios) ──\n\n")
fe <- fixef(m_nb)$cond
se <- sqrt(diag(vcov(m_nb)$cond))
ci_lo <- fe - 1.96 * se
ci_hi <- fe + 1.96 * se

fe_table <- data.frame(
  Estimate     = round(fe, 4),
  SE           = round(se, 4),
  z            = round(fe / se, 3),
  p            = formatC(2 * pnorm(-abs(fe / se)), format = "e", digits = 2),
  exp_Est      = round(exp(fe), 4),
  exp_CI_lower = round(exp(ci_lo), 4),
  exp_CI_upper = round(exp(ci_hi), 4)
)
print(fe_table)
cat("\n")

# ── Random effects ──
cat("── Random effects ──\n\n")
vc <- VarCorr(m_nb)$cond
cat("Site intercept variance:", round(as.numeric(attr(vc$site, "stddev"))^2, 4), "\n")
cat("Site intercept SD:      ", round(as.numeric(attr(vc$site, "stddev")), 4), "\n\n")

# ── Overdispersion check ──
cat("── Overdispersion check ──\n\n")
E2 <- residuals(m_nb, type = "pearson")
disp <- sum(E2^2) / df.residual(m_nb)
cat("Pearson dispersion ratio:", round(disp, 4), "\n")
cat("(1 = no overdispersion; >1.5 suggests concern)\n\n")

# ── AIC / logLik ──
cat("── Model fit ──\n\n")
cat("AIC:   ", round(AIC(m_nb), 1), "\n")
cat("BIC:   ", round(BIC(m_nb), 1), "\n")
cat("logLik:", round(as.numeric(logLik(m_nb)), 1), "\n")

cat("\n================================================================\n")
cat("  End of report\n")
cat("================================================================\n")

sink()



# ── 12. Prediction grid & plot ──────────────────────────────────────────────
d <- whales_long

# --- prediction data (same as your code) ---
bmin <- tapply(d$biomass, d$spp_short, min, na.rm = TRUE)
bmax <- tapply(d$biomass, d$spp_short, max, na.rm = TRUE)
newdat <- do.call(rbind, lapply(levels(d$spp_short), function(sp) {
  data.frame(spp_short = sp,
             biomass   = seq(bmin[[sp]], bmax[[sp]], length.out = 120))
}))

mu_b  <- mean(log1p(d$biomass), na.rm = TRUE)
sd_b  <- sd(log1p(d$biomass),   na.rm = TRUE)
newdat$biomass_z <- (log1p(newdat$biomass) - mu_b) / sd_b

med_year <- median(d$dates, na.rm = TRUE)
newdat$dates_z <- (med_year - mean(d$dates, na.rm = TRUE)) /
  sd(d$dates, na.rm = TRUE)
newdat$site <- NA

p <- predict(m_nb, newdata = newdat, type = "response",
             se.fit = TRUE, re.form = NA)
newdat$fit   <- as.numeric(p$fit)
newdat$lower <- pmax(0, newdat$fit - 1.96 * p$se.fit)
newdat$upper <- newdat$fit + 1.96 * p$se.fit

# --- base R faceted plot ---
# site colours from coredata
idx <- match(tolower(levels(d$site)), tolower(coredata$Site))
site_cols <- setNames(coredata$col[idx], levels(d$site))
spp   <- sort(unique(d$spp_short))
sites <- sort(unique(d$site))
n_spp <- length(spp)

pdf("figures/whalecomp/modelout.pdf",width=9,height=5.5)
par(mfrow = c(2, ceiling(n_spp / 2)), mar = c(2, 2, 2, 1), oma = c(2, 2, 0, 0))

for (s in spp) {
  dd <- d[d$spp_short == s, ]
  xlim <- range(dd$biomass[dd$biomass > 0], na.rm = TRUE)
  ylim <- c(0, max(dd$value, na.rm = TRUE))
  
  plot(1, type = "n", log = "x",
       xlim = xlim, ylim = ylim,
       xlab = "", ylab = "", main = s)
  
  for (si in sites) {
    dds <- dd[dd$site == si, ]
    if (nrow(dds) == 0) next
    
    # prediction for this site
    nd <- data.frame(
      spp_short = s,
      biomass   = seq(xlim[1], xlim[2], length.out = 120),
      site      = si
    )
    nd$biomass_z <- (log1p(nd$biomass) - mu_b) / sd_b
    nd$dates_z   <- (med_year - mean(d$dates, na.rm = TRUE)) /
      sd(d$dates, na.rm = TRUE)
    
    p  <- predict(m_nb, newdata = nd, type = "response",
                  se.fit = TRUE, re.form = NULL)
    nd$fit   <- as.numeric(p$fit)
    nd$lower <- pmax(0, nd$fit - 1.96 * p$se.fit)
    nd$upper <- nd$fit + 1.96 * p$se.fit
    
    ox <- order(nd$biomass)
    
    # confidence ribbon
    col_rgb <- col2rgb(site_cols[si]) / 255
    polygon(c(nd$biomass[ox], rev(nd$biomass[ox])),
            c(nd$lower[ox], rev(nd$upper[ox])),
            col = rgb(col_rgb[1], col_rgb[2], col_rgb[3], 0.15),
            border = NA)
    
    # fitted line
    lines(nd$biomass[ox], nd$fit[ox], col = site_cols[si], lwd = 2)
    
    # raw points
    points(dds$biomass, dds$value,
           pch = 16, col = adjustcolor(site_cols[si], 0.5), cex = 0.8)
  }
}

mtext("Modelled biomass (Tonnes)", side = 1, outer = TRUE, line = 0.5)
mtext("Positive PCR replicates", side = 2, outer = TRUE, line = 0.5)

# legend in last panel
plot.new()
legend("center", legend = names(site_cols), col = site_cols,
       lwd = 2, pch = 16, bty = "n", title = "Site")
dev.off()


# ── 13. Time x whale/metabarcoding plot ──────────────────────────────────────

# Map species short names to biomass dataframes
biomass_list <- list(blue = blue, orca = orca, hump = hump, fin = fin, mink = mink)

spp   <- sort(unique(d$spp_short))
sites <- sort(unique(d$site))
n_spp <- length(spp)

pdf("figures/whalecomp/time_biomass_edna.pdf", width = 10, height = 14)
par(mfrow = c(5,1), mar = c(2, 4, 2, 4), oma = c(2, 0, 0, 0))

for (s in spp) {
  bm <- biomass_list[[s]]
  dd <- d[d$spp_short == s, ]
  
  # ── left axis: modelled biomass trajectory ──
  yr_range <- range(dd$dates, na.rm = TRUE)
  bm_sub   <- bm[bm$year >= yr_range[1] & bm$year <= yr_range[2], ]
  
  plot(bm_sub$year, bm_sub$`0.5`, type = "l", lwd = 2, col = "black",
       xlim = yr_range,
       ylim = range(c(bm_sub$`0.025`, bm_sub$`0.975`)),
       xlab = "", ylab = "", main = s)
  
  polygon(c(bm_sub$year, rev(bm_sub$year)),
          c(bm_sub$`0.025`, rev(bm_sub$`0.975`)),
          col = rgb(0, 0, 0, 0.1), border = NA)
  
  mtext("Modelled biomass", side = 2, line = 2.5, cex = 0.65)
  
  par(new = TRUE)
  ylim2 <- c(0, max(dd$value, na.rm = TRUE))
  plot(1, type = "n", xlim = yr_range, ylim = ylim2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Positive PCR replicates", side = 4, line = 2.5, cex = 0.65)
  
  for (si in sites) {
    dds <- dd[dd$site == si, ]
    if (nrow(dds) == 0) next
    points(dds$dates, dds$value,
           pch = 16, col = adjustcolor(site_cols[si], 0.6), cex = 0.9)
  }
  
  # legend on first panel only, after everything is drawn
  if (s == spp[1]) {
    legend("topright",
           legend = c(names(site_cols), "Biomass (median)", "Biomass (95% CI)"),
           col    = c(site_cols, "black", rgb(0, 0, 0, 0.3)),
           pch    = c(rep(16, length(site_cols)), NA, NA),
           lwd    = c(rep(NA, length(site_cols)), 2, 8),
           lty    = c(rep(NA, length(site_cols)), 1, 1),
           bty    = "n")
  }
}


dev.off()

