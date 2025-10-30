### Whale analysis 


## helper: load, reshape, add pristine year = 1850
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
add_pristine <- function(df, y0 = 1850) {
  r <- df[1, c("0.025","0.5","0.975")]
  r$year <- y0
  out <- rbind(r, df[, c("0.025","0.5","0.975","year")])
  rownames(out) <- NULL
  out
}

## Blue
blue3 <- prep("rawdata/biomassdatafromMilo/Blue_whale__SH.csv")
blue4 <- add_pristine(blue3)

## Orca
orca3 <- prep("rawdata/biomassdatafromMilo/Killer_whale__SH.csv")
orca4 <- add_pristine(orca3)

## Humpback
hump3 <- prep("rawdata/biomassdatafromMilo/Humpback_whale__SH.csv")
hump4 <- add_pristine(hump3)

## Fin
fin3 <- prep("rawdata/biomassdatafromMilo/Fin_whale__SH.csv")
fin4 <- add_pristine(fin3)

## plot using new data
yl <- c(0, max(blue4$`0.975`, orca4$`0.975`, hump4$`0.975`, fin4$`0.975`, na.rm = TRUE))
plot(blue4$year, blue4$`0.5`, type = "l", col = "dodgerblue", ylim = yl, xlab = "Year", ylab = "Biomass")
lines(blue4$year, blue4$`0.975`, lty = 2, col = "dodgerblue")
lines(blue4$year, blue4$`0.025`, lty = 2, col = "dodgerblue")

lines(orca4$year, orca4$`0.5`,   col = "red")
lines(orca4$year, orca4$`0.975`, lty = 2, col = "red")
lines(orca4$year, orca4$`0.025`, lty = 2, col = "red")

lines(hump4$year, hump4$`0.5`,   col = "grey")
lines(hump4$year, hump4$`0.975`, lty = 2, col = "grey")
lines(hump4$year, hump4$`0.025`, lty = 2, col = "grey")

lines(fin4$year, fin4$`0.5`,   col = "purple")
lines(fin4$year, fin4$`0.975`, lty = 2, col = "purple")
lines(fin4$year, fin4$`0.025`, lty = 2, col = "purple")

legend("topright",
       legend = c("Blue", "Orca", "Humpback", "Fin"),
       lty = 1, col = c("dodgerblue","red","grey","purple"), bty = "n")
