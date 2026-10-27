## 0_5: Screen metazoan family assignments against RAMS and build the merged
##      functional-groups file consumed by 1_1_Analysis.R
##
## Endpoint format verified 2026-08-25: the RAMS REST is routed through
## aphia.php with the path carried on a bare "__route__" parameter:
##   https://www.marinespecies.org/rams/aphia.php?p=rest&__route__/AphiaRecordsByName/{name}&like=false&...
## Note "&" (not "?") joins parameters after the route - the URL already
## carries a query string. The clean "/rams/rest/..." path in the RAMS docs'
## own R example 404s and must not be used.
##
## Rules (in order, all mechanical):
##   1. Family in RAMS                      -> keep name (valid_name catches synonyms)
##   2. Not in RAMS, name known to WoRMS    -> climb Order -> Class -> Phylum;
##      first rank with a RAMS record       -> "<Rank> (fam. indet.)"
##   3. Nothing up to Phylum in RAMS        -> register SCOPE gap, not evidence of
##      absence -> retain family name with note (e.g. Xenacoelomorpha/Convolutidae)
##   4. Name unrecognised by WoRMS          -> trace via TYPE GENUS (family names
##      are formed from the type-genus stem) to the accepted WoRMS family, then
##      re-apply rules 1-3 to that family name
##
##
## Outputs:
##   rawdata/family_RAMS_check.csv              (per-family screen results)
##   rawdata/metazoa_functional_groups_RAMS.csv (merged file used for plotting)

library(jsonlite)

rams_base  <- "https://www.marinespecies.org/rams/aphia.php?p=rest&__route__"
worms_base <- "https://www.marinespecies.org/rest"

fetch_json <- function(url) {
  Sys.sleep(0.4)
  tryCatch(fromJSON(url), error = function(e) NULL)
}

aphia_by_name <- function(name, base) {
  sep <- if (grepl("\\?", base)) "&" else "?"   # routed base already has ?p=rest
  url <- sprintf("%s/AphiaRecordsByName/%s%slike=false&marine_only=false",
                 base, URLencode(name, reserved = TRUE), sep)
  rec <- fetch_json(url)
  if (is.null(rec) || length(rec) == 0) return(NULL)
  rec <- as.data.frame(rec)
  acc <- rec[rec$status == "accepted", , drop = FALSE]
  if (nrow(acc) > 0) acc[1, ] else rec[1, ]
}

worms_classification <- function(aphia_id) {
  cl <- fetch_json(sprintf("%s/AphiaClassificationByAphiaID/%s",
                           worms_base, aphia_id))
  out <- character(0)
  while (!is.null(cl)) {
    out[cl$rank] <- cl$scientificname
    cl <- cl$child
    if (is.null(cl$rank)) break
  }
  out
}

in_rams <- function(name) !is.null(aphia_by_name(name, rams_base))

## Rule 4 lookup: NCBI family names not recognised by WoRMS, traced via their
## type genus (mechanical: Obeli-idae -> Obelia; Siphonodentali-idae ->
## Siphonodentalium). Extend here if new unrecognised names appear.
type_genus <- c("Obeliidae"          = "Obelia",
                "Siphonodentaliidae" = "Siphonodentalium")

## Screen one family name through rules 1-3; returns list(name, basis, extras)
screen_family <- function(f) {
  rams_rec <- aphia_by_name(f, rams_base)
  if (!is.null(rams_rec)) {
    return(list(in_RAMS = TRUE, rams_status = rams_rec$status,
                rams_valid_name = rams_rec$valid_name,
                is_marine = rams_rec$isMarine,
                is_freshwater = rams_rec$isFreshwater,
                worms_order = NA, worms_class = NA,
                proposed_name = rams_rec$valid_name,
                basis = if (identical(rams_rec$valid_name, f)) "RAMS: accepted"
                else "RAMS: synonym updated"))
  }
  worms_rec <- aphia_by_name(f, worms_base)
  if (is.null(worms_rec)) {
    return(list(in_RAMS = FALSE, rams_status = NA, rams_valid_name = NA,
                is_marine = NA, is_freshwater = NA,
                worms_order = NA, worms_class = NA,
                proposed_name = NA, basis = "not found in WoRMS"))
  }
  cls <- worms_classification(worms_rec$AphiaID)
  for (rank in c("Order", "Class", "Phylum")) {
    if (!is.na(cls[rank]) && in_rams(cls[rank])) {
      return(list(in_RAMS = FALSE, rams_status = NA, rams_valid_name = NA,
                  is_marine = NA, is_freshwater = NA,
                  worms_order = cls["Order"], worms_class = cls["Class"],
                  proposed_name = paste0(cls[rank], " (fam. indet.)"),
                  basis = paste0("no SO record; nearest RAMS rank = ", rank)))
    }
  }
  ## Rule 3: nothing up to Phylum -> register scope gap, retain the name
  list(in_RAMS = FALSE, rams_status = NA, rams_valid_name = NA,
       is_marine = NA, is_freshwater = NA,
       worms_order = cls["Order"], worms_class = cls["Class"],
       proposed_name = f,
       basis = "retained: RAMS coverage gap (phylum absent from register scope)")
}

## -- run screen ----------------------------------------------------------------
grp  <- read.csv("rawdata/metazoa_functional_groups.csv")
fams <- sort(unique(grp$family))

res <- data.frame(family = fams, in_RAMS = NA, rams_status = NA,
                  rams_valid_name = NA, is_marine = NA, is_freshwater = NA,
                  worms_order = NA, worms_class = NA,
                  proposed_name = NA, basis = NA, stringsAsFactors = FALSE)

for (i in seq_along(fams)) {
  f <- fams[i]
  cat(sprintf("[%2d/%d] %-18s ", i, length(fams), f))
  
  s <- screen_family(f)
  
  ## Rule 4: WoRMS-unrecognised name -> trace via type genus, re-screen
  if (identical(s$basis, "not found in WoRMS")) {
    tg <- type_genus[f]
    if (!is.na(tg)) {
      gen <- aphia_by_name(tg, worms_base)
      if (!is.null(gen) && !is.na(gen$family)) {
        s <- screen_family(gen$family)
        s$basis <- sprintf(
          "NCBI name unrecognised by WoRMS; type genus %s -> %s (WoRMS); %s",
          tg, gen$family, s$basis)
      } else {
        s$proposed_name <- paste(f, "(not found)")
        s$basis <- "type genus not resolvable - inspect manually"
      }
    } else {
      s$proposed_name <- paste(f, "(not found)")
      s$basis <- "not in RAMS; not in WoRMS; no type-genus entry - add to type_genus"
    }
  }
  
  for (col in names(s)) res[[col]][i] <- s[[col]]
  cat("->", res$proposed_name[i], "\n")
}

res$name_changed <- res$proposed_name != res$family

## -- guild assessment ----------------------------------------------------------
## Does the functional group survive the (usually broader) relabelled rank?
## Documented judgments; anything not listed defaults to "holds" with no note.
## FLAGGED = guild no longer secured by taxonomy at the new rank; the plotting
## script (1_1_Analysis.R) greys these dots automatically.
guild_assessment <- list(
  "Barentsiidae" = c("holds",
                     "Entoprocta are uniformly sessile suspension feeders; group holds at phylum rank."),
  "Echinometridae" = c("holds",
                       "Camarodont urchins (incl. Antarctic Sterechinus) are grazers/omnivores; group holds at order rank."),
  "Foveolariidae" = c("holds",
                      "All cheilostome bryozoans are suspension feeders; group holds at order rank."),
  "Obeliidae" = c("holds",
                  "Name change only (NCBI Obeliidae = WoRMS Campanulariidae); same colonial hydroid ecology."),
  "Siphonodentaliidae" = c("holds",
                           "Name change only (Siphonodentalium now placed in Gadilidae); same infaunal scaphopods."),
  "Speleoithonidae" = c("FLAGGED",
                        paste("Cyclopoida spans benthic AND planktonic ecologies (cf. Oithonidae, pelagic in",
                              "this dataset). At order rank the Meiofauna assignment is no longer secured by",
                              "taxonomy; retained provisionally (sediment-derived ASV) and greyed in guild colouring.")),
  "Nephtheidae" = c("holds - with note",
                    paste("Name retained (in RAMS via Gersemia) but NCBI and PR2 disagree on this ASV's",
                          "placement (Nephtheidae vs Swiftia/Plexauridae); any candidate octocoral is a",
                          "sessile suspension feeder, so the group is robust.")),
  "Convolutidae" = c("holds - with note",
                     "Name retained: RAMS does not cover Xenacoelomorpha (register scope gap, not evidence of absence)."),
  "Canthocamptidae" = c("holds - with note",
                        paste("In RAMS with both marine AND freshwater flags; constituent ASV genera span",
                              "freshwater (Bryocamptus) and brackish/marine (Mesochra, Itunella) taxa,",
                              "consistent with mixed meltwater input and marine meiofauna.")),
  "Clavelinidae" = c("holds - with note",
                     "NCBI places Distaplia in Clavelinidae; WoRMS uses Holozoidae - nomenclature only.")
)

grp$display_name <- res$proposed_name[match(grp$family, res$family)]
ga <- guild_assessment[grp$family]
grp$group_check <- vapply(ga, function(x) if (is.null(x)) "holds" else x[1], "")
grp$group_explanation <- vapply(ga, function(x) if (is.null(x)) "" else x[2], "")

merged <- cbind(grp[, c("family", "group", "display_name",
                        "group_check", "group_explanation")],
                res[match(grp$family, res$family),
                    c("in_RAMS", "rams_status", "rams_valid_name",
                      "is_marine", "is_freshwater", "worms_order",
                      "worms_class", "proposed_name", "basis", "name_changed")])

write.csv(res, "rawdata/family_RAMS_check.csv", row.names = FALSE)
write.csv(merged, "rawdata/metazoa_functional_groups_RAMS.csv", row.names = FALSE)

## -- summary & tripwires -------------------------------------------------------
cat("\n---- Summary ----\n")
cat("Families checked:", nrow(res),
    "| in RAMS:", sum(res$in_RAMS, na.rm = TRUE),
    "| names changed:", sum(res$name_changed, na.rm = TRUE),
    "| guild flags:", sum(merged$group_check == "FLAGGED"), "\n\n")
print(res[res$name_changed, c("family", "proposed_name", "basis")],
      row.names = FALSE)

## Tripwires. (1) If these are "absent" the endpoint has broken again - do not
## trust the run; this caught the silent failure in v1. (2) Every family must
## leave with a resolved name. (3) The merged file must cover all input rows.
stopifnot(res$in_RAMS[match(c("Euphausiidae", "Salpidae", "Iphimediidae"),
                            res$family)])
stopifnot(!any(is.na(merged$display_name)),
          !any(grepl("not found|inspect manually", merged$display_name)))
stopifnot(nrow(merged) == nrow(grp))
cat("\nSanity checks passed: RAMS responding; all", nrow(merged),
    "families resolved; merged file written for 1_1_Analysis.R.\n")
