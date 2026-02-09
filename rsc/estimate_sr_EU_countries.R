library(dplyr)
library(tidyr)
library(terra)
library(nanoparquet)

# Set file paths
f_bbuff <- file.path("/lud11/poppman/Europe", "Europe_merged.shp")
f_template <- "/lud11/nobis/Manuel/gbif/data/map_template.tif"
f_taxa <- "/lud11/nobis/Manuel/gbif/data/map_taxa.parquet"
f_cells <- "/lud11/nobis/Manuel/gbif/data/map_cells.parquet"

# Load taxa and frequency table
taxa <- nanoparquet::read_parquet(f_taxa)
freq_tab <- data.frame(nanoparquet::read_parquet(f_cells))
colnames(freq_tab) <- c("row","specID", "cellID", "freq")
freq_tab$cellID <- as.numeric(freq_tab$cellID)

# Load the buffered biomes and rasterise them to get biome IDs per cell
bbuff <- terra::vect(f_bbuff)
bbuff$ID <- 1:nrow(bbuff)
map <- terra::rast(f_template)
bbuffr <- terra::rasterize(bbuff, map, field = "ID")

# Join biome IDs to frequency table
biome_ids <- terra::as.data.frame(bbuffr, xy = FALSE, cells = TRUE)
colnames(biome_ids) <- c("cellID", "biomeID")
freq_tab <- base::merge(freq_tab, biome_ids, by = "cellID")

# Estimate species richness per biome
# Estimate species richness per biome
## Helper function to safely convert to numeric
as.num <- function(x) {
  x <- as.numeric(x)
  if (is.null(x)) {
    return(c(NA))
  }
  if (length(x) == 0) {
    return(c(NA))
  }
  if (is.nan(x)) {
    return(c(NA))
  }
  return(x)
}

## Function to estimate species richness for a given biome ID
get_species_richness <- function(biome_id) {
  sub_tab <- freq_tab[freq_tab$biomeID == biome_id, ]

  if (nrow(sub_tab) == 0) {
    return(
        data.frame(
            ID = biome_id,
            speciesRichnessBa = NA,
            speciesRichnessBaSE = NA,
            speciesRichnessChao1 = NA,
            speciesRichnessChao1SE = NA,
            speciesRichnessACE = NA,
            speciesRichnessACESE = NA,
            speciesRichnessChao2 = NA
            )
        )
    }

  obs_counts <- sub_tab %>%
    dplyr::group_by(specID) %>%
    dplyr::summarise(counts = sum(freq)) %>%
    dplyr::pull(counts)
  
  baway <- tryCatch(
    breakaway::breakaway(obs_counts),
    error = function(e) list(est = NA, se = NA)
  )
  
  chao1 <- vegan::estimateR(obs_counts)
  
  inc_freq <- sub_tab %>%
    dplyr::distinct(cellID, specID) %>%
    dplyr::count(specID, name = "inc")
  Q1 <- sum(inc_freq$inc == 1)
  Q2 <- sum(inc_freq$inc == 2)
  n  <- dplyr::n_distinct(sub_tab$cellID)
  S_obs <- nrow(inc_freq)

  if (Q2 > 0) {
    chao2 <- S_obs + (Q1^2) / (2 * Q2)
  } else {
    chao2 <- S_obs + (Q1 * (Q1 - 1)) / 2
  }
  
  df <- data.frame(
    ID = as.num(biome_id),
    speciesRichnessBa = as.num(baway$est),
    speciesRichnessBaSE = as.num(baway$se),
    speciesRichnessChao1 = as.num(chao1["S.chao1"]),
    speciesRichnessChao1SE = as.num(chao1["se.chao1"]),
    speciesRichnessACE = as.num(chao1["S.ACE"]),
    speciesRichnessACESE = as.num(chao1["se.ACE"]),
    speciesRichnessChao2 = as.num(chao2)
    )
  return(df)
}

biome_ids_unique <- sort(unique(bbuff$ID))

results_df <- do.call(
    rbind,
    lapply(
      X = biome_ids_unique,
      FUN = function(id) {
        tryCatch(
          get_species_richness(id),
          error = function(e) {
            stop(paste("Error in biome ID", id, ":", e$message))
          }
          )
        }
      )
    )

# Save results
terra::merge(bbuff, results_df, by = "ID", all.x = TRUE) %>%
  terra::writeVector(
    filename = "/lud11/poppman/Europe/Europe_species_richness.gpkg",
    overwrite = TRUE
    )
