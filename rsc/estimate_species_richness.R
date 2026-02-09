library(dplyr)
library(tidyr)
library(terra)
library(nanoparquet)

# Set file paths
f_bbuff <- file.path("/lud11/poppman/tmp", "biomes_buff.gpkg")
f_template <- "/lud11/nobis/Manuel/gbif/data/map_template.tif"
f_taxa <- "/lud11/nobis/Manuel/gbif/data/map_taxa.parquet"
f_cells <- "/lud11/nobis/Manuel/gbif/data/map_cells.parquet"

f_out <- file.path("/lud11/poppman/tmp", "biome_species_richness.csv")

# Load taxa and frequency table
taxa <- nanoparquet::read_parquet(f_taxa)
freq_tab <- data.frame(nanoparquet::read_parquet(f_cells))
colnames(freq_tab) <- c("row","specID", "cellID", "freq")
freq_tab$cellID <- as.numeric(freq_tab$cellID)

# Load the buffered biomes and rasterise them to get biome IDs per cell
bbuff <- terra::vect(f_bbuff)
map <- terra::rast(f_template)
bbuffr <- terra::rasterize(bbuff, map, field = "ID")

# Join biome IDs to frequency table
biome_ids <- terra::as.data.frame(bbuffr, xy = FALSE, cells = TRUE)
colnames(biome_ids) <- c("cellID", "biomeID")
freq_tab <- base::merge(freq_tab, biome_ids, by = "cellID")

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
            biomePatchID = biome_id,
            speciesRichnessBa = NA,
            speciesRichnessBaSE = NA,
            speciesRichnessChao1 = NA,
            speciesRichnessChao1SE = NA,
            speciesRichnessACE = NA,
            speciesRichnessACESE = NA,
            speciesRichnessChao2 = NA,
            speciesRichnessChao2SE = NA,
            speciesRichnessjack1 = NA,
            speciesRichnessjack1SE = NA,
            speciesRichnessjack2 = NA,
            speciesRichnessjack2SE = NA,
            speciesRichnessboot = NA,
            speciesRichnessbootSE = NA
            )
        )
    }

  obs_counts <- sub_tab %>%
    dplyr::group_by(specID) %>%
    dplyr::summarise(counts = sum(freq)) %>%
    dplyr::pull(counts)
  
  # Breakaway (abundance-based)
  baway <- tryCatch(
    breakaway::breakaway(obs_counts),
    error = function(e) list(est = NA, se = NA)
  )
  # Chao1 (abundance-based)
  chao1 <- vegan::estimateR(obs_counts)
  
  # Create incidence matrix
  inc_matrix <- tryCatch(
    sub_tab %>%
      dplyr::select(cellID, specID, freq) %>%
      dplyr::mutate(presence = ifelse(freq > 0, 1, 0)) %>%
      tidyr::pivot_wider(
        names_from = specID,
        values_from = presence,
        values_fill = 0
      ) %>%
      dplyr::select(-cellID) %>%
      as.matrix(),
    error = function(e) NULL
  )
  # Chao2 (incidence-based)
  if (length(inc_matrix) == 0) {
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
    specp <- list(
      chao = chao2,
      chao.se = NA,
      jack1 = NA,
      jack1.se = NA,
      jack2 = NA,
      jack2.se = NA,
      boot = NA,
      boot.se = NA
    )
  } else {
    specp <- vegan::specpool(inc_matrix)
  }
  
  df <- data.frame(
    biomePatchID = as.num(biome_id),
    speciesRichnessBa = as.num(baway$est),
    speciesRichnessBaSE = as.num(baway$se),
    speciesRichnessChao1 = as.num(chao1["S.chao1"]),
    speciesRichnessChao1SE = as.num(chao1["se.chao1"]),
    speciesRichnessACE = as.num(chao1["S.ACE"]),
    speciesRichnessACESE = as.num(chao1["se.ACE"]),
    speciesRichnessChao2 = as.num(specp$chao),
    speciesRichnessChao2SE = as.num(specp$chao.se),
    speciesRichnessjack1 = as.num(specp$jack1),
    speciesRichnessjack1SE = as.num(specp$jack1.se),
    speciesRichnessjack2 = as.num(specp$jack2),
    speciesRichnessjack2SE = as.num(specp$jack2.se),
    speciesRichnessboot = as.num(specp$boot),
    speciesRichnessbootSE = as.num(specp$boot.se)
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
write.csv(results_df, file = f_out, row.names = FALSE)