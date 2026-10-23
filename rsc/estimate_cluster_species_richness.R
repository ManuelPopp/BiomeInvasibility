lud11 <- ifelse(Sys.info()["sysname"] == "Windows", "L:", "/lud11")
# install.packages(file.path(lud11, "nobis", "Manuel", "rsdd_0.2.10.tar.gz"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("phyloseq")
# pak::pak("adw96/breakaway")
library("rsdd")
library("dplyr")
library("tidyr")
library("terra")
library("nanoparquet")
library("vegan")
library("breakaway")

# Set file paths
f_bbuff <- file.path(lud11, "poppman/tmp", "biomes_buff.gpkg")
f_template <- file.path(lud11, "nobis/Manuel/gbif/data/map_template.tif")
f_taxa <- file.path(lud11, "nobis/Manuel/gbif-powo_raw_specID-cellID-occN.parquet")
f_native <- file.path(
  lud11, "poppman", "data", "bir", "dat", "lud11", "native_obs.parquet"
)
f_native_rastid <- file.path(
  lud11, "poppman", "data", "bir", "dat", "lud11", "native_obs_rast.parquet"
)

# Settings
C <- 6.7 * 10^(-7)
half_decay_distance <- log(0.5) / (-C)

# Load data
rsdd::dataset("gbif-powo_raw")
taxon_tab <- rsdd::taxa()
f_out <- file.path(lud11, "poppman/tmp", "biome_cluster_species_richness.csv")

# Load frequency table
if (!file.exists(f_native)) {
  ftab <- nanoparquet::read_parquet(f_taxa)
  colnames(ftab) <- c("specID", "cellID", "freq")
  obs_tab <- merge(ftab, rsdd:::ds_cells, by = c("specID", "cellID")) %>%
    dplyr::filter(native == TRUE) %>%
    dplyr::select(specID, cellID, freq)
  nanoparquet::write_parquet(obs_tab, f_native)
  gc()
} else {
  obs_tab <- nanoparquet::read_parquet(f_native)
}

# Load the buffered biomes and rasterise them to get biome IDs per cell
bbuff <- terra::vect(f_bbuff)

if (!file.exists(f_native_rastid)) {
  map <- terra::rast(f_template)
  bbuffr <- terra::rasterize(
    bbuff,
    map,
    field = "ID",
    filename = file.path(
      lud11, "poppman", "data", "bir", "dat", "lud11", "BiomesBuffered.tif"
    ),
    overwrite = TRUE,
    background = NA
  )
  
  # Join biome IDs to frequency table
  biome_ids <- terra::as.data.frame(bbuffr, xy = FALSE, cells = TRUE)
  colnames(biome_ids) <- c("cellID", "biomeID")
  freq_tab <- base::merge(obs_tab, biome_ids, by = "cellID")
  nanoparquet::write_parquet(freq_tab, f_native_rastid)
  gc()
} else {
  freq_tab <- nanoparquet::read_parquet(f_native_rastid)
}

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
get_cluster_species_richness <- function(
    biome_id,
    biome_vect = bbuff, maxdist = half_decay_distance
    ) {
  row_idx <- match(biome_id, biome_vect$ID)
  if (is.na(row_idx)) {
    stop(paste("Biome ID not found:", biome_id))
  }
  buffered_patch <- terra::buffer(biome_vect[row_idx, ], width = maxdist)
  same_biome <- biome_vect$BIOME == buffered_patch$BIOME[1]
  intersects <- terra::relate(buffered_patch, biome_vect, "intersects")[1, ]
  biome_ids <- which(intersects & same_biome)
  sub_tab <- freq_tab[freq_tab$biomeID %in% biome_ids, ]
  
  if (nrow(sub_tab) == 0) {
    return(
      data.frame(
        ID = biome_id,
        speciesRichnessRaw = NA,
        speciesRichnessBa = NA,
        speciesRichnessBaSE = NA
      )
    )
  }
  
  # Counted species
  spec_count <- length(unique(sub_tab$specID))
  
  # Based on total observation counts
  obs_counts <- sub_tab %>%
    dplyr::group_by(specID) %>%
    dplyr::summarise(counts = sum(freq)) %>%
    dplyr::pull(counts)
  
  # Based on occupied cell counts
  cell_counts <- sub_tab %>%
    dplyr::group_by(specID) %>%
    dplyr::summarise(counts = dplyr::n()) %>%
    dplyr::pull(counts)
  
  baway <- tryCatch(
    breakaway::breakaway(obs_counts),
    error = function(e) {
      warning(
        "breakaway failed for biome ID ", biome_id,
        ": ", conditionMessage(e)
      )
      list(est = NA, se = NA)
    }
  )
  
  df <- data.frame(
    ID = as.num(biome_id),
    speciesRichnessRaw = spec_count,
    speciesRichnessBa = as.num(baway$est),
    speciesRichnessBaSE = as.num(baway$se)
  )
  gc()
  return(df)
}

biome_ids_unique <- sort(unique(bbuff$ID))

results_df <- do.call(
  rbind,
  lapply(
    X = biome_ids_unique,
    FUN = function(id) {
      tryCatch(
        get_cluster_species_richness(id),
        error = function(e) {
          stop(paste("Error in biome ID", id, ":", e$message))
        }
      )
    }
  )
)

# Save results
write.csv(results_df, file = f_out, row.names = FALSE)
