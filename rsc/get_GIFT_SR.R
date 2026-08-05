#!/usr/bin/env Rscript

library("sf")
library("GIFT")

# ------------------------------------------------------------------
# Optional command line arguments
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

nativeness_source <- if(length(args) >= 1) {
  args[1]
} else {
  "sinas"
}

biome_def <- if(length(args) >= 2) {
  args[2]
} else {
  "olson"
}

# ------------------------------------------------------------------
# Main directories
# ------------------------------------------------------------------

if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir"
  dir_lud11 <- "L:"
} else {
  dir_main <- "/lud11/poppman/data/bir"
  dir_lud11 <- "/lud11"
}

dir_dat <- file.path(dir_lud11, "poppman", "data", "bir", "dat", "lud11")
dir_imeb <- file.path(dir_dat, "biomes", biome_def, "intermediate_data")
f_bfullinfo <- file.path(dir_imeb, "biomes_full_info.gpkg")

biomes <- sf::read_sf(f_bfullinfo)

biomes$gift_nat_centroid <- NA
biomes$gift_nat_covered <- NA
biomes$gift_nad_centroid <- NA
biomes$gift_nad_covered <- NA


# ------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------

check_gift <- function(patch, overlap) {
  success <- FALSE
  lists <- NULL
  while (!success) {
    success <- tryCatch({
      lists <- GIFT::GIFT_checklists(
        taxon_name = "Tracheophyta",
        floristic_group = "all",
        complete_taxon = FALSE,
        shp = patch,
        overlap = overlap
      )
      TRUE
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (grepl(
        "connection|timeout|server|HTTP|socket|curl",
        msg,
        ignore.case = TRUE
      )) {
        message(
          format(Sys.time()),
          " - Connection error: ",
          msg,
          "\nRetrying in 5 minutes..."
        )
        Sys.sleep(300)
        FALSE
      } else {
        stop(e)
      }
    })
  }
  return(lists)
}

# ------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------

for (i in 1:nrow(biomes)) {
  patch_id <- biomes$ID[i]
  
  patch <- biomes |>
    dplyr::filter(ID == patch_id)
  
  # Check for covered centroids ------------------------------------------------
  lists <- check_gift(patch = patch, overlap = "centroid_inside")
  
  # Check whether the response is empty
  # It is a bit counter intuitive to check.
  if (!is.null(nrow(lists))) { # The responses are in an inconsistent format
    biomes$gift_nat_centroid[i] <- NA
    biomes$gift_nat_covered[i] <- NA
    biomes$gift_nad_centroid[i] <- NA
    biomes$gift_nad_covered[i] <- NA
    next
  }
  
  cl <- dplyr::bind_rows(lists$checklists)
  lists <- NULL
  
  biomes$gift_nat_centroid[i] <- subset(cl, native == 1) %>%
    dplyr::summarise(N = dplyr::n_distinct(work_species)) %>%
    dplyr::pull(N)
  biomes$gift_nad_centroid[i] <- subset(cl, naturalized == 1) %>%
    dplyr::summarise(N = dplyr::n_distinct(work_species)) %>%
    dplyr::pull(N)
  
  # Check for entire covered polygons ------------------------------------------
  lists <- check_gift(patch = patch, overlap = "shape_inside")
  
  if (!is.null(nrow(lists))) {
    biomes$gift_nat_covered[i] <- NA
    biomes$gift_nad_covered[i] <- NA
    next
  }
  
  cl <- dplyr::bind_rows(lists$checklists)
  lists <- NULL
  
  biomes$gift_nat_covered[i] <- subset(cl, native == 1) %>%
    dplyr::summarise(N = dplyr::n_distinct(work_species)) %>%
    dplyr::pull(N)
  biomes$gift_nad_covered[i] <- subset(cl, naturalized == 1) %>%
    dplyr::summarise(N = dplyr::n_distinct(work_species)) %>%
    dplyr::pull(N)
}

sf::st_write(biomes, f_bfullinfo, overwrite = TRUE)
