#!/usr/bin/env Rscript
#>------------------------------------<
##
## Script name: Get species observations
##
## Author: Manuel R. Popp
## Email: manuel.popp@wsl.ch
##
## Date Created: 2025-04-30
##
## ---------------------------
##
## Descripton: Get the observation data from local GBIF download
## Notes: -
##
#>----------------------------------------------------------------------------<|
#> Install/load packages
rm(list = ls())
import <- function(...) {
  #' Import R packages. Install them if necessary.
  #'
  #' @param ... any argument that can be passed to install.packages.
  #' @details The function installs only packages that are missing. Packages
  #' are loaded.
  #' @examples
  #' # Load packages
  #' import("dplyr", "MASS", "terra", dependencies = TRUE)
  #'
  #' @seealso \code{\link[base]{install.packages}}
  #' @export
  args <- list(...)
  packages = args[names(args) == ""]
  kwargs = args[names(args) != ""]
  
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      do.call(install.packages, c(list(package), kwargs))
    }
    require(package, character.only = TRUE)
  }
}

import(
  "dplyr", "devtools",
  dependencies = TRUE
)

if (!require("GBIFhandleR", character.only = TRUE)) {
  devtools::install_github("https://github.com/ManuelPopp/GBIFhandleR")
}

require("GBIFhandleR")

#>----------------------------------------------------------------------------<|
#> Settings
if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "C:/Users/poppman/switchdrive/PhD/prj/bir"
} else {
  dir_main <- "/lud11/poppman/BiomeInvasibility"
}

dir_dat <- file.path(dir_main, "dat", "lud11")
dir_obs <- file.path(dir_dat, "obs")
dir_git_dat <- file.path(dir_main, "dat")

f_species <- file.path(dir_dat, "db", "invasive_db_2013_normalised.csv")
f_database_meta <- file.path(dir_dat, "db", "master_table_speciesID.rds")

# Set up directory for output
dir.create(dir_obs, showWarnings = FALSE)

if (!file.exists(file.path(dir_main, ".gitignore.txt"))) {
  sink(file.path(dir_main, ".gitignore.txt"))
  print("dat/obs")
  sink()
}

#>----------------------------------------------------------------------------<|
#> Functions

#>----------------------------------------------------------------------------<|
#> Read database metadata
db_meta <- readRDS(f_database_meta)

#>----------------------------------------------------------------------------<|
#> Read species list
species_data <- read.csv(f_species, comment.char = "#")
for (i in which(species_data$status != "ACCEPTED")) {
  name_query <- rgbif::name_backbone(species_data$scientificName[i])
  
  if ("usageKey" %in% names(name_query)) {
    usage_key <- name_query$usageKey
    gbif_taxon_data <- rgbif::name_usage(key = usage_key)$data
  } else {
    gbif_taxon_data <- data.frame()
  }
  
  if (nrow(gbif_taxon_data) != 0 & "species" %in% names(gbif_taxon_data)) {
    species_data$canonicalName[i] <- gbif_taxon_data$species
  } else {
    species_data$canonicalName[i] <- NA
  }
}

species_names <- species_data %>%
  dplyr::pull(canonicalName) %>%
  na.omit()

matching_names <- species_names[which(species_names %in% db_meta$species)]

db_species <- db_meta[which(db_meta$species %in% matching_names),]
write.csv(
  db_species,
  file = file.path(dir_git_dat, "matching_species.csv"),
  row.names = FALSE
  )
