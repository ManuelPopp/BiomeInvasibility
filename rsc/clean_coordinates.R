#!/usr/bin/env Rscript
#>------------------------------------<
##
## Script name: Clean coordinates
##
##
## Author: Manuel R. Popp
## Email: manuel.popp@wsl.ch
##
## Date Created: 2024-12-07
##
## ---------------------------
##
## Descripton: Remove bad coordinates from species occurence data
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
  "dplyr", "devtools", "rnaturalearth",
  dependencies = TRUE
)

if (!require("GBIFhandleR", character.only = TRUE)) {
  devtools::install_github("https://github.com/ManuelPopp/GBIFhandleR")
}

require("GBIFhandleR")

#>----------------------------------------------------------------------------<|
#> Settings
if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "C:/Users/poppman/switchdrive/PhD/prj/bir/git/BiomeInvasibility"
  dir_dat <- "C:/Users/poppman/switchdrive/PhD/prj/bir/dat"
  dir_raw <- "L:/nobis/biomes/gbif/gbif-cleaned01-rds"
} else {
  dir_main <- "/lud11/poppman/git/BiomeInvasibility"
  dir_dat <- "/lud11/poppman/data/bir/dat"
  dir_raw <- "/lud11/nobis/biomes/gbif/gbif-cleaned01-rds"
}

dir_lud <- file.path(dir_dat, "lud11")
dir_obs <- file.path(dir_lud, "obs")
dir_cleaned <- file.path(dir_lud, "cln")
dir_git_dat <- "C:/Users/poppman/switchdrive/PhD/prj/bir/git/BiomeInvasibility/dat"

#>----------------------------------------------------------------------------<|
#> Functions
clean_files <- function(
    files, out_dir, out_name = "cleaned_observations", one_file = FALSE,
    quiet = TRUE
    ) {
  # Get reference polygons
  ref_sea <- rnaturalearth::ne_download(
    scale = 110, type = "land", category = "physical", returnclass = "sf"
    )
  ref_urban <- rnaturalearth::ne_download(
    scale = "medium", type = "urban_areas", returnclass = "sf"
    )

  pb <- utils::txtProgressBar(
    min = 0, max = length(files), style = ifelse(quiet, 1, 3)
    )
  cleaned <- lapply(
    X = seq(1:length(files)), FUN = function(i) {
      file <- files[i]
      if (file.exists(file)) {
        df <- readRDS(file) %>%
          as.data.frame()
      } else {
        df <- data.frame()
      }
      
      if (nrow(df) > 0) {
        if (quiet) {
          out <- suppressMessages(
            suppressWarnings(
              capture.output(
                GBIFhandleR::clean_coords(
                  df, ref_urban = ref_urban, ref_sea = ref_sea,
                  species_column = "acceptedScientificName"
                  ),
                file = ifelse(
                  Sys.info()["sysname"] == "Windows", "C:/NUL", "/dev/null"
                  )
                )
            )
          )
        } else {
          out <- GBIFhandleR::clean_coords(
            df, ref_urban = ref_urban, ref_sea = ref_sea
          )
          }
        } else {
          out <- df
        }
      utils::setTxtProgressBar(pb, i)
      return(out)
      }
    )

  names(cleaned) <- sub(".rds", ".csv", basename(files))

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  num_records <- unlist(lapply(X = cleaned, FUN = nrow))
  if (one_file) {
    cleaned_output <- do.call(rbind, cleaned)
    write.csv(
      cleaned_output, file = file.path(out_dir, paste0(out_name, ".csv"))
      )
  } else {
    for (i in 1:length(cleaned)) {
      write.csv(cleaned[[i]], file = file.path(out_dir, names(cleaned)[i]))
    }
  }
  return(num_records)
}

#>----------------------------------------------------------------------------<|
#> Main
species_table <- read.csv(file.path(dir_git_dat, "matching_species.csv")) %>%
  dplyr::mutate(
    path = file.path(
      dir_raw, class, order, family, genus, paste0(species, ".rds")
    )
  )

num_records <- clean_files(
  files = species_table$path,
  out_dir = dir_cleaned
  )