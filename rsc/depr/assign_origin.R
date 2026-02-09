#!/usr/bin/env Rscript
#>------------------------------------<
##
## Script name: Clean coordinates
##
##
## Author: Manuel R. Popp
## Email: manuel.popp@wsl.ch
##
## Date Created: 2025-05-10
##
## ---------------------------
##
## Descripton: Assign species origin
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
  "dplyr", "devtools", "sf", "ggplot2",
  dependencies = TRUE
)

if (!require("GBIFhandleR", character.only = TRUE)) {
  devtools::install_github("https://github.com/ManuelPopp/GBIFhandleR")
}

require("GBIFhandleR")

#>----------------------------------------------------------------------------<|
#> Settings
if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir"
  dir_dat <- "C:/Users/poppman/switchdrive/PhD/prj/bir/dat"
  dir_raw <- "L:/nobis/biomes/gbif/gbif-cleaned01-rds"
} else {
  dir_main <- "/lud11/poppman/git/BiomeInvasibility"
  dir_dat <- "/lud11/poppman/data/bir/dat"
  dir_raw <- "/lud11/nobis/biomes/gbif/gbif-cleaned01-rds"
}

dir_lud <- "L:/poppman/data/bir/dat/lud11"
dir_obs <- file.path(dir_lud, "obs")
dir_cleaned <- "D:/tmp"
f_biomes <- file.path(dir_lud, "shp", "olson_ecoregions", "biomes.shp")

#>----------------------------------------------------------------------------<|
#> Functions
load_and_process <- function(file, biomes) {
  taxon_name <- sub(".csv", "", basename(file))
  observations <- GBIFhandleR::load_observations(file)
  sp_range <- GBIFhandleR::get_range(taxon_name)
  locations <- GBIFhandleR::assign_status(observations, sp_range)
  native_patches <- sf::st_intersects(
    biomes,
    y = locations[locations$status == "native",],
    sparse = TRUE
    )
  exotic_patches <- sf::st_intersects(
    biomes,
    y = locations[locations$status == "exotic",]
    )
  unknown_patches <- sf::st_intersects(
    biomes,
    y = locations[locations$status == "unknown",]
  )
  
  # Plot map
  colours <- c("bisque1", "chartreuse2", "indianred2", "skyblue2")[
    as.numeric(lengths(native_patches) >= 1) +
      2 * as.numeric(lengths(exotic_patches) >= 1) +
      3 * as.numeric(lengths(unknown_patches) >= 1) +
      1
    ]
  
  png(
    filename = file.path(dir_main, "map", paste0(taxon_name, ".png")),
    width = 2048, height = 768
    )
  plot(biomes, col = colours, lty = 0)
  plot(sf::st_geometry(locations), add = TRUE, col = "black", pch = 4)
  dev.off()
  
  return(
    list(
      native_patches = lengths(native_patches),
      exotic_patches = lengths(exotic_patches),
      unknown_patches = lengths(unknown_patches),
      taxon_name = taxon_name
      )
    )
}

#>----------------------------------------------------------------------------<|
#> Load files and assign origin
wwf_biomes <- sf::st_read(f_biomes) %>%
  sf::st_geometry() %>%
  sf::st_make_valid()
  

files <- list.files(dir_cleaned, pattern = ".csv", full.names = TRUE)

data <- lapply(
  X = files,
  FUN = function(x) {
    try(load_and_process(x, biomes = wwf_biomes), silent = TRUE)
    }
)

save(
  data,
  file = file.path("L:/poppman/data/bir/dat/lud11/nat/native_exotic.Rdata")
  )

# Plot
alpha_values <- log(data[[1]]$native_patches + 1)
alpha_native <- alpha_values / (max(alpha_values) + 1e-6)

alpha_values <- log(data[[1]]$exotic_patches + 1)
alpha_exotic <- alpha_values / (max(alpha_values) + 1e-6)

alpha_values <- log(data[[1]]$unknown_patches + 1)
alpha_unknown <- alpha_values / (max(alpha_values) + 1e-6)

plot(wwf_biomes, col = rgb(0, 1, 0, alpha_native), lty = 0)
plot(wwf_biomes, col = rgb(1, 0, 0, alpha_exotic), lty = 0, add = TRUE)
plot(wwf_biomes, col = rgb(0, 0, 1, alpha_unknown), lty = 0, add = TRUE)