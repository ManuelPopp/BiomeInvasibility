#!/usr/bin/env Rscript
#>------------------------------------<
##
## Script name: Download species
## occurence data
##
## Author: Manuel R. Popp
## Email: manuel.popp@wsl.ch
##
## Date Created: 2024-12-07
##
## ---------------------------
##
## Descripton: Download species observations from GBIF
## Notes: -
##
#>----------------------------------------------------------------------------<|
#> Install/load packages
rm(list = ls())

if (!require("GBIFhandleR", character.only = TRUE)) {
  devtools::install_github("https://github.com/ManuelPopp/GBIFhandleR")
}

require("GBIFhandleR")
require("dplyr")

#>----------------------------------------------------------------------------<|
#> Settings
dir_main <- "C:/Users/poppman/switchdrive/PhD/prj/bir/git/BiomeInvasibility"
dir_dat <- file.path(dir_main, "dat")
dir_obs <- file.path(dir_dat, "obs")

f_species <- "species_names.csv"
col_spnames <- "Accepted_SPNAME"

# Set up directory for output
dir.create(dir_obs, showWarnings = FALSE)

if (!file.exists(file.path(dir_main, ".gitignore.txt"))) {
  sink(file.path(dir_main, ".gitignore.txt"))
  print("dat/obs")
  sink()
}

#>----------------------------------------------------------------------------<|
#> Functions
download <- function(species_name) {
  dst <- file.path(dir_obs, paste0(species_name, ".csv"))
  
  if (!file.exists(dst)) {
    observations <- GBIFhandleR::get_observations(species_name)
    write.csv2(dst)
  }
}

#>----------------------------------------------------------------------------<|
#> Read species list
species <- read.csv(file.path(dir_dat, f_species)) %>%
  dplyr::select(col_spnames) %>%
  as.vector()

for (i in 1:length(species)) {
  species_name <- species[i]
  cat(
    "Downloading species", i, "of", length(species), paste0("(", species, ").")
    )
  
  download(species_name)
}
