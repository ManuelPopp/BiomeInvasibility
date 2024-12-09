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
  dir_main <- "C:/Users/poppman/switchdrive/PhD/prj/bir/git/BiomeInvasibility"
  dir_main <- "L:/poppman/BiomeInvasibility"
} else {
  dir_main <- "/lud11/poppman/BiomeInvasibility"
}

dir_dat <- file.path(dir_main, "dat")
dir_obs <- file.path(dir_dat, "obs")
dir_cleaned <- file.path(dir_obs, "cleaned")

#>----------------------------------------------------------------------------<|
#> Functions
clean_files <- function(
    directory, out_dir, out_name = "cleaned_observations", one_file = TRUE
    ) {
  files <- list.files(path = directory, pattern = ".csv", full.names = TRUE)
  cleaned <- lapply(
    X = files, FUN = function(x) {
      df <- read.csv(x, header = TRUE)
      return(GBIFhandleR::clean_coords(df))
      }
    )

  names(cleaned) <- basename(files)

  if (one_file) {
    cleaned_output <- do.call(rbind, cleaned)
    write.csv(
      cleaned_output, file = file.path(out_dir, paste0(out_name, ".csv"))
      )
  } else {
    for (ds in cleaned) {
      write.csv(ds, file = file.path(out_dir, name(ds)))
    }
  }
}

#>----------------------------------------------------------------------------<|
#> Main
clean_files(
  directory = dir_obs,
  out_dir = dir_cleaned
  )
