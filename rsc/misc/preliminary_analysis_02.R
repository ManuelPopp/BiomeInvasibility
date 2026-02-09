# Preliminary analysis based on GloNaf data

if (!require("breakaway")) {
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  remotes::install_github("adw96/breakaway")
}

library(breakaway)
library(dplyr)
library(forcats)
library(terra)


dir_glonaf = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/glonaf"
f_glonaf_names <- file.path(dir_glonaf, "glonaf_taxon_wcvp.csv")
f_glonaf_data <- file.path(dir_glonaf, "glonaf_flora2.csv")
f_glonaf_region_info <- file.path(dir_glonaf, "glonaf_region.csv")
f_glonaf_regions <- file.path(
  dir_glonaf,
  "glonaf_2024_regions/glonaf_2024_regions.shp"
  )

glonaf_names <- read.csv(f_glonaf_names) %>%
  dplyr::mutate(
    taxon = taxa_accepted,
    taxon_rank = wcvp_taxon_rank,
    taxon_wcvp_id = id
    ) %>%
  dplyr::select(taxon_wcvp_id, taxon, taxon_rank)

glonaf_region_info <- read.csv(f_glonaf_region_info) %>%
  dplyr::mutate(region_id = id)

glonaf_data <- read.csv(f_glonaf_data) %>%
  dplyr::mutate(
    status = forcats::fct_recode(status, "Naturalised" = "Naturalized")
    ) %>%
  dplyr::rename(data_id = id) %>%
  dplyr::left_join(glonaf_names, by = "taxon_wcvp_id") %>%
  dplyr::left_join(glonaf_region_info, by = "region_id")


glonaf_regions <- terra::vect(f_glonaf_regions)

glonaf_taxa <- unique(glonaf_data$taxon)

taxon <- glonaf_taxa[1]

plot(glonaf_regions)
df <- glonaf_data[which(glonaf_data$taxon == taxon),]
naturalised_loc_ids <- df[which(df$status == "Naturalised"),]$OBJIDsic
invasive_loc_ids <- df[which(df$status == "Invasive"),]$OBJIDsic
naturalised_regs <- glonaf_regions[which(glonaf_regions$OBJIDsic %in% naturalised_loc_ids),]
invasive_regs <- glonaf_regions[which(glonaf_regions$OBJIDsic %in% invasive_loc_ids),]
terra::plot(naturalised_regs, col = "orange", add = TRUE)
terra::plot(invasive_regs, col = "red", alpha = 0.5, add = TRUE)
