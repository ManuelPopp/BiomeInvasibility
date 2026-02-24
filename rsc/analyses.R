#>=============================================================================<
#>-----------------------------------------------------------------------------<
#> Analysis based on GBIF data of SINAS and GloNAF-listed (invasive) species
#> 
#> Date: 2026-01-27
#> Author: Popp, MR
#> 
#>----------------------------------------------------------------------------->
#<=============================================================================>

print(Sys.time())
library(rsdd)
library(terra)
library(tidyterra)
library(pbapply)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(progress)

#>=============================================================================<
#> Settings
#<=============================================================================>

recompute <- FALSE

biome_names <- c(
  "(Sub)tropical Moist BLF",
  "(Sub)tropical Dry BLF",
  "(Sub)tropical Coniferous Forest",
  "Temperate Mixed Forest",
  "Temperate Coniferous Forest",
  "Taiga",
  "Savanna/Grassland",
  "Temperate Grassland",
  "Flooded Savanna",
  "Montane Grassland",
  "Tundra",
  "Mediterranean",
  "Desert",
  "Mangrove"
)

if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir"
  dir_tmp <- "L:/poppman/tmp"
} else {
  dir_main <- "/lud11/poppman/data/bir"
  dir_tmp <- "/lud11/poppman/tmp"
}

dir_dat <- file.path(dir_main, "dat", "lud11")
dir_stats <- file.path(dir_main, "stats")
f_glm_stats <- file.path(dir_stats, "glm_summaries.txt")
f_mountains <- file.path(
  dir_main, "dat", "lud11",
  "shp", "GMBA_Inventory_v2.0_standard", "GMBA_Inventory_v2.0_standard.shp"
  )
f_bathymetrie <- file.path(
  dir_main, "dat", "lud11",
  "gebco_2025_geotiff", "gebco.gpkg"
  )

#>=============================================================================<
#> Data preparation
#<=============================================================================>

rsdd::dataset("gbif-powo_raw")

f_sinas_places <- file.path(dir_dat, "sinas", "SInAS_Locations")
f_sinas_data <- file.path(dir_dat, "sinas", "SInAS_3.1.1.csv")

sinas_places <- terra::vect(
  f_sinas_places
) %>%
  terra::project("epsg:4326")

sinas_data <- read.table(f_sinas_data, sep = " ", header = TRUE)

taxa <- rsdd::taxa()
tax_avail <- taxa$species

#>-----------------------------------------------------------------------------<
#> Prepare biome map
#>
## Create a buffered version of the biomes to exclude species at the margins or
## in mountains

f_bbuff <- file.path(dir_dat, "biomes", "biomes_buff.gpkg")
if (file.exists(f_bbuff) & !recompute) {
  biomes_buff <- terra::vect(f_bbuff)
} else {
  ## Load original biome file (contiguous biome areas in Olson, 2002)
  biomes <- terra::vect(
    file.path(dir_dat, "biomes/biomes.shp")
  )
  
  biomes$ID <- seq(1:nrow(biomes))
  
  ## Create a buffered mountain vector layer to remove mountains, since montane
  ## species are often untypical or even represent a completely different biome
  mountains_buff <- terra::vect(f_mountains)  %>%
    terra::aggregate() %>%
    terra::makeValid() %>%
    terra::buffer(width = 100)
  
  biomes %>%
    terra::buffer(width = -1000) %>%
    terra::erase(mountains_buff) %>%
    terra::writeVector(
      filename = f_bbuff,
      overwrite = TRUE
    )
  biomes_buff <- terra::vect(f_bbuff)
}

## Create a version of biomes without mountains to get their lowland area
f_berased <- file.path(dir_dat, "biomes", "biomes_wout_mountains.gpkg")
if (file.exists(f_berased) & !recompute) {
  biomes_erased <- terra::vect(f_berased)
} else {
  ## Load original biome file (contiguous biome areas in Olson, 2002)
  biomes <- terra::vect(
    file.path(dir_dat, "biomes", "biomes.shp")
  )
  
  biomes$ID <- seq(1:nrow(biomes))
  
  ## Create a buffered mountain vector layer to remove mountains, since montane
  ## species are often untypical or even represent a completely different biome
  biomes %>%
    terra::erase(terra::vect(f_mountains)) %>%
    terra::writeVector(
      filename = f_berased,
      overwrite = TRUE
    )
  biomes_erased <- terra::vect(f_berased)
}

f_bfullinfo <- file.path(dir_dat, "biomes", "biomes_full_info.gpkg")
if (file.exists(f_bfullinfo) & !recompute) {
  biomes <- terra::vect(f_bfullinfo)
} else {
  biomes <- terra::vect(
    file.path(dir_dat, "biomes", "biomes.shp")
  )
  biomes$ID <- seq(1:nrow(biomes))
  biomes$total_area <- terra::expanse(biomes)
  biomes$lowland_area <- data.frame(
    ID = biomes_erased$ID,
    Area = terra::expanse(biomes_erased)
    ) %>%
    dplyr::right_join(
      as.data.frame(biomes), by = "ID"
      ) %>%
    dplyr::pull(Area)
  
  biomes$clusterIDbathy <- terra::vect(f_bathymetrie) %>%
    terra::extract(terra::centroids(biomes, inside = TRUE)) %>%
    dplyr::pull(ID)
  
  terra::writeVector(
    biomes,
    filename = f_bfullinfo,
    overwrite = TRUE
  )
}

# Get biome patch connectivity
if (!"dECA" %in% names(biomes) | recompute) {
  for (biome_id in unique(biomes$BIOME)) {
    cmd <- paste("Rscript get_dECA.R", biome_id)
    system(cmd)
  }
  
  dat <- do.call(
    rbind,
    lapply(
      X = list.files(
        file.path(dir_tmp, "dECA"),
        pattern = ".csv",
        full.names = TRUE
        ),
      FUN = read.csv
    )
  ) %>%
    dplyr::rename(clusterIDmaxdist = cluster)
  
  if (any(duplicated(dat$ID))) {
    stop("Duplicates in dECA file biome patch IDs found.")
  }
  
  biomes <- terra::merge(
    terra::vect(f_bfullinfo),
    dat,
    by = "ID",
    all.x = TRUE
  )
  
  terra::writeVector(
    biomes,
    filename = f_bfullinfo,
    overwrite = TRUE
  )
}

#>-----------------------------------------------------------------------------<
#> Add species richness estimates
f_est_div <- file.path(dir_tmp, "biome_species_richness.csv")
est_div <- read.csv(f_est_div) %>%
  dplyr::rename(ID = biomePatchID)

biomes <- terra::merge(biomes, est_div, by = "ID", all.x = TRUE)
biomes$speciesRichnessBa[which(biomes$speciesRichnessBa > 1e5)] <- NA

vals <- biomes$speciesRichnessChao2
vals_log <- log10(vals)

breaks <- pretty(vals_log, n = 12)
cols <- viridis::viridis(length(breaks) - 1)

col_idx <- cut(vals_log, breaks = breaks, include.lowest = TRUE)
plot_cols <- cols[col_idx]
plot_cols[is.na(vals)] <- "grey50"

png(
  filename = file.path(dir_main, "fig", "SpeciesRichness.png"),
  width = 700, height = 500
  )
plot(
  biomes, col = plot_cols, border = NA,
  main = "Estimated richness of vascular plant species"
  )

legend(
  "left",
  fill = cols,
  legend = signif(10 ^ breaks[-1], 3),
  title = "Species richness"
)
dev.off()


#>-----------------------------------------------------------------------------<
#> Assign species to biome patches
assign_patches <- function(taxon) {
  ## Ensure to use a valid taxon name. Since alternative names ans subspecies
  ## were aggregated in rsdd, we use species-level names. This may introduce
  ## some errror if multiple subspecies are native or invasive within the same
  ## region and the wrong one is selected, but this will likely result in minor
  ## noise.
  taxon_species_lvl <- paste(
    strsplit(taxon, split = " ")[[1]][c(1, 2)],
    collapse = " "
    )
  if (taxon %in% tax_avail) {
    rsdd_taxon <- taxon
  } else if (taxon_species_lvl %in% tax_avail) {
    rsdd_taxon <- taxon_species_lvl
  } else {
    return(data.frame())
  }
  
  ds <- sinas_data[sinas_data$taxon == taxon, ]
  native_loc_ids <- ds[ds$establishmentMeans == "native", ]$location
  introduced_loc_ids <- ds[ds$establishmentMeans == "introduced", ]$location
  
  if (length(native_loc_ids) == 0 | length(introduced_loc_ids) == 0) {
    return(data.frame())
  }
  
  native <- sinas_places[sinas_places$Location %in% native_loc_ids, ]
  introduced <- sinas_places[sinas_places$Location %in% introduced_loc_ids, ]
  
  if (terra::is.empty(native) | terra::is.empty(introduced)) {
    return(data.frame())
  }
  
  obs <- rsdd::get_taxon(
    rsdd_taxon, status = "all", # Here, we use all observations to employ the SINaS native/Exotic definitions
    format = "centroids"
    ) %>%
    terra::vect(geom = c("x", "y"), crs = "epsg:4326")
  
  native_obs <- terra::intersect(obs, native)
  introduced_obs <- terra::intersect(obs, introduced)
  
  if (length(native_obs) < 5 | length(introduced_obs) < 5) {
    return(data.frame())
  }
  
  # Find patches of interest
  native_patches <- biomes_buff[
    which(
      terra::is.related(
        biomes_buff, native_obs, relation = "covers"
      )
    ),
  ]
  if (length(native_patches) == 0) {
    return(data.frame())
  }
  
  # Obtain number of observations and area
  native_patches$obs_count <- terra::relate(
    native_patches, native_obs, "contains"
  ) %>%
    rowSums() %>%
    as.vector()
  
  # Find patches of interest
  introduced_patches <- biomes_buff[
    which(
      terra::is.related(
        biomes_buff, introduced_obs, relation = "covers"
      )
    ),
  ]
  if (length(introduced_patches) == 0) {
    return(data.frame())
  }
  
  # Obtain number of observations and area
  introduced_patches$obs_count <- terra::relate(
    introduced_patches, introduced_obs, "contains"
  ) %>%
    rowSums() %>%
    as.vector()
  
  # Create data frames
  df_origin <- data.frame(
    Species = rep(taxon, length(native_patches$ID)),
    Biome = native_patches$BIOME,
    PatchID = native_patches$ID,
    Count = native_patches$obs_count,
    Status = rep("Donor", length(native_patches$ID))
  )
  
  df_destination <- data.frame(
    Species = rep(taxon, length(introduced_patches$ID)),
    Biome = introduced_patches$BIOME,
    PatchID = introduced_patches$ID,
    Count = introduced_patches$obs_count,
    Status = rep("Receiver", length(introduced_patches$ID))
  )
  
  df <- rbind(df_origin, df_destination)
  orig_biomes <- df_origin %>%
    dplyr::group_by(Biome) %>%
    dplyr::summarise(MaxCount = max(Count)) %>%
    dplyr::filter(MaxCount == max(MaxCount)) %>%
    dplyr::pull(Biome)
  
  df$MainOrigin <- ifelse(length(orig_biomes) == 1, orig_biomes, NA)
  
  return(df)
}


specs <- unique(sinas_data$taxon)
#specs <- sample(specs, size = 20000)

f_out <- file.path(dir_dat, "df_species_patches.csv")
if (file.exists(f_out) & !recompute) {
  df_species_patches <- read.csv(f_out)
} else {
  Sys.time()
  df_species_patches <- do.call(
    rbind,
    pbapply::pblapply(
      X = specs,
      FUN = function(spec) {
        tryCatch(
          assign_patches(spec),
          error = function(e) {
            stop(paste("Error for species", spec, ":", e$message))
          }
        )
      }
    )
  )
  
  write.csv(df_species_patches, file = f_out, row.names = FALSE)
}

merged <- terra::merge(
  dplyr::rename(df_species_patches, ID = PatchID),
  as.data.frame(biomes),
  by = "ID",
  all.x = TRUE
  ) %>%
  dplyr::mutate(
    clusterID = paste(
      as.character(Biome),
      as.character(clusterIDmaxdist),
      as.character(clusterIDbathy),
      sep = "-"
    )
  ) %>%
  dplyr::mutate(
    total_area_km2 = total_area / 1e6,
    lowland_area_km2 = lowland_area / 1e6
  )

merged$Biome <- factor(merged$Biome)
merged$clusterID <- factor(merged$clusterID)

#>=============================================================================<
#> Plot maps and boxplots
#<=============================================================================>

# Plot maps of donor and receiver regions of introduced species
map_stats <- merged %>%
  dplyr::group_by(ID, Status) %>%
  dplyr::summarise(
    clusterID = dplyr::first(clusterID),
    biome = dplyr::first(Biome),
    total_area = dplyr::first(total_area),
    lowland_area = dplyr::first(lowland_area),
    species_richness = dplyr::first(speciesRichnessChao2), # Select Chao2 (best estimate)
    species_count = dplyr::n_distinct(Species, na.rm = TRUE)
  ) %>%
  dplyr::group_by(clusterID) %>%
  dplyr::mutate(
    max_cluster_richness = if(
      all(is.na(species_richness))
      ) NA else max(
        species_richness, na.rm = TRUE
        )
    ) %>%
  dplyr::ungroup()


remerge <- terra::merge(biomes, map_stats, by = "ID", all.x = TRUE)

polygons_sf <- sf::st_as_sf(remerge)
gg_dr <- ggplot2::ggplot(polygons_sf[which(!is.na(polygons_sf$Status)),]) +
  ggplot2::geom_sf(data = biomes, colour = "black", fill = NA) +
  ggplot2::geom_sf(aes(fill = species_count)) +
  ggplot2::facet_wrap(~ Status, ncol = 2) +
  ggplot2::scale_fill_viridis_c(
    name = "Species count",
    option = "viridis",
    na.value = "grey90"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "right",
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(dir_main, "fig", "SInAS_donors_and_receivers.png"),
  plot = gg_dr,
  width = 13, height = 3
  )

#-------------------------------------------------------------------------------
# Plot boxplots
boxplot_data <- merged %>%
  dplyr::select(
    Species, Status, lowland_area_km2, speciesRichnessChao2
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::summarise(
    max_lowland_km2 = if(
      all(is.na(lowland_area_km2))
    ) NA else max(
      lowland_area_km2, na.rm = TRUE
    ),
    max_richness = if(
      all(is.na(speciesRichnessChao2))
    ) NA else max(
      speciesRichnessChao2, na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(max_lowland_1e3km2 = max_lowland_km2 / 1e3) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(
    cols = c(max_lowland_1e3km2, max_richness),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::filter(
    is.finite(value)
  ) %>%
  dplyr::group_by(metric, Species) %>%
  dplyr::filter(
    all(c("Donor", "Receiver") %in% Status)
  ) %>%
  dplyr::ungroup()

# Calculate statistics
library(rstatix)
## Calculate statistics with log-scale aware positioning
metric_max_values <- boxplot_data %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    max_value = max(value, na.rm = TRUE),
    log_max_value = log10(max(value, na.rm = TRUE))
  )

## Wilcoxon test
test_res <- boxplot_data |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_test(
    value ~ Status,
    paired = TRUE,
    alternative = "greater"
  ) |>
  rstatix::adjust_pvalue(method = "BH")

## Effect size
eff_res <- boxplot_data |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_effsize(
    value ~ Status,
    paired = TRUE
  )

## Create stats dataframe - using log-space coordinates
stats <- test_res %>%
  dplyr::left_join(eff_res %>% select(metric, effsize), by = "metric") %>%
  dplyr::left_join(metric_max_values, by = "metric") %>%
  dplyr::mutate(
    label = paste0(
      "W = ", statistic,
      "\n", "p = ", signif(p, 3),
      "\n", "r = ", round(effsize, 2)
    ),
    y.position = max_value,
    .y. = "value"
  ) %>%
  dplyr::select(
    metric, group1, group2, statistic, p, effsize, label, y.position, .y.
  ) %>%
  tidyr::pivot_longer(
    cols = c(group1, group2),
    names_to = "Group",
    values_to = "Status"
  ) %>%
  dplyr::mutate(value = y.position)

## Simple boxplots by Status (facetted by metric)
ggplot2::ggplot(boxplot_data, aes(x = Status, y = value, fill = Status)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    ~ metric, scales = "free_y",
    labeller = as_labeller(
      c(
        "max_lowland_1e3km2" = "Maximum Lowland Area",
        "max_richness" = "Maximum Species Richness (Chao2)"
      )
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c("Donor" = "#2E86AB", "Receiver" = "#A23B72")
  ) +
  # Boxplots are shown on a logarithmic scale to accommodate the large dynamic range of values.
  #ggplot2::scale_y_log10() +
  ggplot2::coord_transform(y = "log10") +
  ggplot2::labs(
    title = NULL,
    y = expression("Area (" * 10^3 * " km"^2 * ") / Estimated species count"),
    x = "Biome patch status"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_rect(fill = "lightgray"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::geom_text(
    data = stats,
    ggplot2::aes(x = 1.5, label = label),
    vjust = 7
    ) +
  ggplot2::scale_y_continuous(
    breaks = c(1e-2, 1, 1e2, 1e4),
    labels = scales::label_scientific()
  )

#>=============================================================================<
#> Add sPlot data
#<=============================================================================>
# sPlot Open


#>=============================================================================<
#> Statistical analyses
#<=============================================================================>

cat(
  ">--------------------Statistical model outputs--------------------<",
  file = f_glm_stats, append = FALSE
)

stats_received <- merged %>%
  dplyr::filter(Status == "Receiver") %>%
  dplyr::group_by(ID, Biome, Status) %>%
  dplyr::summarise(
    clusterID = dplyr::first(clusterID),
    biome = dplyr::first(Biome),
    total_area = dplyr::first(total_area),
    lowland_area = dplyr::first(lowland_area),
    species_richness = dplyr::first(speciesRichnessChao2), # Select Chao2 (best estimate)
    dECA = dplyr::first(dECA),
    localECA = dplyr::first(localECA),
    clusterECA = dplyr::first(clusterECA),
    introduced_species_count = dplyr::n_distinct(Species, na.rm = TRUE)
  )

stats_received_long <- stats_received %>%
  tidyr::pivot_longer(
    cols = c(total_area, lowland_area, species_richness),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::ungroup()

ggplot2::ggplot(
  data = stats_received_long,
  ggplot2::aes(x = log(value), y = log(introduced_species_count), colour = Biome)
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::facet_wrap(.~metric, scales = "free_x")

library(MASS)
mod <- stats::glm(
  introduced_species_count / species_richness ~ species_richness * lowland_area * dECA + Biome,
  data = stats_received
  )
summary(mod)
