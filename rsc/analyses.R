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
library("rsdd")
library("terra")
library("tidyterra")
library("ggplot2")
library("viridis")
library("dplyr")
library("tidyr")
library("lme4")
library("glmmTMB")
library("boot")
library("future.apply")
library("progressr")

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

f_sinas_places <- file.path(dir_dat, "sinas", "SInAS_Locations", "SInAS_Locations.shp")
f_sinas_data <- file.path(dir_dat, "sinas", "SInAS_3.1.1.csv")

sinas_places <- terra::vect(
  f_sinas_places
) %>%
  terra::project("epsg:4326")

sinas_places <- sf::st_read(
  f_sinas_places
) %>%
  sf::st_transform(4326) %>%
  terra::vect()

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
est_div <- read.csv(f_est_div)

biomes <- terra::merge(biomes, est_div, by = "ID", all.x = TRUE)
biomes$speciesRichnessBa[which(biomes$speciesRichnessBa > 1e5)] <- NA

vals <- biomes$speciesRichnessBa
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

# just to check how many names can be matched somehow
# check_avail <- function(taxon) {
#   taxon_species_lvl <- paste(
#     strsplit(taxon, split = " ")[[1]][c(1, 2)],
#     collapse = " "
#   )
#   if (taxon %in% tax_avail) {
#     return(1)
#   } else if (taxon_species_lvl %in% tax_avail) {
#     return(1)
#   } else {
#     return(0)
#   }
# }
# n_avail <- sum(unlist(lapply(X = specs, FUN = check_avail)))


f_out <- file.path(dir_dat, "df_species_patches.csv")
if (file.exists(f_out) & !recompute) {
  df_species_patches <- read.csv(f_out)
} else {
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  future::plan(future::multicore)
  Sys.time()
  progressr::with_progress({
    p <- progressr::progressor(along = specs)
    df_species_patches <- dplyr::bind_rows(
      future.apply::future_lapply(
        X = specs,
        FUN = function(spec) {
          tryCatch(
            {
              res <- assign_patches(spec)
              p()
              res
            },
            error = function(e) {
              message(paste("Error for species", spec, ":", e$message))
              NULL
            }
          )
        }
      )
    )
  })
  
  write.csv(df_species_patches, file = f_out, row.names = FALSE)
}

#>=============================================================================<
#> Add GIST status for the most invasive species
#<=============================================================================>
gisd_codes <- c("MV", "MR", "MO", "MN", "MC")
gisd_levels <- c("Massive", "Major", "Moderate", "Minor", "Minimal concern")
gisd <- read.csv(
  file.path(dir_dat, "db", "gisd_2026-02-25", "converted.csv")
) %>%
  dplyr::filter(EICAT.Category %in% gisd_codes) %>%
  dplyr::mutate(
    Countries.of.impact = dplyr::recode(
      Countries.of.impact, "Ethiopa" = "Ethiopia"
    )
  )

worst <- read.csv(
  file.path(dir_dat, "db", "gisd_2026-02-25", "100_of_the_Worlds_Worst.csv")
  )$Species

check_gist <- function(taxon, biome_patch_id) {
  if (taxon %in% gisd$Species) {
    locations <- sinas_places[
      which(
        terra::is.related(
          sinas_places,
          biomes[which(biomes$ID == biome_patch_id),],
          relation = "intersects"
        )
      ),
    ]$Location
    
    categories <- gisd[
      gisd$Countries.of.impact %in% locations & gisd$Species == taxon,
    ]$EICAT.Category
    if (length(categories) < 1) {
      return(NA)
    }
    if (length(categories) == 1) {
      categories
    }
    return(min(match(categories, gisd_codes)))
  } else {
    return(NA)
  }
}


df_species_patches <- df_species_patches %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gisd_impact_level = check_gist(Species, PatchID)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    gist_impact_label = ifelse(
      is.na(gisd_impact_level),
      "Not evaluated",
      gisd_levels[gisd_impact_level]
    )
  ) %>%
  dplyr::mutate(
    gist_impact_label = factor(
      gist_impact_label, levels = c(gisd_levels, "Not evaluated")
      ),
    gist_100_worst = Species %in% worst
  )

#>=============================================================================<
#> Merge with biome patch information
#<=============================================================================>
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

# Plot
ggdf <- merged |>
  dplyr::filter(Status == "Receiver") |>
  dplyr::add_count(gist_impact_label, name = "N") |>
  dplyr::mutate(
    gist_label_N = paste0(gist_impact_label, "\nN = ", N),
    gist_label_N = factor(
      gist_label_N,
      levels = paste0(
        levels(gist_impact_label),
        "\nN = ",
        tapply(N, gist_impact_label, unique)
        )
      )
    )

gg_gist <- ggplot2::ggplot(
  data = ggdf,
  ggplot2::aes(x = gist_label_N, y = lowland_area_km2, fill = gist_impact_label)
) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab("GIST impact class") +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(
    values = c("firebrick3", "coral3", "lightsalmon3", "darkgoldenrod3", "slateblue3")
    ) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_y_log10(
    name = expression("Area in km"^2~"(log scaled)"),
    labels = scales::label_number()
    )

ggplot2::ggsave(
  filename = file.path(dir_main, "fig", "BoxplotGISDclasses.svg"),
  plot = gg_gist,
  width = 7, height = 5
)

gg_gist_Ba <- ggplot2::ggplot(
  data = ggdf,
  ggplot2::aes(x = gist_label_N, y = speciesRichnessBa, fill = gist_impact_label)
) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab("GIST impact class") +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(
    values = c("firebrick3", "coral3", "lightsalmon3", "darkgoldenrod3", "slateblue3")
  ) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_y_log10(
    name = "Estimated species richness (log scaled)",
    labels = scales::label_number()
  )

ggplot2::ggsave(
  filename = file.path(dir_main, "fig", "BoxplotGISDclassesBa.svg"),
  plot = gg_gist_Ba,
  width = 7, height = 5
)

gg_gist_worst <- ggplot2::ggplot(
  data = ggdf,
  ggplot2::aes(x = gist_100_worst, y = log(speciesRichnessBa), fill = gist_100_worst)
) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab("GIST assessment") +
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(
    labels = c(
      "TRUE" = "GIST 100 worst",
      "FALSE" = "Other"
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c("firebrick3", "slateblue3")
  ) +
  ggplot2::theme(legend.position = "none")


ggplot2::ggsave(
  filename = file.path(dir_main, "fig", "BoxplotGISD100worst.svg"),
  plot = gg_gist_worst,
  width = 5, height = 5
)

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
    species_richness = dplyr::first(speciesRichnessBa),
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
    Species, Status, lowland_area_km2, speciesRichnessBa
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::summarise(
    max_lowland_km2 = if(
      all(is.na(lowland_area_km2))
    ) NA else max(
      lowland_area_km2, na.rm = TRUE
    ),
    max_richness = if(
      all(is.na(speciesRichnessBa))
    ) NA else max(
      speciesRichnessBa, na.rm = TRUE
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
      "\n", "effect size = ", round(effsize, 2)
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
gg_area <- ggplot2::ggplot(boxplot_data, aes(x = Status, y = value, fill = Status)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    ~ metric, scales = "free_y",
    labeller = as_labeller(
      c(
        "max_lowland_1e3km2" = "Maximum Lowland Area",
        "max_richness" = "Maximum Species Richness"
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
    y = expression("Area (" * 10^3 * " km"^2 * ")"),
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
    vjust = 4.5
    ) +
  ggplot2::scale_y_continuous(
    breaks = c(1e-2, 1, 1e2, 1e4),
    labels = scales::label_scientific(),
    sec.axis = sec_axis(
      ~ ., name = "Estimated species count", breaks = NULL, labels = NULL
      )
  )

ggplot2::ggsave(
  filename = file.path(dir_main, "fig", "MaxAreaRichnessBoxplot.svg"),
  plot = gg_area,
  width = 10, height = 4
)

#>=============================================================================<
#> Add sPlot data
#<=============================================================================>
# sPlot Open


#>=============================================================================<
#> Statistical analyses
#<=============================================================================>

df_donor <- merged %>%
  dplyr::mutate(
    PatchID = ID,
    Area = lowland_area,
    SpeciesRichness = round(speciesRichnessBa, 0)
    ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::filter(
    Status == "Donor",
    Count == max(Count),
    SpeciesRichness > 50
    ) %>%
  dplyr::group_by(Biome, PatchID, Area, SpeciesRichness, Status) %>%
  dplyr::summarise(
    N = dplyr::n()
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(Biome, PatchID, Area, SpeciesRichness),
    names_from = Status,
    values_from = N
  ) %>%
  tidyr::replace_na(replace = list(Donor = 0, Receiver = 0)) %>%
  dplyr::rename(
    DonorOf = Donor
  )

df_within$failures <- pmax(df_within$SpeciesRichness - df_within$DonorOf, 0) 

## GLM that includes biome
modInvasiveBreeder_glm <- stats::glm(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + factor(Biome),
  data = df_within,
  family = stats::quasibinomial()
)

summary(modInvasiveBreeder_glm)
coef_tab_glm <- summary(modInvasiveBreeder_glm)$coefficients

## GLMM that controls for biome
cor(df_donor$Area, df_donor$SpeciesRichness, use = "complete.obs")
modInvasiveBreeder_glmm <- lme4::glmer(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + (1 | Biome),
  data = df_within,
  family = binomial
)

summary(modInvasiveBreeder_glmm)
anova(modInvasiveBreeder_glmm)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  sum(rp^2) / rdf
}

overdisp_fun(modInvasiveBreeder_glmm)

## Betabinomial GLMM
modInvasiveBreeder_bb <- glmmTMB::glmmTMB(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + (1 | Biome),
  data = df_within,
  family = glmmTMB::betabinomial()
)

summary(modInvasiveBreeder_bb)


library(DHARMa)
res <- DHARMa::simulateResiduals(modInvasiveBreeder_bb)
DHARMa::testDispersion(res)

DHARMa::plotResiduals(res)
DHARMa::testZeroInflation(res)







df_within <- merged %>%
  dplyr::mutate(
    PatchID = ID,
    Area = lowland_area,
    SpeciesRichness = round(speciesRichnessBa, 0)
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::filter(
    #Biome == MainOrigin,
    Count > 10,
    SpeciesRichness > 50
  ) %>%
  dplyr::group_by(Biome, Species, clusterID, Status) %>%
  dplyr::summarise(
    Count = sum(Count),
    clusterRichness = max(speciesRichnessBa, na.rm = TRUE),
    localECA = max(localECA, na.rm = TRUE),
    clusterECA = max(clusterECA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Biome, Species, Status) %>%
  dplyr::filter(
    Count == max(Count, na.rm = TRUE)
  ) %>%
  dplyr::filter(# Addiditonal filter, since in some cases, multiple clusters have max(Count)
    clusterRichness == max(clusterRichness)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Biome, Species) %>%
  dplyr::filter(
    !duplicated(clusterID),
    dplyr::n_distinct(Status) == 2
    ) %>%
  dplyr::ungroup()

ggplot2::ggplot(
  data = df_within,
  ggplot2::aes(y = log(clusterRichness), colour = Status)
  ) +
  ggplot2::geom_boxplot()

df_wide_richness <- df_within %>%
  dplyr::select(Biome, Species, Status, clusterRichness) %>%
  tidyr::pivot_wider(
    names_from = Status,
    values_from = clusterRichness
  ) %>%
  dplyr::filter(!is.na(Donor) & !is.na(Receiver)) %>%
  dplyr::mutate(diff = Donor - Receiver)

# Wilcoxon signed-rank test (paired)
wilcox.test(
  df_wide_richness$Donor,
  df_wide_richness$Receiver,
  paired = TRUE,
  alternative = "greater"
)

# Bootstrap effect size (median difference)
boot_median <- boot::boot(
  data = df_wide_richness$diff,
  statistic = function(x, i) median(x[i]),
  R = 10000
)
boot::boot.ci(boot_median, type = "perc")

df_wide_richness <- df_within %>%
  dplyr::select(Biome, Species, Status, clusterRichness) %>%
  tidyr::pivot_wider(
    names_from = Status,
    values_from = clusterRichness
  ) %>%
  dplyr::filter(!is.na(Donor) & !is.na(Receiver)) %>%
  dplyr::mutate(
    diff = Donor - Receiver,
    ratio = Receiver / Donor
    )

# Wilcoxon signed-rank test (paired)
wilcox.test(
  df_wide_richness$Donor,
  df_wide_richness$Receiver,
  paired = TRUE,
  alternative = "greater"
)

# Bootstrap effect size (median difference)
boot_median <- boot::boot(
  data = df_wide_richness$diff,
  statistic = function(x, i) median(x[i]),
  R = 10000
)
boot::boot.ci(boot_median, type = "perc")

# Same for cluster ECA
df_wide_ECA <- df_within %>%
  dplyr::select(Biome, Species, Status, clusterECA) %>%
  tidyr::pivot_wider(
    names_from = Status,
    values_from = clusterECA
  ) %>%
  dplyr::filter(!is.na(Donor) & !is.na(Receiver)) %>%
  dplyr::mutate(diff = Donor - Receiver)

# Wilcoxon signed-rank test (paired)
wilcox.test(
  df_wide_ECA$Donor,
  df_wide_ECA$Receiver,
  paired = TRUE,
  alternative = "greater"
)

# Bootstrap effect size (median difference)
boot_median <- boot::boot(
  data = df_wide_ECA$diff,
  statistic = function(x, i) median(x[i]),
  R = 10000
)
boot::boot.ci(boot_median, type = "perc")
