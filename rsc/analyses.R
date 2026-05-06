#>=============================================================================<
#>-----------------------------------------------------------------------------<
#> Analysis based on GBIF data of SINAS or GloNAF-listed (invasive) species
#> 
#> Date: 2026-01-27
#> Author: Popp, MR
#> 
#>----------------------------------------------------------------------------->
#<=============================================================================>

print(Sys.time())
library("rsdd")
library("ggplot2")
library("viridis")
library("dplyr")
library("tidyr")
library("terra")
library("tidyterra")
library("rgbif")
library("lme4")
library("mgcv")
library("glmmTMB")
library("boot")
library("coin")
library("rstatix")
library("performance")
library("MuMIn") # Alternative to performance
library("DHARMa")
library("collapse")
library("future.apply")
library("progressr")


#>=============================================================================<
#> Settings
#<=============================================================================>

recompute <- FALSE
set.seed(161)

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
  "Mangrove",
  rep(NA, 83),
  "Lake",
  "Rock and Ice"
)

nativeness_source <- "sinas"
biome_def <- "olson"

# Main directories
if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir"
  dir_lud11 <- "L:"
} else {
  dir_main <- "/lud11/poppman/data/bir"
  dir_lud11 <- "/lud11"
}

dir_dat <- file.path(dir_lud11, "poppman", "data", "bir", "dat", "lud11")
dir_stats <- file.path(dir_main, "stats")
dir_imeb <- file.path(dir_dat, "biomes", biome_def, "intermediate_data")  # Only biome-level data
dir_imed <- file.path(dir_dat, "biomes", biome_def, "intermediate_data", nativeness_source)  # Patch-level data
dir_fig <- file.path(dir_main, "fig", nativeness_source)

if (!dir.exists(dir_fig)) {
  dir.create(dir_fig, recursive = TRUE)
}

# Invasive species data base files
f_sinas_places <- file.path(dir_dat, "sinas", "SInAS_locations_3.1.gpkg")
f_sinas_data <- file.path(dir_dat, "sinas", "SInAS_3.1.1.csv")
f_sinas_taxon_match <- file.path(dir_dat, "sinas", "taxonomy_matching.csv")

f_glonaf_places <- file.path(dir_dat, "glonaf", "glonaf_2024_regions", "glonaf_2024_regions.shp")
f_glonaf_data <- file.path(dir_dat, "glonaf", "glonaf_flora2.csv")
f_glonaf_taxa <- file.path(dir_dat, "glonaf", "glonaf_taxon_wcvp.csv")

# Original biomes file
f_biomes <- file.path(dir_dat, "biomes", biome_def, "biomes.shp")

# Intermediate data files
f_bbuff <- file.path(dir_imeb, "biomes_buff.gpkg")
f_berased <- file.path(dir_imeb, "biomes_wout_mountains.gpkg") # Just w/out mountains for area calculation
f_bfullinfo <- file.path(dir_imeb, "biomes_full_info.gpkg") # Original biomes with added information
f_out <- file.path(dir_imed, "df_species_patches.rds") # Output of invasive species assignment to patches

# Miscancellous data files
f_mountains <- file.path(
  dir_dat,
  "shp", "GMBA_Inventory_v2.0_standard", "GMBA_Inventory_v2.0_standard.shp"
  )
f_bathymetrie <- file.path(
  dir_dat,
  "gebco_2025_geotiff", "gebco.gpkg"
  )

f_est_div <- file.path( # Species richness estimates
  dir_imeb, "biome_species_richness.csv"
  )

f_env_dat <- file.path( # Environmental data
  dir_imeb, "biomes_env.rds"
)

f_spl_eft <- file.path( # Total GBIF observation counts
  dir_imeb, "sampling_effort.csv"
)

f_hum_mod <- file.path(
  dir_imeb, "human_modification.csv"
)

f_climate_overlap <- file.path(dir_imeb, "climate_overlap.csv")
f_climate_restime <- file.path(dir_imeb, "clim_res_time.csv")

f_gisd <- file.path(dir_dat, "db", "gisd_2026-02-25", "converted.csv")
f_100_worst <- file.path(
  dir_dat, "db", "gisd_2026-02-25", "100_of_the_Worlds_Worst.csv"
  )

f_sPlot <- file.path(
  #dir_lud11, "poppman/data/bir/dat/lud11/sPlotOpen/sPlotOpen.RData"
  "D:/sPlot_closed", "extracted_data.RData"
)

f_gadm <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/db", "gadm_countries.gpkg"
)

# Environmental raster data
f_envstack <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/environment", "environmentPC.tif"
)
f_tasmean <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/environment", "tasmean_clim.tif"
  )
f_prec <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/environment", "pr_clim.tif"
)
f_ghm <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/environment", "gHM_resampled_CHELSA.tif"
)
f_npp <- file.path(
  dir_lud11, "poppman/data/bir/dat/lud11/environment", "npp_mean.tif"
)

# Create missing directories
if (!dir.exists(dir_imed)) {
  dir.create(dir_imed, showWarnings = FALSE, recursive = TRUE)
}


#>=============================================================================<
#> Functions
#<=============================================================================>

## Create a model formula from a vector of predictors
make_formula <- function(
    predictors, k = 3, response = "Invaded", biome = TRUE
    ) {
  smooth_terms <- sprintf("s(%s, k = %s)", predictors, k)
  rhs <- paste(smooth_terms, collapse = " + ")
  if (biome) {
    rhs <- paste(rhs, "+ Biome")
  }
  as.formula(paste(response, "~", rhs))
}

## Fail-save computation of Bhattacharyya distances
save_bhattacharyya <- function(mu1, mu2, Sigma1, Sigma2) {
  fail <- function(reason) {
    list(value = NA_real_, reason = reason)
  }
  if (is.null(mu1) || is.null(mu2) || is.null(Sigma1) || is.null(Sigma2)) {
    return(fail("null_input"))
  }
  
  mu1 <- as.numeric(mu1)
  mu2 <- as.numeric(mu2)
  Sigma1 <- tryCatch(as.matrix(Sigma1), error = function(e) NULL)
  Sigma2 <- tryCatch(as.matrix(Sigma2), error = function(e) NULL)
  
  if (is.null(Sigma1) || is.null(Sigma2)) {
    return(fail("matrix_coercion_failed"))
  }
  
  # Dimension checks
  if (length(mu1) != length(mu2)) {
    return(fail("mean_dim_mismatch"))
  }
  
  if (
    nrow(Sigma1) != ncol(Sigma1) ||
    nrow(Sigma2) != ncol(Sigma2) ||
    nrow(Sigma1) != length(mu1) ||
    nrow(Sigma2) != length(mu2)
  ) {
    return(fail("cov_dim_mismatch"))
  }
  
  if (
    any(!is.finite(mu1)) ||
    any(!is.finite(mu2)) ||
    any(!is.finite(Sigma1)) ||
    any(!is.finite(Sigma2))
  ) {
    return(fail("non_finite_values"))
  }
  
  if (det(Sigma1) < 1e-12 || det(Sigma2) < 1e-12) {
    return(fail("singular_covariance"))
  }
  
  val <- tryCatch(
    fpc::bhattacharyya.dist(mu1, mu2, Sigma1, Sigma2),
    error = function(e) {NA_real_}
  )
  
  if (is.na(val)) {
    return(fail("numerical_failure"))
  }
  
  if (!is.finite(val)) {
    return(fail("infinite_result"))
  }
  
  return(list(value = val, reason = "ok"))
}


#>=============================================================================<
#> Data preparation
#<=============================================================================>

## GBIF
rsdd::dataset("gbif-powo_raw")
taxa <- rsdd::taxa()
tax_avail <- taxa$species

## SINaS
sinas_places <- terra::vect(
  f_sinas_places
)

sinas_data <- read.table(f_sinas_data, sep = " ", header = TRUE)

## GloNaF
glonaf_places <- terra::vect(
  f_glonaf_places
)

glonaf_taxa <- read.csv(f_glonaf_taxa, header = TRUE)

if (!"taxon_gbif" %in% names(glonaf_taxa) | all(is.na(glonaf_taxa$taxon_gbif))) {
  glonaf_taxa$taxon_gbif <- NA
  failures <- c()
  pb <- progress::progress_bar$new(total = nrow(glonaf_taxa))
  for (i in 1:nrow(glonaf_taxa)) {
    taxon_i <- glonaf_taxa$taxa_and_authors_accepted[i]
    if (is.na(taxon_i)) {
      taxon_i <- glonaf_taxa$taxon_corrected[i]
    }
    gbif_match <- rgbif::name_backbone(
      name = taxon_i,
      rank = glonaf_taxa$wcvp_taxon_rank[i]
    )
    exact <- dplyr::filter(gbif_match, matchType == "EXACT")
    if (nrow(exact) == 1) {
      glonaf_taxa$taxon_gbif[i] <- exact$canonicalName[1]
      cat("Matched")
    } else {
      failures <- c(failures, i)
      #cat("\nFailed for taxon", glonaf_taxa$taxa_and_authors_accepted[i])
    }
    pb$tick()
  }
  pb <- progress::progress_bar$new(total = length(failures))
  for (i in failures) {
    taxon_i <- glonaf_taxa$taxa_and_authors_accepted[i]
    if (is.na(taxon_i)) {
      taxon_i <- glonaf_taxa$taxon_corrected[i]
    }
    gbif_match <- rgbif::name_backbone(
      name = taxon_i,
      rank = "Species"
    )
    exact <- dplyr::filter(gbif_match, matchType == "EXACT")
    if (nrow(exact) == 1) {
      glonaf_taxa$taxon_gbif[i] <- exact$canonicalName[1]
      cat("Matched")
    } else {
      #cat("\nFailed for taxon", glonaf_taxa$taxa_and_authors_accepted[i])
    }
    pb$tick()
  }
  rm(pb)
  write.csv(glonaf_taxa, file = f_glonaf_taxa, row.names = FALSE)
}

glonaf_taxa <- glonaf_taxa %>%
  dplyr::filter(!is.na(taxon_gbif))

glonaf_data <- read.csv(f_glonaf_data, header = TRUE) %>%
  dplyr::mutate(
    taxon = glonaf_taxa[match(id, glonaf_taxa$id), "taxon_gbif"],
    wcvp_taxon = glonaf_taxa[match(id, glonaf_taxa$id), "taxon_corrected"],
    wcvp_taxon_short = stringr::str_extract(
      glonaf_taxa[match(id, glonaf_taxa$id), "taxon_corrected"],
      "^\\S+\\s+\\S+"
      ),
    rank = glonaf_taxa[match(id, glonaf_taxa$id), "wcvp_taxon_rank"]
  ) #%>%
  #dplyr::filter(rank == "Species") %>%  # Exclude subspecies

env_stack <- terra::rast(f_envstack)
tasmean <- terra::rast(f_tasmean)
prec <- terra::rast(f_prec)
ghm <- terra::rast(f_ghm)
npp <- terra::rast(f_npp)

#>-----------------------------------------------------------------------------<
#> Prepare biome map
#>
## Create a buffered version of the biomes to exclude species at the margins or
## in mountains
if (file.exists(f_bbuff) & !recompute) {
  biomes_buff <- terra::vect(f_bbuff)
} else {
  ## Load original biome file (contiguous biome areas in Olson, 2002)
  biomes_original <- terra::vect(f_biomes)
  
  biomes_original$ID <- seq(1:nrow(biomes_original))
  
  ## Create a buffered mountain vector layer to remove mountains, since montane
  ## species are often untypical or even represent a completely different biome
  mountains_buff <- terra::vect(f_mountains)  %>%
    terra::aggregate() %>%
    terra::makeValid() %>%
    terra::buffer(width = 100)
  
  biomes_original %>%
    terra::buffer(width = -1000) %>%
    terra::erase(mountains_buff) %>%
    terra::writeVector(
      filename = f_bbuff,
      overwrite = TRUE
    )
  biomes_buff <- terra::vect(f_bbuff)
}

## Create a version of biomes without mountains to get their lowland area
if (file.exists(f_berased) & !recompute) {
  biomes_erased <- terra::vect(f_berased)
} else {
  ## Load original biome file
  biomes_original <- terra::vect(f_biomes)
  
  biomes_original$ID <- seq(1:nrow(biomes_original))
  
  ## Create a buffered mountain vector layer to remove mountains, since montane
  ## species are often untypical or even represent a completely different biome
  biomes_original %>%
    terra::erase(terra::vect(f_mountains)) %>%
    terra::writeVector(
      filename = f_berased,
      overwrite = TRUE
    )
  biomes_erased <- terra::vect(f_berased)
}


if (file.exists(f_bfullinfo) & !recompute) {
  biomes <- terra::vect(f_bfullinfo)
} else {
  biomes <- terra::vect(f_biomes)
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
## Note: Connectedness is calculated as the distance-weighted area of all other
## patches of the same biome with an exponential decay function. A combined
## measure of "local ECA" could be simply the sum of patch area and
## distance-weighted neighbour patch area: focalECA = connectedness + area
if (!"dECA" %in% names(biomes) | recompute) {
  for (biome_id in unique(biomes$BIOME)) {
    cmd <- paste("Rscript get_dECA.R", biome_id)
    system(cmd)
  }
  
  dat <- do.call(
    rbind,
    lapply(
      X = list.files(
        file.path(dir_imeb, "dECA"),
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
  
  terra::merge(
    terra::vect(f_bfullinfo)[, 1:8],
    dat,
    by = "ID",
    all.x = TRUE
  ) %>%
    dplyr::mutate(
      focalECA = max(total_area, 0, na.rm = TRUE) + connectedness
    ) %>%
    terra::writeVector(
      filename = f_bfullinfo,
      overwrite = TRUE
      )
}

#>-----------------------------------------------------------------------------<
#> Add species richness estimates, climate stability, sampling effort, etc
est_div <- read.csv(f_est_div)
spl_eff <- read.csv(f_spl_eft)
hum_mod <- read.csv(f_hum_mod)

clim_stab <- read.csv(f_climate_overlap) %>%
  dplyr::rename(
    min_clim_similarity = pminolap,
    clim_integral = pcumolap,
    bottleneck_clim_integral = pcumbnolap,
    min_clim_cluster_similarity = cminolap
  )

clim_rest <- read.csv(f_climate_restime) %>%
  dplyr::rename(
    climate_velocity_kmpa = v, # mean climate velocity (km/yr)
    circle_diam_km = d, # diameter of the equivalent circle (km)
    climate_restime_a = resTim # residence time (years) as the ratio D/vel
  ) %>%
  dplyr::select(ID, climate_velocity_kmpa, climate_restime_a)

biomes <- terra::vect(f_bfullinfo) %>%
  dplyr::mutate(
    clusterID = paste(
      as.character(BIOME),
      as.character(clusterIDmaxdist),
      as.character(clusterIDbathy),
      sep = "-"
    )
  ) %>%
  terra::merge(est_div, by = "ID", all.x = TRUE) %>%
  terra::merge(spl_eff, by = "ID", all.x = TRUE) %>%
  terra::merge(hum_mod, by = "ID", all.x = TRUE) %>%
  terra::merge(clim_stab, by = "ID", all.x = TRUE) %>%
  terra::merge(clim_rest, by = "ID", all.x = TRUE) %>%
  dplyr::mutate(
    log_total_area = log(total_area),
    log_SpeciesRichness = log(speciesRichnessBa),
    log_Connectedness = log(connectedness),
    log_sampling_effort = log(sampling_effort),
    biome_name = biome_names[as.numeric(BIOME)]
  )

#>-----------------------------------------------------------------------------<
#> Correct species richness estimates for sampling bias

biomes$speciesRichnessBa[which(biomes$speciesRichnessBa > 1e5)] <- NA

df_fit <- as.data.frame(biomes) %>%
  dplyr::filter(log_SpeciesRichness > log(100))

m1 <- mgcv::gam(
  log_SpeciesRichness ~ s(log_sampling_effort, k = 3),
  data = df_fit
)
mgcv::gam.check(m1)

effort_fit <- stats::predict(m1)
max_sampling_effect <- stats::predict(
  m1,
  newdata = df_fit %>%
    dplyr::slice_min(
      abs(
        log_sampling_effort - quantile(log_sampling_effort, 0.90, na.rm = TRUE)
      ),
      n = 1
    )
)
biomes$correctedSR <- NA
biomes$correctedSR[which(biomes$log_SpeciesRichness > log(100))] <- data.frame(
  modelled = df_fit$log_SpeciesRichness - effort_fit + max_sampling_effect[1],
  estimated = df_fit$log_SpeciesRichness
) %>%
  dplyr::mutate(
    maximum = pmax(modelled, estimated, na.rm = TRUE)
  ) %>%
  dplyr::pull(maximum)

#>-----------------------------------------------------------------------------<
#> Plot species richness estimates

for (i in 1:2) {
  if (i == 1) {
    vals <- biomes$speciesRichnessBa
    fn <- "EstimatedSpeciesRichness.png"
  } else {
    vals <- biomes$correctedSR
    fn <- "CorrectedSpeciesRichness.png"
  }
  vals_log <- log10(vals)
  
  breaks <- pretty(vals_log, n = 12)
  cols <- viridis::viridis(length(breaks) - 1)
  
  col_idx <- cut(vals_log, breaks = breaks, include.lowest = TRUE)
  plot_cols <- cols[col_idx]
  plot_cols[is.na(vals)] <- "grey50"
  
  png(
    filename = file.path(dir_fig, fn),
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
}

#>-----------------------------------------------------------------------------<
#> Assign species to biome patches

assign_patches <- function(
    taxon, nativeness_source = "powo", load_spatial_data = FALSE
    ) {
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
  
  if (nativeness_source == "sinas") {
    ds <- sinas_data[which(sinas_data$taxon == taxon), ]
    native_loc_ids <- ds[ds$establishmentMeans == "native", ]$locationID
    introduced_loc_ids <- ds[ds$establishmentMeans == "introduced", ]$locationID
    
    if (length(native_loc_ids) == 0 | length(introduced_loc_ids) == 0) {
      return(data.frame())
    }
    
    native <- sinas_places[which(sinas_places$locationID %in% native_loc_ids), ]
    introduced <- sinas_places[which(sinas_places$locationID %in% introduced_loc_ids), ]
    
    if (terra::is.empty(native) | terra::is.empty(introduced)) {
      return(data.frame())
    }
    
    obs <- rsdd::get_taxon(
      rsdd_taxon, status = c("native", "non-native", "N/A"), # Here, we use all observations to employ the SINaS native/Exotic definitions
      format = "centroids"
    ) %>%
      terra::vect(geom = c("x", "y"), crs = "epsg:4326")
    
    native_obs <- terra::intersect(obs, native)
    introduced_obs <- terra::intersect(obs, introduced)
  } else if (nativeness_source == "glonaf") {
    ds <- glonaf_data[which(glonaf_data$taxon == taxon), ]
    introduced_loc_ids <- ds$region_id
    
    introduced <- glonaf_places[which(glonaf_places$OBJIDsic %in% introduced_loc_ids), ]
    
    native_obs <- rsdd::get_taxon(
      rsdd_taxon, status = "native",
      format = "centroids"
    ) %>%
      terra::vect(geom = c("x", "y"), crs = "epsg:4326") %>%
      terra::mask(introduced, inverse = TRUE)
    
    if (terra::is.empty(native_obs)) {
      return(data.frame())
    }
    
    introduced_obs <- rsdd::get_taxon(
      rsdd_taxon, status = c("non-native", "N/A"),
      format = "centroids"
    ) %>%
      terra::vect(geom = c("x", "y"), crs = "epsg:4326") %>%
      terra::intersect(introduced)
    
    if (terra::is.empty(introduced_obs)) {
      return(data.frame())
    }
  } else {
    native_obs <- rsdd::get_taxon(
      rsdd_taxon, status = "native",
      format = "centroids"
    ) %>%
      terra::vect(geom = c("x", "y"), crs = "epsg:4326")
    introduced_obs <- rsdd::get_taxon(
      rsdd_taxon, status = c("non-native", "N/A"),
      format = "centroids"
    ) %>%
      terra::vect(geom = c("x", "y"), crs = "epsg:4326")
  }
  
  if (length(native_obs) < 5 | length(introduced_obs) < 5) {
    return(data.frame())
  }
  
  # Get buffered biomes
  if (load_spatial_data) {
    biomes_buff <- terra::vect(f_bbuff)
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
  
  # Get climate and productivity variables
  if (load_spatial_data) {
    env_stack <- terra::rast(f_envstack)
    tasmean <- terra::rast(f_tasmean)
    prec <- terra::rast(f_prec)
    ghm <- terra::rast(f_ghm)
    npp <- terra::rast(f_npp)
  }
  native_tasmean <- terra::extract(tasmean, native_obs, ID = FALSE)[, 1]
  native_psum <- terra::extract(prec, native_obs, ID = FALSE)[, 1]
  native_ghm <- terra::extract(ghm, native_obs, ID = FALSE)[, 1]
  native_npp <- terra::extract(npp, native_obs, ID = FALSE)[, 1]
  native_env <- terra::extract(env_stack, native_obs, ID = FALSE)
  intro_tasmean <- terra::extract(tasmean, introduced_obs, ID = FALSE)[, 1]
  intro_psum <- terra::extract(prec, introduced_obs, ID = FALSE)[, 1]
  intro_ghm <- terra::extract(ghm, introduced_obs, ID = FALSE)[, 1]
  intro_npp <- terra::extract(npp, introduced_obs, ID = FALSE)[, 1]
  intro_env <- terra::extract(env_stack, introduced_obs, ID = FALSE)
  
  if ("swe" %in% names(native_env)) { # Fill zeros for areas w/out snowfall
    native_env$swe <- tidyr::replace_na(native_env$swe, 0)
    intro_env$swe <- tidyr::replace_na(intro_env$swe, 0)
  }
  
  # Allow only complete observations for covariance matrix
  complete_native_env <- native_env[complete.cases(native_env), , drop = FALSE]
  if (nrow(complete_native_env) >= 2) {
    Sigma_native <- cov(complete_native_env)
  } else {
    Sigma_native <- NA
  }
  centre_native <- colMeans(native_env, na.rm = TRUE)
  
  complete_intro_env <- intro_env[complete.cases(intro_env), , drop = FALSE]
  if (nrow(complete_intro_env) >= 2) {
    Sigma_introduced <- cov(complete_intro_env)
  } else {
    Sigma_introduced <- NA
  }
  centre_introduced <- colMeans(intro_env, na.rm = TRUE)
  
  # Create data frames
  df_origin <- data.frame(
    Species = rep(taxon, length(native_patches$ID)),
    Biome = native_patches$BIOME,
    PatchID = native_patches$ID,
    Count = native_patches$obs_count,
    Status = rep("Donor", length(native_patches$ID)),
    tasmean = mean(native_tasmean, na.rm = TRUE),
    tasmean_sd = sd(native_tasmean, na.rm = TRUE),
    prec = mean(native_psum, na.rm = TRUE),
    prec_sd = sd(native_psum, na.rm = TRUE),
    ghm = mean(native_ghm, na.rm = TRUE),
    npp = mean(native_npp, na.rm = TRUE)
  )
  df_origin$centre <- rep(list(centre_native), nrow(df_origin))
  df_origin$Sigma <- rep(list(Sigma_native), nrow(df_origin))
  
  df_destination <- data.frame(
    Species = rep(taxon, length(introduced_patches$ID)),
    Biome = introduced_patches$BIOME,
    PatchID = introduced_patches$ID,
    Count = introduced_patches$obs_count,
    Status = rep("Receiver", length(introduced_patches$ID)),
    tasmean = mean(intro_tasmean, na.rm = TRUE),
    tasmean_sd = sd(intro_tasmean, na.rm = TRUE),
    prec = mean(intro_psum, na.rm = TRUE),
    prec_sd = sd(intro_psum, na.rm = TRUE),
    ghm = mean(intro_ghm, na.rm = TRUE),
    npp = mean(intro_npp, na.rm = TRUE)
  )
  df_destination$centre <- rep(list(centre_introduced), nrow(df_destination))
  df_destination$Sigma <- rep(list(Sigma_introduced), nrow(df_destination))
  
  df <- rbind(df_origin, df_destination)
  orig_biomes <- df_origin %>%
    dplyr::group_by(Biome) %>%
    dplyr::summarise(MaxCount = max(Count)) %>%
    dplyr::filter(MaxCount == max(MaxCount)) %>%
    dplyr::pull(Biome)
  
  df$MainOrigin <- ifelse(length(orig_biomes) == 1, orig_biomes, NA)
  
  return(df)
}


if (nativeness_source == "sinas") {
  specs <- unique(sinas_data$taxon)
} else if (nativeness_source == "glonaf") {
  specs <- unique(glonaf_data$taxon)
} else {
  specs <- tax_avail
}


# just to check how many names can be matched
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


if (file.exists(f_out) & !recompute) {
  df_species_patches <- readRDS(f_out)
} else {
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  if (Sys.info()["sysname"] == "Windows") {
    future::plan(future::multisession, workers = parallel::detectCores() - 1)
    lsd <- TRUE # Load spatial data which is not shared between sessions
    stop("Not implemented error: Windows is garbage and this won't run.")
  } else {
    future::plan(future::multicore)
    lsd <- FALSE
  }
  Sys.time()
  progressr::with_progress({
    p <- progressr::progressor(along = specs)
    df_species_patches <- dplyr::bind_rows(
      future.apply::future_lapply(
        X = specs,
        FUN = function(spec) {
          tryCatch(
            {
              res <- assign_patches(
                spec,
                nativeness_source = nativeness_source,
                load_spatial_data = lsd
                )
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
  
  saveRDS(df_species_patches, file = f_out)
}


#>=============================================================================<
#> Add GIST status for the most invasive species
#<=============================================================================>

gisd_codes <- c("MV", "MR", "MO", "MN", "MC")
gisd_levels <- c("Massive", "Major", "Moderate", "Minor", "Minimal concern")
gisd <- read.csv(f_gisd) %>%
  dplyr::filter(EICAT.Category %in% gisd_codes) %>%
  dplyr::mutate(
    Countries.of.impact = dplyr::recode(
      Countries.of.impact, "Ethiopa" = "Ethiopia"
    )
  )

worst <- read.csv(f_100_worst)$Species

gadm <- terra::vect(f_gadm)

check_gist <- function(taxon, biome_patch_id) {
  if (taxon %in% gisd$Species) {
    countries <- gadm[
      which(
        terra::is.related(
          gadm,
          biomes[which(biomes$ID == biome_patch_id),],
          relation = "intersects"
        )
      ),
    ]$COUNTRY
    
    categories <- gisd[
      gisd$Countries.of.impact %in% countries & gisd$Species == taxon,
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

f_out_gist <- sub("df_species_patches.rds", "df_species_patches_gist.rds", f_out)
if (!file.exists(f_out_gist)) {
  df_species_patches %>%
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
      gist_100_worst = factor(Species %in% worst)
    ) %>%
    saveRDS(
      file = f_out_gist
    )
}

df_species_patches <- readRDS(f_out_gist)


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
    total_area_km2 = total_area / 1e6,
    lowland_area_km2 = lowland_area / 1e6
  )

merged$BiomeID <- as.numeric(merged$Biome)
merged$Biome <- factor(biome_names[merged$BiomeID], levels = biome_names)
merged$clusterID <- factor(merged$clusterID)

# Plot
ggdf <- merged |>
  dplyr::filter(Status == "Receiver") |>
  dplyr::add_count(gist_impact_label, name = "N") |>
  dplyr::mutate(
    gist_100_worst = factor(gist_100_worst),
    gist_label_N = paste0(gist_impact_label, "\nN = ", N),
    gist_label_N = factor(
      gist_label_N,
      levels = unique(gist_label_N)[
        match(
          substr(c(gisd_levels, "Not evaluated"), 1, 5),
          substr(unique(gist_label_N), 1, 5)
          )
      ]
    )
    )|>
  dplyr::filter(BiomeID <= 13 & BiomeID != 10)

gg_gist <- ggplot2::ggplot(
  data = ggdf,
  ggplot2::aes(x = gist_label_N, y = lowland_area_km2, fill = gist_impact_label)
) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab("GIST impact class") +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(
    values = c(
      "firebrick3", "coral3", "lightsalmon3", "darkgoldenrod3", "slateblue3"
      )
    ) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_y_log10(
    name = expression("Area in km"^2~"(log scaled)"),
    labels = scales::label_number()
    )

ggplot2::ggsave(
  filename = file.path(dir_fig, "BoxplotGISDclasses.svg"),
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
  filename = file.path(dir_fig, "BoxplotGISDclassesBa.svg"),
  plot = gg_gist_Ba,
  width = 7, height = 5
)

gg_gist_worst <- ggplot2::ggplot(
  data = ggdf %>% dplyr::filter(Status == "Receiver"),
  ggplot2::aes(x = gist_100_worst, y = speciesRichnessBa, fill = gist_100_worst)
) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab("GIST assessment") +
  ggplot2::ylab("Estimated species richness") +
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(
    labels = c(
      "TRUE" = "GIST 100 worst",
      "FALSE" = "Other"
    )
  ) +
  ggplot2::scale_y_continuous(transform = "log10") +
  ggplot2::scale_fill_manual(
    values = c("firebrick3", "slateblue3")
  ) +
  ggplot2::theme(legend.position = "none")

# Direct comparison of group mean for 100 worst and other invasives:
# Effect size is extremely tiny
t.test(speciesRichnessBa ~ gist_100_worst, data = ggdf)
effectsize::cohens_d(
  speciesRichnessBa ~ gist_100_worst,
  data = ggdf,
  alternative = "greater"
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "BoxplotGISD100worst.svg"),
  plot = gg_gist_worst,
  width = 5, height = 5
)


# coin::wilcox_test(
#   speciesRichnessBa ~ gist_100_worst | Biome,
#   data = ggdf,
#   distribution = "approximate"
#   )
# 
# mod <- lme4::lmer(
#   log(speciesRichnessBa) ~ gist_100_worst + (1 | Biome),
#   data = ggdf
#   )
# summary(mod)
# sim_res <- DHARMa::simulateResiduals(mod)
# plot(sim_res)
# 
# DHARMa::testUniformity(sim_res)
# DHARMa::testDispersion(sim_res)
# DHARMa::testOutliers(sim_res)

#confint(mod, parm = "gist_100_worstTRUE", method = "boot")
#                            2.5 %    97.5 %
#   gist_100_worstTRUE -0.4279766 -0.3223543


#>=============================================================================<
#> Plot maps and boxplots
#<=============================================================================>

# Plot maps of donor and receiver regions of introduced species
map_stats <- merged %>%
  dplyr::group_by(ID, Status) %>%
  dplyr::summarise(
    clusterID = dplyr::first(na.omit(clusterID)),
    BiomeID = dplyr::first(na.omit(BiomeID)),
    Biome = dplyr::first(na.omit(Biome)),
    total_area = dplyr::first(na.omit(total_area)),
    lowland_area = dplyr::first(na.omit(lowland_area)),
    species_richness = dplyr::first(na.omit(speciesRichnessBa)),
    species_count = dplyr::n_distinct(Species, na.rm = TRUE),
    .groups = "drop"
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


remerge <- biomes[rep(1:nrow(biomes), each = 2), ]
values(remerge) <- biomes %>%
  terra::as.data.frame() %>%
  dplyr::slice(rep(1:n(), each = 2)) %>%
  dplyr::mutate(Status = rep(c("Donor", "Receiver"), times = nrow(biomes))) %>%
  dplyr::left_join(
    map_stats[, -which(names(map_stats) %in% names(biomes)[-1])],
    by = c("ID", "Status")
  ) %>%
  dplyr::mutate(
    BiomeID = BIOME,
    Biome = factor(biome_names[BiomeID], levels = biome_names),
    species_count = replace_na(species_count, 0)
    )

polygons_sf <- sf::st_as_sf(remerge)
gg_dr <- ggplot2::ggplot(polygons_sf[which(!is.na(polygons_sf$Status)),]) +
  ggplot2::geom_sf(data = biomes, colour = "black", fill = NA) +
  ggplot2::geom_sf(ggplot2::aes(fill = species_count)) +
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
  filename = file.path(dir_fig, "SInAS_donors_and_receivers.png"),
  plot = gg_dr,
  width = 13, height = 3
  )


#-------------------------------------------------------------------------------
# Plot boxplots
boxplot_data <- merged %>%
  dplyr::select(
    Species, Status, lowland_area_km2, speciesRichnessBa, focalECA
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
    max_focalECA_m2 = if(
      all(is.na(focalECA))
    ) NA else max(
      focalECA, na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    max_lowland_1e3km2 = max_lowland_km2 / 1e3,
    max_focalECA = max_focalECA_m2 / 1e9
    ) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(
    cols = c(max_lowland_1e3km2, max_richness, max_focalECA),
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
## Calculate statistics with log-scale aware positioning
metric_max_values <- boxplot_data %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    max_value = max(value, na.rm = TRUE),
    log_max_value = log10(max(value, na.rm = TRUE)),
    place_low = median(value, na.rm = TRUE) < (max(value, na.rm = TRUE) / 2)
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
    metric, group1, group2, statistic, p, effsize, label, y.position, .y.,
    place_low
  ) %>%
  tidyr::pivot_longer(
    cols = c(group1, group2),
    names_to = "Group",
    values_to = "Status"
  ) %>%
  dplyr::mutate(
    value = y.position,
    place_low = place_low
    )

## Simple boxplots by Status (facetted by metric)
gg_area <- ggplot2::ggplot(boxplot_data, aes(x = Status, y = value, fill = Status)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    ~ metric, scales = "free_y",
    labeller = as_labeller(
      c(
        "max_lowland_1e3km2" = "Maximum Lowland Area",
        "max_richness" = "Maximum Species Richness",
        "max_focalECA" = "Maximum Connected Area"
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
    ggplot2::aes(
      x = 1.5, label = label, vjust = ifelse(place_low, 4.5, 2)
      )
    ) +
  ggplot2::scale_y_continuous(
    breaks = c(1e-2, 1, 1e2, 1e4),
    labels = scales::label_scientific(),
    sec.axis = sec_axis(
      ~ ., name = "Estimated species count", breaks = NULL, labels = NULL
      )
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "MaxAreaRichnessBoxplot.svg"),
  plot = gg_area,
  width = 10, height = 4
)


#>=============================================================================<
#> Add sPlot data
#<=============================================================================>

# sPlot Open¨
if (FALSE) {
e <- new.env()
load(f_sPlot, envir = e)
sPlot <- as.list(e)

sPlotData <- sPlot$vegetation %>%
  dplyr::rename(Species = Resolved_name) %>%
  dplyr::mutate(
    Invader = Species %in% unique(sinas_data$taxon) # TODO: CHeck for each potential invader its actual status according to SiNAS
  ) %>%
  dplyr::left_join(
    y = sPlot$header,
    by = "PlotObservationID"
  )

sPlotVect <- sPlotData %>%
  dplyr::group_by(
    PlotObservationID, Latitude, Longitude, Releve_area, Cover_bare_soil,
    Elevation, Naturalness, Invader
    ) %>%
  dplyr::summarise(
    Relative_cover = sum(Relative_cover, na.rm = TRUE),
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Status = ifelse(Invader, "invasive", "native")) %>%
  tidyr::pivot_wider(
    id_cols = c(
      PlotObservationID, Latitude, Longitude, Releve_area, Cover_bare_soil,
      Elevation, Naturalness
      ),
    names_from = Status,
    values_from = c(Relative_cover, N)
    ) %>%
  tidyr::replace_na(
    replace = list(
      "Relative_cover_native" = 0, "Relative_cover_invasive" = 0,
      "N_native" = 0, "N_invasive" = 0
      )
    ) %>%
  dplyr::mutate(N_species = N_native + N_invasive) %>%
  dplyr::filter(
    !is.na(Releve_area)
  ) %>%
  dplyr::mutate(
    BiomePatchID = terra::extract(biomes, .) %>%
      dplyr::pull(ID)
  ) %>%
  terra::vect(geom = c("Longitude", "Latitude"), crs = "epsg:4326")

df_sPlotBiome <- sPlotVect %>%
  dplyr::filter(!is.na(BiomePatchID), Naturalness == "Natural") %>%
  dplyr::group_by(BiomePatchID) %>%
  dplyr::summarise(
    Native_cover = stats::weighted.mean(
      Relative_cover_native,
      w = Releve_area,
      na.rm = TRUE
      ),
    Invasive_cover = stats::weighted.mean(
      Relative_cover_invasive,
      w = Releve_area,
      na.rm = TRUE
    ),
    Bare_soil_cover = stats::weighted.mean(
      Cover_bare_soil,
      w = Releve_area,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    biomes %>%
      as.data.frame() %>%
      dplyr::mutate(BiomePatchID = ID),
    by = "BiomePatchID"
  ) %>%
  dplyr::mutate(
    Biome = factor(BIOME, levels = 1:length(biome_names), labels = biome_names)
  ) %>%
  dplyr::filter(
    as.numeric(Biome) < 14,
    speciesRichnessBa > 30
    ) %>%
  dplyr::mutate(
    log_speciesRichnessBa = log(speciesRichnessBa),
    log_sampling_effort = log(sampling_effort),
    log_Connectedness = log(connectedness),
    log_focalECA = log(focalECA)
  ) %>%
  as.data.frame()

# Fit a GLMM to estimate the impact of our variables
## Assumptions violated:
# n <- nrow(df_sPlotBiome)
# df_sPlotBiome$Invasive_cov <- (df_sPlotBiome$Invasive_cover * (n - 1) + 0.5) / n
# 
# mod <- glmmTMB::glmmTMB(
#   Invasive_cov ~ clim_integral + log(speciesRichnessBa) +
#     (1 | Biome),
#   family = beta_family(),
#   data = df_sPlotBiome
#   )
# 
# summary(mod)
# performance::r2_nakagawa(mod)
# 
# sim <- DHARMa::simulateResiduals(mod)
# plot(sim)

# GAM per biome
bt <- table(df_sPlotBiome$Biome)

par(mfrow = c(2, 2))
for (b in names(bt[bt > 10])) {
  sdf <- df_sPlotBiome[which(df_sPlotBiome$Biome == b), ]
  mod <- mgcv::gam(
    Invasive_cover ~
      s(log_speciesRichnessBa, k = 3) +
      s(bottleneck_clim_integral, k = 3) +
      s(log_focalECA, k = 3) +
      s(log_sampling_effort, k = 3),
    family = quasibinomial(link = "logit"),
    data = sdf
  )
  sm <- summary(mod)
  cat(b, sm$dev.expl, "(GAM)\n")
  
  vars <- c(
    "log_speciesRichnessBa", "bottleneck_clim_integral", "log_sampling_effort",
    "log_focalECA"
    )
  
  for (v in vars) {
    x <- sdf[[v]]
    if (all(!is.finite(x))) next
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng))) next
    grid <- data.frame(
      log_speciesRichnessBa = mean(sdf$log_speciesRichnessBa, na.rm = TRUE),
      bottleneck_clim_integral = mean(sdf$bottleneck_clim_integral, na.rm = TRUE),
      log_sampling_effort = mean(sdf$log_sampling_effort, na.rm = TRUE),
      log_focalECA = mean(sdf$log_focalECA, na.rm = TRUE)
    )
    
    xseq <- seq(rng[1], rng[2], length.out = 100)
    grid <- grid[rep(1, 100), ]
    grid[[v]] <- xseq
    
    pred <- predict(mod, newdata = grid, type = "response")
    
    plot(xseq, pred,
         type = "l",
         main = paste(b, v),
         xlab = v,
         ylab = "Invasive_cover")
  }
}
par(mfrow = c(1, 1))

} # End of sPlot exclusion
#>=============================================================================<
#> Statistical analyses
#<=============================================================================>

env_dat <- readRDS(f_env_dat)
env_df <- tibble::tibble(
  ID = sapply(env_dat, function(x) x$ID),
  patchMu = lapply(env_dat, function(x) x$mu),
  patchSigma = lapply(env_dat, function(x) x$Sigma)
)
env_df$ID <- unlist(env_df$ID)

merged_env <- merge(merged, env_df, by = "ID", all.x = TRUE) %>%
  dplyr::rename(
    climateStability = climate_restime_a,
    climateVelocity = climate_velocity_kmpa
  )

# Run estimate_max_bhattacharyya.R to estimate the maximum environmental distance
# for sampling potentially invaded plots
max_lobBhat <- 1.34

# Create sample data frame
future::plan(multisession, workers = parallel::detectCores() - 1)

load_tab <- TRUE
if (load_tab) {
  load("C:/Users/poppman/Desktop/df_invasion.Rsave")
} else {
  df_invasion <- dplyr::bind_rows(
    future.apply::future_lapply(
      X = unique(merged_env$Species),
      FUN = function(species) {
        donor <- merged_env %>%
          dplyr::filter(
            Species == species,
            !is.na(speciesRichnessBa),
            Status == "Donor"
            )
        if (nrow(donor) == 0) return(NULL)
        
        donor_stats <- donor %>%
          dplyr::summarise(
            modalBiome = collapse::fmode(Biome),
            max_total_area_km2 = max(total_area_km2, na.rm = TRUE),
            max_lowland_area_km2 = max(lowland_area_km2, na.rm = TRUE),
            max_focalECA = max(focalECA, na.rm = TRUE),
            max_connectedness = max(connectedness, na.rm = TRUE),
            maxRichness = max(speciesRichnessBa, na.rm = TRUE),
            maxCorrectedRichness = max(correctedSR, na.rm = TRUE),
            max_climateStability = max(climateStability, na.rm = TRUE),
            min_climateVelocity = min(climateVelocity, na.rm = TRUE)
          )
        
        native_environment <- donor %>%
          dplyr::filter(lowland_area_km2 == donor_stats$max_lowland_area_km2) %>%
          dplyr::slice_head(n = 1) %>%
          dplyr::select(centre, Sigma)
        
        if (
          nrow(native_environment) == 0 ||
          is.null(native_environment$centre[[1]]) ||
          is.null(native_environment$Sigma[[1]])
        ) {
          return(NULL)
        }
        
        receiver <- merged_env %>%
          dplyr::filter(
            Species == species,
            !is.na(speciesRichnessBa),
            Status == "Receiver"
            ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            Invaded = 1,
            tmp = list(
                save_bhattacharyya(
                mu1 = patchMu,
                mu2 = native_environment$centre[[1]],
                Sigma1 = patchSigma,
                Sigma2 = native_environment$Sigma[[1]]
                )
            ),
            Bhattacharyya = as.numeric(tmp$value),
            Bhattacharyya_reason = tmp$reason
          )
        
        id_has_species <- merged_env %>%
          dplyr::filter(Species == species) %>%
          dplyr::pull(ID) %>%
          unique()
        
        background <- merged_env %>%
          dplyr::filter(
            !ID %in% id_has_species,
            !is.na(speciesRichnessBa)
            ) %>%
          dplyr::distinct(ID, .keep_all = TRUE) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            Status = "absent",
            Invaded = 0,
            tmp = list(
              save_bhattacharyya(
                mu1 = patchMu,
                mu2 = native_environment$centre[[1]],
                Sigma1 = patchSigma,
                Sigma2 = native_environment$Sigma[[1]]
              )
            ),
            Bhattacharyya = as.numeric(tmp$value),
            Bhattacharyya_reason = tmp$reason
          ) %>%
          dplyr::filter(
            Bhattacharyya_reason == "ok",
            Bhattacharyya < exp(max_lobBhat)
          )
          
        # Generate output data.frame
        out <- dplyr::bind_rows(
          receiver,
          background
        ) %>%
          dplyr::select(
            ID, clusterID,
            Species, Biome, Count,
            Status, Invaded,
            total_area_km2, lowland_area_km2, focalECA,
            speciesRichnessBa, sampling_effort, correctedSR,
            Bhattacharyya, Bhattacharyya_reason,
            climateStability, climateVelocity
            ) %>%
          dplyr::mutate(
            donor_modalBiome = donor_stats$modalBiome,
            donor_max_total_area_km2 = donor_stats$max_total_area_km2,
            donor_max_lowland_area_km2 = donor_stats$max_lowland_area_km2,
            donor_max_focalECA = donor_stats$max_focalECA,
            donor_max_connectedness = donor_stats$max_connectedness,
            donor_maxRichness = donor_stats$maxRichness,
            donor_maxCorrectedRichness = donor_stats$maxCorrectedRichness,
            donor_maxClimateStability = donor_stats$max_climateStability,
            donor_minclimateVelocity = donor_stats$min_climateVelocity
          )
        
        return(out)
      }
    )
  ) %>%
    dplyr::mutate(
      areaRatio = lowland_area_km2 / donor_max_lowland_area_km2,
      connAreaRatio = focalECA / donor_max_focalECA,
      richnessRatio = speciesRichnessBa / donor_maxRichness,
      correctedRichnessRatio = correctedSR / donor_maxCorrectedRichness,
      relClimStability = climateStability / donor_maxClimateStability,
      relClimVelocity = climateVelocity / donor_minclimateVelocity,
      logAreaRatio = log(areaRatio),
      logConnAreaRatio = log(connAreaRatio),
      logRichnessRatio = log(richnessRatio),
      logCorrectedRichnessRatio = log(correctedRichnessRatio),
      logBhattacharyya = log(Bhattacharyya),
      logSpeciesRichnessBa = log(speciesRichnessBa),
      logRelClimStability = log(relClimStability),
      logRelClimVelocity = log(relClimVelocity),
      logSamplingEffort = log(sampling_effort)
    ) %>%
    as.data.frame()
  
  save(df_invasion, file = "C:/Users/poppman/Desktop/df_invasion.Rsave")
}

## Debugging/analytics
reason_summary <- table(df_invasion$Bhattacharyya_reason)
cat(
  "Bhattacharyya distances computed for",
  sum(df_invasion$Bhattacharyya_reason == "ok"), "of", nrow(df_invasion),
  "observations (i.e.,",
  round(100 * mean(df_invasion$Bhattacharyya_reason == "ok"), 1), "%).\n",
  "Failure breakdown:\n"
)
print(reason_summary)

par(mfrow = c(2, 2))
hist(df_invasion$logBhattacharyya, main = "Bhattacharyya (log transformed)")
hist(df_invasion$logAreaRatio, main = "Area ratio (log transformed)")
hist(df_invasion$logRichnessRatio, main = "Spechies richness ratio (log transformed)")
hist(df_invasion$logCorrectedRichnessRatio, main = "Spechies richness ratio (corrected, log transformed)")
hist(df_invasion$logSamplingEffort, main = "Sampling effort (log transformed)")
hist(df_invasion$logRelClimStability, main = "Relative climate stability (log transformed)")
hist(df_invasion$logRelClimVelocity, main = "Relative climate velocity (log transformed)")
par(mfrow = c(1, 1))

# Filter data
potential_predictors <- c(
  "logConnAreaRatio", "logSpeciesRichnessBa", "logCorrectedRichnessRatio",
  "climateStability"
)
df_invasion_filtered <- df_invasion[ # Filter out incomplete cases
  stats::complete.cases(
    df_invasion[, potential_predictors]
    ),
  ] %>%
  dplyr::filter( # Filter out non-finite values
    dplyr::if_all(
      dplyr::all_of(potential_predictors),
      ~ is.finite(.)
    )
  ) %>%
  dplyr::mutate(
    Specieslvl = stringr::word(Species, 1, 2)
  ) %>%
  dplyr::filter(
    as.numeric(Biome) <= 12
  ) %>%
  dplyr::distinct( # Drop unclear subspecies
    dplyr::across(-Species),
    .keep_all = TRUE
  ) %>%
  dplyr::group_by(Species, Invaded) %>% # Limit to 1000 pseudo-absences per species
  dplyr::slice_min(order_by = Bhattacharyya, n = 1000) %>%
  dplyr::group_by(Species) %>% # Demand a min of 5 pseudo-absences per species
  dplyr::filter(sum(Invaded == 0, na.rm = TRUE) >= 5) %>%
  dplyr::ungroup()
  

num_samples <- df_invasion_filtered %>%
  dplyr::group_by(Species, Invaded) %>%
  dplyr::summarise(
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Invaded) %>%
  dplyr::summarise(
    min = min(N, na.rm = TRUE),
    mean = mean(N, na.rm = TRUE),
    max = max(N, na.rm = TRUE),
    median = median(N, na.rm = TRUE),
    sum = sum(N, na.rm = TRUE)
  )

mean_samples <- num_samples %>%
  dplyr::filter(Invaded == 1) %>%
  dplyr::pull(mean) %>%
  as.integer()

# df_mod <- df_invasion_filtered %>%
#   dplyr::group_by(Species) %>%
#   dplyr::group_modify(
#     ~ {
#       invaded <- dplyr::filter(.x, Invaded == 1)
#       background <- dplyr::filter(.x, Invaded == 0)
#       background <- background %>%
#         dplyr::slice_sample(n = min(nrow(background), mean_samples))
#       dplyr::bind_rows(invaded, background)
#     }
#   ) %>%
#   dplyr::ungroup()

# if (length(which(df_mod$Invaded == 0)) != length(which(df_mod$Invaded == 1))) {
#   warning("Unbalanced samples. Use weights in models.")
# }


# Fit models and compute explained deviance by predictor
## Full model
w <- num_samples %>%
  dplyr::select(Invaded, sum) %>%
  tidyr::pivot_wider(names_from = Invaded, values_from = sum) %>%
  dplyr::transmute(w = `0` / `1`) %>%
  dplyr::pull(w)

df_mod <- df_invasion_filtered %>%
  dplyr::mutate(weight = ifelse(Invaded == 1, w, 1))

predictors <- c(
  "logRichnessRatio",
  "logConnAreaRatio",
  "logRelClimVelocity",
  "logSamplingEffort"
)

frml_full <- make_formula(predictors, biome = TRUE)
mod_gam <- mgcv::gam(
  frml_full,
  data = df_mod,
  method = "REML",
  family = stats::binomial("logit"),
  weights = weight
)

summary(mod_gam)
adjD2_full <- ecospat::ecospat.adj.D2.glm(mod_gam)

mgcv::gam.check(mod_gam)
mgcv::concurvity(mod_gam)

# mod_interact <- mgcv::gam(
#   Invaded ~ s(logRichnessRatio, k = 3) + 
#     te(climateStability, logConnAreaRatio, k = 3) +
#     s(logSamplingEffort, k = 3) + 
#     Biome,
#   data = df_mod,
#   method = "REML",
#   family = stats::binomial("logit"),
#   weights = weight
# )
# summary(mod_interact)
# ecospat::ecospat.adj.D2.glm(mod_interact)

df_pred_d2 <- data.frame()
pb <- progress::progress_bar$new(total = length(predictors) + 1)
for (i in 1:(length(predictors) + 1)) {
  if (i > length(predictors)) {
      p <- "Biome"
      frml_p <- as.formula("Invaded ~ Biome")
    } else {
      p <- predictors[i]
      frml_p <- make_formula(p, biome = FALSE)
    }
  
  df_rand <- df_mod %>%
    dplyr::mutate(
      "{p}" := sample(.data[[p]])
    )
  
  mod_r <- mgcv::gam(
    frml_full,
    data = df_rand,
    method = "REML",
    family = stats::binomial("logit"),
    weights = weight
    )
  adjD2_other <- ecospat::ecospat.adj.D2.glm(mod_r)
  
  mod_p <- mgcv::gam(
    frml_p,
    data = df_rand,
    family = stats::binomial("logit"),
    weights = weight
    )
  adjD2_single <- ecospat::ecospat.adj.D2.glm(mod_p)
  
  df_pred_d2 <- rbind(
    df_pred_d2,
    data.frame(
      Predictor = p,
      Delta_adjD2 = adjD2_full - adjD2_other,
      single_adjD2 = adjD2_single
    )
  )
  pb$tick()
}
rm(pb)

df_pred_d2 <- rbind(
  df_pred_d2 %>%
    dplyr::mutate(
      Predictor = sapply(
        Predictor,
        function(x) trimws(sub("log", "", gsub("([A-Z])", " \\1", x)))
      )
    ),
  data.frame(
    Predictor = c("Shared", "Unexplained"),
    Delta_adjD2 = c(adjD2_full - sum(df_pred_d2$Delta_adjD2), 1 - adjD2_full),
    single_adjD2 = c(NA, NA)
    )
)

gg_pie <- ggplot2::ggplot(
  data = df_pred_d2,
  ggplot2::aes(x = "", y = Delta_adjD2, fill = Predictor)
  ) +
  ggplot2::geom_bar(stat = "identity", width = 1) +
  ggplot2::coord_polar("y", start = 0) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank()
    ) +
  ggplot2::ggtitle("Explained Deviance") +
  ggplot2::scale_fill_manual(
    values = c(
      rgb(0, 102/255, 102/255),
      "red",
      "green",
      rgb(70/255, 100/255, 170/255),
      "yellow",
      "orange",
      #"grey25",
      "grey75"
    )
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "DevExplPie.svg"),
  plot = gg_pie,
  width = 5, height = 5
  )

mods_gam <- list()
d2s_gam <- c()
N <- 50
pb <- progress::progress_bar$new(total = N)
for (i in 1:N) {
  set.seed(i)
  df_ss <- df_mod %>%
    dplyr::group_by(Specieslvl, Status) %>%
    dplyr::sample_n(size = 1)
  
  mods_gam[[i]] <- mgcv::gam(
    frml_full,
    data = df_ss,
    family = stats::binomial("logit")
  )
  d2s_gam <- c(
    d2s_gam,
    ecospat::ecospat.adj.D2.glm(mods_gam[[i]])
  )
  pb$tick()
}
rm(pb)

hist(d2s_gam)


# Fit models to subsets of the data and plot response shapes
means <- colMeans(df_mod[, predictors], na.rm = TRUE)

plot_dfs <- lapply(predictors, function(v) {
  x_seq <- seq(
    min(df_mod[[v]], na.rm = TRUE),
    max(df_mod[[v]], na.rm = TRUE),
    length.out = 500
  )
  
  newdata <- as.data.frame(as.list(means))[rep(1, 500), ]
  newdata[[v]] <- x_seq
  newdata$Biome <- df_mod$Biome[1]
  
  # predictions: matrix (rows = x, cols = models)
  pred_mat <- sapply(
    X = mods_gam,
    FUN = function(m) {predict(m, newdata = newdata, type = "response")}
    )
  
  # summary stats
  data.frame(
    x = x_seq,
    mean = rowMeans(pred_mat, na.rm = TRUE),
    lo = apply(pred_mat, 1, quantile, probs = 0.1, na.rm = TRUE),
    hi = apply(pred_mat, 1, quantile, probs = 0.9, na.rm = TRUE),
    variable = v
  )
})

df_plot <- do.call(rbind, plot_dfs) %>%
  dplyr::mutate(
    Predictor = sapply(
      variable,
      function(x) trimws(sub("log", "", gsub("([A-Z])", " \\1", x)))
      )
  )

gg_resp <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = mean)) +
  ggplot2::geom_line(colour = rgb(0, 102/255, 102/255)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lo, ymax = hi),
    alpha = 0.2, fill = rgb(0, 102/255, 102/255)
    ) +
  ggplot2::labs(
    x = "Predictor value (log scaled)",
    y = "Predicted invasion probability",
  ) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(. ~ Predictor, scales = "free_x")

ggplot2::ggsave(
  filename = file.path(dir_fig, "ResponseShapes.svg"),
  plot = gg_resp,
  width = 7, height = 5
)



df_donor <- merged %>%
  dplyr::mutate(
    PatchID = ID,
    Area = lowland_area,
    OlsonArea = total_area,
    SpeciesRichness = replace_na(round(speciesRichnessBa, 0), 0)
    ) %>%
  dplyr::filter(
    Area > 0,
    Status == "Donor",
    SpeciesRichness > 50
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::filter(
    Count == max(Count)
    ) %>%
  dplyr::group_by(
    Biome, PatchID, Area, OlsonArea, dECA, focalECA, connectedness,
    SpeciesRichness, Status
    ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(
      Biome, PatchID, Area, OlsonArea, dECA, focalECA, connectedness,
      SpeciesRichness
      ),
    names_from = Status,
    values_from = N
  ) %>%
  tidyr::replace_na(replace = list(Donor = 0, Receiver = 0)) %>%
  dplyr::rename(
    DonorOf = Donor
  ) %>%
  dplyr::mutate(
    PropInvasive = DonorOf / SpeciesRichness,
    failures = pmax(SpeciesRichness - DonorOf, 0),
    successRate = DonorOf / SpeciesRichness,
    successProb = successRate / (1 - successRate)
  )

plot(log(successProb) ~ log(SpeciesRichness), data = df_donor)

## GLM that includes biome
modInvasiveBreeder_glm <- stats::glm(
  cbind(DonorOf, failures) ~ log(focalECA) + log(SpeciesRichness) + factor(Biome),
  data = df_donor,
  family = stats::quasibinomial()
)

summary(modInvasiveBreeder_glm)
coef_tab_glm <- summary(modInvasiveBreeder_glm)$coefficients

## GLMM that controls for biome
cor(df_donor$Area, df_donor$SpeciesRichness, use = "complete.obs")
modInvasiveBreeder_glmm <- lme4::glmer(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + (1 | Biome),
  data = df_donor,
  family = binomial
)

summary(modInvasiveBreeder_glmm)
performance::r2(modInvasiveBreeder_glmm)
anova(modInvasiveBreeder_glmm)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  sum(rp^2) / rdf
}

overdisp_fun(modInvasiveBreeder_glmm)

## Betabinomial GLMM
modInvasiveBreeder_bb <- glmmTMB::glmmTMB(
  cbind(DonorOf, failures) ~ log(Area) + log(dECA) + log(SpeciesRichness) + (1 | Biome),
  data = df_donor,
  family = glmmTMB::betabinomial()
)

summary(modInvasiveBreeder_bb)
performance::r2(modInvasiveBreeder_bb)
MuMIn::r.squaredGLMM(modInvasiveBreeder_bb)
res <- DHARMa::simulateResiduals(modInvasiveBreeder_bb)
DHARMa::testDispersion(res)
DHARMa::plotResiduals(res)
DHARMa::testZeroInflation(res)

df_donor_copy <- df_donor
r2 <- c()
for (i in 1:100) {
  df_donor_copy$DonorOf <- sample(df_donor$DonorOf)
  df_donor_copy$failures <- pmax(df_donor_copy$SpeciesRichness - df_donor_copy$DonorOf, 0)
  mod <- glmmTMB::glmmTMB(
    cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + (1 | Biome),
    data = df_donor_copy,
    family = glmmTMB::betabinomial()
  )
  r2 <- c(r2, MuMIn::r.squaredGLMM(mod)[1])
}
max(r2)
hist(r2)
abline(v = MuMIn::r.squaredGLMM(modInvasiveBreeder_bb)[1])


modInvasiveBreeder_gam <- mgcv::gam(
  cbind(DonorOf, failures) ~ s(log(SpeciesRichness), k = 3) + Biome,
  data = df_donor,
  method = "REML",
  family = stats::binomial("logit")
  )

plot(modInvasiveBreeder_gam, pages = 1, shade = TRUE)
summary(modInvasiveBreeder_gam)
ecospat::ecospat.adj.D2.glm(modInvasiveBreeder_gam)


# sequence across observed range
newdata <- expand.grid(
  SpeciesRichness = seq(min(df_donor$SpeciesRichness), max(df_donor$SpeciesRichness), length.out = 100),
  Biome = levels(df_donor$Biome)
)
newdata$log_SpeciesRichness <- log(newdata$SpeciesRichness)

pred <- predict(modInvasiveBreeder_gam, newdata, type = "response", se.fit = TRUE)
newdata$fit <- pred$fit
newdata$lower <- pred$fit - 1.96 * pred$se.fit
newdata$upper <- pred$fit + 1.96 * pred$se.fit

ggplot(newdata, aes(x = SpeciesRichness, y = fit, color = Biome, fill = Biome)) +
  geom_line(colour = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  labs(x = "Species Richness", y = "Predicted probability") +
  theme_minimal()





library("mgcViz")
check_model <- function(mod) {
  # simulate residuals for DHARMa
  sim <- DHARMa::simulateResiduals(mod, plot = FALSE)
  
  # run tests
  disp_test <- DHARMa::testDispersion(sim)
  outlier_test <- DHARMa::testOutliers(sim)
  uni_test <- DHARMa::testUniformity(sim)
  
  # count minor/major issues
  minor <- sum(c(disp_test$p.value, uni_test$p.value) < 0.05)
  major <- sum(c(outlier_test$p.value) < 0.05)
  
  if (minor == 0 && major == 0) {
    return("OK")
  } else {
    paste0(
      if (major > 0) paste0(major, " major issue", if (major > 1) "s" else "") else "",
      if (major > 0 & minor > 0) ", " else "",
      if (minor > 0) paste0(minor, " minor issue", if (minor > 1) "s" else "") else ""
    )
  }
}

pdf(file.path(dir_fig, "biome_response_curves.pdf"), width = 10, height = 5)
par(
  mfrow = c(1, 2),
  mar = c(4, 4, 3, 1),   # inner margins
  oma = c(0, 0, 4, 0)    # outer margins (TOP increased)
)
for (biome in sort(unique(df_donor$Biome))) {
  
  df_sset <- df_donor[df_donor$Biome == biome, ]
  n_patches <- dplyr::n_distinct(df_sset$PatchID)
  
  if (n_patches < 15) {
    next
  }
  
  # ----- GLM -----
  mod_glm <- stats::glm(
    cbind(DonorOf, failures) ~ log(SpeciesRichness),
    data = df_sset,
    family = stats::quasibinomial()
  )
  
  dev_glm <- 1 - mod_glm$deviance / mod_glm$null.deviance
  
  # ----- GAM -----
  mod_gam <- mgcv::gam(
    cbind(DonorOf, failures) ~ s(log(SpeciesRichness), k = 3),
    data = df_sset,
    method = "REML",
    family = stats::binomial("logit")
  )
  
  dev_gam <- ecospat::ecospat.adj.D2.glm(mod_gam)
  
  # ----- Prediction grid -----
  x_seq <- seq(
    min(df_sset$SpeciesRichness, na.rm = TRUE),
    max(df_sset$SpeciesRichness, na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- data.frame(
    SpeciesRichness = x_seq
  )
  
  # GLM prediction
  pred_glm <- predict(
    mod_glm,
    newdata = newdata,
    type = "response"
  )
  
  # GAM prediction
  pred_gam <- predict(
    mod_gam,
    newdata = newdata,
    type = "response"
  )
  
  # Check models
  mod_glm_check <- stats::glm(
    cbind(DonorOf, failures) ~ log(SpeciesRichness),
    data = df_sset,
    family = stats::binomial()
  )
  
  diag_glm <- check_model(mod_glm_check)
  diag_gam <- check_model(mod_gam)
  
  # ----- Plot -----
  ids <- unique(df_sset$PatchID)
  
  # plot base map
  plot(biomes, border = "grey75")
  
  # highlight subset polygons
  subset_polys <- biomes[biomes$ID %in% ids, ]
  plot(subset_polys, col = "red", add = TRUE)
  
  # compute centroids
  cent <- terra::centroids(subset_polys)
  coords <- terra::crds(cent)
  
  # define an offset for arrows (adjust depending on map scale)
  offset <- max(terra::ext(biomes)[2:3] - terra::ext(biomes)[1:2]) * 0.02
  
  # draw arrows from offset positions pointing to centroids
  arrows(
    x0 = coords[,1] + offset,
    y0 = coords[,2] + offset,
    x1 = coords[,1],
    y1 = coords[,2],
    length = 0.1,
    col = "blue",
    lwd = 2.5
  )
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # GLM plot
  plot(
    df_sset$SpeciesRichness,
    df_sset$DonorOf / (df_sset$DonorOf + df_sset$failures),
    pch = 16, cex = log(df_sset$OlsonArea) / log(max(df_sset$OlsonArea)),
    xlab = "Species Richness",
    ylab = "Probability",
    main = paste0("GLM\nD² = ", round(dev_glm, 3), "\nModel checks:", diag_glm),
    log = "x"
  )
  lines(x_seq, pred_glm, lwd = 2)
  
  # GAM plot
  plot(
    df_sset$SpeciesRichness,
    df_sset$DonorOf / (df_sset$DonorOf + df_sset$failures),
    pch = 16, cex = log(df_sset$OlsonArea) / log(max(df_sset$OlsonArea)),
    xlab = "Species Richness",
    ylab = "Probability",
    main = paste0("GAM\nD² = ", round(dev_gam, 3), "\nModel checks:", diag_gam),
    log = "x"
  )
  lines(x_seq, pred_gam, lwd = 2)
  
  # page title
  mtext(
    paste0("Biome: ", biome, "   |   N patches = ", n_patches),
    outer = TRUE,
    line = 1,
    cex = 1.2
  )
}

dev.off()





for (i in 1:14) {
  biome <- biome_names[i]
  df_sset <- df_donor[df_donor$Biome == biome, ]
  n_patches <- dplyr::n_distinct(df_sset$PatchID)
  
  if (n_patches < 15) {
    next
  }
  
  dev_glms <- c()
  dev_gams <- c()
  glms <- list()
  gams <- list()
  for (i in 1:30) {
    # Number of bins along SpeciesRichness
    n_bins <- 10
    max_per_bin <- 3
    
    df_sset_sample <- df_sset %>%
      mutate(bin = cut(SpeciesRichness, breaks = n_bins, include.lowest = TRUE)) %>%
      group_by(bin) %>%
      group_modify(~ slice_sample(.x, n = min(nrow(.x), max_per_bin)), .keep = FALSE) %>%
      ungroup() %>%
      select(-bin)
    
    # Check distribution
    table(cut(df_sset_sample$SpeciesRichness, breaks = n_bins))
    
    # ----- GLM -----
    mod_glm <- stats::glm(
      cbind(DonorOf, failures) ~ log(SpeciesRichness),
      data = df_sset_sample,
      family = stats::quasibinomial()
    )
    
    dev_glm <- 1 - mod_glm$deviance / mod_glm$null.deviance
    dev_glms <- c(dev_glms, dev_glm)
    glms[[i]] <- mod_glm
    
    # ----- GAM -----
    mod_gam <- mgcv::gam(
      cbind(DonorOf, failures) ~ s(log(SpeciesRichness), k = 3),
      data = df_sset_sample,
      method = "REML",
      family = stats::binomial("logit")
    )
    
    dev_gam <- ecospat::ecospat.adj.D2.glm(mod_gam)
    dev_gams <- c(dev_gams, dev_gam)
    gams[[i]] <- mod_gam
  }
  
  hist(dev_glms)
  hist(dev_gams)
  
  dev_glms[dev_glms > 1] <- 0
  dev_gams[dev_gams > 1] <- 0
  
  x_seq <- seq(
    min(df_sset$SpeciesRichness, na.rm = TRUE),
    max(df_sset$SpeciesRichness, na.rm = TRUE),
    length.out = 100
  )
  
  newdata <- data.frame(
    SpeciesRichness = x_seq
  )
  
  # GLM plot
  grDevices::svg(
    filename = file.path(
      "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/model_tests",
      paste0("GLM_", biome, ".svg")
    ), width = 7, height = 5
  )
  plot(
    df_sset$SpeciesRichness,
    df_sset$DonorOf / (df_sset$DonorOf + df_sset$failures),
    pch = 16, cex = log(df_sset$OlsonArea) / log(max(df_sset$OlsonArea)),
    xlab = "Species Richness",
    ylab = "Probability",
    main = paste0(biome, ": GLM\nmean D² = ", round(mean(dev_glms), 2)),
    log = "x"
  )
  
  # GLM prediction
  for (i in 1:length(glms)) {
    pred_glm <- predict(
      glms[[i]],
      newdata = newdata,
      type = "response"
    )
    lines(x_seq, pred_glm, lwd = 2, col = rgb(1, 0, 0, pmax(dev_glms[i], 0)))
  }
  dev.off()
  
  # GAM plot
  grDevices::svg(
    filename = file.path(
      "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/model_tests",
      paste0("GAM_", biome, ".svg")
      ), width = 7, height = 5
  )
  plot(
    df_sset$SpeciesRichness,
    df_sset$DonorOf / (df_sset$DonorOf + df_sset$failures),
    pch = 16, cex = log(df_sset$OlsonArea) / log(max(df_sset$OlsonArea)),
    xlab = "Species Richness",
    ylab = "Probability",
    main = paste0(biome, ": GAM\nD² = ", round(mean(dev_gams), 2)),
    log = "x"
  )
  
  # GAM prediction
  for (i in 1:length(gams)) {
    pred_gam <- predict(
      gams[[i]],
      newdata = newdata,
      type = "response"
    )
    lines(x_seq, pred_gam, lwd = 2, col = rgb(0, 1, 0, pmax(dev_gams[i], 0)))
  }
  dev.off()
}



















df_donor$SpeciesRichness <- sample(df_donor$SpeciesRichness)
df_donor$failures <- pmax(df_donor$SpeciesRichness - df_donor$DonorOf, 0) 



df_within <- merged %>%
  dplyr::mutate(
    PatchID = ID,
    Area = lowland_area,
    SpeciesRichness = round(speciesRichnessBa, 0)
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::filter(
    Count > 10,
    SpeciesRichness > 50
  ) %>%
  dplyr::group_by(BiomeID, Biome, Species, clusterID, Status) %>%
  dplyr::summarise(
    Count = sum(Count),
    clusterRichness = max(speciesRichnessBa, na.rm = TRUE),
    focalECA = max(focalECA, na.rm = TRUE),
    connectedness = max(connectedness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(BiomeID, Biome, Species, Status) %>%
  dplyr::filter(
    Count == max(Count, na.rm = TRUE)
  ) %>%
  dplyr::filter(# Addiditonal filter, since in some cases, multiple clusters have max(Count)
    clusterRichness == max(clusterRichness)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BiomeID, Biome, Species) %>%
  dplyr::filter(!duplicated(clusterID)) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 2) %>%
  dplyr::ungroup()

ggplot2::ggplot(
  data = df_within,
  ggplot2::aes(y = log(clusterRichness), colour = Status)
  ) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(.~Biome) +
  ggplot2::theme_bw()

df_wide_richness <- df_within %>%
  dplyr::select(Biome, Species, Status, clusterRichness) %>%
  tidyr::pivot_wider(
    names_from = Status,
    values_from = clusterRichness
  ) %>%
  dplyr::filter(!is.na(Donor) & !is.na(Receiver)) %>%
  dplyr::mutate(
    diff = Donor - Receiver,
    ratio = (Donor - Receiver) / (Donor + Receiver)
    )

mn <- mean(df_wide_richness$ratio)
md <- median(df_wide_richness$ratio)

hist(
  df_wide_richness$ratio,
  breaks = 30,
  col = rgb(0, 102/255, 102/255, 0.25),
  border = "grey50",
  main = "Richness index (Donor - Receiver) / (Donor + Receiver)"
  )
abline(v = mn, col = "forestgreen", lwd = 2, lty = 2)
abline(v = md, col = "steelblue", lwd = 2, lty = 3)
legend(
  "topright",
  legend = c("Mean", "Median"),
  col = c("forestgreen", "steelblue"),
  lty = c(2, 3), lwd = 2
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
  dplyr::select(Biome, Species, Status, connectedness) %>%
  tidyr::pivot_wider(
    names_from = Status,
    values_from = connectedness
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




















# Boxplots of cluster richness of source/sink patch for each species facetted
# by biome (only entries where source and sink biome are the same type)

## Calculate statistics with log-scale aware positioning
df_within_bp <- df_within %>%
  dplyr::mutate(value = clusterRichness) %>%
  dplyr::filter(BiomeID <= 12) %>%
  base::droplevels()

metric_max_values <- df_within_bp %>%
  dplyr::group_by(Biome) %>%
  dplyr::summarise(
    max_value = max(value, na.rm = TRUE),
    log_max_value = log10(max(value, na.rm = TRUE))
  )

## Wilcoxon test
test_res <- df_within_bp |>
  dplyr::group_by(Biome) |>
  rstatix::wilcox_test(
    value ~ Status,
    paired = TRUE,
    alternative = "greater"
  ) |>
  rstatix::adjust_pvalue(method = "BH")

## Effect size
eff_res <- df_within_bp |>
  dplyr::group_by(Biome) |>
  rstatix::wilcox_effsize(
    value ~ Status,
    paired = TRUE
  )

## Create stats dataframe - using log-space coordinates
stats <- test_res %>%
  dplyr::left_join(eff_res %>% select(Biome, effsize), by = "Biome") %>%
  dplyr::left_join(metric_max_values, by = "Biome") %>%
  dplyr::mutate(
    label = paste0(
      "W = ", statistic,
      "\n", "p = ", signif(p, 3),
      "\n", "effect size = ", round(effsize, 2),
      "\nN = ", n1
    ),
    y.position = max_value,
    .y. = "value"
  ) %>%
  dplyr::select(
    Biome, group1, group2, statistic, p, effsize, label, y.position, .y.
  ) %>%
  tidyr::pivot_longer(
    cols = c(group1, group2),
    names_to = "Group",
    values_to = "Status"
  ) %>%
  dplyr::mutate(value = y.position)

## Simple boxplots by Status (facetted by metric)
gg_cluster_by_biome <- ggplot2::ggplot(
  data = df_within_bp,
  ggplot2::aes(x = Status, y = value, fill = Status)
  ) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    .~ Biome, scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c("Donor" = "#2E86AB", "Receiver" = "#A23B72")
  ) +
  ggplot2::coord_transform(y = "log10") +
  ggplot2::labs(
    title = NULL,
    y = "Maximum species richness of the biome patch cluster", # expression("Area (" * 10^3 * " km"^2 * ")"),
    x = "Biome patch cluster status"
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
    vjust = 3
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(1e-2, 1, 1e2, 1e4),
    labels = scales::label_scientific(),
    sec.axis = sec_axis(
      ~ ., name = "Estimated species count", breaks = NULL, labels = NULL
    )
  )
