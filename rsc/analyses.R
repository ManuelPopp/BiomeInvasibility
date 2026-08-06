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
library("stringr")
library("boot")
library("rsdd")
library("ggplot2")
library("patchwork")
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
#library("coin")
library("lmodel2")
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

nativeness_source <- "sinas" # "powo"
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
dir_imed <- file.path(dir_imeb, nativeness_source)  # Patch-level data
dir_fig <- file.path(dir_main, "fig", nativeness_source)
dir_tab <- "C:/Users/poppman/Dropbox/Apps/Overleaf/BiomeInvasibility/tab"

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

f_powo_places <- file.path(dir_dat, "powo", "tdwg_level3.geojson")
f_powo_taxa <- file.path(dir_dat, "powo", "wcvp_names.csv")
f_powo_taxa_dist <- file.path(dir_dat, "powo", "wcvp_distribution.csv")

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

f_road_dns <- file.path(
  dir_imeb, "road_density.csv"
)

f_gdp <- file.path(
  dir_imeb, "gdp_1990to2020.csv"
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

## Ensure fallback operator exists
`%||%` <- function(x, y) if (is.null(x)) y else x

## Create a model formula from a vector of predictors
make_formula <- function(
    predictors, k = 3, response = "Invaded", biome = TRUE, interact = NA
    ) {
  if (k > 0) {
    smooth_terms <- sprintf("s(%s, k = %s)", predictors, k)
    rhs <- paste(smooth_terms, collapse = " + ")
    if (!is.na(interact)) {
      warning("Interactions not implemented for smooth terms.")
    }
  } else {
    if (is.na(interact)) {
      rhs <- paste(predictors, collapse = " + ")
    } else {
      if (is.integer(interact)) {
        interact <- c(interact)
      }
      if (!is.numeric(interact)) {
        stop("Parameter interact must be a vector of integers.")
      }
      operator <- c(rep("+", length(predictors) - 1), "")
      operator[interact] <- "*"
      rhs <- paste(as.vector(t(cbind(predictors, operator))), collapse = "")
    }
  }
  if (biome) {
    rhs <- paste(rhs, "+ Biome")
  }
  as.formula(paste(response, "~", rhs))
}

## Check dispersion
dispersion <- function(mod) {
  return(sum(residuals(mod, type = "pearson")^2) / mod$df.residual)
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

# Check whether a species is considered native at a given location
sinas_status <- function(taxon, sinas_loc_ids, sinas_data) {
  if (!(taxon %in% sinas_data$taxon)) {
    return(rep(NA_character_, length(sinas_loc_ids)))
  }
  
  ds <- sinas_data[which(sinas_data$taxon == taxon), ]
  native_loc_ids <- ds[ds$establishmentMeans == "native", ]$locationID
  introduced_loc_ids <- ds[ds$establishmentMeans == "introduced", ]$locationID
  
  if (
    length(native_loc_ids) == 0 &&
    length(introduced_loc_ids) == 0
  ) {
    return(rep("unknown", length(sinas_loc_ids)))
  }
  
  is_native <- sinas_loc_ids %in% native_loc_ids
  is_introduced <- sinas_loc_ids %in% introduced_loc_ids
  out <- ifelse(
    is_native & !is_introduced,
    "native",
    ifelse(
      is_introduced & !is_native,
      "introduced",
      "unknown"
    )
  )
  return(out)
}


boot_conf_int <- function(sdf, frml, fam, coef_names) {
  mod <- stats::glm(
    frml,
    family = fam,
    data = sdf
  )
  
  estimates <- stats::coef(mod)
  coef_names <- names(estimates)
  
  boot_coef <- function(dat, indices) {
    d <- dat[indices, ]
    stats::coef(
      stats::glm(
        frml,
        family = fam,
        data = d
      )
    )
  }
  
  boot_res <- boot::boot(
    data = sdf,
    statistic = boot_coef,
    R = 1000
  )
  
  boot_mat <- boot_res$t
  
  sign_stability <- sapply(
    seq_along(estimates),
    function(i) {
      mean(
        sign(boot_mat[, i]) == sign(estimates[i]),
        na.rm = TRUE
      )
    }
  )
  
  data.frame(
    predictor = coef_names,
    estimate = estimates,
    
    ci90_low  = apply(boot_mat, 2, quantile, probs = 0.05,  na.rm = TRUE),
    ci90_high = apply(boot_mat, 2, quantile, probs = 0.95,  na.rm = TRUE),
    
    ci95_low  = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    ci95_high = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    
    ci99_low  = apply(boot_mat, 2, quantile, probs = 0.005, na.rm = TRUE),
    ci99_high = apply(boot_mat, 2, quantile, probs = 0.995, na.rm = TRUE),
    
    sign_stability = round(100 * sign_stability, 1),
    
    row.names = NULL,
    check.names = FALSE
  )
}


#>=============================================================================<
#> Data preparation
#<=============================================================================>

## GBIF
rsdd::dataset("gbif-powo_raw")
taxa <- rsdd::taxa()
tax_avail <- taxa$species

## SINaS
sinas_places <- terra::vect(f_sinas_places)
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

if (!"gift_nat_centroid" %in% names(biomes) | recompute) {
  system(paste("Rscript get_GIFT_SR.R", nativeness_source, biome_def))
  biomes <- terra::vect(f_bfullinfo)
}

#>-----------------------------------------------------------------------------<
#> Add species richness estimates, climate stability, sampling effort, etc
est_div <- read.csv(f_est_div)
spl_eff <- read.csv(f_spl_eft)
hum_mod <- read.csv(f_hum_mod)
road_dns <- read.csv(f_road_dns)
gdp <- read.csv(f_gdp)

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
  terra::merge(road_dns, by = "ID", all.x = TRUE) %>%
  terra::merge(gdp, by = "ID", all.x = TRUE) %>%
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

biomes$logCorrectedSR <- NA
biomes$logCorrectedSR[which(biomes$log_SpeciesRichness > log(100))] <- data.frame(
  corrected = df_fit$log_SpeciesRichness - effort_fit + max_sampling_effect[1],
  estimated = df_fit$log_SpeciesRichness
) %>%
  dplyr::mutate(
    corrected = ifelse(
      # Very high estimates are probably approx. accurate while very high corrected
      # values are probably only accurate, it the estimate itself is also high.
      exp(corrected) > 3 * exp(estimated) | exp(corrected) > 3.2e+4,
      NA,
      corrected
      ),
    estimated = ifelse(
      exp(estimated) > 6e+4, # Amazon rainforest = some 50 k, higher is unrealistic
      NA,
      estimated
    ),
    maximum = pmax(corrected, estimated, na.rm = TRUE),
    average = log(
      rowMeans(
        cbind(exp(corrected), exp(estimated)),
        na.rm = TRUE
      )
    ),
    use_average = ifelse(
      is.na(corrected),
      FALSE,
      exp(corrected) > 2 * exp(estimated) | exp(estimated) > 2 * exp(corrected)
      ),
    use_value = ifelse(use_average, average, maximum)
  ) %>%
  dplyr::pull(use_value)


#>-----------------------------------------------------------------------------<
#> Plot species richness estimates

gg_sr_df <- as.data.frame(biomes) |>
  dplyr::mutate(
    logSpeciesRichnessBa = log(speciesRichnessBa),
    logGIFTcentroid = log(gift_nat_centroid),
    logGIFTcovered = log(gift_nat_covered)
  ) |>
  dplyr::select(
    ID,
    logSpeciesRichnessBa,
    logCorrectedSR,
    logGIFTcentroid,
    logGIFTcovered
    ) |>
  dplyr::select(-ID) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ replace(.x, !is.finite(.x), NA_real_)
    )
  )

var_labels <- c(
  logSpeciesRichnessBa = "Breakaway (log)",
  logCorrectedSR = "Corrected Breakaway (log)",
  logGIFTcentroid = "GIFT centroid-based (log)",
  logGIFTcovered = "GIFT fully covered (log)"
)

panel_fun <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  fit <- lmodel2::lmodel2(
    y ~ x,
    data = data.frame(x = x, y = y),
    nperm = 0
  )
  ma <- fit$regression.results[
    fit$regression.results$Method == "MA",
  ]
  slope <- ma$Slope
  intercept <- ma$Intercept
  
  grid <- seq(min(x), max(x), length.out = 100)
  ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_point(colour = rgb(0, 102/255, 102/255, 0.5), size = 1) +
    ggplot2::geom_abline(
      slope = slope,
      intercept = intercept,
      colour = "royalblue4"
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "grey40"
    )
}

gg_sr_comp <- GGally::ggpairs(
  gg_sr_df,
  upper = list(continuous = "cor"),
  lower = list(continuous = panel_fun),
  labeller = ggplot2::as_labeller(var_labels)
) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = file.path(dir_fig, "SR_estimate_comparison.svg"),
  plot = gg_sr_comp, width = 7, height = 7
  )

cat(
  "\nBreakaway:", length(which(is.finite(gg_sr_df$logSpeciesRichnessBa))),
  "\nBreakaway corrected:", length(which(is.finite(gg_sr_df$logCorrectedSR))),
  "\nGIFT centroid:", length(which(is.finite(gg_sr_df$logGIFTcentroid))),
  "\nGIFT covered:", length(which(is.finite(gg_sr_df$logGIFTcovered)))
)


for (i in 1:2) {
  if (i == 1) {
    vals <- biomes$speciesRichnessBa
    fn <- "EstimatedSpeciesRichness.png"
    mn <- "Estimated richness of vascular plant species"
    scale_max <- max(vals, na.rm = TRUE)
  } else {
    vals <- exp(biomes$logCorrectedSR)
    #vals[which(vals > scale_max)] <- scale_max
    fn <- "CorrectedSpeciesRichness.png"
    mn <- "Corrected richness of vascular plant species"
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
    biomes, col = plot_cols, border = rgb(0.5, 0.5, 0.5, 0.5), lwd = 0.2,
    main = mn
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
  ggplot2::facet_wrap(~ Status, nrow = 2) +
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
  width = 7, height = 6 # 13 and 3 for two-column version
  )


#-------------------------------------------------------------------------------
# Plot boxplots
metrics <- c(
  "max_lowland_1e3km2",
  "max_focalECA",
  #"max_richness",
  #"max_grichness",
  "max_crichness"
)

boxplot_data <- merged %>%
  dplyr::select(
    Species, Status, lowland_area_km2, speciesRichnessBa, focalECA,
    logCorrectedSR, gift_nat_centroid
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
    max_crichness = if(
      all(is.na(exp(logCorrectedSR)))
    ) NA else max(
      logCorrectedSR, na.rm = TRUE
    ),
    max_grichness = if(
      all(is.na(gift_nat_centroid))
    ) NA else max(
      gift_nat_centroid, na.rm = TRUE
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
    cols = c(
      max_lowland_1e3km2,
      max_richness, max_crichness, max_grichness,
      max_focalECA
      ),
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
  dplyr::ungroup() %>%
  dplyr::filter(metric %in% metrics) %>%
  dplyr::mutate(
    metric = factor(metric)
  )

# Calculate statistics
## Calculate statistics with log-scale aware positioning
metric_max_values <- boxplot_data %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    max_value = max(value, na.rm = TRUE),
    log_max_value = log10(max(value, na.rm = TRUE)),
    place_low = median(log10(value), na.rm = TRUE) > (
      max(log10(value), na.rm = TRUE) - (
        max(log10(value), na.rm = TRUE) - min(log10(value), na.rm = TRUE)
        ) / 2
      )
  )

## Wilcoxon test
test_res <- boxplot_data |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_test(
    value ~ Status,
    paired = TRUE,
    alternative = "greater",
    p.adjust.method = "BH"
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
      "\n", "p = ", signif(p.adj, 3),
      "\n", "effect size = ", sprintf("%.2f", effsize)
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

# Results
# Metric  p-value     effectsize Test stat (W) N
# SR      5.24e-118   0.232      13345744      10529
# C. SR   1.34e-228   0.302      15143531      10523
# GIFT SR 2.63e-118   0.220      11311265      10383
# ECA     7.65e-240   0.330      19587786      10529
# Area    1.86e- 16   0.106      15901928      10510

## Simple boxplots by Status (facetted by metric)
gg_area <- ggplot2::ggplot(
  boxplot_data,
  ggplot2::aes(x = Status, y = value, fill = Status)
  ) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    ~ metric, scales = "free_y", nrow = 1,
    labeller = as_labeller(
      c(
        "max_lowland_1e3km2" = "Maximum Lowland Area",
        "max_richness" = "Maximum Species Richness",
        "max_focalECA" = "Maximum Connected Area",
        "max_crichness" = "Maximum Species Richness"
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
    y = expression("Area in " * 10^3 * " km"^2 * " or estimated number or species"),
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
      x = 1.5, label = label, vjust = ifelse(place_low, 4.5, 1.5)
      )
    ) +
  ggplot2::scale_y_continuous(
    breaks = scales::breaks_log(n = 5),
    labels = scales::label_scientific()#,
    # sec.axis = sec_axis(
    #   ~ ., name = "Estimated species count", breaks = NULL, labels = NULL
    #   )
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "MaxAreaRichnessBoxplot.svg"),
  plot = gg_area,
  width = 3.3 * length(metrics), height = 4
)


#>=============================================================================<
#> Analyses based on sPlot data
#<=============================================================================>

# sPlot
f_splotdata <- file.path(dir_imed, "splot.RData")
if (file.exists(f_splotdata)) {
  load(f_splotdata)
} else {
  if (nativeness_source == "sinas") {
    listed_taxa  <- unique(sinas_data$taxon)
  } else if (nativeness_source == "powo") {
    taxa <- data.table::fread(f_powo_taxa, header = TRUE, sep = "|") %>%
      dplyr::select(
        plant_name_id, species, taxon_name, taxon_rank, taxon_status,
        accepted_plant_name_id
      )
    locs <- data.table::fread(f_powo_taxa_dist, header = TRUE, sep = "|") %>%
      dplyr::select(
        plant_name_id, area_code_l3, introduced, extinct, location_doubtful
      )
    powo <- dplyr::left_join(
      locs,
      dplyr::filter(taxa, taxon_status == "Accepted"),
      by = "plant_name_id"
      ) %>%
      dplyr::filter(
        taxon_rank %in% c("Species", "Subspecies")
      ) %>%
      dplyr::mutate(
        locationCode = factor(area_code_l3),
        locationID = as.numeric(locationCode)
      )
    listed_taxa <- unique(powo$taxon_name)
  }
  
  e <- new.env()
  load(f_sPlot, envir = e)
  sPlot <- as.list(e)
  rm(e)
  gc()
  
  sPlot$header <- sPlot$header %>%
    dplyr::select(
      PlotObservationID,
      Longitude, Latitude,
      Releve_area, Cover_scale, Naturalness, Grassland, Shrubland,
      Forest, Wetland,
      Cover_total, Cover_forbs, Cover_bare_soil, Cover_bare_rock, Cover_open_water,
      Altitude
    )
  
  sPlotData <- sPlot$vegetation %>%
    dplyr::rename(Species = Resolved_name) %>%
    dplyr::select(
      PlotObservationID, Species, Taxon_group, Layer, Abundance_scale, Abundance,
      Rel_Abund_Layer, Rel_Abund_Plot
    ) %>%
    dplyr::filter(
      Taxon_group %in% c("Angiosperm", "Gymnosperm", "Vascular plant")
    ) %>%
    dplyr::left_join(
      y = sPlot$header,
      by = "PlotObservationID"
    ) %>%
    dplyr::filter(
      Releve_area >= 1,
      Releve_area <= 10000,
      is.finite(Longitude),
      is.finite(Latitude)
    ) %>%
    dplyr::mutate(
      Taxon_group = factor(Taxon_group),
      Cover_scale = factor(Cover_scale),
      Invader = Species %in% listed_taxa
    )
  
  rm(sPlot)
  gc()
  
  if (nativeness_source == "powo") {
    sPlotSpecies <- sPlotData %>%
      dplyr::mutate(
        Species = stringr::str_replace(
          Species,
          "^([A-Z][a-z-]+) ([a-z-]+) subsp\\. \\2$",
          "\\1 \\2"
        )
      ) %>%
      dplyr::group_by(Species) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
        dplyr::filter(taxa, taxon_status == "Accepted"),
        by = c("Species" = "taxon_name")
      )
    
    accepted_taxa <- taxa %>%
      dplyr::filter(taxon_status == "Accepted") %>%
      dplyr::select(
        plant_name_id,
        taxon_name
      )
    
    synonyms <- sPlotSpecies %>%
      dplyr::filter(is.na(plant_name_id)) %>%
      dplyr::select(Species) %>%
      dplyr::left_join(
        dplyr::filter(taxa, taxon_status == "Synonym") %>%
          dplyr::select(
            taxon_name,
            accepted_plant_name_id
          ),
        by = c("Species" = "taxon_name")
      ) %>%
      dplyr::left_join(
        accepted_taxa,
        by = c("accepted_plant_name_id" = "plant_name_id")
      ) %>%
      dplyr::transmute(
        Species,
        plant_name_id = accepted_plant_name_id,
        taxon_name
      ) %>%
      dplyr::group_by(Species) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    sPlotSpecies <- sPlotSpecies %>%
      dplyr::filter(!is.na(plant_name_id)) %>%
      dplyr::transmute(
        Species,
        plant_name_id,
        taxon_name = Species
      ) %>%
      dplyr::select(Species, plant_name_id, taxon_name) %>%
      dplyr::bind_rows(synonyms) %>%
      dplyr::distinct(Species, .keep_all = TRUE) %>%
      dplyr::rename(
        Species_sPlot = Species,
        Species = taxon_name
      )
    
    sPlotData <- sPlotData %>%
      dplyr::rename(Species_sPlot = Species) %>%
      dplyr::mutate(
        Species_sPlot = stringr::str_replace(
          Species_sPlot,
          "^([A-Z][a-z-]+) ([a-z-]+) subsp\\. \\2$",
          "\\1 \\2"
        )
      ) %>%
      dplyr::left_join(sPlotSpecies, by = "Species_sPlot")
  }
  
  # Get SINaS location ID for each sPlot plot
  # Check if sinas_places contains multiple levels (e.g. continents + countries)
  #overlapping_places <- terra::relate(sinas_places, relation = "overlaps")
  
  plot_locs <- sPlotData %>%
    dplyr::group_by(PlotObservationID) %>%
    dplyr::summarise(
      Longitude = dplyr::first(Longitude),
      Latitude = dplyr::first(Latitude),
      has_invaders = any(Invader)
    ) %>%
    dplyr::filter(has_invaders) %>%
    dplyr::select(PlotObservationID, Longitude, Latitude) %>%
    terra::vect(geom = c("Longitude", "Latitude"), crs = "epsg:4326")
  
  # Rasterise places for faster value extraction
  if (nativeness_source == "sinas") {
    r_loc <- terra::rasterize(
      sinas_places,
      terra::rast(f_envstack, lyrs = 1),
      field = "locationID",
      touches = TRUE
    )
  } else if (nativeness_source == "powo") {
    powo_places <- terra::vect(f_powo_places) %>%
      dplyr::mutate(
        locationCode = factor(LEVEL3_COD),
        locationID = as.numeric(locationCode)
        ) %>%
      dplyr::select(locationCode, locationID)
    
    r_loc <- terra::rasterize(
      powo_places,
      terra::rast(f_envstack, lyrs = 1),
      field = "locationID",
      touches = TRUE
    )
  }
  
  # Extract SINaS location IDs
  loc_ids <- terra::extract(r_loc, plot_locs) %>%
    dplyr::mutate(PlotObservationID = plot_locs$PlotObservationID)
  
  # Lookup table
  plot_to_loc <- setNames(
    loc_ids$locationID,
    loc_ids$PlotObservationID
  )
  
  # Sinas taxa: 41187
  # sPlot taxa: 100222
  # Shared taxa: 12763 before filtering, 11613 after filtering
  sPlotData$Status <- NA_character_
  
  # Create look-up tables for native/introduced location IDs
  if (nativeness_source == "sinas") {
    invader_species <- unique(sPlotData$Species[which(sPlotData$Invader)])
    taxa_small <- sinas_data %>%
      dplyr::filter(taxon %in% invader_species) %>%
      dplyr::select(taxon, locationID, establishmentMeans)
  } else if (nativeness_source == "powo") {
    invader_species <- unique(
      sPlotData$Species[which(!is.na(sPlotData$plant_name_id))]
      )
    taxa_small <- powo %>%
      dplyr::filter(taxon_name %in% invader_species) %>%
      dplyr::mutate(
        taxon = taxon_name,
        establishmentMeans = ifelse(
          introduced == 1, "introduced", "native"
        )
      ) %>%
      dplyr::select(taxon, locationID, establishmentMeans)
  }
  
  native_map <- split(
    taxa_small$locationID[taxa_small$establishmentMeans == "native"],
    taxa_small$taxon[taxa_small$establishmentMeans == "native"]
  )
  
  intro_map <- split(
    taxa_small$locationID[taxa_small$establishmentMeans == "introduced"],
    taxa_small$taxon[taxa_small$establishmentMeans == "introduced"]
  )
  
  native_map <- setNames(
    lapply(invader_species, function(sp) {native_map[[sp]] %||% integer(0)}),
    invader_species
  )
  
  intro_map <- setNames(
    lapply(invader_species, function(sp) {intro_map[[sp]] %||% integer(0)}),
    invader_species
  )
  
  # Get native/introduced status of invader species observations
  sPlotDataSmall <- dplyr::filter(sPlotData, Invader) %>%
    dplyr::select(PlotObservationID, Species) %>%
    dplyr::mutate(
      sinas_loc_id = unname(plot_to_loc[as.character(PlotObservationID)])
    )
  
  is_native <- mapply(
    FUN = function(t, loc) {loc %in% native_map[[t]]},
    t = sPlotDataSmall$Species,
    loc = sPlotDataSmall$sinas_loc_id,
    SIMPLIFY = TRUE
  )
  
  is_intro <- mapply(
    FUN = function(t, loc) {loc %in% intro_map[[t]]},
    t = sPlotDataSmall$Species,
    loc = sPlotDataSmall$sinas_loc_id,
    SIMPLIFY = TRUE
  )
  
  sPlotDataSmall$Status[is_native & !is_intro] <- "native"
  sPlotDataSmall$Status[is_intro & !is_native] <- "introduced"
  sPlotDataSmall$Status[is_native & is_intro] <- "unknown"
  sPlotDataSmall$Status[!is_native & !is_intro] <- "unknown"
  
  sPlotData$Status[sPlotData$Invader] <- sPlotDataSmall$Status
  
  save(sPlotData, file = f_splotdata)
}


f_splotvect <- file.path(dir_imed, "splot.gpkg")
if (file.exists(f_splotvect)) {
  sPlotVect <- terra::vect(f_splotvect)
} else {
  # Allowed abundance scales for harmonisation of relative abundance
  allowed_abundance_scales <- c("Cover", "Biomass")
  
  sPlotVect <- sPlotData %>%
    dplyr::filter(Abundance_scale %in% allowed_abundance_scales) %>%
    dplyr::rename(Elevation = Altitude) %>%
    dplyr::group_by(
      PlotObservationID, Latitude, Longitude, Elevation, Releve_area,
      Cover_bare_soil, Cover_bare_rock, Cover_open_water, Cover_forbs,
      Naturalness, Grassland, Shrubland, Forest, Wetland, Layer, Status
    ) %>%
    dplyr::summarise(
      Total_abundance = sum(Rel_Abund_Layer, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = c(
        PlotObservationID, Latitude, Longitude, Elevation, Releve_area,
        Cover_bare_soil, Cover_bare_rock, Cover_open_water, Cover_forbs,
        Naturalness, Grassland, Shrubland, Forest, Wetland,
        Elevation, Naturalness, Layer
      ),
      names_from = Status,
      values_from = c(Total_abundance, N)
    ) %>%
    tidyr::replace_na(
      replace = list(
        "Total_abundance_native" = 0, "Total_abundance_introduced" = 0,
        "N_native" = 0, "N_introduced" = 0
      )
    ) %>%
    dplyr::mutate(N_species = N_native + N_introduced) %>%
    dplyr::filter(
      !is.na(Releve_area)
    ) %>%
    dplyr::group_by(dplyr::across(-Layer)) %>%
    dplyr::summarise(
      Total_abundance_native = max(Total_abundance_native, na.rm = TRUE),
      Total_abundance_introduced = max(Total_abundance_introduced, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    terra::vect(geom = c("Longitude", "Latitude"), crs = "epsg:4326") %>%
    dplyr::mutate(
      Relative_abundance_native = Total_abundance_native / (
        Total_abundance_introduced + Total_abundance_native
      ),
      Relative_abundance_introduced = Total_abundance_introduced / (
        Total_abundance_introduced + Total_abundance_native
      ),
      BiomePatchID = terra::extract(biomes, .) %>%
        dplyr::pull(ID),
      HumanModification = terra::extract(ghm, ., method = "bilinear") %>%
        dplyr::pull(gHM)
    )
  terra::writeVector(sPlotVect, filename = f_splotvect, overwrite = TRUE)
}


df_sPlotBiome <- sPlotVect %>%
  dplyr::filter(
    !is.na(BiomePatchID),
    #Naturalness == "Natural" # Current sPlot version: Mostly NA, naturalness encoded as integers
    ) %>%
  dplyr::group_by(# Summarise plots within a biome patch by relative area
    BiomePatchID, Grassland, Shrubland, Forest, Wetland
    ) %>%
  dplyr::summarise(
    Native_abundance = stats::weighted.mean(
      Relative_abundance_native,
      w = Releve_area,
      na.rm = TRUE
      ),
    Introduced_abundance = stats::weighted.mean(
      Relative_abundance_introduced,
      w = Releve_area,
      na.rm = TRUE
      ),
    Bare_soil_cover = stats::weighted.mean(
      Cover_bare_soil,
      w = Releve_area,
      na.rm = TRUE
    ),
    HumanModification = stats::weighted.mean(
      HumanModification,
      w = Releve_area,
      na.rm = TRUE
      ),
    Sum_area = sum(Releve_area, na.rm = TRUE),
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
    as.numeric(Biome) <= 12,
    speciesRichnessBa > 30
    ) %>%
  dplyr::mutate(# Remove veg. types that do not match expected climax veg.
    Biome = base::droplevels(Biome),
    BiomeOK = BIOME %in% c(
      1:6, 12
      ) & Forest | BIOME %in% 7:11 & Grassland | BIOME %in% c(
        9, 11, 14
        ) & Wetland | BIOME %in% c(8, 9, 12, 13) & Shrubland
  ) %>%
  dplyr::filter(BiomeOK) %>%
  dplyr::group_by(# Summarise again, to average over allowed vegetation types
    dplyr::across(-c(Grassland, Shrubland, Forest, Wetland))
    ) %>%
  dplyr::summarise(
    Native_abundance = stats::weighted.mean(
      Native_abundance,
      w = Sum_area,
      na.rm = TRUE
      ),
    Introduced_abundance = stats::weighted.mean(
      Introduced_abundance,
      w = Sum_area,
      na.rm = TRUE
      ),
    Bare_soil_cover = stats::weighted.mean(
      Bare_soil_cover,
      w = Sum_area,
      na.rm = TRUE
      ),
    HumanModification = stats::weighted.mean(
      HumanModification,
      w = Sum_area,
      na.rm = TRUE
      ),
    .groups = "drop"
    ) %>%
  dplyr::mutate(
    log_introducedAbundance = log1p(Introduced_abundance),
    log_speciesRichnessBa = log(speciesRichnessBa),
    log_sampling_effort = log(sampling_effort),
    log_Connectedness = log(connectedness),
    log_focalECA = log(focalECA)
  ) %>%
  as.data.frame()


bi_tab <- table(df_sPlotBiome$Biome)
gg_abundance_list <- list()
for (bi in names(bi_tab)[bi_tab > 10]) {
  summary_tab <- df_sPlotBiome %>%
    dplyr::filter(
      Biome == bi
    ) %>%
    dplyr::summarise(
      median_invasion = median(
        Introduced_abundance[Introduced_abundance > 0],
        na.rm = TRUE
      ),
      n_none = sum(Introduced_abundance == 0, na.rm = TRUE),
      n_low = sum(
        Introduced_abundance > 0 &
          Introduced_abundance < median_invasion,
        na.rm = TRUE
      ),
      n_high = sum(
        Introduced_abundance >= median_invasion,
        na.rm = TRUE
      )
    )
  
  lvls <- c(
    paste0("none (N=", summary_tab$n_none, ")"),
    paste0("low (N=", summary_tab$n_low, ")"),
    paste0("high (N=", summary_tab$n_high, ")")
  )
  
  
  df_biome <- df_sPlotBiome %>%
    dplyr::filter(
      Biome == bi,
      !is.na(Introduced_abundance)
    ) %>%
    dplyr::mutate(
      Invasion_level = factor(
        dplyr::case_when(
          Introduced_abundance == 0 ~ lvls[1],
          Introduced_abundance < summary_tab$median_invasion ~ lvls[2],
          Introduced_abundance >= summary_tab$median_invasion ~ lvls[3]
        ),
        levels = lvls,
        ordered = TRUE
      )
    )
  
  gg_abundance_list[[bi]] <- ggplot2::ggplot(
    data = df_biome,
    ggplot2::aes(x = Invasion_level, y = focalECA, fill = Invasion_level)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(bi) +
    ggplot2::xlab("Introduced invasive species abundance") +
    ggplot2::ylab("Connected area (log scaled)")
}


# Fit a GLMM to estimate the impact of our variables
## Assumptions violated:
# n <- nrow(df_sPlotBiome)
# df_sPlotBiome$Introduced_abund <- (df_sPlotBiome$Introduced_abundance * (n - 1) + 0.5) / n
# 
# mod <- glmmTMB::glmmTMB(
#   Introduced_abund ~ HumanModification + log_speciesRichnessBa + log_focalECA +
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

vars <- c(
  "log_focalECA",
  "HumanModification",
  "climateVelocityRank"
)

fam <- stats::gaussian()

frml <- make_formula(
  predictors = vars,
  response = "Introduced_abundance_logit_boundary_corrected",
  k = 0,
  biome = FALSE
  )

# if(!dir.exists(file.path(dir_fig, "sPlotResponse"))) {
#   dir.create(file.path(dir_fig, "sPlotResponse"), recursive = TRUE)
# }
subtables_tex = character(0)
binames <- names(bi_tab[bi_tab >= 15])
for (bi in binames) {
  biname_clean <- gsub("[^[:alpha:]]", "", bi, perl = TRUE)
  
  sdf <- df_sPlotBiome %>%
    dplyr::filter(Biome == bi) %>%
    dplyr::mutate(
      climateVelocityRank = qnorm(
        rank(climate_velocity_kmpa) / (length(climate_velocity_kmpa) + 1)
      )
    )
  cor(sdf[vars], use = "pairwise.complete.obs")
  
  n <- sum(!is.na(sdf$Introduced_abundance))
  sdf$Introduced_abundance_logit_boundary_corrected <- qlogis(
    (sdf$Introduced_abundance * (n - 1) + 0.5) / n
  )
  mod_full <- stats::glm(
    frml,
    family = fam,
    data = sdf
  )
  
  ci <- boot_conf_int(sdf = sdf, frml = frml, fam = fam)
  
  d2_adj_full <- ecospat::ecospat.adj.D2.glm(mod_full)
  cat("\n\n", bi, "D2:", d2_adj_full, "(GLM)\n")
  coefs <- summary(mod_full)$coefficients
  
  par(mfrow = c(2, 2), mar = c(2, 2, 2, 0))
  plot(mod_full)
  par(mfrow = c(1, 1))
  
  # Check variable importance
  stats::drop1(mod_full, test = "Chisq")
  df_list_vars <- list()
  for (i in 1:length(vars)) {
    var <- vars[i]
    frml_red <- make_formula(
      predictors = vars[-i],
      response = "Introduced_abundance_logit_boundary_corrected",
      k = 0,
      biome = FALSE
    )
    mod_red <- glm(
      frml_red,
      family = fam,
      data = sdf
    )
    delta_aic <- AIC(mod_red) - AIC(mod_full)
    tabanova <- anova(mod_red, mod_full, test = "Chisq")
    dev <- tabanova$Deviance[2]
    p <- tabanova[2, 5]
    df <- tabanova$Df[2]
    df_list_vars[[i]] <- data.frame(
      Variable = var, df = df, dev = dev, p_Chisq = p, delta_aic = delta_aic
      )
  }
  
  df_vars <- ci %>%
    dplyr::rename(
      Variable = predictor
    ) %>%
    dplyr::left_join(
      do.call(rbind, df_list_vars),
      by = "Variable"
    )  %>%
    dplyr::relocate(Variable)
  
  print(df_vars)
  
  df_out <- df_vars %>%
    dplyr::filter(Variable != "(Intercept)") %>%
    dplyr::mutate(
      Variable = dplyr::recode(
        Variable,
        "log_focalECA" = "Connected area",
        "HumanModification" = "Human modification",
        "climateVelocityRank" = "Climate velocity"
      ),
      estimate_ci = sprintf(
        "%.2f [%.2f, %.2f]",
        estimate,
        ci95_low,
        ci95_high
      ),
      sign_stability = sprintf("%.1f", sign_stability),
      p_Chisq = sprintf("%.3f", p_Chisq),
      delta_aic = sprintf("%.2f", delta_aic)
    ) %>%
    dplyr::select(
      Variable,
      estimate_ci,
      sign_stability,
      p_Chisq,
      delta_aic
    )
  
  
  rows <- apply(
    df_out,
    1,
    function(x) paste(x, collapse = " & ")
  )
  
  rows <- paste0(rows, " \\\\")
  
  subtables_tex <- c(
    subtables_tex,
    "\\begin{subtable}[t]{0.95\\textwidth}",
    "\\centering",
    sprintf(
      "\\caption{%s}",
      paste0(bi, " (\\(D^2\\) = ", sprintf("%.3f", d2_adj_full), ")")
      ),
    "\\begin{tabular}{p{4cm}cS[table-format=3.1]S[table-format=1.3]S[table-format=2.2]}",
    "\\toprule",
    "{Variable} & {Estimate [95\\% CI]} & {Sign stability (\\%)} & {$p$} & {$\\Delta$AIC} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{subtable}",
    ifelse(bi != binames[length(binames)], "\\par\\medskip", "")
  )
}

writeLines(
  subtables_tex,
  file.path(dir_tab, paste("glm", nativeness_source, "subtab.tex", sep = "_"))
)


#>=============================================================================<
#> Biome-patch-level analyses on invasion probability
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
future::plan(multisession, workers = parallel::detectCores() - 2)

f_invasion <- file.path(dir_imed, "df_invasion.Rsave")
if (file.exists(f_invasion)) {
  load(f_invasion)
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
            maxCorrectedRichness = max(logCorrectedSR, na.rm = TRUE),
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
            speciesRichnessBa, sampling_effort, logCorrectedSR,
            gHM, RoadDensity, gdp1990to2020,
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
      correctedRichnessRatio = logCorrectedSR / donor_maxCorrectedRichness,
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
      logSamplingEffort = log(sampling_effort),
      logGDP = log(gdp1990to2020),
      logRoadDensity = log(RoadDensity)
    ) %>%
    as.data.frame()
  
  save(df_invasion, file = f_invasion)
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
hist(df_invasion$logGDP, main = "GDP (log transformed)")
hist(df_invasion$logRoadDensity, main = "Road density (log transformed)")
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
  dplyr::mutate(
    climateVelocityRank = qnorm(
      rank(climateVelocity) / (length(climateVelocity) + 1)
      ),
    weight = ifelse(Invaded == 1, w, 1),
    )

predictors <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "climateVelocityRank",
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
  gc()
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
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank()
    ) +
  #ggplot2::ggtitle("Explained Deviance") +
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
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "DevExplPie.svg"),
  plot = gg_pie,
  width = 5, height = 5.5
  )


# Fit model per biome, estimate explained deviance contributions
# and plot response shapes

#
# Maybe worth a try: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/ginla.html
#
stop("With 100 repetitions per variable, this now takes forever. Consider using the existing plots.")

response_labeller <- function(x) {
  is_log <- grepl("^log", x)
  x <- gsub("^log", "", x)
  x <- gsub("([a-z])([A-Z])", "\\1 \\2", x)
  x <- gsub("Conn", "Connected", x)
  x <- gsub("Rank", "(rank-normalised)", x)
  x <- trimws(x)
  x <- paste0(
    toupper(substr(x, 1, 1)),
    substr(x, 2, nchar(x))
  )
  x <- ifelse(
    is_log,
    paste0(x, " (log-scaled)"),
    x
  )
  x
}

frml_bi <- make_formula(predictors, biome = FALSE)
dir_plots <- file.path(dir_fig, "invasion_prob_by_biome")
if(!dir.exists(dir_plots)) {
  dir.create(dir_plots, showWarnings = FALSE)
}
plot_list <- list()
for (bi in sort(unique(df_mod$Biome))) {
  w <- subset(df_mod, Biome == bi) %>%
    dplyr::group_by(Species, Invaded) %>%
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::group_by(Invaded) %>%
    dplyr::summarise(
      sum = sum(N, na.rm = TRUE)
    ) %>%
    tidyr::pivot_wider(names_from = Invaded, values_from = sum) %>%
    dplyr::transmute(w = `0` / `1`) %>%
    dplyr::pull(w)
  
  df_bi <- subset(df_mod, Biome == bi) %>%
    dplyr::mutate(
      climateVelocityRank = qnorm(
        rank(climateVelocity) / (length(climateVelocity) + 1)
      ),
      weight = ifelse(Invaded == 1, w, 1),
    )
  
  mod_bi <- mgcv::gam(
    frml_bi,
    data = df_bi,
    method = "REML",
    family = stats::binomial("logit"),
    weights = weight
  )
  adjD2_bi <- ecospat::ecospat.adj.D2.glm(mod_bi)
  
  # Get delta D2 for each predictor
  delta_adjD2_p <- c()
  for (p in predictors) {
    adjD2_rand_c <- c()
    for (i in 1:100) {
      set.seed(i)
      df_rand <- df_bi %>%
        dplyr::mutate(
          "{p}" := sample(.data[[p]])
        )
      
      mod_rand <- mgcv::gam(
        frml_bi,
        data = df_rand,
        method = "REML",
        family = stats::binomial("logit"),
        weights = weight
      )
      adjD2_rand_c <- c(adjD2_rand_c, ecospat::ecospat.adj.D2.glm(mod_rand))
    }
    adjD2_rand <- mean(adjD2_rand_c, na.rm = TRUE)
    delta_adjD2_p <- c(delta_adjD2_p, adjD2_bi - adjD2_rand)
  }
  names(delta_adjD2_p) <- predictors
  
  # Plot smooth responses
  pred_dfs <- list()
  for (p in predictors) {
    nd <- data.frame(
      matrix(
        ncol = length(predictors),
        nrow = 200
      )
    )
    
    names(nd) <- predictors
    
    # Set all variables to median
    for (v in predictors) {
      nd[[v]] <- median(df_bi[[v]], na.rm = TRUE)
    }
    
    nd[[p]] <- seq(
      min(df_bi[[p]], na.rm = TRUE),
      max(df_bi[[p]], na.rm = TRUE),
      length.out = 200
    )
    
    pred <- predict(
      mod_bi,
      newdata = nd,
      type = "lpmatrix"#"link",
      #se.fit = TRUE
    )
    b <- coef(mod_bi)
    V <- vcov(mod_bi)
    sim_b <- MASS::mvrnorm(n = 2000, mu = b, Sigma = V)
    fit_sim <- plogis(sim_b %*% t(pred))
    fit_mean <- colMeans(fit_sim)
    fit_low <- apply(fit_sim, 2, quantile, 0.025)
    fit_high <- apply(fit_sim, 2, quantile, 0.975)
    
    d <- data.frame(
      predictor = p,
      x = nd[[p]],
      fit = fit_mean,
      lower = fit_low,
      upper = fit_high,
      deltaD2 = delta_adjD2_p[p]
    )
    
    pred_dfs[[p]] <- d
  }
  
  plot_df <- dplyr::bind_rows(pred_dfs)
  
  # Plot response shapes
  pred_plots <- list()
  for (p in predictors) {
    d <- plot_df %>%
      dplyr::filter(predictor == p)
    
    delta_d2 <- unique(d$deltaD2)
    
    # Response curve
    p_curve <- ggplot2::ggplot(
      data = d,
      ggplot2::aes(x = x, y = fit)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = lower,
          ymax = upper
        ),
        fill = "grey70",
        alpha = 0.4
      ) +
      ggplot2::geom_line(
        linewidth = 1
      ) +
      ggplot2::labs(
        title = response_labeller(p),
        x = NULL,
        y = NULL
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 10,
          face = "bold"
        ),
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        plot.margin = grid::unit(
          c(2, 2, 2, 2),
          "pt"
        )
      )
    
    # Explained deviance bar
    bar_df <- data.frame(
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = pmax(delta_d2 / adjD2_bi, 0)
    )
    
    p_bar <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = 0,
          xmax = 1,
          ymin = 0,
          ymax = 1
        ),
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ) +
      ggplot2::geom_rect(
        data = bar_df,
        ggplot2::aes(
          xmin = 0,
          xmax = 1,
          ymin = 0,
          ymax = fill
        ),
        fill = "black"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = grid::unit(
          c(2, 2, 2, 2),
          "pt"
        )
      )
    pred_plots[[p]] <- p_curve + p_bar +
      patchwork::plot_layout(
        widths = c(12, 1)
      )
  }
  # --- Combine all predictors vertically ---
  p_final <- patchwork::wrap_plots(
    pred_plots,
    ncol = 1
  ) +
    patchwork::plot_annotation(
      title = paste0(bi, " | adj. D² = ", round(adjD2_bi, 3))
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 8,
        face = "bold"
      )
    )
  
  plot_list[[bi]] <- p_final
  bi_idx <- which(levels(df_mod$Biome) == bi)
  ggplot2::ggsave(
    filename = file.path(dir_plots, paste0("Response_", bi_idx, ".svg")),
    plot = p_final, width = 3, height = 7
    )
  ggplot2::ggsave(
    filename = file.path(dir_plots, paste0("Response_", bi_idx, ".pdf")),
    plot = p_final, width = 3, height = 7
  )
}

plot_list[[1]]


#>=============================================================================<
#> Invasive species breeder model
#>=============================================================================<


# Open questions: Do we also need patches that donate zero invasives to meet
# distribution assumptions?
#
# Should we assign each invasive species to the patch with most observations,
# or should we count each invasive for each patch where it has been observed?
# The former one might solve the "issue" of zero observations.
#
# Why would sampling intensity explain anything here? (It is a bit logical if we assign by "most counts")

df_donor <- merged_env %>%
  dplyr::filter(
    as.numeric(Biome) <= 12,
    as.numeric(Biome) != 10,
    total_area > 100,
    #Status == "Donor",
    logCorrectedSR > log(100)
  ) %>%
  dplyr::mutate(
    Biome = base::droplevels(Biome),
    PatchID = ID,
    correctedSR = exp(logCorrectedSR),
    logConnArea = log(focalECA),
    SpeciesRichness = tidyr::replace_na(round(correctedSR, 0), 0)
    ) %>%
  # dplyr::group_by(Species, Status) %>%
  # dplyr::filter( # Assign the origin of each invasive species to the patch with its highest frequency
  #   Count == max(Count)
  #   ) %>%
  dplyr::group_by(
    Biome, PatchID, log_total_area, focalECA, connectedness,
    SpeciesRichness, Status, climateVelocity, climateStability,
    log_sampling_effort
    ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(
      Biome, PatchID, log_total_area, focalECA, connectedness,
      SpeciesRichness, climateVelocity, climateStability,
      log_sampling_effort
      ),
    names_from = Status,
    values_from = N
  ) %>%
  tidyr::replace_na(replace = list(Donor = 0, Receiver = 0)) %>%
  dplyr::rename(
    DonorOf = Donor
  ) %>%
  dplyr::filter(
    DonorOf < SpeciesRichness
    ) %>%
  dplyr::mutate(
    log_focalECA = log(focalECA),
    PropInvasive = DonorOf / SpeciesRichness,
    failures = pmax(SpeciesRichness - DonorOf, 0),
    successRate = DonorOf / SpeciesRichness,
    successProb = successRate / (1 - successRate)
  )


## GLM by biome -> variance exceeds the expected under binomial sampling process
d2s <- c()
for (bi in levels(df_donor$Biome)) {
  df_curr <- df_donor %>%
    dplyr::filter(Biome == bi) %>%
    dplyr::mutate(
      climateVelocityRank = qnorm(
        rank(climateVelocity) / (length(climateVelocity) + 1)
      )
    )
  
  modInvasiveBreeder_glm <- stats::glm(
    cbind(DonorOf, failures) ~ log_focalECA + climateVelocityRank + log_sampling_effort,
    data = df_curr,
    family = stats::quasibinomial()
  )
  d2 <- ecospat::ecospat.adj.D2.glm(modInvasiveBreeder_glm)
  summary(modInvasiveBreeder_glm)
  coef_tab_glm <- summary(modInvasiveBreeder_glm)$coefficients
  d2s <- c(d2s, d2)
  cat("\n", bi, "D2:", d2, "\n")
  print(coef_tab_glm)
  par(mfrow = c(2, 2))
  plot(modInvasiveBreeder_glm)
  par(mfrow = c(1, 1))
  
  cat("\nDispersion:", dispersion(modInvasiveBreeder_glm), "\n\n")
}


## Betabinomial GLMM -> Model assumptions still violated
modInvasiveBreeder_bb <- glmmTMB::glmmTMB(
  cbind(DonorOf, failures) ~ log(focalECA) + log_sampling_effort + (1|Biome),
  data = df_donor,
  family = glmmTMB::betabinomial()
)

summary(modInvasiveBreeder_bb)
performance::check_singularity(modInvasiveBreeder_bb)
performance::r2(modInvasiveBreeder_bb)
MuMIn::r.squaredGLMM(modInvasiveBreeder_bb)
res <- DHARMa::simulateResiduals(modInvasiveBreeder_bb)
DHARMa::testDispersion(res)
DHARMa::plotResiduals(res)
DHARMa::testZeroInflation(res)


## GAM
predictors <- c("log_focalECA", "climateVelocityRank", "log_sampling_effort")
frml <- make_formula(
  predictors = predictors, response = "cbind(DonorOf, failures)", biome = FALSE
)
d2s <- c()
for (bi in levels(df_donor$Biome)) {
  df_curr <- df_donor %>%
    dplyr::filter(Biome == bi) %>%
    dplyr::mutate(
      climateVelocityRank = qnorm(
        rank(climateVelocity) / (length(climateVelocity) + 1)
      )
    )
  
  modInvasiveBreeder_gam <- mgcv::gam(
    frml,
    family = stats::quasibinomial(),
    data = df_curr
  )
  summary(modInvasiveBreeder_gam)
  
  # Plot some diagnostics
  par(mfrow = c(2, 2))
  gam.check(modInvasiveBreeder_gam, cex = 5)
  mtext(bi, side = 3, line = 3)
  par(mfrow = c(1, 1))
  d2 <- ecospat::ecospat.adj.D2.glm(modInvasiveBreeder_gam)
  par(mfrow = c(2, 2))
  plot(modInvasiveBreeder_gam, residuals = TRUE, cex = 5, main = bi)
  par(mfrow = c(1, 1))
  # modInvasiveBreeder_gamm <- glmmTMB::glmmTMB(
  #   frml,
  #   family = glmmTMB::betabinomial(),
  #   data = df_curr
  # )
  # par(mfrow = c(2, 2))
  # plot(modInvasiveBreeder_gamm)
  # par(mfrow = c(1, 1))
  # d2 <- ecospat::ecospat.adj.D2.glm(modInvasiveBreeder_gamm)
  # summary(modInvasiveBreeder_gamm)
  # res <- DHARMa::simulateResiduals(modInvasiveBreeder_gamm)
  # plot(res)
  # DHARMa::testDispersion(res)
  # summary(modInvasiveBreeder_gamm)
  # MuMIn::r.squaredGLMM(modInvasiveBreeder_gamm)
  d2s <- c(d2s, d2)
  cat("\n", bi, "D2:", d2, "\n")
}
