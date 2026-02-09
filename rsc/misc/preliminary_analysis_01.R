# Preliminary analysis based on sPlot data of SINAS listed (invasive) species
# to assess introduced species cover fractions

library(rsdd)
library(terra)
library(tidyterra)
library(pbapply)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

recompute <- TRUE

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

f_stats <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/stats/glm_summaries.txt"
cat(
  ">--------------------Statistical model outputs--------------------<",
  file = f_stats, append = FALSE
)

# Label ruderal species and treat them differently
# Check with Gengchen the biome assignment workflow
rsdd::dataset(2)

f_sinas_places <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/sinas/SInAS_Locations"
f_sinas_data <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/sinas/SInAS_3.1.1.csv"

sinas_places <- terra::vect(
  f_sinas_places
) %>%
  terra::project("epsg:4326")

sinas_data <- read.table(f_sinas_data, sep = " ", header = TRUE)

biomes <- terra::vect(
  "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/biomes/biomes_treediv.gpkg"
)
biomes$ID <- seq(1:nrow(biomes))

# sPlot Open
e <- new.env()
load(
  "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/sPlotOpen/3474_76_sPlotOpen.RData",
  envir = e
)
sPlot <- as.list(e)

sPlotData <- sPlot$DT2.oa %>%
  dplyr::mutate(
    Invader = Species %in% unique(sinas_data$taxon)
  ) %>%
  dplyr::left_join(
    y = sPlot$header.oa,
    by = "PlotObservationID"
  )

## Status levels
## Unknown: Known invader but unknown status at the site
## Assumed native: Species not known as invader and, thus, assumed to be native
## Native: Native here, but introduced elsewhere according to db
## Introduced: Introduced according to db

assign_status <- function(row) {
  if (length(unique(row$Species)) != 1) {
    stop("Only single species data frames are allowed.")
  }
  
  if (!row$Invader[1]) {
    row$Status <- "Assumed native"
    return(row)
  }
  
  ds <- sinas_data[sinas_data$taxon == row$Species[1], ]
  if (nrow(ds) == 0) {
    row$Status <- "Assumed native"
    return(row)
  }
  
  native_loc_ids <- ds[ds$establishmentMeans == "native", ]$location
  introduced_loc_ids <- ds[ds$establishmentMeans == "introduced", ]$location
  native <- sinas_places[sinas_places$Location %in% native_loc_ids, ]
  introduced <- sinas_places[sinas_places$Location %in% introduced_loc_ids, ]
  
  if (terra::is.empty(native) | terra::is.empty(introduced)) {
    row$Status <- "Unknown"
    return(row)
  }
  
  location <- terra::vect(
    as.matrix(row[, c("Longitude", "Latitude")]),
    type = "points",
    crs = "EPSG:4326"
  )
  
  row$Status <- NA
  row$Status[
    which(terra::is.related(location, native, relation = "within"))
  ] <- "Native"
  
  row$Status[
    which(terra::is.related(location, introduced, relation = "within"))
  ] <- "Introduced"
  
  return(row)
}

sPlot_split <- split(sPlotData, sPlotData$Species)

f_dat <- "C:/Users/poppman/Desktop/preliminary_analysis_01.Rdata"

if (file.exists(f_dat)) {
  load(f_dat)
} else {
  print(Sys.time())
  n_workers <- 10
  species_split <- split(
    sPlot_split, rep(1:n_workers, length.out = length(sPlot_split))
  )
  
  cl <- parallel::makeCluster(n_workers)
  parallel::clusterExport(
    cl,
    varlist = c("sinas_data", "f_sinas_places", "assign_status")
  )
  result_chunks <- parallel::parLapply(
    cl, species_split,
    function(chunk) {
      library(terra)
      library(dplyr)
      sinas_places <- terra::vect(f_sinas_places) |> terra::project("EPSG:4326")
      do.call(rbind, lapply(chunk, assign_status))
    }
  )
  
  sPlotData_processed <- do.call(rbind, result_chunks)
  stopCluster(cl)
  
  status_list <- lapply(
    sPlot_split,
    FUN = assign_status
  )
  
  save(sPlotData_processed, file = f_dat)
}

sdat <- sPlotData_processed %>%
  dplyr::filter(!is.na(Releve_area)) %>%
  dplyr::filter(Releve_area < 500) %>%
  dplyr::group_by(
    PlotObservationID, Status, Biome, Latitude, Longitude, Elevation, Releve_area
    ) %>%
  dplyr::summarise(
    Species_count = dplyr::n(),
    Abundance = sum(Original_abundance),
    Abundance_scale = first(Abundance_scale),
    Cover = sum(Relative_cover)
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(PlotObservationID, Biome, Latitude, Longitude, Elevation, Releve_area),
    names_from = Status,
    values_from = c(Species_count, Cover)
  ) %>%
  dplyr::mutate(
    Species_count_probnative = sum(
      `Species_count_Assumed native`, Species_count_Native,
      Species_count_Unknown, 0
    ),
    Introduced_fraction_n = sum(
      Species_count_Introduced, 0,
      na.rm = TRUE
      ) / sum(
        `Species_count_Assumed native`, Species_count_Native,
        Species_count_Unknown, Species_count_Introduced, Cover_NA, 0,
        na.rm = TRUE
        ),
    Introduced_fraction_cover = sum(
      Cover_Introduced, 0,
      na.rm = TRUE
    ) / sum(
      `Cover_Assumed native`, Cover_Native, Cover_Unknown, Cover_Introduced,
      Cover_NA, 0,
      na.rm = TRUE
    )
    )

biomes$Area <- terra::expanse(biomes, unit = "m")
biomes$Olson_biome <- biome_names[biomes$BIOME]

svect <- terra::vect(
  sdat, geom = c("Longitude", "Latitude"), crs = "EPSG:4326"
  ) %>%
  terra::intersect(biomes)

sdata <- as.data.frame(svect) %>%
  dplyr::filter(
    !is.na(Species_count_probnative),
    !is.na(Introduced_fraction_cover),
    !is.na(Releve_area),
    Releve_area > 0,
    Species_count_probnative > 0,
    !is.na(Olson_biome)
  ) %>%
  dplyr::group_by(
    Olson_biome
  ) %>%
  dplyr::filter(
    dplyr::n() > 3
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Biome = factor(Biome),
    BIOME = factor(BIOME)
  ) %>%
  dplyr::filter(
    Olson_biome != "Flooded Savanna" # Few samples & unable to fit SAR model
  )

# Use a species-area-relationship model to correct for sampling plot area
sar_models <- do.call(
  rbind,
  lapply(
    X = unique(sdata$Olson_biome),
    FUN = function(biome) {
      df <- sdata[which(sdata$Olson_biome == biome),]
      mod <- lm(log(Species_count_probnative) ~ log(Releve_area), data = df)
      coefs <- as.vector(coefficients(mod))
      rsqadj <- summary(mod)$adj.r.squared
      pvals <- as.vector(summary(mod)$coefficients[, 4])
      out <- data.frame(
        Biome = biome,
        Intercept = coefs[1],
        Slope = coefs[2],
        adjR2 = rsqadj,
        p_intercept = pvals[1],
        p_slope = pvals[2]
        )
      return(out)
    }
  )
)

# Join model coefficients into the data.frame and apply a correction
A0 <- 100

sdata <- sdata %>%
  dplyr::left_join(
    sar_models %>% dplyr::select(Biome, Intercept, Slope),
    by = c("Olson_biome" = "Biome")
  ) %>%
  dplyr::mutate(
    Native_count_corrected = Species_count_probnative * 
      exp(Slope * (log(A0) - log(Releve_area)))
  )


test_biotic_resistance <- sdata %>%
  dplyr::group_by(Olson_biome) %>%
  dplyr::filter(dplyr::n() >= 10) %>%
  dplyr::summarise(
    n = dplyr::n(),
    spearman = list(
      cor.test(
        Native_count_corrected, Introduced_fraction_cover,
        method = "spearman",
        use = "complete.obs"
        )
      ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    rho = spearman$estimate,
    p_value = spearman$p.value
  ) %>%
  select(Olson_biome, n, rho, p_value) %>%
  dplyr::mutate(
    signif = ifelse(p_value > 0.05, "", ifelse(p_value > 0.01, "*", "**"))
  )

print(test_biotic_resistance)

#
test_biome_resistance <- sdata %>%
  dplyr::filter(
    dplyr::if_all(
      c(species_richness, Introduced_fraction_cover),
      ~ !is.na(.x)
    )
  ) %>%
  dplyr::group_by(Olson_biome) %>%
  dplyr::filter(dplyr::n() >= 10) %>%
  dplyr::summarise(
    n = dplyr::n(),
    spearman = list(
      cor.test(
        species_richness, Introduced_fraction_cover,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    rho = spearman$estimate,
    p_value = spearman$p.value
  ) %>%
  select(Olson_biome, n, rho, p_value) %>%
  dplyr::mutate(
    signif = ifelse(p_value > 0.05, "", ifelse(p_value > 0.01, "*", "**"))
  )

print(test_biome_resistance)

library(clipr)

write_clip(test_biotic_resistance)









specs <- unique(sinas_data$taxon)
#specs <- sample(specs, size = 20000)
specs <- specs[which(specs %in% trees$species)] # Use only invasive tree species

f_out <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/df_area.csv"
if (file.exists(f_out) & !recompute) {
  df_area <- read.csv(f_out)
} else {
  df_area <- do.call(
    rbind,
    pbapply::pblapply(
      X = specs,
      FUN = assign_patches
    )
  )
  
  write.csv(df_area, file = f_out, row.names = FALSE)
}

stats <- df_area %>%
  dplyr::group_by(Species, Biome, Status) %>%
  dplyr::summarise(
    MainOrigin = median(MainOrigin),
    MaxArea = max(Area),
    MaxRichness = max(SpeciesRichness),
    MedianArea = median(Area),
    MedianRichness = median(SpeciesRichness),
    MeanArea = mean(Area),
    MeanRichness = mean(SpeciesRichness)
  ) %>%
  dplyr::mutate(
    BiomeIsMain = Biome == MainOrigin
  )

withinInvasions = stats %>%
  dplyr::filter(BiomeIsMain)


# For each biome patch, get the number of invasive species as well as the number
# of native species invasive elsewhere
df_overall <- df_area %>%
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
  dplyr::mutate(
    FractionInvasive = Donor / SpeciesRichness,
    FractionInvaders = Receiver / SpeciesRichness
  ) %>%
  dplyr::rename(
    DonorOf = Donor,
    ReceiverOf = Receiver
  )

df_within <- df_area %>%
  dplyr::group_by(Biome, PatchID, Area, SpeciesRichness, Status) %>%
  dplyr::filter(Biome == MainOrigin) %>%
  dplyr::summarise(
    N = dplyr::n()
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(Biome, PatchID, Area, SpeciesRichness),
    names_from = Status,
    values_from = N
  ) %>%
  tidyr::replace_na(replace = list(Donor = 0, Receiver = 0)) %>%
  dplyr::mutate(
    FractionInvasive = Donor / SpeciesRichness,
    FractionInvaders = Receiver / SpeciesRichness
  ) %>%
  dplyr::rename(
    DonorOf = Donor,
    ReceiverOf = Receiver
  )

ggplot2::ggplot(
  data = df_within,
  ggplot2::aes(x = log(Area), y = log(FractionInvaders), colour = factor(Biome))
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", se = FALSE)

ggplot2::ggplot(
  data = df_within,
  ggplot2::aes(x = log(Area), y = log(FractionInvasive), colour = factor(Biome))
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", se = FALSE)


# Statistical models
## Receiver models
modReceiver <- MASS::glm.nb(
  ReceiverOf ~ log(Area) + log(SpeciesRichness) + factor(Biome),
  data = df_within
)

summary(modReceiver)

cat(
  "\n\nModel of overall probability of a patch to receive an invasive species\n",
  file = f_stats, append = TRUE
)
sink(f_stats, append = TRUE)
print(summary(modReceiver))
sink()

for (biome in unique(df_within$Biome)) {
  modInvasiveReceiverCurr <- MASS::glm.nb(
    ReceiverOf ~ log(Area) + log(SpeciesRichness),
    data = df_within[which(df_within$Biome == biome),]
  )
  
  cat(
    "\nSingle biome model: ", paste0(biome_names[biome], "\n"),
    file = f_stats, append = TRUE
  )
  sink(f_stats, append = TRUE)
  print(summary(modInvasiveReceiverCurr))
  sink()
}


## Donor models
modDonor <- MASS::glm.nb(
  DonorOf ~ log(Area) + log(SpeciesRichness) + factor(Biome),
  data = df_within
)

summary(modDonor)

df_within$failures <- pmax(df_within$SpeciesRichness - df_within$DonorOf, 0) 
modInvasiveBreeder <- stats::glm(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + factor(Biome),
  data = df_within,
  family = stats::quasibinomial()
)

summary(modInvasiveBreeder)

cat(
  "\n\nModel of overall probability of a patch to breed a successful invader\n",
  file = f_stats, append = TRUE
)
sink(f_stats, append = TRUE)
print(summary(modInvasiveBreeder))
sink()

for (biome in unique(df_within$Biome)) {
  modInvasiveBreederCurr <- stats::glm(
    cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness),
    data = df_within[which(df_within$Biome == biome),],
    family = stats::quasibinomial()
  )
  
  cat(
    "\nSingle biome model: ", paste0(biome_names[biome], "\n"),
    file = f_stats, append = TRUE
  )
  sink(f_stats, append = TRUE)
  print(summary(modInvasiveBreederCurr))
  sink()
}
