# Preliminary analysis based on GBIF data of SINAS listed (invasive) species

print(Sys.time())
library(rsdd)
library(terra)
library(tidyterra)
library(pbapply)
library(ggplot2)
library(dplyr)
library(tidyr)
library(progress)

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

f_stats <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/stats/glm_summaries.txt"
cat(
  ">--------------------Statistical model outputs--------------------<",
  file = f_stats, append = FALSE
)

rsdd::dataset(2)

f_invasive_tree_db <- "L:/poppman/data/bir/dat/lud11/db/invasive_tree_db/invasive_db_2013_normalised.csv"
trees <- read.table(f_invasive_tree_db, header = TRUE, sep = ",")

f_sinas_places <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/sinas/SInAS_Locations"
f_sinas_data <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/sinas/SInAS_3.1.1.csv"

sinas_places <- terra::vect(
  f_sinas_places
) %>%
  terra::project("epsg:4326")

sinas_data <- read.table(f_sinas_data, sep = " ", header = TRUE)

taxa <- rsdd::taxa()
tax_avail <- taxa$species

biomes <- terra::vect(
  "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/biomes/biomes_treediv.gpkg"
)
biomes$ID <- seq(1:nrow(biomes))

f_bbuff <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/dat/biomes/biomes_buff.gpkg"
if (file.exists(f_bbuff)) {
  biomes_buff <- terra::vect(f_bbuff)
} else {
  biomes_buff <- biomes %>%
    terra::buffer(width = -1000) %>%
    terra::writeVector(
      filename = f_bbuff,
      overwrite = TRUE
    )
}

biomes <- biomes[which(biomes$ID %in% biomes_buff$ID),]

assign_patches <- function(taxon) {
  ds <- sinas_data[sinas_data$taxon == taxon, ]
  native_loc_ids <- ds[ds$establishmentMeans == "native", ]$location
  introduced_loc_ids <- ds[ds$establishmentMeans == "introduced", ]$location
  native <- sinas_places[sinas_places$Location %in% native_loc_ids, ]
  introduced <- sinas_places[sinas_places$Location %in% introduced_loc_ids, ]
  
  if (terra::is.empty(native) | terra::is.empty(introduced)) {
    return(data.frame())
  }
  
  rst <- rsdd::get_taxon(taxon, format = "SpatRaster")
  
  native_obs <- terra::mask(
    rst, mask = native
  ) %>%
    terra::as.points()
  
  introduced_obs <- terra::mask(
    rst, introduced
  ) %>%
    terra::as.points()
  
  if (length(native_obs) < 10 | length(introduced_obs) < 10) {
    return(data.frame())
  }
  
  # Find patches of interest
  native_patches <- biomes[
    which(
      terra::is.related(
        biomes_buff, native_obs, relation = "covers"
        )
      ),
  ]
  
  # Obtain number of observations and area
  native_patches$obs_count <- terra::relate(
    native_patches, native_obs, "contains"
    ) %>%
    rowSums() %>%
    as.vector()
  
  native_patches$area <- terra::expanse(native_patches, unit = "m")
  
  # Find patches of interest
  introduced_patches <- biomes[
    which(
      terra::is.related(
        biomes_buff, introduced_obs, relation = "covers"
        )
      ),
  ]
  
  # Obtain number of observations and area
  introduced_patches$obs_count <- terra::relate(
    introduced_patches, introduced_obs, "contains"
    ) %>%
    rowSums() %>%
    as.vector()
  
  introduced_patches$area <- terra::expanse(introduced_patches, unit = "m")
  
  # Create data frames
  df_origin <- data.frame(
    Species = rep(taxon, length(native_patches$area)),
    Biome = native_patches$BIOME,
    PatchID = native_patches$ID,
    Area = native_patches$area,
    SpeciesRichness = native_patches$species_richness,
    Count = native_patches$obs_count,
    Status = rep("Donor", length(native_patches$area))
  )
  
  df_destination <- data.frame(
    Species = rep(taxon, length(introduced_patches$area)),
    Biome = introduced_patches$BIOME,
    PatchID = introduced_patches$ID,
    Area = introduced_patches$area,
    SpeciesRichness = introduced_patches$species_richness,
    Count = introduced_patches$obs_count,
    Status = rep("Receiver", length(introduced_patches$area))
  )
  
  df <- rbind(df_origin, df_destination)
  df$MainOrigin <- df_origin %>%
    dplyr::group_by(Biome) %>%
    dplyr::summarise(MaxCount = max(Count)) %>%
    dplyr::filter(MaxCount == max(MaxCount)) %>%
    dplyr::pull(Biome)
  
  df$LargestSource <- df_origin %>%
    dplyr::filter(Count == max(Count)) %>%
    dplyr::pull(PatchID)
  
  df$LargestSink <- df_destination %>%
    dplyr::filter(Count == max(Count)) %>%
    dplyr::pull(PatchID)
  
  df$RelFrequency = df$Count / df$Area
  
  return(df)
}

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

library(ggalluvial)

df_long <- df_area %>%
  dplyr::count(MainOrigin, Biome, name = "freq") %>%
  dplyr::mutate(
    Origin = biome_names[MainOrigin],
    Invaded = biome_names[Biome]
  )

gg_sankey <- ggplot2::ggplot(
  df_long,
  ggplot2::aes(
    axis1 = Origin,
    axis2 = Invaded,
    y = freq
  )
) +
  ggalluvial::geom_alluvium(
    ggplot2::aes(fill = Origin),
    discern = TRUE
    ) +
  ggalluvial::geom_stratum() +
  ggplot2::geom_label(
    stat = "stratum",
    ggplot2::aes(label = ggplot2::after_stat(stratum)),
    fill = "white",
    label.size = 0.25
    ) +
  ggplot2::scale_x_discrete(limits = c("Origin", "Invaded")) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_void()+
  ggplot2::theme(
    legend.position = "none"
  )

ggplot2::ggsave(
  filename = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/Sankey.svg",
  plot = gg_sankey,
  width = 10,
  height = 10
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

df_receiver <- data.frame()
for (biome in unique(df_within$Biome)) {
  modInvasiveReceiverCurr <- MASS::glm.nb(
    ReceiverOf ~ log(Area) + log(SpeciesRichness),
    data = df_within[which(df_within$Biome == biome),]
  )
  df <- as.data.frame(summary(modInvasiveReceiverCurr)$coefficients)
  df_receiver <- rbind(
    df_receiver,
    data.frame(
      Biome = biome,
      Area = df$Estimate[2],
      SpeciesRichness = df$Estimate[3],
      pArea = df[2, 4],
      pSpeciesRichness = df[3, 4]
    )
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

df_breeder <- data.frame()
for (biome in unique(df_within$Biome)) {
  modInvasiveBreederCurr <- stats::glm(
    cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness),
    data = df_within[which(df_within$Biome == biome),],
    family = stats::quasibinomial()
  )
  df <- as.data.frame(summary(modInvasiveBreederCurr)$coefficients)
  df_breeder <- rbind(
    df_breeder,
    data.frame(
      Biome = biome,
      Area = df$Estimate[2],
      SpeciesRichness = df$Estimate[3],
      pArea = df[2, 4],
      pSpeciesRichness = df[3, 4]
      )
  )
  cat(
    "\nSingle biome model: ", paste0(biome_names[biome], "\n"),
    file = f_stats, append = TRUE
    )
  sink(f_stats, append = TRUE)
  print(summary(modInvasiveBreederCurr))
  sink()
}

