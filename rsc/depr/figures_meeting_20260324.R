library("terra")
library("tidyterra")
library("rsdd")
library("leaflet")
library("leaflet.extras2")
library("progress")
rsdd::dataset("gbif-powo_0p5")
taxa <- rsdd::taxa()
tax_avail <- taxa$species

nativeness_source <- "sinas"
biome_def <- "olson"

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
  "Lake",
  "Rock and Ice"
)

if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir"
  dir_lud11 <- "L:"
} else {
  dir_main <- "/lud11/poppman/data/bir"
  dir_lud11 <- "/lud11"
}

dir_dat <- file.path(dir_lud11, "poppman", "data", "bir", "dat", "lud11")
dir_imed <- file.path(
  dir_dat, "biomes", biome_def, "intermediate_data", nativeness_source
  )

f_est_div <- file.path( # Species richness estimates
  dir_imed, "biome_species_richness.csv"
)

# Invasive species data base files
f_sinas_places <- file.path(dir_dat, "sinas", "SInAS_locations_3.1.gpkg")
f_sinas_data <- file.path(dir_dat, "sinas", "SInAS_3.1.1.csv")

# Biomes file
f_bfullinfo <- file.path(dir_imed, "biomes_full_info.gpkg")

biomes <- terra::vect(f_bfullinfo)

sinas_places <- terra::vect(
  f_sinas_places
)

sinas_data <- read.table(f_sinas_data, sep = " ", header = TRUE)

taxa <- unique(sinas_data$taxon)
pb <- progress_bar$new(
  format = "[:bar] :percent (:current/:total) eta: :eta",
  total = length(taxa),
  clear = FALSE,
  width = 60
)

if (!file.exists("D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map/introduced_richness.tif")) {
  native_richness <- 0
  introduced_richness <- 0
  for (taxon in taxa) {
    pb$tick()
    taxon_species_lvl <- paste(
      strsplit(taxon, split = " ")[[1]][c(1, 2)],
      collapse = " "
    )
    if (taxon %in% tax_avail) {
      rsdd_taxon <- taxon
    } else if (taxon_species_lvl %in% tax_avail) {
      rsdd_taxon <- taxon_species_lvl
    } else {
      next
    }
    
    native_obs <- rsdd::get_taxon(
      rsdd_taxon, status = "native",
      format = "spatraster"
    )
    if (!is.null(native_obs)) {
      native_richness <- native_richness + native_obs
    }
    
    introduced_obs <- rsdd::get_taxon(
      rsdd_taxon, status = c("non-native", "N/A"),
      format = "spatraster"
    )
    if (!is.null(introduced_obs)) {
      introduced_richness <- introduced_richness + introduced_obs
    }
  }
  
  native_masked <- terra::mask(native_richness, biomes)
  
  terra::writeRaster(
    native_masked,
    filename = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map/native_richness.tif",
    overwrite = TRUE
  )
  
  introduced_masked <- terra::mask(introduced_richness, biomes)
  terra::writeRaster(
    introduced_masked,
    filename = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map/introduced_richness.tif",
    overwrite = TRUE
  )
}

valid_extent <- terra::ext(-179, 179, -85, 85)
native_masked <- terra::rast(
  "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map/native_richness.tif"
  ) %>%
  terra::project("epsg:4326") %>%
  #terra::aggregate(fact = 3, fun = "max") %>%
  terra::crop(valid_extent)
introduced_masked <- terra::rast(
  "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map/introduced_richness.tif"
  ) %>%
  terra::project("epsg:4326") %>%
  #terra::aggregate(fact = 3, fun = "max") %>%
  terra::crop(valid_extent)

names(native_masked) <- "native"
names(introduced_masked) <- "introduced"
native_masked <- native_masked[["native"]]
introduced_masked <- introduced_masked[["introduced"]]

est_div <- read.csv(f_est_div)

biomes <- terra::vect(f_bfullinfo) %>%
  terra::merge(est_div, by = "ID", all.x = TRUE)
biomes$speciesRichnessBa[which(biomes$speciesRichnessBa > 1e5)] <- NA

biomes$BIOME[which(biomes$BIOME == 98)] <- 15
biomes$BIOME[which(biomes$BIOME == 99)] <- 16

biomes_sf <- biomes %>%
  dplyr::mutate(
    species_richness = round(speciesRichnessBa),
    BiomeID = BIOME,
    Biome = factor(
      biome_names[BIOME],
      levels = biome_names
    ),
    total_area_km2 = total_area / 1e6,
    lowland_area_km2 = lowland_area / 1e6
    ) %>%
  dplyr::select(
    ID, BiomeID, Biome, total_area_km2, lowland_area_km2, dECA, species_richness
    ) %>%
  sf::st_as_sf()

biome_colours <- c(
  "#098742", "#c7b839", "#9cce4e", "#1f7762", "#087186", "#81c4a1", "#ffa42c",
  "#ffd338", "#66d0c3", "#cea675", "#bddf98", "#ff2e17",
  "grey75", "pink", "steelblue1", "grey95"
)[unique(rev(biomes_sf$BiomeID))]#[match(1:16, unique(biomes_sf$BiomeID))]

biome_colours_named <- setNames(biome_colours, biome_names)
pal_biome <- colorFactor(
  palette = biome_colours_named,
  domain = levels(biomes_sf$Biome)
)
pal_nat <- colorNumeric(
  "Greens",
  domain = terra::minmax(native_masked),
  na.color = "transparent"
)
pal_int <- colorNumeric(
  "Reds",
  domain = terra::minmax(introduced_masked),
  na.color = "transparent"
  )

setwd("D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/map")
leaflet::leaflet() %>%
  leaflet::addTiles() %>%
  # Native raster
  leaflet::addRasterImage(
    native_masked,
    colors = pal_nat,
    #opacity = 0.7,
    group = "Native",
    options = leaflet::gridOptions(maxBytes = 50 * 1024 * 1024),
    project = TRUE
  ) %>%
  # Introduced raster
  leaflet::addRasterImage(
    introduced_masked,
    colors = pal_int,
    #opacity = 0.7,
    group = "Introduced",
    options = leaflet::gridOptions(maxBytes = 50 * 1024 * 1024),
    project = TRUE
  ) %>%
  # Biomes polygons
  leaflet::addPolygons(
    data = biomes_sf,
    fillColor = ~pal_biome(Biome),
    fillOpacity = 0.6,
    color = "black",
    weight = 0.1,
    popup = ~paste0(
        paste0(
          "<b>ID:</b> ", ID, "<br>",
          "<b>Biome:</b> ", Biome, "<br>",
          "<b>Area:</b> ", round(total_area_km2, 0), " km² <br>",
          "<b>Est. spec. richness:</b> ", species_richness
          )
        ),
    labelOptions = labelOptions(
      direction = "auto",
      sticky = TRUE,
      opacity = 0.8
    ),
    group = "Biomes"
  ) %>%
  leaflet::addLayersControl(
    overlayGroups = c("Native", "Introduced", "Biomes"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  ) %>%
  leaflet::addLegend(
    "topleft",
    pal = pal_biome,
    values = biomes_sf$Biome,
    title = "Biome",
    group = "Biomes"
  ) %>%
  htmlwidgets::saveWidget(
    file = "map.html",
    selfcontained = TRUE
  )
