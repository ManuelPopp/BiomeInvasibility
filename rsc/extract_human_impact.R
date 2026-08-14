#>=============================================================================<
#>-----------------------------------------------------------------------------<
#> Extract patch-level information from raster data
#<=============================================================================>

print(Sys.time())
library("dplyr")
library("terra")
library("tidyterra")


#>=============================================================================<
#> Settings
#<=============================================================================>
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
dir_imeb <- file.path(dir_dat, "biomes", biome_def, "intermediate_data")

f_bfullinfo <- file.path(dir_imeb, "biomes_full_info.gpkg")
f_road_density <- file.path(dir_dat, "grip_road_density", "GRIP4_tp1plus2.tif")
f_gdp_rst <- file.path(
  dir_dat, "gdp_data", "gdp_1990to2020_mean.tif"
  )

f_road_dns <- file.path(
  dir_imeb, "road_density.csv"
)

f_gdp <- file.path(
  dir_imeb, "gdp_1990to2020.csv"
)

#>=============================================================================<
#> Extract
#<=============================================================================>
road_density <- terra::rast(f_road_density) %>%
  terra::project("epsg:6933", method = "bilinear")
names(road_density) <- "RoadDensity"

gdp <- terra::rast(f_gdp_rst) %>%
  terra::project("epsg:6933")
names(gdp) <- "gdp1990to2020"

biomes <- terra::vect(f_bfullinfo) %>%
  terra::project("epsg:6933", method = "bilinear")

terra::extract(
  road_density,
  biomes,
  fun = mean,
  na.rm = TRUE
) %>%
  dplyr::mutate(
    RoadDensity = base::replace(RoadDensity, is.na(RoadDensity), 0)
    ) %>%
  utils::write.csv(
    file = f_road_dns,
    row.names = FALSE
  )

terra::extract(
  gdp,
  biomes,
  fun = mean,
  na.rm = TRUE
) %>%
  dplyr::mutate(
    gdp1990to2020 = base::replace(gdp1990to2020, is.na(gdp1990to2020), 0)
  ) %>%
  utils::write.csv(
    file = f_gdp,
    row.names = FALSE
  )