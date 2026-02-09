#' Combine polygons from
#' World Geographical Scheme for Recording Plant Distributions (WGSRPD)
#'
#' Data source:
#' https://github.com/tdwg/wgsrpd/tree/master/geojson

require("sf")
require("terra")
require("tidyterra")
bbox_sf <- sf::st_as_sfc(
  sf::st_bbox(
    c(xmin = -180, ymin = -90, xmax = 180, ymax = 90), crs = 4326)
  )

lvl1 <- sf::st_read("C:/Users/poppman/Downloads/level1.gpkg") %>%
  dplyr::rename(
    code = LEVEL1_COD,
    name = LEVEL1_NAM
  ) %>%
  dplyr::mutate(
    code = name,
    level = 1
  ) %>%
  dplyr::select(
    code, name, level
  ) %>%
  sf::st_intersection(
    bbox_sf
  )

lvl2 <- sf::st_read("C:/Users/poppman/Downloads/level2.gpkg") %>%
  dplyr::rename(
    code = LEVEL2_COD,
    name = LEVEL2_NAM,
    l1c = LEVEL1_COD,
    l1n = LEVEL1_NAM
  ) %>%
  dplyr::mutate(
    code = name,
    level = 2
  ) %>%
  dplyr::select(
    code, name, level
  ) %>%
  sf::st_intersection(
    bbox_sf
  )

lvl3 <- sf::st_read("C:/Users/poppman/Downloads/level3.gpkg") %>%
  dplyr::rename(
    code = LEVEL3_COD,
    name = LEVEL3_NAM,
    l1c = LEVEL1_COD,
    l2c = LEVEL2_COD
  ) %>%
  dplyr::mutate(
    level = 3
  ) %>%
  dplyr::select(
    code, name, level
  ) %>%
  sf::st_intersection(
    bbox_sf
  )

lvl4 <- sf::st_read("C:/Users/poppman/Downloads/level4.gpkg") %>%
  dplyr::rename(
    code = Level4_cod,
    name = Level_4_Na,
    l1c = Level1_cod,
    l2c = Level2_cod,
    l3c = Level3_cod
  ) %>%
  dplyr::mutate(
    level = 4
  ) %>%
  dplyr::select(
    code, name, level
  ) %>%
  sf::st_intersection(
    bbox_sf
  )

tdwg <- rbind(lvl1, lvl2, lvl3, lvl4)
  
setwd("D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/git/GBIFhandleR")
usethis::use_data(tdwg, overwrite = TRUE)
