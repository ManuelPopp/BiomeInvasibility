#!/usr/bin/env Rscript
#>------------------------------------<
##
## Script name: Download species
## occurence data
##
## Author: Manuel R. Popp
## Email: manuel.popp@wsl.ch
##
## Date Created: 2024-12-07
##
## ---------------------------
##
## Descripton: Download species observations from GBIF
## Notes: -
##
#>----------------------------------------------------------------------------<|
#> Install/load packages
rm(list = ls())
import <- function(...) {
  #' Import R packages. Install them if necessary.
  #'
  #' @param ... any argument that can be passed to install.packages.
  #' @details The function installs only packages that are missing. Packages
  #' are loaded.
  #' @examples
  #' # Load packages
  #' import("dplyr", "MASS", "terra", dependencies = TRUE)
  #'
  #' @seealso \code{\link[base]{install.packages}}
  #' @export
  args <- list(...)
  packages = args[names(args) == ""]
  kwargs = args[names(args) != ""]
  
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      do.call(install.packages, c(list(package), kwargs))
    }
    require(package, character.only = TRUE)
  }
}

import(
  "dplyr", "devtools", "stringdist", "terra", "tidyterra", "sf",
  dependencies = TRUE
)

if (!require("GBIFhandleR", character.only = TRUE)) {
  devtools::install_github("https://github.com/ManuelPopp/GBIFhandleR")
}

require("GBIFhandleR")

#>----------------------------------------------------------------------------<|
#> Settings
if (Sys.info()["sysname"] == "Windows") {
  dir_main <- "C:/Users/poppman/switchdrive/PhD/prj/bir"
} else {
  dir_main <- "/lud11/poppman/BiomeInvasibility"
}

dir_dat <- file.path(dir_main, "dat")
dir_out <- file.path(dir_main, "out")
dir_lud <- file.path(dir_dat, "lud11")
dir_cln <- file.path(dir_lud, "cln")

f_species <- file.path(dir_lud, "db", "invasive_db_2013.csv")
f_biomes <- file.path(dir_lud, "shp", "olson_ecoregions", "wwf_terr_ecos.shp")
f_countries <- file.path(
  dir_lud, "shp", "world-administrative-boundaries.geojson"
  )

f_areas <- file.path(dir_out, "Areas.csv")

# Set up directory for output
dir.create(dir_cln, showWarnings = FALSE)
dir.create(dir_out, showWarnings = FALSE)

if (!file.exists(file.path(dir_main, ".gitignore.txt"))) {
  sink(file.path(dir_main, ".gitignore.txt"))
  print("dat/obs")
  sink()
}

place_names_abbr <- c(
  "NAm" = "Northern America",
  "Af" = "Africa",
  "ME" = "Middle East",
  "IO" = "Indian Ocean",
  "As" = "Asia",
  "PI" = "Pacific Islands",
  "Au" = "Australia",
  "In" = "India",
  "Papua New Guinea" = "Papua New Guinea",
  "Malesia" = "Malaysia",
  "CAm" = "Central America",
  "SAm" = "South America",
  "Venezuela" = "Venezuela",
  "SE Nam" = "Namibia",
  "SAf" = "South Africa",
  "Tasmania" = "Australia and New Zealand",
  "NSW" = "Australia",
  "Florida" = "United States of America",
  "Car" = "Caribbean",
  "China" = "China",
  "Jap" = "Japan",
  "Eu" = "Europe",
  "W As" = "Western Asia",
  "Japan" = "Japan",
  "E NAm" = "Eastern Asia",
  "Taiwan" = "Taiwan",
  "SE As" = "South-Eastern Asia",
  "Iran" = "Iran (Islamic Republic of)",
  "Brazil" = "Brazil",
  "W Af" = "Western Africa",
  "Philipp" = "Philippines",
  "N Au" = "Northern Australia",
  "S Eu" = "Southern Europe",
  "E As" = "Eastern Asia",
  "Russia" = "Russian Federation",
  "NZ" = "New Zealand",
  "S As" = "Southern Asia",
  "northern Mexico" = "Mexico",
  "Madagascar" = "Madagascar",
  "India" = "India",
  "Mexico" = "Mexico",
  "Sam" = "Samoa",
  "N Af" = "Northern Africa",
  "Peru" = "Peru",
  "SE Brazil" = "Brazil",
  "argentina" = "Argentina",
  "old intro" = "Old Introductions",
  "SE US" = "United States of America",
  "Sri Lanka" = "Sri Lanka",
  "AI" = "Atlantic Islands",
  "Canary Is" = "Canary Islands",
  "Panama" = "Panama",
  "Bangladesh" = "Bangladesh",
  "Korea" = "Republic of Korea",
  "Anaman" = "Andaman Islands",
  "NG" = "Papua New Guinea",
  "Amapa" = "Brazil",
  "Maranhao" = "Brazil",
  "Para" = "Brazil",
  "Nam" = "Namibia",
  "Ethiopia" = "Ethiopia",
  "Europe" = "Europe",
  "N As" = "Northern Asia",
  "E Siberia" = "Russian Federation",
  "E Af" = "Eastern Africa",
  "Iberian" = "Portugal",
  "Queensland" = "Australia",
  "Canary" = "Canary Islands",
  "SW Eur" = "Southern Europe",
  "SAmer" = "South America",
  "Philippines" = "Philippines",
  "hybrid Eu" = "Europe",
  "New Guinea" = "Papua New Guinea",
  "Palau" = "Palau",
  "Papua" = "Papua New Guinea",
  "Asia" = "Asia",
  "Vietnam" = "Vietnam",
  "Solomon Is" = "Solomon Islands",
  "Belize" = "Belize",
  "Guatemala" = "Guatemala",
  "Turkey" = "Turkey",
  "bot gard" = "Botanical Gardens",
  "Fiji" = "Fiji",
  "Paraguay" = "Paraguay",
  "Argentina" = "Argentina",
  "Comoros" = "Comoros",
  "Malaya" = "Malaysia",
  "New Caledonia and Loyalty Islands" = "New Caledonia",
  "Cam" = "Cameroon",
  "SW Au" = "Southwestern Australia",
  "Canary Islands" = "Canary Islands",
  "Persian Gulf" = "Western Asia",
  "up to Lena River" = "Russian Federation",
  "Java" = "Indonesia",
  "Sumatra" = "Indonesia",
  "N Mexico" = "Mexico",
  "Guadalupe" = "Guadeloupe",
  "Spain" = "Spain",
  "Hainan" = "China",
  "S Japan" = "Japan",
  "SE China" = "China",
  "SW Eu" = "Southern Europe",
  "SW As" = "Southwestern Asia",
  "N Amer" = "Northern America",
  "GB" = "U.K. of Great Britain and Northern Ireland",
  "Eur" = "Europe",
  "England" = "U.K. of Great Britain and Northern Ireland",
  "France" = "France",
  "Ireland" = "Ireland",
  "Iran" = "Iran (Islamic Republic of)",
  "W Eu" = "Western Europe",
  "W" = "West",
  "centr Eu" = "Central Europe",
  "Philopp" = "Philippines",
  "N" = "North",
  "Am" = "Americas",
  "Britain" = "U.K. of Great Britain and Northern Ireland",
  "E Au" = "Eastern Australia",
  "W NAm" = "Western Northern America",
  "Germany" = "Germany",
  "Hawaii" = "United States of America",
  "New Caledonia" = "New Caledonia",
  "N SAm" = "Northern South America",
  "Bolivia" = "Bolivia",
  "W Au" = "Western Australia",
  "cult Eu" = "Europe",
  "SE Eu" = "Southeastern Europe",
  "Syria" = "Syrian Arab Republic",
  "Lesser Antiles" = "Caribbean",
  "Timor" = "Timor-Leste",
  "PO" = "Pacific Ocean",
  "Mauritius" = "Mauritius",
  "S SAm" = "Southern South America",
  "NAf" = "Northern Africa",
  "Baja CA" = "Mexico",
  "NAm" = "Northern America",
  "Af" = "Africa",
  "Au" = "Australia"
)

#>----------------------------------------------------------------------------<|
#> Functions
download <- function(species_name) {
  dst <- file.path(dir_cln, paste0(species_name, ".csv"))
  
  if (!file.exists(dst)) {
    observations <- GBIFhandleR::get_observations(
      species_name,
      basisOfRecord = c("OBSERVATION", "HUMAN_OBSERVATION")
    )
    write.csv(observations$data, file = dst, row.names = FALSE)
  }
}

get_origin <- function(
    countries = NULL, regions = NULL, continents = NULL,
    admin_boundaries = f_countries
    ) {
  if (is.null(countries)) {
    countries <- c()
  }
  if (is.null(regions)) {
    regions <- c()
  }
  if (is.null(continents)) {
    continents <- c()
  }
  
  bounds <- terra::vect(f_countries) %>%
    dplyr::filter(
      name %in% countries | region %in% regions | continent %in% continents
    ) %>%
    terra::aggregate()
  
  return(bounds)
}

translate_origin <- function(strings, places_df, place_names_abbr) {
  places_df_names <- sort(as.character(na.omit(places_df$name)))
  names(places_df_names) <- places_df_names
  names_abbr <- c(
    place_names_abbr, places_df_names
  )
  cleaned_strings <- sub("India(old intro.)", "", strings)
  cleaned_strings <- sub(
    "ME",
    paste0(
      "Bahrain,Cyprus,Egypt,Iran,Iraq,Israel,Jordan,Kuwait,Lebanon,Oman,",
      "Qatar,Saudi Arabia,Syria,Turkey,United Arab Emirates,Yemen"
      ),
    strings
    )
  cleaned_strings <- gsub("[^[:alpha:][:space:]]", ",", cleaned_strings)
  cleaned_strings <- gsub("\\s+", " ", cleaned_strings)
  cleaned_strings <- strsplit(cleaned_strings, ",")
  cleaned_strings <- lapply(cleaned_strings, FUN = trimws)
  cleaned_strings <- lapply(
    cleaned_strings, FUN = function(x){return(x[which(x != "")])}
    )
  place_names <- as.character(names_abbr[unlist(cleaned_strings)])
  
  return(place_names)
}

place_by_type <- function(cleaned_strings, places_df) {
  out_list = list(
    "countries" = c(),
    "regions" = c(),
    "continents" = c()
  )
  
  for (string in cleaned_strings) {
    idx <- which(places_df$name == string)
    if (length(idx) < 1) {
      stop(paste0("\nNo match for place name ", string))
    }
    if (places_df$type[idx] == "country") {
      out_list[["countries"]] <- c(out_list[["countries"]], string)
    } else if (places_df$type[idx] == "region") {
      out_list[["regions"]] <- c(out_list[["regions"]], string)
    } else {
      out_list[["continents"]] <- c(out_list[["continents"]], string)
    }
  }
  return(out_list)
}

get_region_wrapper <- function(
    strings, places_df, place_names_abbr, admin_boundaries = f_countries
) {
  cleaned <- translate_origin(strings, places_df, place_names_abbr)
  args <- place_by_type(cleaned, places_df)
  origin <- get_origin(
    countries = args[["countries"]], regions = args[["regions"]],
    continents = args[["continents"]],
    admin_boundaries = f_countries
  )
  return(origin)
}

split_sp_name <- function(x) {
  parts <- strsplit(x, " ")[[1]]
  sp_name <- paste(parts[which(parts != "")][c(1, 2)], collapse = " ")
  sp_author <- paste(parts[which(parts != "")][-c(1, 2)], collapse = " ")
  return(c(sp_name, sp_author))
}

escape_special_chars <- function(string) {
  gsub("([\\(){}\\[\\]])", "\\\\\\1", string, perl = TRUE)
  #gsub("([\\(){}\\[\\]])", "", string, perl = TRUE)
}

#>----------------------------------------------------------------------------<|
#> Load biome and other data
# biomes <- sf::st_read(f_biomes, layer = "wwf_terr_ecos") %>%
#   sf::st_cast("POLYGON") %>%
#   dplyr::mutate(BIOME = factor(BIOME)) %>%
#   sf::st_make_valid() %>%
#   dplyr::group_by(BIOME) %>%
#   dplyr::summarise(geometry = sf::st_union(geometry)) %>%
#   sf::st_cast("POLYGON") %>%
#   terra::vect()

biomes <- terra::vect(f_biomes) %>%
  terra::aggregate(by = "BIOME") %>%
  terra::union() %>%
  terra::disagg()

countries <- terra::vect(f_countries)
places <- data.frame(
  name = c(
    sort(unique(countries$name)), sort(unique(countries$continent)),
    sort(unique(countries$region))
    ),
  type = c(
    rep("country", length(unique(countries$name))),
    rep("continent", length(unique(countries$continent))),
    rep("region", length(unique(na.omit(countries$region))))
  )
)

#>----------------------------------------------------------------------------<|
#> Read species list
species_data <- read.csv(f_species) %>%
  dplyr::mutate(Species_name = paste(Genus, Species)) %>%
  dplyr::filter(!Origin %in% c("Eu (bot gard)"))

species <- species_data %>%
  dplyr::pull("Species_name")

observation_data <- list.files(dir_cln, pattern = ".csv", full.names = TRUE)

# Loop through species
no_observations <- c()# Vector to fill with species names where observations are empty

row <- paste("Species", "Native_m2", "Invaded_m2", collapse = ",")
cat(row, "\n", file = f_areas, append = FALSE)

for (file in observation_data) {
  species <- sub(".csv", "", basename(file))
  # Create SpatVector from observations
  species_observations <- read.csv(file = file) %>%
    dplyr::select(
      scientificName, continent, decimalLongitude, decimalLatitude
      ) %>%
    terra::vect(
      geom = c("decimalLongitude", "decimalLatitude"), crs = "epsg:4326"
    )
  
  if (nrow(species_observations) < 1) {
    print(paste("No observations for", species))
    next
  }
  
  # Get native range
  species_name <- species_observations$scientificName[1]
  # The database is... shitty. Found a trailing whitespace in the first entry,
  # thus, we need to try fuzzy matching.
  idx <- stringdist::amatch(
    species_name, species_data$Species_name, maxDist = Inf
  )
  
  origin_codes <- species_data$Origin[idx]
  
  origin <- get_region_wrapper(
    strings = origin_codes, places_df = places,
    place_names_abbr = place_names_abbr, admin_boundaries = f_countries
  ) %>%
    terra::buffer(width = 10000) # Accuracy of some locations has errors > 1 km
  
  species_observations$native <- terra::relate(
    species_observations, origin, relation = "intersects"
  )
  
  native_range_obs <- species_observations %>%
    dplyr::filter(native) %>%
    terra::aggregate()
  
  exotic_range_obs <- species_observations %>%
    dplyr::filter(!native) %>%
    terra::aggregate()
  
  sp_name <- escape_special_chars(trimws(split_sp_name(species_name)[1]))
  sp_author <- escape_special_chars(trimws(split_sp_name(species_name)[2]))
  suffix <- escape_special_chars(
    trimws(paste("Origin:", paste(origin_codes, collapse = ", ")))
    )
  
  title_expr <- bquote(
    italic(.(sp_name)) ~ .(sp_author) ~ .(suffix)
  )
  plot(
    countries,
    col = "yellow3",
    main = paste(sp_name, sp_author, suffix)
    )
  plot(origin, add = TRUE, col = "orange3")
  plot(native_range_obs, add = TRUE, col = "green", cex = 0.3)
  plot(exotic_range_obs, add = TRUE, col = "red", cex = 0.3)
  
  native_biomes <- biomes[
    which(terra::relate(biomes, native_range_obs, relation = "intersects")),
  ]
  
  invaded_biomes <- biomes[
    which(terra::relate(biomes, exotic_range_obs, relation = "intersects")),
  ]
  
  plot(biomes, main = paste(sp_name, sp_author, suffix))
  plot(native_biomes, col = "orange3", add = TRUE)
  plot(invaded_biomes, col = "yellow3", add = TRUE)
  plot(native_range_obs, add = TRUE, col = "green", cex = 0.3)
  plot(exotic_range_obs, add = TRUE, col = "red", cex = 0.3)
  
  area_native <- sum(terra::expanse(native_biomes, unit = "m"))
  area_invaded <- terra::expanse(invaded_biomes, unit = "m")
  
  for (a in area_invaded) {
    row <- paste(species_name, area_native, area_invaded, collapse = ",")
    cat(row, "\n", file = f_areas, append = TRUE)
  }
  
  flush.console()
}
