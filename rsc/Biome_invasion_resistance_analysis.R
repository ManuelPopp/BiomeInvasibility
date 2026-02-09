#------------------------------------------------------------------------------|
# Import packages
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

import("httr", "jsonlite", "sf", "dplyr")

#------------------------------------------------------------------------------|
# Settings
dir_prj <- "C:/Users/poppman/switchdrive/PhD/prj/bir"
dir_dat_large <- "L:/poppman/shared/biome_invasion/dat"

f_olson_biomes <- file.path(dir_dat_large, "Olson_Biomes_TEOW")
olson_biomes <- sf::st_read(f_olson_biomes, layer = "wwf_terr_ecos") %>%
  sf::st_cast("POLYGON") %>%
  dplyr::mutate(BIOME = factor(BIOME)) %>%
  sf::st_make_valid() %>%
  dplyr::group_by(BIOME) %>%
  dplyr::summarise(geometry = sf::st_union(geometry))

plot(olson_biomes[, c("BIOME", "geometry")])


species_name <- "Quercus robur"  # Replace with your species of interest
species_name_encoded <- URLencode(species_name)
# Build the search URL
search_url <- paste0("https://powo.science.kew.org/api/2/search?q=", species_name_encoded)

# Make the API request
response <- httr::GET(search_url)

# Check the response status
if (httr::status_code(response) != 200) {
  stop("API request failed with status: ", httr::status_code(response))
}

# Parse the JSON content
content <- httr::content(response, as = "text", encoding = "UTF-8")
json_data <- jsonlite::fromJSON(content)
if (json_data$totalResults == 0) {
  stop("No results found for species: ", species_name)
}

# Build the taxon URL
taxon_url <- paste0("https://powo.science.kew.org/api/2/", json_data$results$url[1])

# Make the API request
taxon_response <- httr::GET(taxon_url)

# Check the response status
if (httr::status_code(taxon_response) != 200) {
  stop("Failed to retrieve taxon data with status: ", httr::status_code(taxon_response))
}

# Parse the JSON content
taxon_content <- httr::content(taxon_response, as = "text", encoding = "UTF-8")
taxon_json <- jsonlite::fromJSON(taxon_content)
taxon_json$locations
