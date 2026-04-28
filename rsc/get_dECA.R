library(dplyr)
library(terra)
library(sf)
library(tidyterra)
library(igraph)

args <- commandArgs(trailingOnly = TRUE)
biome_id <- as.numeric(args[1])

if (Sys.info()["sysname"] == "Windows") {
  lud11 <- "L:"
} else {
  lud11 <- "/lud11"
}

biome_def <- "olson"

dir_dat <- file.path(lud11, "poppman", "data", "bir", "dat", "lud11")
dir_imeb <- file.path(dir_dat, "biomes", biome_def, "intermediate_data")

if (biome_def == "olson") {
  biomes <- terra::vect(
    file.path(lud11, "poppman/data/bir/dat/lud11/shp/olson_ecoregions/biomes.shp")
  )
} else {
  stop("Biome input file not defined.")
}

biomes$ID <- 1:nrow(biomes)
biomes_sub <- dplyr::filter(biomes, BIOME == biome_id)

dist_mat <- biomes_sub %>%
  sf::st_as_sf() %>%
  sf::st_distance() %>%
  as.matrix() %>%
  units::drop_units()

# Get areas
areas <- terra::expanse(biomes_sub, unit = "m")

# Define connectivity (exponential decay)
alpha <- 0.67 / 1e6
p_matrix <- exp(-dist_mat * alpha)

# Calculate equivalent connected area (ECA) and dECA.
ECA_full <- sqrt(sum(outer(areas, areas) * p_matrix))

get_deca_i <- function(i) {
  areas_sub <- areas[-i]
  p_sub <- p_matrix[-i, -i]
  ECA_sub <- sqrt(sum(outer(areas_sub, areas_sub) * p_sub))
  
  return(ECA_full - ECA_sub)
}

dECA <- sapply(seq_len(length(areas)), get_deca_i)

# Compute connectedness of each patch (distance-weighted area of neighbour patches)
patch_connectedness <- sapply(
  seq_along(areas),
  function(i) {
    sum(areas[-i] * exp(-dist_mat[i, -i] * alpha))
    }
  )

# Compute connectedness with alternative values for alpha
alpha_low <- 0.5 / 1e6
patch_connectedness_low <- sapply(
  seq_along(areas),
  function(i) {
    sum(areas[-i] * exp(-dist_mat[i, -i] * alpha_low))
  }
)

alpha_high <- 1 / 1e6
patch_connectedness_high <- sapply(
  seq_along(areas),
  function(i) {
    sum(areas[-i] * exp(-dist_mat[i, -i] * alpha_high))
  }
)

alpha_max <- 1.3 / 1e6
patch_connectedness_highest <- sapply(
  seq_along(areas),
  function(i) {
    sum(areas[-i] * exp(-dist_mat[i, -i] * alpha_max))
  }
)

# Assign clusters based on a distance threshold
f <- 0.5
max_dist <- log(f) / (-alpha)
max_dist <- round(max_dist / 1e6, 0) * 1e6
adj <- dist_mat <= max_dist
diag(adj) <- FALSE

g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
clusters <- igraph::components(g)$membership

# Compute cluster-wise ECA
clECA <- numeric(length(clusters))
for (cl in unique(clusters)) {
  indices <- which(clusters == cl)
  areas_cl <- areas[indices]
  p_cl <- p_matrix[indices, indices]
  ECA_cl <- sqrt(sum(outer(areas_cl, areas_cl) * p_cl))
  clECA[indices] <- ECA_cl
}

dECA_df <- data.frame(
  ID = biomes_sub$ID,
  dECA = dECA,
  connectedness = patch_connectedness,
  connectedness_low_a = patch_connectedness_low,
  connectedness_high_a = patch_connectedness_high,
  connectedness_highest_a = patch_connectedness_highest,
  clusterECA = clECA,
  cluster = clusters
)

write.csv(
  dECA_df,
  file = file.path(
    dir_imeb, "dECA", paste0("biome", biome_id, ".csv")
    ),
  row.names = FALSE
)