library(dplyr)
library(terra)
library(tidyterra)
library(sf)
library(igraph)

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

if (Sys.info()["sysname"] == "Windows") {
  lud11 <- "L:"
} else {
  lud11 <- "/lud11"
}

dir_fig <- "C:/Users/poppman/Dropbox/Apps/Overleaf/BiomeInvasibility/suppl_files"

biomes <- terra::vect(
  file.path(lud11, "poppman/data/bir/dat/lud11/shp/olson_ecoregions/biomes.shp")
)
biomes$ID <- 1:nrow(biomes)

df <- data.frame()
for (biome_id in sort(unique(biomes$BIOME))) {
  biomes_sub <- dplyr::filter(biomes, BIOME == biome_id)
  
  dist_mat <- biomes_sub %>%
    sf::st_as_sf() %>%
    sf::st_distance() %>%
    as.matrix() %>%
    units::drop_units()
  
  # Get areas
  areas <- terra::expanse(biomes_sub)
  
  # Define connectivity (exponential decay)
  alpha <- 0.67 / 1e6
  p_matrix <- exp(-dist_mat * alpha)
  
  # Assign clusters based on a distance threshold
  for (f in c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)) {
    max_dist <- log(f) / (-alpha)
    adj <- dist_mat <= max_dist
    diag(adj) <- FALSE
    
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    clusters <- igraph::components(g)$membership
    df <- rbind(
      df,
      data.frame(
        biome_id = biome_id,
        f = f,
        num_clusters = max(clusters),
        num_patches = nrow(biomes_sub)
      )
    )
  }
}

dfp <- df %>%
  dplyr::filter(biome_id <= 14) %>%
  dplyr::mutate(
    biome = biome_names[biome_id]
  )

dfp$Biome <- factor(dfp$biome, levels = biome_names)

gg <- ggplot2::ggplot(
  data = dfp,
  ggplot2::aes(
    x = f,
    y = num_clusters,
    colour = Biome,
    linetype = Biome,
    shape = Biome
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::xlab("Threshold connectivity fraction") +
  ggplot2::ylab("Number of resulting clusters") +
  ggplot2::theme_bw() +
  ggplot2::scale_linetype_manual(values = 1:nlevels(dfp$Biome)) +
  ggplot2::scale_shape_manual(values = 1:nlevels(dfp$Biome))

ggplot2::ggsave(
  filename = file.path(dir_fig, "Threshold_selection.pdf"),
  plot = gg,
  width = 8,
  height = 5
)
