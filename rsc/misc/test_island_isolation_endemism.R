# Test whether island endemism increases with isolation to potentially establish
# a relationship to use as an estimator for a reasonable decay function
# parameterisation for ECA metrics

library("GIFT")
library("dplyr")
library("terra")
library("rnaturalearth")
library("rnaturalearthdata")

# Download countries at 110m resolution
countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Dissolve countries into continents
continents <- countries %>%
  dplyr::group_by(continent) %>%
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
  sf::st_cast("POLYGON") %>%
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry)) / 1e6) %>%
  dplyr::filter(area_km2 > 5e6) %>%
  sf::st_union()

# Get degree of endemism for islands
meta <- GIFT::GIFT_lists()
island_ids <- meta$entity_ID[
  meta$entity_class %in% c("Island", "Island Group")
  ]

richness_df <- GIFT::GIFT_richness(taxon_name = "Tracheophyta") %>%
  dplyr::filter(entity_ID %in% island_ids) %>%
  dplyr::mutate(frac_endemic = endemic_min / native) %>%
  dplyr::filter(
    !is.na(frac_endemic)
    )

island_shapes <- GIFT::GIFT_shapes(entity_ID = richness_df$entity_ID) %>%
  dplyr::arrange(match(entity_ID, richness_df$entity_ID))

# Add island distance to mainland
island_shapes$dist_to_mainland <- sf::st_distance(island_shapes, continents)
island_shapes$area <- as.numeric(sf::st_area(island_shapes))

df <- merge(island_shapes, richness_df, by = "entity_ID") %>%
  sf::st_drop_geometry()

df$area_quantile <- cut(
  df$area,
  breaks = quantile(df$area, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

ggplot2::ggplot(
  data = df[which(df$frac_endemic > 0) & df$area > median(df$area),],
  ggplot2::aes(x = as.numeric(dist_to_mainland), y = frac_endemic, colour = factor(area_quantile))
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm")


mod <- glm(frac_endemic ~ dist_to_mainland + area, data = df)
summary(mod)
plot(frac_endemic ~ log(area), data = df[which(df$frac_endemic > 0),])
