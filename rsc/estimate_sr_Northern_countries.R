library("dplyr")
library("tidyr")
library("terra")
library("nanoparquet")
library("GIFT")
library("ggplot2")
library("nortest")
library("rnaturalearth")

# Set file paths
f_template <- "/lud11/nobis/Manuel/gbif/data/map_template.tif"
f_taxa <- "/lud11/nobis/Manuel/gbif/data/map_taxa.parquet"
f_cells <- "/lud11/nobis/Manuel/gbif/data/map_cells.parquet"
dir_out <- "/lud11/poppman/tmp"

# Load taxa and frequency table
taxa <- nanoparquet::read_parquet(f_taxa)
freq_tab <- data.frame(nanoparquet::read_parquet(f_cells))
colnames(freq_tab) <- c("row","specID", "cellID", "freq")
freq_tab$cellID <- as.numeric(freq_tab$cellID)

# Load GIFT countries and rasterise them to get country IDs per cell
gift_shapes <- GIFT::GIFT_shapes()
countries <- gift_shapes %>%
  dplyr::filter(
    entity_type == "Botanical Country",
    y_max > 30, # Avoid most places with high uncertainties
    x_min > -180 # Remove broken shapes
  ) %>%
  terra::vect()

countries$ID <- 1:nrow(countries)
map <- terra::rast(f_template)
countries_rst <- terra::rasterize(countries, map, field = "ID")

# Join country IDs to frequency table
country_ids <- terra::as.data.frame(countries_rst, xy = FALSE, cells = TRUE)
colnames(country_ids) <- c("cellID", "countryID")
freq_tab <- base::merge(freq_tab, country_ids, by = "cellID")

# Estimate species richness per country
# Estimate species richness per country
## Helper function to safely convert to numeric
as.num <- function(x) {
  x <- as.numeric(x)
  if (is.null(x)) {
    return(c(NA))
  }
  if (length(x) == 0) {
    return(c(NA))
  }
  if (is.nan(x)) {
    return(c(NA))
  }
  return(x)
}

## Function to estimate species richness for a given country ID
get_species_richness <- function(country_id) {
  sub_tab <- freq_tab[freq_tab$countryID == country_id, ]

  if (nrow(sub_tab) == 0) {
    return(
        data.frame(
            ID = country_id,
            speciesRichnessBa = NA,
            speciesRichnessBaSE = NA,
            speciesRichnessChao1 = NA,
            speciesRichnessChao1SE = NA,
            speciesRichnessACE = NA,
            speciesRichnessACESE = NA,
            speciesRichnessChao2 = NA
            )
        )
    }

  obs_counts <- sub_tab %>%
    dplyr::group_by(specID) %>%
    dplyr::summarise(counts = sum(freq)) %>%
    dplyr::pull(counts)
  
  baway <- tryCatch(
    breakaway::breakaway(obs_counts),
    error = function(e) list(est = NA, se = NA)
  )
  
  chao1 <- vegan::estimateR(obs_counts)
  
  inc_freq <- sub_tab %>%
    dplyr::distinct(cellID, specID) %>%
    dplyr::count(specID, name = "inc")
  Q1 <- sum(inc_freq$inc == 1)
  Q2 <- sum(inc_freq$inc == 2)
  n  <- dplyr::n_distinct(sub_tab$cellID)
  S_obs <- nrow(inc_freq)

  if (Q2 > 0) {
    chao2 <- S_obs + (Q1^2) / (2 * Q2)
  } else {
    chao2 <- S_obs + (Q1 * (Q1 - 1)) / 2
  }
  
  df <- data.frame(
    ID = as.num(country_id),
    speciesRichnessBa = as.num(baway$est),
    speciesRichnessBaSE = as.num(baway$se),
    speciesRichnessChao1 = as.num(chao1["S.chao1"]),
    speciesRichnessChao1SE = as.num(chao1["se.chao1"]),
    speciesRichnessACE = as.num(chao1["S.ACE"]),
    speciesRichnessACESE = as.num(chao1["se.ACE"]),
    speciesRichnessChao2 = as.num(chao2)
    )
  return(df)
}

country_ids_unique <- sort(unique(countries$ID))

results_df <- do.call(
    rbind,
    lapply(
      X = country_ids_unique,
      FUN = function(id) {
        tryCatch(
          get_species_richness(id),
          error = function(e) {
            stop(paste("Error in country ID", id, ":", e$message))
          }
          )
        }
      )
    )

# Add GIFT species counts
countries$GIFT <- NA
for (i in 1:nrow(countries)) {
  taxa <- GIFT::GIFT_checklists(
    taxon_name = "Tracheophyta",
    complete_taxon = FALSE,
    floristic_group = "native",
    shp = countries[i,] %>% sf::st_as_sf(),
    overlap = "centroid_inside",
    list_set_only = FALSE
  )
  countries$GIFT[i] <- length(unique(taxa$checklists$work_species))
}

save(results_df, file = file.path(dir_out, "species_richness_Northern.RData"))

countries_df <- as.data.frame(countries)
save(
  countries_df,
  file = file.path(dir_out, "GIFT_richness_Northern.RData")
)

spec_richness <- terra::merge(countries, results_df, by = "ID", all.x = TRUE)

# Test corellation between Chao2 estimated species richness and GIFT
corr <- lm(speciesRichnessChao2 ~ GIFT, data = spec_richness)
res <- residuals(corr)
spec_richness$residuals <- res
nortest::ad.test(res)

resmod <- lm(abs(residuals) ~ GIFT, data = spec_richness)
summary(resmod)

cor.test(
  spec_richness$speciesRichnessChao2,
  spec_richness$GIFT,
  method = "spearman"
)

# Add information on continents
continents <- rnaturalearth::ne_countries(
  scale = "medium", returnclass = "sf"
) %>%
  dplyr::select(continent) %>%
  terra::vect()

spec_richness$Continent <- terra::extract(
  continents,
  terra::centroids(spec_richness)
) %>%
  dplyr::pull(continent)

# Reshape into a long format for plotting
df_long <- spec_richness %>%
  as.data.frame() %>%
  dplyr::select(
    GIFT, starts_with("speciesRichness"), area, geo_entity, Continent
  ) %>%
  dplyr::select(-ends_with("SE")) %>%
  tidyr::pivot_longer(
    cols = -c(GIFT, area, geo_entity, Continent),
    names_to = "Method",
    values_to = "Estimate"
  ) %>%
  dplyr::mutate(
    Method = sub(pattern = "speciesRichness", replacement = "", Method)
  ) %>%
  dplyr::mutate(
    Method = sub(pattern = "Ba", replacement = "Breakaway", Method)
  ) %>%
  dplyr::rename(Area = area)

# Compute R2 per method
r2_tab <- df_long %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(GIFT), !is.na(Estimate)) %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    r2 = summary(lm(Estimate ~ GIFT))$r.squared,
    .groups = "drop"
  )

r2_tab_ctnt <- df_long %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(GIFT), !is.na(Estimate), !is.na(Continent)) %>%
  dplyr::group_by(Method, Continent) %>%
  dplyr::summarise(
    r2 = summary(lm(Estimate ~ GIFT))$r.squared,
    .groups = "drop"
  )

df_long <- df_long %>%
  dplyr::group_by(Method) %>%
  dplyr::mutate(
    Residual = resid(lm(Estimate ~ GIFT))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Label = ifelse(
      abs(Residual) > 0.5 * Estimate & Estimate > 1000,
      geo_entity,
      NA
    )
  )

lims <- range(c(df_long$GIFT, df_long$Estimate), na.rm = TRUE)

gg <- ggplot2::ggplot(
  df_long,
  ggplot2::aes(x = GIFT, y = Estimate, colour = Continent)
  ) +
  ggplot2::geom_point(alpha = 0.6, ggplot2::aes(size = Area)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::geom_smooth(
    method = "lm", se = FALSE, colour = rgb(0, 102, 102, maxColorValue = 255)
  ) +
  ggplot2::facet_wrap(~ Method, scales = "fixed") +
  ggplot2::geom_text(
    data = r2_tab,
    ggplot2::aes(
      x = -Inf,
      y = Inf,
      label = paste0("R² = ", round(r2, 2))
    ),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_fixed(ratio = 1, xlim = lims, ylim = lims) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "GIFT species richness",
    y = "Estimated species richness"
  )

ggplot2::ggsave(
  filename = file.path(dir_out, "species_richness_Northern_countries.png"),
  plot = gg,
  width = 10,
  height = 7
)

gg_label <- gg +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = Label),
    colour = "black", size = 2, hjust = 0, vjust = 0.5
  )

ggplot2::ggsave(
  filename = file.path(
    dir_out,
    "species_richness_Northern_countries_label.png"
  ),
  plot = gg_label,
  width = 10,
  height = 7
)

gg_res <- ggplot2::ggplot(
  df_long,
  ggplot2::aes(x = Estimate, y = Residual, colour = Continent)
  ) +
  ggplot2::geom_point(alpha = 0.6, ggplot2::aes(size = Area)) +
  ggplot2::geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  ggplot2::geom_smooth(
    method = "lm", se = FALSE, colour = rgb(0, 102, 102, maxColorValue = 255)
  ) +
  ggplot2::facet_wrap(~ Method, scales = "fixed") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "Estimated species richness",
    y = "Residuals"
  )

ggplot2::ggsave(
  filename = file.path(
    dir_out, "species_richness_Northern_countries_residuals.png"
    ),
  plot = gg_res,
  width = 10,
  height = 7
)
