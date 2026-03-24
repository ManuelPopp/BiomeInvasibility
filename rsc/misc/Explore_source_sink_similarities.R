df_species_patches
names(merged)
df <- merged %>%
  dplyr::group_by(Species, ID) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 1) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::summarise(
    tasmean = mean(tasmean, na.rm = TRUE),
    tasmean_sd = mean(tasmean_sd, na.rm = TRUE),
    prec = mean(prec, na.rm = TRUE),
    prec_sd = mean(prec_sd, na.rm = TRUE),
    npp = mean(npp, na.rm = TRUE),
    ghm = mean(ghm, na.rm = TRUE),
    mean_species_richness = mean(speciesRichnessBa, na.rm = TRUE),
    max_species_richness = max(speciesRichnessBa, na.rm = TRUE),
    total_area_km2 = sum(total_area_km2, na.rm = TRUE),
    mean_total_area_km2 = mean(total_area_km2, na.rm = TRUE),
    max_total_area_km2 = max(total_area_km2, na.rm = TRUE),
    lowland_area_km2 = sum(lowland_area_km2, na.rm = TRUE),
    mean_lowland_area_km2 = mean(lowland_area_km2, na.rm = TRUE),
    max_lowland_area_km2 = max(lowland_area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = Species,
    names_from = Status,
    values_from = c(
      tasmean, tasmean_sd, prec, prec_sd, npp, ghm,
      mean_species_richness, max_species_richness,
      total_area_km2, mean_total_area_km2, max_total_area_km2,
      lowland_area_km2, mean_lowland_area_km2, max_lowland_area_km2
      )
  )

svg(
  filename = file.path(dir_main, "fig", "Predictor_overview.svg"),
  width = 9, height = 6
  )
par(mfrow = c(2, 3))
plot(
  tasmean_Receiver ~ tasmean_Donor,
  data = df,
  main = "Mean annual temperature",
  xlab = "Invaded range", ylab = "Native range",
  col = rgb(201/255, 74/255, 3/255, 0.25), pch = 16
  )
abline(a = 0, b = 1)
plot(
  prec_Receiver ~ prec_Donor,
  data = df,
  main = "Total annual precipitation",
  xlab = "Invaded range", ylab = "Native range",
  col = rgb(70/255, 100/255, 170/255, 0.25), pch = 16
  )
abline(a = 0, b = 1)
plot(
  npp_Receiver ~ npp_Donor,
  data = df,
  main = "NPP (MODIS, mean 2002–2018)",
  xlab = "Invaded range", ylab = "Native range",
  col = rgb(0/255, 201/255, 25/255, 0.25), pch = 16
  )
abline(a = 0, b = 1)
plot(
  sqrt(mean_species_richness_Receiver) ~ sqrt(mean_species_richness_Donor),
  data = df,
  main = "Mean species richness (samples = species)",
  xlab = expression(sqrt("Donor patches")), ylab = expression(sqrt("Receiver patches")),
  col = rgb(0/255, 102/255, 102/255, 0.25), pch = 16
)
abline(a = 0, b = 1)
plot(
  sqrt(max_species_richness_Receiver) ~ sqrt(max_species_richness_Donor),
  data = df,
  main = "Maximum species richness (samples = species)",
  xlab = expression(sqrt("Donor patches")), ylab = expression(sqrt("Receiver patches")),
  col = rgb(0/255, 102/255, 102/255, 0.25), pch = 16
)
abline(a = 0, b = 1)
plot(
  ghm_Receiver ~ ghm_Donor,
  data = df,
  main = "Human modification index",
  xlab = "Invaded range", ylab = "Native range",
  col = rgb(181/255, 145/255, 2/255, 0.25), pch = 16
)
abline(a = 0, b = 1)
dev.off()

df2 <- merged %>%
  dplyr::group_by(Species, ID) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 1) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::summarise(Biome = collapse::fmode(Biome), .groups = "drop") %>%
  dplyr::mutate(Biome_num = as.numeric(Biome)) %>%
  tidyr::pivot_wider(
    id_cols = Species,
    names_from = Status,
    values_from = c(Biome, Biome_num)
  ) %>%
  dplyr::count(Biome_Donor, Biome_Receiver, name = "freq") %>%
  dplyr::mutate(
    Origin = Biome_Donor,
    Invaded = Biome_Receiver
  )

df3 <- merged %>%
  dplyr::group_by(Species, ID) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 1) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(dplyr::n_distinct(Status) == 2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Species, Status, Biome) %>%
  dplyr::summarise(sum_count = sum(Count, na.rm = TRUE)) %>%
  dplyr::filter(sum_count == max(sum_count)) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Biome_num = as.numeric(Biome)) %>%
  tidyr::pivot_wider(
    id_cols = Species,
    names_from = Status,
    values_from = c(Biome, Biome_num)
  ) %>%
  dplyr::count(Biome_Donor, Biome_Receiver, name = "freq") %>%
  dplyr::mutate(
    Origin = Biome_Donor,
    Invaded = Biome_Receiver
  )

gg_sankey <- ggplot2::ggplot(
  df3,
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
  ggalluvial::stat_stratum(
    ggplot2::aes(label = ggplot2::after_stat(stratum)),
    geom = "label",
    fill = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_x_discrete(limits = c("Origin", "Invaded")) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "none"
  )

ggplot2::ggsave(
  filename = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/Sankey.svg",
  plot = gg_sankey,
  width = 10,
  height = 15
)














# 1. Summarize donor values per species
donor_summary <- merged %>%
  dplyr::filter(Status == "Donor") %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    tasmean_donor = mean(tasmean, na.rm = TRUE),
    tasmean_sd_donor = mean(tasmean_sd, na.rm = TRUE),
    prec_donor = mean(prec, na.rm = TRUE),
    prec_sd_donor = mean(prec_sd, na.rm = TRUE),
    npp_donor = mean(npp, na.rm = TRUE),
    ghm_donor = mean(ghm, na.rm = TRUE),
    mean_species_richness_donor = mean(speciesRichnessBa, na.rm = TRUE),
    max_species_richness_donor = max(speciesRichnessBa, na.rm = TRUE),
    total_area_km2_donor = sum(total_area_km2, na.rm = TRUE),
    lowland_area_km2_donor = sum(lowland_area_km2, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Create data frame for model
use_vars <- c(
  "tasmean", "prec", "ghm", "speciesRichnessBa", "lowland_area_km2",
  "tasmean_donor", "prec_donor", "ghm_donor", "max_species_richness_donor", "lowland_area_km2_donor"
  )
df_receivers <- merged %>%
  dplyr::left_join(donor_summary, by = "Species") %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(use_vars), ~ !is.na(.))) %>%
  dplyr::filter(lowland_area_km2 > 1, lowland_area_km2_donor > 1) %>%
  dplyr::filter(Status == "Receiver",) %>%
  dplyr::mutate(
    delta_tas = tasmean - tasmean_donor,
    delta_prec = prec - prec_donor,
    delta_ghm = ghm - ghm_donor,
    ratio_species_richness = speciesRichnessBa / max_species_richness_donor,
    ratio_area = lowland_area_km2 / lowland_area_km2_donor,
    received = TRUE
  ) %>%
  dplyr::filter()

df_null <- merged %>%
  dplyr::left_join(donor_summary, by = "Species") %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(use_vars), ~ !is.na(.))) %>%
  dplyr::filter(lowland_area_km2 > 1, lowland_area_km2_donor > 1) %>%
  dplyr::mutate(
    tasmean = sample(tasmean),
    prec = sample(prec),
    ghm = sample(ghm),
    speciesRichnessBa = sample(speciesRichnessBa),
    lowland_area_km2 = sample(lowland_area_km2)
    ) %>%
  dplyr::filter(Status == "Receiver") %>%
  dplyr::mutate(
    delta_tas = tasmean - tasmean_donor,
    delta_prec = prec - prec_donor,
    delta_ghm = ghm - ghm_donor,
    ratio_species_richness = speciesRichnessBa / max_species_richness_donor,
    ratio_area = lowland_area_km2 / lowland_area_km2_donor,
    received = FALSE
  )

df_test <- rbind(df_receivers, df_null)
df_clean <- df_test[complete.cases(df_test[, use_vars]), ]

model <- lme4::glmer(
  received ~ delta_tas + delta_prec + log(ratio_species_richness) + log(ratio_area) + (1 | Biome),
  data = df_clean,
  family = stats::binomial()
  )

summary(model)
r2 <- performance::r2(model)
r2


mod_richness <- lmer(
  log(SpeciesRichness) ~ log(Area) + Productivity + Disturbance + (1 | Biome),
  data = df
)

mod_donor <- glmmTMB(
  cbind(DonorOf, failures) ~ log(Area) + log(SpeciesRichness) + Productivity + (1 | Biome),
  family = betabinomial(),
  data = df
)

mod_receiver <- glmmTMB(
  InvasiveCount ~ log(SpeciesRichness) + Disturbance + HumanModification + (1 | Biome),
  family = poisson(),
  data = df
)

psem(mod_richness, mod_donor, mod_receiver)