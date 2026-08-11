potential_predictors_h <- c(
  "logGDP",
  "logRoadDensity"
)

df_invasion_filtered_h <- df_invasion[ # Filter out incomplete cases
  stats::complete.cases(
    df_invasion[, potential_predictors_h]
  ),
] %>%
  dplyr::filter( # Filter out non-finite values
    dplyr::if_all(
      dplyr::all_of(potential_predictors_h),
      ~ is.finite(.)
    )
  ) %>%
  dplyr::mutate(
    Specieslvl = stringr::word(Species, 1, 2)
  ) %>%
  dplyr::filter(
    as.numeric(Biome) <= 12
  ) %>%
  dplyr::distinct( # Drop unclear subspecies
    dplyr::across(-Species),
    .keep_all = TRUE
  ) %>%
  dplyr::group_by(Species, Invaded) %>% # Limit to 1000 pseudo-absences per species
  dplyr::slice_min(order_by = Bhattacharyya, n = 1000) %>%
  dplyr::group_by(Species) %>% # Demand a min of 5 pseudo-absences per species
  dplyr::filter(sum(Invaded == 0, na.rm = TRUE) >= 5) %>%
  dplyr::ungroup()

num_samples_h <- df_invasion_filtered_h %>%
  dplyr::group_by(Species, Invaded) %>%
  dplyr::summarise(
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Invaded) %>%
  dplyr::summarise(
    min = min(N, na.rm = TRUE),
    mean = mean(N, na.rm = TRUE),
    max = max(N, na.rm = TRUE),
    median = median(N, na.rm = TRUE),
    sum = sum(N, na.rm = TRUE)
  )

mean_samples_h <- num_samples_h %>%
  dplyr::filter(Invaded == 1) %>%
  dplyr::pull(mean) %>%
  as.integer()

w_h <- num_samples_h %>%
  dplyr::select(Invaded, sum) %>%
  tidyr::pivot_wider(names_from = Invaded, values_from = sum) %>%
  dplyr::transmute(w = `0` / `1`) %>%
  dplyr::pull(w)

df_mod_h <- df_invasion_filtered_h %>%
  dplyr::mutate(
    weight = ifelse(Invaded == 1, w_h, 1),
  )

# Fit human impact model
for (p in potential_predictors_h) {
  frml_h <- make_formula(p, biome = TRUE)
  mod_gam_h <- mgcv::gam(
    frml_h,
    data = df_mod_h,
    method = "REML",
    family = stats::binomial("logit"),
    weights = weight
  )
  cat(
    "\nPredictor:", p,
    "\nadj. D2:", ecospat::ecospat.adj.D2.glm(mod_gam_h),
    "\n"
  )
}
