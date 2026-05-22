df_sPlotBiome <- sPlotVect %>%
  dplyr::filter(
    !is.na(BiomePatchID),
    #Naturalness == "Natural" # Current sPlot version: Mostly NA, naturalness encoded as integers
  ) %>%
  dplyr::group_by(BiomePatchID, Grassland, Shrubland, Forest, Wetland) %>%
  dplyr::summarise(
    Native_abundande = stats::weighted.mean(
      Relative_abundance_native,
      w = Releve_area,
      na.rm = TRUE
    ),
    Introduced_abundance = stats::weighted.mean(
      Relative_abundance_introduced,
      w = Releve_area,
      na.rm = TRUE
    ),
    Bare_soil_cover = stats::weighted.mean(
      Cover_bare_soil,
      w = Releve_area,
      na.rm = TRUE
    ),
    HumanModification = stats::weighted.mean(
      HumanModification,
      w = Releve_area,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    biomes %>%
      as.data.frame() %>%
      dplyr::mutate(BiomePatchID = ID),
    by = "BiomePatchID"
  ) %>%
  dplyr::mutate(
    Biome = factor(BIOME, levels = 1:length(biome_names), labels = biome_names)
  ) %>%
  dplyr::filter(
    as.numeric(Biome) <= 12,
    speciesRichnessBa > 30
  ) %>%
  dplyr::mutate(
    log_introducedAbundance = log1p(Introduced_abundance),
    log_speciesRichnessBa = log(speciesRichnessBa),
    log_sampling_effort = log(sampling_effort),
    log_Connectedness = log(connectedness),
    log_focalECA = log(focalECA),
    Biome = base::droplevels(Biome),
    BiomeOK = BIOME %in% c(1:6, 12) & Forest | BIOME %in% 7:11 & Grassland | BIOME %in% c(9, 11, 14) & Wetland | BIOME %in% c(8, 9, 12, 13) & Shrubland
  ) %>%
  as.data.frame() %>%
  dplyr::filter(
    BiomeOK
  )