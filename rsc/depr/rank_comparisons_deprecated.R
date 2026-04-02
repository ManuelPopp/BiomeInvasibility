f_b <- "L:/poppman/data/bff/dat/lud11/biomes/olson_ecoregions/wwf_terr_ecos.shp"
olson <- vect(f_b)

ctrs <- biomes_original %>%
  terra::centroids(inside = TRUE)
extr <- terra::extract(olson, ctrs) %>%
  dplyr::pull(ECO_NAME)

biomes_original$Name <- extr
biomes$Name <- extr

df_within <- merged %>%
  dplyr::mutate(
    PatchID = ID,
    Area = lowland_area,
    SpeciesRichness = round(speciesRichnessBa, 0)
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::filter(
    Count > 10,
    SpeciesRichness > 50
  ) %>%
  dplyr::group_by(BiomeID, Biome, clusterID, Status) %>%
  dplyr::summarise(
    Count = dplyr::n_distinct(Species),
    clusterRichness = max(speciesRichnessBa, na.rm = TRUE),
    localECA = max(localECA, na.rm = TRUE),
    clusterECA = max(clusterECA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(BiomeID, Biome, clusterID, clusterRichness, localECA, clusterECA),
    names_from = Status,
    values_from = Count
  )%>%
  dplyr::mutate(
    Donor = tidyr::replace_na(Donor, 0),
    Receiver = tidyr::replace_na(Receiver, 0),
    jDonor = jitter(tidyr::replace_na(Donor, 0), factor = 1e-5),
    jReceiver = jitter(tidyr::replace_na(Receiver, 0), factor = 1e-5)
  )

df <- df_within %>%
  dplyr::mutate(
    fracDonated = Donor / clusterRichness,
    fracReceived = Receiver / clusterRichness,
    rankDonor = dplyr::min_rank(jDonor),
    rankReceiver = dplyr::min_rank(jReceiver),
    data = "Observed"
  )

simulations <- list()
for (i in 1:500) {
  simulations[[i]] <- df_within %>%
    dplyr::mutate(
      jDonor = sample(jitter(tidyr::replace_na(Donor, 0), factor = 1e-5)),
      jReceiver = jitter(tidyr::replace_na(Receiver, 0), factor = 1e-5),
      fracDonated = (Donor + 1e-3) / clusterRichness,
      fracReceived = (Receiver + 1e-3) / clusterRichness,
      rankDonor = dplyr::min_rank(fracDonated),
      rankReceiver = dplyr::min_rank(fracReceived),
      data = "Observed"
    )
}

df2 <- do.call(
  rbind,
  simulations
)

library("MethComp")

reg <- MethComp::Deming(x = df$rankDonor, y = df$rankReceiver)
reg2 <- MethComp::Deming(x = df2$rankDonor, y = df2$rankReceiver)

ggplot2::ggplot(data = df, ggplot2::aes(x = rankReceiver, y = rankDonor)) +
  ggplot2::geom_point(colour = rgb(0, 102/255, 102/255), shape = 16, alpha = 0.2) +
  ggplot2::geom_abline(slope = reg[2], intercept = reg[1], colour = "lightblue", size = 2) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black", size = 1.2) +
  ggplot2::theme_bw() +
  ggplot2::xlab("Patch cluster rank as source") +
  ggplot2::ylab("Patch cluster rank as sink")

ids <- df %>%
  dplyr::filter(rankReceiver < 30 & rankDonor > 50) %>%
  dplyr::pull(clusterID) %>%
  as.character()

df[which(df$clusterID == "4-6-24691"),]

outliers <- biomes %>%
  dplyr::mutate(
    clusterID = paste(
      as.character(BIOME),
      as.character(clusterIDmaxdist),
      as.character(clusterIDbathy),
      sep = "-"
    )
    ) %>%
  dplyr::filter(
    clusterID %in% ids
  ) %>%
  dplyr::select(Name)

plot(outliers, col = "red")

dfb <- df_within %>%
  dplyr::group_by(BiomeID, Biome, clusterID, Status, clusterRichness, localECA, clusterECA) %>%
  dplyr::summarise(
    Count = dplyr::n_distinct(Species),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(BiomeID, Biome, clusterID, clusterRichness, localECA, clusterECA),
    names_from = Status,
    values_from = Count
  ) %>%
  dplyr::mutate(
    Donor = tidyr::replace_na(Donor, 0),
    Receiver = tidyr::replace_na(Receiver, 0),
    fracDonated = Donor + 1e-3 / clusterRichness,
    fracReceived = Receiver + 1e-3 / clusterRichness
  ) %>%
  dplyr::group_by(BiomeID) %>%
  dplyr::mutate(
    Biome = dplyr::first(Biome),
    clusterID = dplyr::first(clusterID),
    clusterRichness = dplyr::first(clusterRichness),
    localECA = dplyr::first(localECA),
    clusterECA = dplyr::first(clusterECA),
    rankDonor = dplyr::min_rank(fracDonated),
    rankReceiver = dplyr::min_rank(fracReceived)
  )

ggplot2::ggplot(data = dfb, ggplot2::aes(x = rankDonor, y = rankReceiver, colour = Biome)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggplot2::theme_bw() +
  ggplot2::xlab("Patch cluster rank as source") +
  ggplot2::ylab("Patch cluster rank as sink")
