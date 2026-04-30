future::plan(multisession, workers = parallel::detectCores() - 1)

df_invasion <- dplyr::bind_rows(
  future.apply::future_lapply(
    X = unique(merged_env$Species),
    FUN = function(species) {
      donor <- merged_env %>%
        dplyr::filter(
          Species == species,
          !is.na(speciesRichnessBa),
          Status == "Donor"
        )
      if (nrow(donor) == 0) return(NULL)
      
      donor_stats <- donor %>%
        dplyr::summarise(
          modalBiome = collapse::fmode(Biome),
          max_total_area_km2 = max(total_area_km2, na.rm = TRUE),
          max_lowland_area_km2 = max(lowland_area_km2, na.rm = TRUE),
          max_focalECA = max(focalECA, na.rm = TRUE),
          max_connectedness = max(connectedness, na.rm = TRUE),
          maxRichness = max(speciesRichnessBa, na.rm = TRUE),
          max_climateStability = max(climateStability, na.rm = TRUE)
        )
      
      native_environment <- donor %>%
        dplyr::filter(lowland_area_km2 == donor_stats$max_lowland_area_km2) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::select(centre, Sigma)
      
      if (
        nrow(native_environment) == 0 ||
        is.null(native_environment$centre[[1]]) ||
        is.null(native_environment$Sigma[[1]])
        ) {
        return(NULL)
        }
      
      receiver <- merged_env %>%
        dplyr::filter(
          Species == species,
          !is.na(speciesRichnessBa),
          Status == "Receiver"
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Invaded = 1,
          tmp = list(
            save_bhattacharyya(
              mu1 = patchMu,
              mu2 = native_environment$centre[[1]],
              Sigma1 = patchSigma,
              Sigma2 = native_environment$Sigma[[1]]
            )
          ),
          Bhattacharyya = as.numeric(tmp$value),
          Bhattacharyya_reason = tmp$reason
        )
      
      id_has_species <- merged_env %>%
        dplyr::filter(Species == species) %>%
        dplyr::pull(ID) %>%
        unique()
      
      background <- merged_env %>%
        dplyr::filter(
          !ID %in% id_has_species,
          !is.na(speciesRichnessBa)
        ) %>%
        dplyr::distinct(ID, .keep_all = TRUE) %>%
        dplyr::slice_sample(n = 100) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Status = "absent",
          Invaded = 0,
          tmp = list(
            save_bhattacharyya(
              mu1 = patchMu,
              mu2 = native_environment$centre[[1]],
              Sigma1 = patchSigma,
              Sigma2 = native_environment$Sigma[[1]]
            )
          ),
          Bhattacharyya = as.numeric(tmp$value),
          Bhattacharyya_reason = tmp$reason
        )
      
      # Generate output data.frame
      out <- dplyr::bind_rows(
        receiver,
        background
      ) %>%
        dplyr::select(
          ID, clusterID,
          Species, Biome, Count,
          Status, Invaded,
          total_area_km2, lowland_area_km2, focalECA,
          speciesRichnessBa, sampling_effort,
          Bhattacharyya, Bhattacharyya_reason,
          climateStability
        ) %>%
        dplyr::mutate(
          donor_modalBiome = donor_stats$modalBiome,
          donor_max_total_area_km2 = donor_stats$max_total_area_km2,
          donor_max_lowland_area_km2 = donor_stats$max_lowland_area_km2,
          donor_max_focalECA = donor_stats$max_focalECA,
          donor_max_connectedness = donor_stats$max_connectedness,
          donor_maxRichness = donor_stats$maxRichness,
          donor_maxClimateStability = max(climateStability, na.rm = TRUE)
        )
      
      return(out)
    }
  )
) %>%
  dplyr::mutate(
    areaRatio = lowland_area_km2 / donor_max_lowland_area_km2,
    connAreaRatio = focalECA / donor_max_focalECA,
    richnessRatio = speciesRichnessBa / donor_maxRichness,
    logAreaRatio = log(areaRatio),
    logConnAreaRatio = log(connAreaRatio),
    logRichnessRatio = log(richnessRatio),
    logBhattacharyya = log(Bhattacharyya),
    logSpeciesRichnessBa = log(speciesRichnessBa),
    logSamplingEffort = log(sampling_effort)
  ) %>%
  as.data.frame()

## Debugging/analytics
reason_summary <- table(df_invasion$Bhattacharyya_reason)
cat(
  "Bhattacharyya distances computed for",
  sum(df_invasion$Bhattacharyya_reason == "ok"), "of", nrow(df_invasion),
  "observations (i.e.,",
  round(100 * mean(df_invasion$Bhattacharyya_reason == "ok"), 1), "%).\n",
  "Failure breakdown:\n"
)
print(reason_summary)

par(mfrow = c(2, 2))
hist(df_invasion$logBhattacharyya, main = "Bhattacharyya (log transformed)")
hist(df_invasion$logAreaRatio, main = "Area ratio (log transformed)")
hist(df_invasion$logRichnessRatio, main = "Spechies richness ratio (log transformed)")
hist(df_invasion$logSamplingEffort, main = "Sampling effort (log transformed)")
par(mfrow = c(1, 1))
hist(log(df_invasion$climateStability), main = "Climate stability")

df_mod <- df_invasion[stats::complete.cases(df_invasion),] %>%
  dplyr::mutate(
    Specieslvl = stringr::word(Species, 1, 2)
  ) %>%
  dplyr::filter(
    as.numeric(Biome) <= 14
  )

# Fit models and compute explained deviance by predictor
## Full model
predictors <- c(
  "logBhattacharyya",
  "logRichnessRatio",
  "logConnAreaRatio",
  "logSamplingEffort",
  "climateStability"
)

frml_full <- make_formula(predictors, biome = TRUE)

mod_gam <- mgcv::gam(
  frml_full,
  data = df_mod,
  #method = "REML",
  family = stats::binomial("logit")
)

summary(mod_gam)
ecospat::ecospat.adj.D2.glm(mod_gam)

mgcv::gam.check(mod_gam)
mgcv::concurvity(mod_gam)


mods_gam <- list()
d2s_gam <- c()
N <- 50
pb <- progress::progress_bar$new(total = N)
for (i in 1:N) {
  set.seed(i)
  df_ss <- df_mod %>%
    dplyr::group_by(Specieslvl, Status) %>%
    dplyr::sample_n(size = 1)
  
  mods_gam[[i]] <- mgcv::gam(
    frml_full,
    data = df_ss,
    family = stats::binomial("logit")
  )
  d2s_gam <- c(
    d2s_gam,
    ecospat::ecospat.adj.D2.glm(mods_gam[[i]])
  )
  pb$tick()
}
rm(pb)


# Fit models to subsets of the data and plot response shapes
means <- colMeans(df_mod[, predictors], na.rm = TRUE)

plot_dfs <- lapply(predictors, function(v) {
  x_seq <- seq(
    min(df_mod[[v]], na.rm = TRUE),
    max(df_mod[[v]], na.rm = TRUE),
    length.out = 500
  )
  
  newdata <- as.data.frame(as.list(means))[rep(1, 500), ]
  newdata[[v]] <- x_seq
  newdata$Biome <- df_mod$Biome[1]
  
  # predictions: matrix (rows = x, cols = models)
  pred_mat <- sapply(
    X = mods_gam,
    FUN = function(m) {predict(m, newdata = newdata, type = "response")}
  )
  
  # summary stats
  data.frame(
    x = x_seq,
    mean = rowMeans(pred_mat, na.rm = TRUE),
    lo = apply(pred_mat, 1, quantile, probs = 0.1, na.rm = TRUE),
    hi = apply(pred_mat, 1, quantile, probs = 0.9, na.rm = TRUE),
    variable = v
  )
})

df_plot <- do.call(rbind, plot_dfs) %>%
  dplyr::mutate(
    Predictor = sapply(
      variable,
      function(x) trimws(sub("log", "", gsub("([A-Z])", " \\1", x)))
    )
  )


# Calculate maximum Bhattacharyya distance to limit selected habitats
xseq <- seq(
  min(df_mod$logBhattacharyya, na.rm = TRUE),
  max(df_mod$logBhattacharyya, na.rm = TRUE),
  length.out = 500
)

newdata <- data.frame(
  logBhattacharyya = xseq,
  logRichnessRatio = mean(df_mod$logRichnessRatio, na.rm = TRUE),
  logConnAreaRatio = mean(df_mod$logConnAreaRatio, na.rm = TRUE),
  logSamplingEffort = mean(df_mod$logSamplingEffort, na.rm = TRUE),
  climateStability = mean(df_mod$climateStability, na.rm = TRUE),
  Biome = dplyr::first(df_mod$Biome)
)

df_pred <- newdata
df_pred$yhat <- stats::predict(mod_gam, newdata = newdata, type = "response")
dy <- diff(df_pred$yhat) / diff(df_pred$logBhattacharyya)
i_min <- which.min(dy)
x_star <- xseq[i_min]

max_lobBhat <- round(x_star, 2)
# max_lobBhat = 1.34

gg_resp_bhat <- ggplot2::ggplot(
  dplyr::filter(df_plot, Predictor == "Bhattacharyya"),
  ggplot2::aes(x = x, y = mean)
  ) +
  ggplot2::geom_line(colour = rgb(0, 102/255, 102/255)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lo, ymax = hi),
    alpha = 0.2, fill = rgb(0, 102/255, 102/255)
  ) +
  ggplot2::labs(
    x = "Predictor value (log scaled)",
    y = "Predicted invasion probability",
  ) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(. ~ Predictor, scales = "free_x") +
  ggplot2::geom_vline(xintercept = max_bhat)

ggplot2::ggsave(
  filename = file.path(dir_fig, "ResponseShapes_Bhattacharyya_included.svg"),
  plot = gg_resp_bhat,
  width = 7, height = 5
)