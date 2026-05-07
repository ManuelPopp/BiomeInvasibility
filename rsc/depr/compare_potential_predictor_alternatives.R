predictors_0 <- c(
  "logRichnessRatio",
  "logConnAreaRatio",
  "climateStability",
  "logSamplingEffort"
)

predictors_1 <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "climateStability",
  "logSamplingEffort"
)

predictors_2 <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "logRelClimStability",
  "logSamplingEffort"
)

predictors_3 <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "logRelClimVelocity", # Has NA
  "logSamplingEffort"
)

predictors_4 <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "climateVelocityNorm1",
  "logSamplingEffort"
)

predictors_5 <- c(
  "logCorrectedRichnessRatio",
  "logConnAreaRatio",
  "climateVelocityNorm2",
  "logSamplingEffort"
)

predictors_6 <- c(
  "logCorrectedRichnessRatio",
  "logAreaRatio",
  "climateVelocityNorm2",
  "logSamplingEffort"
)

predictors_7 <- c(
  "logCorrectedRichnessRatio",
  "logAreaRatio",
  "logConnAreaRatio",
  "logSamplingEffort"
)

df_mod$climateVelocityNorm1 <- qnorm(rank(df_mod$climateVelocity) / (length(df_mod$climateVelocity) + 1))
df_mod$climateVelocityNorm2 <- log1p(pmax(df_mod$climateVelocity, -1 + 1e-8))
hist(df_mod$climateVelocityNorm2)

pred_set <- 1
gams <- list()
d2s <- c()
concs <- list()
df_pred_d2s <- list()
plots <- list()
for (preds in list(predictors_0, predictors_1, predictors_2, predictors_3, predictors_4, predictors_5)) {
  cat("\nPredictors:", preds)
  df <- df_mod %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(preds),
        ~ is.finite(.)
      )
    )
  cat("\nLength of data set:", nrow(df))
  frml_full <- make_formula(preds, biome = TRUE)
  mod_gam <- mgcv::gam(
    frml_full,
    data = df,
    method = "REML",
    family = stats::binomial("logit"),
    weights = weight
  )
  
  adjD2_full <- ecospat::ecospat.adj.D2.glm(mod_gam)
  conc <- mgcv::concurvity(mod_gam)
  d2s <- c(d2s, adjD2_full)
  gams[[pred_set]] <- mod_gam
  concs[[pred_set]] <- conc
  
  df_pred_d2 <- data.frame()
  pb <- progress::progress_bar$new(total = length(preds) + 1)
  for (i in 1:(length(preds) + 1)) {
    if (i > length(preds)) {
      p <- "Biome"
      frml_p <- as.formula("Invaded ~ Biome")
    } else {
      p <- preds[i]
      frml_p <- make_formula(p, biome = FALSE)
    }
    
    df_rand <- df %>%
      dplyr::mutate(
        "{p}" := sample(.data[[p]])
      )
    
    mod_r <- mgcv::gam(
      make_formula(preds, biome = TRUE),
      data = df_rand,
      method = "REML",
      family = stats::binomial("logit"),
      weights = weight
    )
    adjD2_other <- ecospat::ecospat.adj.D2.glm(mod_r)
    
    mod_p <- mgcv::gam(
      frml_p,
      data = df_rand,
      family = stats::binomial("logit"),
      weights = weight
    )
    adjD2_single <- ecospat::ecospat.adj.D2.glm(mod_p)
    
    df_pred_d2 <- rbind(
      df_pred_d2,
      data.frame(
        Predictor = p,
        Delta_adjD2 = adjD2_full - adjD2_other,
        single_adjD2 = adjD2_single
      )
    )
    pb$tick()
  }
  rm(pb)
  
  df_pred_d2 <- rbind(
    df_pred_d2 %>%
      dplyr::mutate(
        Predictor = sapply(
          Predictor,
          function(x) trimws(sub("log", "", gsub("([A-Z])", " \\1", x)))
        )
      ),
    data.frame(
      Predictor = c("Shared", "Unexplained"),
      Delta_adjD2 = c(adjD2_full - sum(df_pred_d2$Delta_adjD2), 1 - adjD2_full),
      single_adjD2 = c(NA, NA)
    )
  )
  df_pred_d2s[[pred_set]] <- df_pred_d2
  gg_pie <- ggplot2::ggplot(
    data = df_pred_d2,
    ggplot2::aes(x = "", y = Delta_adjD2, fill = Predictor)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Explained Deviance") +
    ggplot2::scale_fill_manual(
      values = c(
        rgb(0, 102/255, 102/255),
        "red",
        "green",
        rgb(70/255, 100/255, 170/255),
        "yellow",
        "orange",
        #"grey25",
        "grey75"
      )
    )
  ggplot2::ggsave(filename = paste0("C:/Users/poppman/Downloads/Pie_", pred_set, ".png"), plot = gg_pie)
  pred_set <- pred_set + 1
}
