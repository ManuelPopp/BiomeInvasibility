f_local <- "L:/poppman/tmp/df_eval.Rsave"
f_cluster <- "/lud11/poppman/tmp/df_eval.Rsave"

if (!file.exists(f_local)) {
  library("breakaway")
  
  set.seed(1)
  n_plots <- 1000
  
  # --- 1. Species pool ---
  n_species <- 150000
  sp <- paste0("sp_", 1:n_species)
  
  # uneven species prevalence (realistic SAD)
  prob <- 1 / (1:n_species)^1.2
  prob <- prob / sum(prob)
  
  # species-specific detectability
  detectability <- rbeta(n_species, 0.5, 5)
  names(detectability) <- sp
  
  # true richness per plot
  S_true <- round(rlnorm(n_plots, meanlog = 9, sdlog = 0.75))
  S_true[S_true > n_species] <- n_species
  
  # sampling effort per plot (drives bias)
  effort <- rlnorm(n_plots, meanlog = 5, sdlog = 2)
  
  # --- 3. Generate assemblages + observations ---
  results <- lapply(
    X = seq_len(n_plots),
    FUN = function(i) {
      # --- true assemblage ---
      present <- sample(sp, size = S_true[i], prob = prob, replace = FALSE)
      
      # --- number of observations (sample size varies) ---
      n_obs <- max(1, rpois(1, lambda = effort[i]))
      
      # --- biased sampling probabilities ---
      abund <- rgamma(length(present), shape = 0.5, rate = 1)
      p <- detectability[present] * abund
      p <- p / sum(p)
      
      # --- draw observations (with replacement) ---
      obs <- sample(present, size = n_obs, replace = TRUE, prob = p)
      
      # --- convert to counts per species ---
      counts <- table(obs)
      
      # breakaway expects integer counts vector (no zeros)
      obs_counts <- as.numeric(counts)
      
      # --- richness estimates ---
      est <- tryCatch(
        breakaway::breakaway(obs_counts),
        error = function(e) NA
      )
      
      est_richness <- if (is.list(est)) est$estimate else NA
      
      # --- observed richness ---
      obs_richness <- length(counts)
      
      data.frame(
        plot = i,
        trueSR = S_true[i],
        observedSR = obs_richness,
        estimatedSR = est_richness,
        sampling_effort = effort[i],
        n_obs = n_obs
      )
    }
  )
  
  # --- 4. Combine results ---
  df_eval <- do.call(rbind, results)
  
  # --- 5. Quick diagnostics ---
  summary(df_eval)
  
  # correlation checks
  cor(df_eval$trueSR, df_eval$observedSR, use = "complete.obs")
  cor(df_eval$trueSR, df_eval$estimatedSR, use = "complete.obs")
  cor(df_eval$sampling_effort, df_eval$observedSR, use = "complete.obs")
  cor(df_eval$sampling_effort, df_eval$estimatedSR, use = "complete.obs")
  save(df_eval, file = f_cluster)
} else {
  load("L:/poppman/tmp/df_eval.Rsave")
}

plot(log(estimatedSR) ~ log(sampling_effort), data = df_eval)


# Try to correct
df_biomes <- df_eval %>%
  dplyr::mutate(
    log_SpeciesRichness = log(estimatedSR),
    log_sampling_effort = log(sampling_effort)
  )

df_fitSR <- as.data.frame(df_biomes) %>%
  dplyr::filter(log_SpeciesRichness > log(100))

m1SR <- mgcv::gam(
  log_SpeciesRichness ~ s(log_sampling_effort, k = 3),
  data = df_fitSR
)
mgcv::gam.check(m1SR)

effort_fitSR <- stats::predict(m1SR)
max_sampling_effectSR <- stats::predict(
  m1SR,
  newdata = df_fitSR %>%
    dplyr::slice_min(
      abs(
        log_sampling_effort - quantile(log_sampling_effort, 0.90, na.rm = TRUE)
      ),
      n = 1
    )
)
df_biomes$correctedSR <- NA
df_biomes$correctedSR[
  which(df_biomes$log_SpeciesRichness > log(100))
  ] <- data.frame(
  corrected = df_fitSR$log_SpeciesRichness - effort_fitSR + max_sampling_effectSR[1],
  estimated = df_fitSR$log_SpeciesRichness
) %>%
  dplyr::mutate(
    maximum = pmax(corrected, estimated, na.rm = TRUE)
  ) %>%
  dplyr::pull(maximum)


lm_nocorrect <- lm(log_SpeciesRichness ~ log(trueSR), data = df_biomes)
r2_nocorrect <- round(summary(lm_nocorrect)$r.squared, 2)
lm_corrected <- lm(correctedSR ~ log(trueSR), data = df_biomes)
r2_corrected <- round(summary(lm_corrected)$r.squared, 2)

pdf(
  file = "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/species_richness_est_correction.pdf",
  width = 10, height = 5
  )
par(
  mfrow = c(1, 2),
  mar = c(3, 3, 2, 0.25),
  mgp = c(1.6, 0.5, 0),
  tcl = -0.2
)
plot(
  log_SpeciesRichness ~ log(trueSR), data = df_biomes,
  main = bquote("R"^2==.(r2_nocorrect)),
  xlab = "True species richness (log-scaled)",
  ylab = "Estimated species richness (log-scaled)",
  pch = 16, col = rgb(0, 70/255, 130/255, 0.5)
  )
abline(
  a = coefficients(lm_nocorrect)[1], b = coefficients(lm_nocorrect)[2],
  col = "orange", lwd = 2
  )

plot(
  correctedSR ~ log(trueSR), data = df_biomes,
  main = bquote("R"^2==.(r2_corrected)),
  xlab = "True species richness (log-scaled)",
  ylab = "Corrected estimate (log-scaled)",
  pch = 16, col = rgb(0, 102/255, 102/255, 0.5)
  )
abline(
  a = coefficients(lm_corrected)[1], b = coefficients(lm_corrected)[2],
  col = "orange", lwd = 2
)

dev.off()
