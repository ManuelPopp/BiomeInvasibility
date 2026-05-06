# Prerequisites: Generate biomes SpatVector object in analyses.R

resp <- "correctedSR"
all_vars <- c(
  "log_total_area",
  "log_Connectedness",
  "log_sampling_effort",
  "climate_restime_a",
  "climate_velocity_kmpa"
)

vars <- c(
  "log_total_area",
  "log_Connectedness",
  "climate_restime_a" # or climate_velocity_kmpa
)

mod_dat <- as.data.frame(biomes) %>%
  dplyr::filter(
    as.numeric(BIOME) <= 12,
    as.numeric(BIOME) != 10
  ) %>%
  dplyr::mutate(
    Biome = factor(BIOME)
  ) %>%
  dplyr::select(dplyr::all_of(c(resp, vars, "Biome"))) %>%
  dplyr::filter(
    stats::complete.cases(dplyr::across(dplyr::all_of(c(vars, resp))))
  )%>%
  dplyr::filter(dplyr::if_all(everything(), is.finite)) %>%
  dplyr::filter(correctedSR > log(100)) %>%
  as.data.frame()

# Approach 1: Model on corrected estimate
frml <- make_formula(predictors = vars, response = "correctedSR", biome = FALSE)
bis <- sort(unique(mod_dat$Biome))
d2s <- c()
for (bi in bis) {
  mod_sr <- mgcv::gam(
    frml,
    data = mod_dat[which(mod_dat$Biome == bi),]
  )
  adjD2 <- ecospat::ecospat.adj.D2.glm(mod_sr)
  d2s <- c(d2s, adjD2)
  
  # Plot
  bi_num <- as.numeric(bi)
  bi_label <- biome_names[bi_num]
  dat_bi <- mod_dat[mod_dat$Biome == bi, ]
  # ---- compact layout ----
  par(
    mfrow = c(2, 2),
    mar = c(2.2, 2.2, 1.2, 0.5),   # bottom, left, top, right
    mgp = c(1.3, 0.3, 0),          # axis title, labels, line
    tcl = -0.2,                    # short ticks
    oma = c(0, 0, 2, 0)            # outer margin for biome title
  )
  
  for (v in vars) {
    
    x_seq <- seq(
      min(dat_bi[[v]], na.rm = TRUE),
      max(dat_bi[[v]], na.rm = TRUE),
      length.out = 100
    )
    
    newdat <- data.frame(
      log_total_area = mean(dat_bi$log_total_area, na.rm = TRUE),
      log_Connectedness = mean(dat_bi$log_Connectedness, na.rm = TRUE),
      climate_restime_a = mean(dat_bi$climate_restime_a, na.rm = TRUE),
      climate_velocity_kmpa = mean(dat_bi$climate_velocity_kmpa, na.rm = TRUE)
    )
    
    newdat <- newdat[rep(1, 100), ]
    newdat[[v]] <- x_seq
    
    pr <- predict(mod_sr, newdata = newdat, type = "terms", se.fit = TRUE)
    term_name <- paste0("s(", v, ")")
    
    fit <- pr$fit[, term_name]
    se  <- pr$se.fit[, term_name]
    
    # ---- plot ----
    plot(
      x_seq, fit, type = "l", lwd = 1.5,
      xlab = "", ylab = "",
      axes = FALSE
    )
    
    # confidence band
    polygon(
      c(x_seq, rev(x_seq)),
      c(fit + 2 * se, rev(fit - 2 * se)),
      col = adjustcolor("grey80", alpha.f = 0.6),
      border = NA
    )
    
    lines(x_seq, fit, lwd = 1.5)
    
    axis(1)
    axis(2)
    
    box()
    
    mtext(v, side = 1, line = 1.2, cex = 0.8)
    mtext("Effect", side = 2, line = 1.2, cex = 0.8)
  }
  
  # empty 4th panel
  plot.new()
  
  # ---- biome title once per panel ----
  mtext(bi_label, outer = TRUE, line = 0.5, cex = 1.1, font = 2)
}

res1 <- data.frame(Biome = bis, D2 = d2s)
round(d2s * 100)
