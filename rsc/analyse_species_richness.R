# Prerequisites: Generate biomes SpatVector object in analyses.R

f_npp <- file.path(
  dir_lud11, "poppman", "data", "bir", "dat", "npp", "npp_max_by_biome.csv"
  )

resp <- "correctedSR"
all_vars <- c(
  "log_total_area",
  "log_Connectedness",
  "log_sampling_effort",
  "climate_restime_a",
  "climate_velocity_kmpa",
  "npp"
)

vars <- c(
  "npp",
  "log_total_area",
  "log_Connectedness",
  "climate_restime_a" # or climate_velocity_kmpa
)

npp <- read.csv(f_npp) %>%
  dplyr::select(ID, max) %>%
  dplyr::rename(npp = max)

mod_dat <- as.data.frame(biomes) %>%
  dplyr::left_join(npp, by = "ID") %>%
  dplyr::filter(
    as.numeric(BIOME) <= 12,
    as.numeric(BIOME) != 10
  ) %>%
  dplyr::mutate(
    Biome = factor(BIOME),
    correctedSR = round(exp(logCorrectedSR)),
    algorithmSR = round(speciesRichnessBa),
    meanSR = round((correctedSR + algorithmSR) / 2)
  ) %>%
  dplyr::select(
    dplyr::all_of(
      c(resp, vars, "Biome", "algorithmSR", "correctedSR", "meanSR")
      )
    ) %>%
  dplyr::filter(
    stats::complete.cases(dplyr::across(dplyr::all_of(c(vars, resp))))
  ) %>%
  dplyr::filter(dplyr::if_all(everything(), is.finite)) %>%
  dplyr::filter(correctedSR > 100) %>%
  as.data.frame()

# Approach 1: Model on corrected estimate

if (!dir.exists(file.path(dir_fig, "species_richness"))) {
  dir.create(file.path(dir_fig, "species_richness"))
}

frml <- make_formula(predictors = vars, response = resp, biome = FALSE)
bis <- sort(unique(mod_dat$Biome))
d2s <- c()
for (bi in bis) {
  delta_adjD2_p <- c()
  mod_sr <- mgcv::gam(
    frml,
    method = "REML",
    family = mgcv::nb(),
    data = mod_dat[which(mod_dat$Biome == bi),]
  )
  
  # frm_glm <- make_formula(
  #   predictors = vars, k = 0, response = resp, biome = FALSE, interact = 2
  # )
  # mod_sr <- MASS::glm.nb(
  #   frm_glm,
  #   data = mod_dat[which(mod_dat$Biome == bi),]
  # )
  
  
  cat("\nBiome:", bi, "Dispersion:", dispersion(mod_sr))
  adjD2 <- ecospat::ecospat.adj.D2.glm(mod_sr)
  d2s <- c(d2s, adjD2)
  
  # Plot
  bi_num <- as.numeric(bi)
  bi_label <- biome_names[bi_num]
  dat_bi <- mod_dat[mod_dat$Biome == bi, ]
  biname_clean <- gsub("[^[:alpha:]]", "", bi_label, perl = TRUE)
  
  svg(
    filename = file.path(dir_fig, "species_richness", paste0(biname_clean, "_SR.svg")),
    width = 7, height = 5
    )
  # ---- compact layout ----
  par(
    mfrow = c(2, 2),
    mar = c(2.2, 2.2, 1.2, 2.0),
    mgp = c(1.3, 0.3, 0),
    tcl = -0.2,
    oma = c(0, 0, 2, 0)
  )
  
  for (v in vars) {
    dat_rand <- dat_bi
    adjD2_rand_c <- c()
    for (i in 1:100) {
      dat_rand[[v]] <- sample(dat_rand[[v]])
      
      mod_rand <- mgcv::gam(
        frml,
        method = "REML",
        family = mgcv::nb(),
        data = dat_rand
      )
      
      adjD2_rand_c <- c(adjD2_rand_c, ecospat::ecospat.adj.D2.glm(mod_rand))
    }
    adjD2_rand <- mean(adjD2_rand_c, na.rm = TRUE)
    delta_adjD2_p[v] <- adjD2 - adjD2_rand
    
    x_seq <- seq(
      min(dat_bi[[v]], na.rm = TRUE),
      max(dat_bi[[v]], na.rm = TRUE),
      length.out = 100
    )
    
    newdat <- data.frame(
      npp = mean(dat_bi$npp, na.rm = TRUE),
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
    # ----- deviance bar -----
    fill_col <- ifelse(delta_adjD2_p[v] >= 0, "black", "red")
    frac <- abs(delta_adjD2_p[v] / adjD2) %>%
      pmin(1, na.rm = TRUE)
    
    usr <- par("usr")
    
    xleft <- usr[2] + 0.04 * (usr[2] - usr[1])
    xright <- usr[2] + 0.10 * (usr[2] - usr[1])
    
    ybottom <- usr[3]
    ytop <- usr[4]
    
    # outline
    rect(
      xleft,
      ybottom,
      xright,
      ytop,
      xpd = NA
    )
    
    # filled fraction
    rect(
      xleft,
      ybottom,
      xright,
      ybottom + frac * (ytop - ybottom),
      col = fill_col,
      border = NA,
      xpd = NA
    )
    
    mtext(v, side = 1, line = 1.2, cex = 0.8)
  }
  
  # ---- biome title once per panel ----
  mtext(
    bquote(
      .(bi_label) ~ "|" ~ N==.(nrow(dat_bi)) ~ "|" ~ D[adj.]^2 == .(
        round(adjD2, 3)
        )
      ),
    outer = TRUE, line = -0.5, cex = 1.1, font = 2
    )
  mtext("Effect", outer = TRUE, side = 2, line = -1, cex = 0.8)
  # empty 4th panel
  #plot.new()
  dev.off()
}

res1 <- data.frame(Biome = bis, D2 = d2s)
round(d2s * 100)

