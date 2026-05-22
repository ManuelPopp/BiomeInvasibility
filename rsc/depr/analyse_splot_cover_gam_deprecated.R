vars <- c(
  "log_focalECA",
  "HumanModification",
  "climateVelocityRank"
)

fam <- mgcv::betar(link = "logit")

frml <- make_formula(
  predictors = vars,
  response = "Introduced_abundance",
  k = 3,
  biome = FALSE
)

if(!dir.exists(file.path(dir_fig, "sPlotResponse"))) {
  dir.create(file.path(dir_fig, "sPlotResponse"))
}
for (bi in names(bi_tab[bi_tab >= 15])) {
  sdf <- df_sPlotBiome %>%
    dplyr::filter(Biome == bi) %>%
    dplyr::mutate(
      climateVelocityRank = qnorm(
        rank(climate_velocity_kmpa) / (length(climate_velocity_kmpa) + 1)
      )
    )
  cor(sdf[vars], use = "pairwise.complete.obs")
  
  mod <- mgcv::gam(
    frml,
    family = fam,
    method = "REML",
    data = sdf
  )
  
  d2_adj_full <- ecospat::ecospat.adj.D2.glm(mod)
  cat(bi, d2_adj_full, "(GAM)\n")
  rm(mod)
  gc()
  
  plot_list <- list()
  for (v in vars) {
    x <- sdf[[v]]
    # Fit response curve
    if (all(!is.finite(x))) next
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng))) next
    grid <- as.data.frame(
      lapply(
        vars,
        function(var) {
          mean(sdf[[var]], na.rm = TRUE)
        }
      )
    )
    names(grid) <- vars
    
    xseq <- seq(
      rng[1],
      rng[2],
      length.out = 1000
    )
    
    grid <- grid[rep(1, 1000), ]
    grid[[v]] <- xseq
    
    pr <- predict(
      mod,
      newdata = grid,
      type = "link",
      se.fit = TRUE
    )
    fit <- plogis(pr$fit)
    lwr <- plogis(pr$fit - 1.96 * pr$se.fit)
    upr <- plogis(pr$fit + 1.96 * pr$se.fit)
    
    # Randomise focal predictor
    sdf_rand <- sdf
    d2_adj_c <- c()
    for (i in 1:30) {
      sdf_rand[[v]] <- sample(sdf_rand[[v]])
      
      mod_rand <- mgcv::gam(
        frml,
        family = fam,
        method = "REML",
        data = sdf_rand
      )
      d2_adj_c <- c(d2_adj_c, ecospat::ecospat.adj.D2.glm(mod_rand))
    }
    d2_adj_c <- pmax(pmin(d2_adj_c, 1), 0)
    d2_adj_rand <- mean(d2_adj_c)
    
    delta_d2 <- d2_adj_full - d2_adj_rand
    pred <- predict(mod, newdata = grid, type = "response")
    
    plot_list[[length(plot_list) + 1]] <- data.frame(
      Biome = bi,
      Variable = v,
      x = xseq,
      fit = fit,
      lwr = lwr,
      upr = upr,
      delta_d2 = delta_d2
    )
  }
  
  df_plot <- do.call(rbind, plot_list)
  delta_tab <- df_plot |>
    dplyr::distinct(Variable, delta_d2) |>
    dplyr::mutate(
      label = paste0(
        Variable,
        " ~ Delta*D[adj.]^2==",
        format(round(delta_d2, 3), nsmall = 3)
      )
    )
  
  gg <- ggplot2::ggplot(
    data = df_plot,
    ggplot2::aes(x = x, y = fit)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(
      ~ Variable,
      scales = "free_x",
      labeller = ggplot2::labeller(
        Variable = ggplot2::as_labeller(
          setNames(delta_tab$label, delta_tab$Variable),
          default = ggplot2::label_parsed
        )
      )
    ) +
    ggplot2::labs(x = NULL, y = "Invasive cover") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(
      bquote(
        .(bi) ~ "|"
        ~ D[adj.]^2 == .(round(d2_adj_full, 3)) ~ "|"
        ~ N == .(nrow(sdf))
      )
    )
  
  biname_clean <- gsub("[^[:alpha:]]", "", bi, perl = TRUE)
  ggplot2::ggsave(
    filename = file.path(dir_fig, "sPlotResponse", paste0(biname_clean, ".svg")),
    plot = gg, width = 7, height = 5
  )
}