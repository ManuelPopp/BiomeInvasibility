# connectedness
# connectedness_low_a
# connectedness_high_a
# connectedness_highest_a

boxplot_data <- merged %>%
  dplyr::mutate(
    focalECA = max(total_area, 0, na.rm = TRUE) + connectedness,
    focalECA_low_a = max(total_area, 0, na.rm = TRUE) + connectedness_low_a,
    focalECA_high_a = max(total_area, 0, na.rm = TRUE) + connectedness_high_a,
    focalECA_highest_a = max(total_area, 0, na.rm = TRUE) + connectedness_highest_a,
  ) %>%
  dplyr::select(
    Species, Status, focalECA, focalECA_low_a, focalECA_high_a, focalECA_highest_a
  ) %>%
  dplyr::group_by(Species, Status) %>%
  dplyr::summarise(
    max_focalECA_m2 = if(
      all(is.na(focalECA))
    ) NA else max(
      focalECA, na.rm = TRUE
    ),
    max_focalECA_low_a_m2 = if(
      all(is.na(focalECA_low_a))
    ) NA else max(
      focalECA_low_a, na.rm = TRUE
    ),
    max_focalECA_high_a_m2 = if(
      all(is.na(focalECA_high_a))
    ) NA else max(
      focalECA_high_a, na.rm = TRUE
    ),
    max_focalECA_highest_a_m2 = if(
      all(is.na(focalECA_highest_a))
    ) NA else max(
      focalECA_highest_a, na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    max_focalECA = max_focalECA_m2 / 1e9,
    max_focalECA_low_a = max_focalECA_low_a_m2 / 1e9,
    max_focalECA_high_a = max_focalECA_high_a_m2 / 1e9,
    max_focalECA_highest_a = max_focalECA_highest_a_m2 / 1e9
  ) %>%
  dplyr::group_by(Species) %>%
  dplyr::filter(n() == 2) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(
    cols = c(max_focalECA, max_focalECA_low_a, max_focalECA_high_a, max_focalECA_highest_a),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::filter(
    is.finite(value)
  ) %>%
  dplyr::group_by(metric, Species) %>%
  dplyr::filter(
    all(c("Donor", "Receiver") %in% Status)
  ) %>%
  dplyr::ungroup()

# Calculate statistics
## Calculate statistics with log-scale aware positioning
metric_max_values <- boxplot_data %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    max_value = max(value, na.rm = TRUE),
    log_max_value = log10(max(value, na.rm = TRUE)),
    place_low = median(value, na.rm = TRUE) < (max(value, na.rm = TRUE) / 2)
  )

## Wilcoxon test
test_res <- boxplot_data |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_test(
    value ~ Status,
    paired = TRUE,
    alternative = "greater"
  ) |>
  rstatix::adjust_pvalue(method = "BH")

## Effect size
eff_res <- boxplot_data |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_effsize(
    value ~ Status,
    paired = TRUE
  )

## Create stats dataframe - using log-space coordinates
stats <- test_res %>%
  dplyr::left_join(eff_res %>% select(metric, effsize), by = "metric") %>%
  dplyr::left_join(metric_max_values, by = "metric") %>%
  dplyr::mutate(
    label = paste0(
      "W = ", statistic,
      "\n", "p = ", signif(p, 3),
      "\n", "effect size = ", round(effsize, 2)
    ),
    y.position = max_value,
    .y. = "value"
  ) %>%
  dplyr::select(
    metric, group1, group2, statistic, p, effsize, label, y.position, .y.,
    place_low
  ) %>%
  tidyr::pivot_longer(
    cols = c(group1, group2),
    names_to = "Group",
    values_to = "Status"
  ) %>%
  dplyr::mutate(
    value = y.position,
    place_low = place_low
  )

## Simple boxplots by Status (facetted by metric)
gg_area_comp <- ggplot2::ggplot(boxplot_data, aes(x = Status, y = value, fill = Status)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  ggplot2::facet_wrap(
    ~ metric, scales = "free_y",
    labeller = as_labeller(
      c(
        "max_focalECA_low_a" = "alpha == 5 %*% 10^{-7}",
        "max_focalECA" = "alpha == 6.7 %*% 10^{-7}",
        "max_focalECA_high_a" = "alpha == 1 %*% 10^{-6}",
        "max_focalECA_highest_a" = "alpha == 1.3 %*% 10^{-6}"
      ),
      label_parsed
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c("Donor" = "#2E86AB", "Receiver" = "#A23B72")
  ) +
  # Boxplots are shown on a logarithmic scale to accommodate the large dynamic range of values.
  #ggplot2::scale_y_log10() +
  ggplot2::coord_transform(y = "log10") +
  ggplot2::labs(
    title = NULL,
    y = expression("Area (" * 10^3 * " km"^2 * ")"),
    x = "Biome patch status"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_rect(fill = "lightgray"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::geom_text(
    data = stats,
    ggplot2::aes(
      x = 1.5, label = label, vjust = ifelse(place_low, 4.5, 2)
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(1e-2, 1, 1e2, 1e4),
    labels = scales::label_scientific(),
    sec.axis = sec_axis(
      ~ ., name = "Estimated species count", breaks = NULL, labels = NULL
    )
  )

ggplot2::ggsave(
  filename = file.path(dir_fig, "MaxAreaRichnessBoxplot_alpha_sensitivity.pdf"),
  plot = gg_area_comp,
  width = 10, height = 10
)
