plot_biomes <- polygons_sf[which(!is.na(polygons_sf$speciesRichnessBa)),] %>%
  dplyr::mutate(
    srph = log(speciesRichnessBa / lowland_area * 1e4)
    )


gg_srlog <- ggplot2::ggplot(polygons_sf[which(!is.na(polygons_sf$speciesRichnessBa)),]) +
  ggplot2::geom_sf(data = biomes, colour = "black", fill = NA) +
  ggplot2::geom_sf(ggplot2::aes(fill = log(speciesRichnessBa))) +
  ggplot2::scale_fill_viridis_c(
    name = "Species count\n(log)",
    option = "viridis",
    na.value = "grey90"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "right",
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )


gg_srph <- ggplot2::ggplot(plot_biomes) +
  ggplot2::geom_sf(data = plot_biomes, colour = "black", fill = NA) +
  ggplot2::geom_sf(ggplot2::aes(fill = srph)) +
  ggplot2::scale_fill_viridis_c(
    name = "Species count\nper ha (log)",
    option = "viridis",
    na.value = "grey90"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "right",
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )
require("gridExtra")
gg_sr <- gridExtra::grid.arrange(gg_srlog, gg_srph, ncol = 2)
ggplot2::ggsave(
  filename = file.path(dir_fig, "est_species_counts.png"),
  plot = gg_sr,
  width = 13, height = 3
)
