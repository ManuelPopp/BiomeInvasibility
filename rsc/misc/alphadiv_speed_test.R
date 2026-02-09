# Compare the speed of different algorithms to estimate alpha diversity from
# opportunistic sampling
install.packages("BiocManager")
BiocManager::install("phyloseq")
remotes::install_github("adw96/breakaway")

library(microbenchmark)
library(iNEXT)
library(vegan)
library(breakaway)

set.seed(42)

# Parameters
n_polygons <- 100
n_species <- seq(50, 20000, 200)
total_species <- 530000
max_obs_per_species <- 5000

# Create polygon names
polygons <- paste0("Poly", 1:n_polygons)

# Initialize a list to hold abundances per polygon
abund_list <- vector("list", length = n_polygons)
names(abund_list) <- polygons

for (i in 1:n_polygons) {
  n_spp_in_poly <- n_species[i]
  species_ids <- sample(1:total_species, n_spp_in_poly)
  abundances_raw <- rlnorm(n_spp_in_poly)
  abundances <- round(abundances_raw/max(abundances_raw)*max_obs_per_species)
  abund_list[[i]] <- abundances
}

results <- data.frame()
for (i in 1:5) {
  r <- microbenchmark::microbenchmark(
    iNEXT::iNEXT(
      x = abund_list[[i]],
      q = 0,
      datatype = "abundance"
    ),
    vegan::estimateR(abund_list[[i]]),
    breakaway::breakaway(abund_list[[i]]),
    times = 3
  )
  df <- data.frame(
    algorithm = rep(c("iNEXT", "estimateR", "breakaway"), 3),
    value = c(summary(r)$min, summary(r)$mean, summary(r)$max)
  )
  df$run <- i
  results <- rbind(results, df)
}

results$n_samples <- as.numeric(results$run * 200 - 150)

ggplot2::ggplot(
  data = results,
  ggplot2::aes(
    x = n_samples, y = value * 1e-9, colour = algorithm, shape = algorithm
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth() +
  ggplot2::theme_bw()

mtimes <- c()

for (i in 1:50) {
  r <- microbenchmark::microbenchmark(
    breakaway::breakaway(abund_list[[i]]),
    times = 3
  )
  mtimes <- c(mtimes, summary(r)$mean)
}

results_brwy <- data.frame(time = mtimes)
results_brwy$n_samples <- n_species[1:50]

ggplot2::ggplot(
  data = results_brwy,
  ggplot2::aes(
    x = n_samples, y = time * 1e-9
  )
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::theme_bw()

estimate_richness <- function(abundances) {
  n_samples <- 1:sum(abundances)/1000
  n_sampled <- unlist(
    lapply(
      X = n_samples,
      FUN = function(n) {
        species <- sample(
          1:length(abundances),
          size = n,
          prob = abundances,
          replace = TRUE
          )
        return(length(unique(species)))
      }
    )
  )
  plot(x = n_samples, y = n_sampled)
}