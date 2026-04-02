library("terra")
library("corrplot")
library("dplyr")
library("ggplot2")

dir_env <- "L:/poppman/data/bir/dat/lud11/environment"
f_env_out <- "L:/poppman/data/bir/dat/lud11/environment/environment.tif"
f_env_pc <- "L:/poppman/data/bir/dat/lud11/environment/environmentPC.tif"
f_plt_corr <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/environment_corr.svg"
f_plt_local <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/prj/bir/fig/environment_pc.svg"
f_plt_online <- "C:/Users/poppman/Dropbox/Apps/Overleaf/BiomeInvasibility/suppl_files/environment_pc.pdf"

# Load biomes as land mask
f_bfullinfo <- file.path(
  "L:/poppman/data/bir/dat/lud11/biomes/olson/intermediate_data",
  "biomes_full_info.gpkg"
  )
biomes <- terra::vect(f_bfullinfo)

# Load environmental variables
chelsa_files <- list.files(
  dir_env,
  pattern = "CHELSA_",
  full.names = TRUE
) %>%
  subset(endsWith(., suffix = ".tif"))

vars = sub(".*CHELSA_(.*?)_1981.*", "\\1", chelsa_files)

env_stack <- terra::rast(chelsa_files)
names(env_stack) <- vars
pts <- terra::project(biomes, "epsg:6933") %>%
  terra::spatSample(size = 100000, method = "random") %>%
  terra::project("epsg:4326")
plot(biomes, col = "lightgreen")
plot(pts, col = "darkblue", cex = 0.2, add = TRUE)

values <- terra::extract(env_stack, pts)

# Plot correlations
v <- values[, -1]
lookup = setNames(
  c("T mean", "T seasonality", "T max", "T min", "P seasonality", "P warmest quart", "GDD5", "GSL", "GSP", "NGD5", "SWB", "Snow mass"),
  c("bio01", "bio04", "bio05", "bio06", "bio15", "bio18", "gdd5", "gsl", "gsp", "ngd5", "swb", "swe")
  )
names(v) <- lookup[vars]
M <- cor(v, use = "complete.obs")
corrplot::corrplot(M)

D <- 1 - abs(M)
hc <- hclust(as.dist(D), method = "average")
plot(hc, xlab = "", ylab = "", main = "Variable clustering")

svg(filename = f_plt_corr, width = 5, height = 5)
corrplot::corrplot(M, order = "AOE", hclust.method = "average")
dev.off()

# Select a subset of variables
v_select <- c(
  "bio06", # T min
  "bio15", # Precipitation seasonality
  "bio18", # Precipitation of the warmest 3 months
  "gdd5", # Growing degree days above 5 degrees
  "swb"#, # Site water balance
  #"swe" # Snow water equivalent # Zero-inflated distribution plus Philipp said the variable has high uncertainties
)
f_select <- chelsa_files[which(vars %in% v_select)]
var_names = sub(".*CHELSA_(.*?)_1981.*", "\\1", f_select)
r <- terra::rast(f_select) %>%
  terra::mask(biomes)
names(r) <- var_names

v_select <- values[, which(names(values) %in% var_names)]

par(mfrow = c(ceiling(sqrt(length(var_names))), ceiling(sqrt(length(var_names)))))
lapply(
  X = var_names,
  FUN = function(x) {hist(v_select[, x], main = x)}
)

scale_fun <- function(x, name) {
  if (name == "bio06") {
    (x - min(x, na.rm = TRUE))^1.5
  } else if (name == "bio15") {
    sqrt(x)
  } else if (name == "bio18") {
    log(x + 1, base = 1.5)
  } else if (name == "gdd5") {
    sqrt(x)
  } else if (name == "swb") {
    (x - min(x, na.rm = TRUE))^2
  } else {
    x
  }
}

par(mfrow = c(ceiling(sqrt(length(var_names))), ceiling(sqrt(length(var_names)))))
for (i in var_names) {
  hist(scale_fun(x = v_select[, i], name = i), main = i)
}

# TODO: Re-write the scaling functions to make them work properly with raster data
scale_fun_rst = function(r) {
  nms = names(r)
  # Precompute minima
  mins = global(r, "min", na.rm = TRUE)[,1]
  names(mins) = nms
  out = lapply(nms, function(nm) {
    x = r[[nm]]
    if (nm == "bio06") {
      (x - mins[nm])^1.5
    } else if (nm == "bio15") {
      sqrt(x)
    } else if (nm == "bio18") {
      log(x + 1) / log(1.5)
    } else if (nm == "gdd5") {
      sqrt(x)
    } else if (nm == "swb") {
      (x - mins[nm])^2
    } else if (nm == "swe") {
      log(x + 1)
    } else {
      x
    }
  })
  
  rast(out)
}

r_scaled = scale_fun_rst(r)

hist(r_scaled)

terra::writeRaster(
  r_scaled,
  filename = f_env_out,
  overwrite = TRUE
)


pca <- terra::prcomp(
  r_scaled,
  center = TRUE,
  scale. = TRUE,
  maxcell = 50000000
)

plotdf <- t(summary(pca)$importance) %>%
  as.data.frame() %>%
  dplyr::mutate(n_comp = dplyr::row_number()) %>%
  dplyr::rename_with(.fn = function(x) {gsub(" ", "_", x)})

gg_pca <- ggplot2::ggplot(data = plotdf, ggplot2::aes(x = n_comp)) +
  ggplot2::geom_bar(
    ggplot2::aes(y = Cumulative_Proportion),
    stat = "identity", fill = grDevices::rgb(0, 102/255, 102/255, 0.5)
    ) +
  ggplot2::geom_line(ggplot2::aes(y = Proportion_of_Variance), linewidth = 1.25) +
  ggplot2::geom_point(ggplot2::aes(y = Proportion_of_Variance), size = 2) +
  ggplot2::xlab("Number of principal components") +
  ggplot2::ylab("(Cumulative) proportion of variance") +
  ggplot2::theme_bw()

ggplot2::ggsave(filename = f_plt_local, plot = gg_pca, width = 7, height = 5)
ggplot2::ggsave(filename = f_plt_online, plot = gg_pca, width = 7, height = 5)


k <- 3

r_pca <- terra::predict(r_scaled, pca)
r_pc <- r_pca[[1:k]]

# "Buffer" by one pixel (1 km) into NA areas to avoid NAs resulting from minor
# location inaccuracies
r_filled <- terra::ifel(
  is.na(r_pc),
  terra::focal(r_pc, w = 3, fun = function(x) x[which(!is.na(x))[1]]),
  r_pc
)

r_filled %>%
  terra::writeRaster(
    filename = f_env_pc,
    overwrite = TRUE
  )

hist(r_pca)
