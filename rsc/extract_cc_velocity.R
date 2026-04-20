
library("dplyr")
library("raster")
library("sp")
library("terra")
library("tidyterra")
library("VoCC")

if (Sys.info()["sysname"] == "Windows") {
  lud11 <- "L:"
} else {
  lud11 <- "/lud11"
}

dir_tas <- file.path(
  lud11,
  "poppman/data/bir/dat/lud11/environment/Bioclim/bio01_proj"
)

dir_imd <- file.path(
  lud11,
  "poppman/data/bir/dat/lud11/environment",
  "cc_vel"
  )

dir_tiles <- file.path(dir_imd, "sg_tiles")
if (!dir.exists(dir_tiles)) {
  dir.create(dir_tiles)
}

if (!dir.exists(dir_imd)) {
  dir.create(dir_imd)
}

f_tr <- file.path(dir_imd, "tr.tif")
f_sg <- file.path(dir_imd, "sg.tif")
f_v <- file.path(dir_imd, "v.tif")

ncol_tiles <- 12
nrow_tiles <- 6
overlap <- 9

# Load biomes
biomes_4326 <- terra::vect(
  file.path(lud11, "poppman/data/bir/dat/lud11/shp/olson_ecoregions/biomes.shp")
)
biomes_4326$ID <- 1:nrow(biomes_4326)

cat("\nTranslating biomes to sp object")
biomes <- biomes_4326 %>%
  as("Spatial")

files <- list.files(dir_tas, pattern = ".tif$", full.names = TRUE)

get_year <- function(x) {
  parts <- stringr::str_split(basename(x), "_")
  year <- as.numeric(parts[[1]][3]) * 1000
  return(year)
}

df <- data.frame(file_path = files) %>%
  dplyr::mutate(
    year = unlist(lapply(X = file_path, FUN = get_year))
  )

df <- df[order(df$year),]

cat("\nLoading raster stack...")
rstack <- raster::stack(df$file_path)

if (!file.exists(f_tr)) {
  cat("\nComputing temporal trend...")
  fit_fun <- function(y) {
    if (all(is.na(y))) return(NA)
    x <- seq_along(y)
    x_mean <- mean(x)
    y_mean <- mean(y)
    sum((x - x_mean) * (y - y_mean)) / sum((x - x_mean)^2)
  }
  
  tr <- terra::app(rstack, fun = fit_fun)
  
  #tr <- VoCC::tempTrend(rstack, th = 10)
  raster::writeRaster(
    tr,
    filename = f_tr,
    format = "GTiff",
    overwrite = TRUE
    )
  cat("\nWrote data to:", f_tr)
} else {
  tr <- raster::raster(f_tr)
}

if (!file.exists(f_sg)) {
  cat("\nComputing spatial gradients...")
  # base raster geometry
  r0 <- rstack[[1]]
  
  # tile boundaries (in cell indices)
  col_breaks <- round(seq(1, ncol(r0) + 1, length.out = ncol_tiles + 1))
  row_breaks <- round(seq(1, nrow(r0) + 1, length.out = nrow_tiles + 1))
  
  # output raster (empty)
  sg_out <- raster::raster(r0)
  
  for (i in seq_len(ncol_tiles)) {
    for (j in seq_len(nrow_tiles)) {
      cat("\nTile:", i, j)
      
      tile_id <- paste0("tile_", i, "_", j)
      f_tile <- file.path(dir_tiles, paste0(tile_id, ".tif"))
      
      if (file.exists(f_tile)) {
        cat(" -> already exists, skipping")
        next
      }
      
      # define tile with overlap
      col_min <- max(1, col_breaks[i] - overlap)
      col_max <- min(ncol(r0), col_breaks[i + 1] - 1 + overlap)
      
      row_min <- max(1, row_breaks[j] - overlap)
      row_max <- min(nrow(r0), row_breaks[j + 1] - 1 + overlap)
      
      e <- raster::extent(
        raster::xFromCol(r0, col_min),
        raster::xFromCol(r0, col_max),
        raster::yFromRow(r0, row_max),
        raster::yFromRow(r0, row_min)
      )
      
      # crop stack
      rst_tile <- raster::crop(rstack, e)
      if (raster::cellStats(rst_tile, stat = "sum", na.rm = TRUE) == 0) {
        cat(" -> skipped (too few valid cells)")
        next
      }
      
      # compute spatial gradient
      sg_tile <- VoCC::spatGrad(rst_tile, th = -Inf, projected = TRUE)
      
      # remove overlap (inner extent only)
      inner_col_min <- col_breaks[i]
      inner_col_max <- col_breaks[i + 1] - 1
      
      inner_row_min <- row_breaks[j]
      inner_row_max <- row_breaks[j + 1] - 1
      
      e_inner <- raster::extent(
        r0, inner_row_min, inner_row_max, inner_col_min, inner_col_max
        )
      
      sg_inner <- raster::crop(sg_tile, e_inner)
      
      raster::writeRaster(
        sg_inner,
        filename = f_tile,
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }
  
  cat("\nMerging tiles...")
  tile_files <- list.files(dir_tiles, pattern = ".tif$", full.names = TRUE)
  
  if (length(tile_files) == 0) {
    stop("No tiles were created.")
  }
  
  sg <- do.call(
    raster::mosaic,
    c(lapply(tile_files, raster::raster), fun = mean, na.rm = TRUE)
    )
  
  raster::writeRaster(
    sg,
    filename = f_sg,
    format = "GTiff",
    overwrite = TRUE
  )
  cat("\nWrote data to:", f_sg)
} else {
  sg <- raster::raster(f_sg)
}

if (!file.exists(f_v)) {
  v <- VoCC::gVoCC(tr, sg)
  raster::writeRaster(
    v,
    filename = f_v,
    format = "GTiff",
    overwrite = TRUE
  )
  cat("\nWrote data to:", f_v)
} else {
  v <- raster::raster(f_v)
}

cat("\nComputing gradient-based climate velocity...")
for (i in 1:nrow(biomes)) {
  biome <- biomes[i,]
  vel <- raster::crop(
    v[[1]],
    raster::extent(
      sp::spTransform(biome, raster::crs(v))
    )
  )
  
  # ensure single layer
  if (raster::nlayers(vel) > 1) {
    vel <- vel[[1]]
  }
  vel <- raster::raster(vel)
  
  cat("\nTransforming and aggregating raster to reduce matrix size...")
  vel_4326 <- raster::projectRaster(vel, crs = "+proj=longlat +datum=WGS84")
  
  if (!raster::hasValues(vel_4326)) {
    next
  }
  
  cat("\nComputing climate residence time within polygons...")
  res_time <- VoCC::resTime(biome, vel_4326, areapg = NA)
  
  cat("\nWriting output...")
  write.csv(
    res_time,
    file = file.path(
      lud11,
      "poppman/data/bir/dat/lud11/", "biomes", "olson", "intermediate_data",
      "clim_res_time.csv"
      ),
    row.names = FALSE,
    append = TRUE
    )
}
cat("\nFinished.")
