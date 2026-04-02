library("dplyr")
library("nanoparquet")
library("terra")

f_gbif <- "L:/nobis/Manuel/gbif-powo_raw_specID-cellID-occN.parquet"
f_rst <- "L:/nobis/Manuel/gbif/data/map_template.tif"

df_gbif <- nanoparquet::read_parquet(f_gbif)

gbif_records <- df_gbif %>%
  dplyr::group_by(cellID) %>%
  dplyr::summarise(N = sum(occN))
head(gbif_records)

template <- terra::rast(f_rst)
template[gbif_records$cellID] <- gbif_records$N

terra::writeRaster(
  template,
  filename = "L:/poppman/data/bir/dat/lud11/GBIFobsCount.tif",
  overwrite = TRUE
)
