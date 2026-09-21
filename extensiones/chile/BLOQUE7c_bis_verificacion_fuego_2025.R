library(terra)
library(dplyr)
library(readr)

VAR_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile/variables"
FIRE_DIR <- file.path(VAR_DIR, "fire_annual")
dir.create(FIRE_DIR, showWarnings = FALSE, recursive = TRUE)

anos <- 2018:2025
resumen <- list()

# AOI de Jubaea chilensis
aoi <- ext(-72.5, -70.0, -35.5, -29.7)

for (y in anos) {
  url <- sprintf("https://storage.googleapis.com/mapbiomas-public/initiatives/chile/fire/collection1/annual_burned_v1/mapbiomas_fire_chile_col1_annual_burned_%d.tif", y)
  dst <- file.path(FIRE_DIR, sprintf("annual_%d.tif", y))
  if (!file.exists(dst)) {
    cat("Descargando", y, "...\n")
    download.file(url, dst, mode = "wb", quiet = TRUE)
  }
  r <- rast(dst)
  r_aoi <- crop(r, aoi)
  # área quemada = píxeles con valor > 0 (0 o NA es no quemado)
  vals <- values(r_aoi)[,1]
  n_quemados <- sum(!is.na(vals) & vals > 0)
  area_km2 <- n_quemados * mean(values(cellSize(r_aoi, unit = "km")), na.rm = TRUE)
  resumen[[as.character(y)]] <- tibble(anio = y, km2_quemados_AOI = round(area_km2, 1),
                                        n_pixeles = n_quemados)
  cat(sprintf("  %d: %d px quemados (%.1f km2) en la AOI\n", y, n_quemados, area_km2))
}

df <- bind_rows(resumen)
cat("\n=== Resumen area quemada anual dentro de la AOI de Jubaea chilensis ===\n")
print(df)

media_2018_2024 <- mean(df$km2_quemados_AOI[df$anio %in% 2018:2024])
cat(sprintf("\nMedia 2018-2024: %.1f km2/anio\n", media_2018_2024))
cat(sprintf("2025:            %.1f km2\n", df$km2_quemados_AOI[df$anio == 2025]))
cat(sprintf("Ratio 2025/media histórica: %.2f\n", df$km2_quemados_AOI[df$anio == 2025] / media_2018_2024))

write_csv(df, "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs/mapbiomas/Col2_2024/verificacion_fuego_anual.csv")
