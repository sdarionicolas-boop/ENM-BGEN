# =============================================================================
#   BLOQUE 8 REDO: recalcula el ranking Ipc de partidos bonaerenses usando
#   los rasters EMmean (BLOQUE4_REDO_EMmean.R) en vez de los oficiales
#   (EMwmeanByTSS), para ver si el fix del sesgo de ponderación cambia el
#   ranking que la tesis usa para priorizar campañas de recolección.
#
#   Misma lógica que BLOQUE8_partidos_bonaerenses.R, apuntando a
#   outputs_emmean/ en vez de outputs/. No descarga GADM de nuevo (usa el
#   cache en variables/). No genera el mapa leaflet, solo las tablas.
# =============================================================================

library(terra)
library(geodata)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs_emmean")
RANK_DIR   <- file.path(OUT_DIR, "ranking_partidos_emmean")
dir.create(RANK_DIR, showWarnings = FALSE, recursive = TRUE)

gadm_arg <- gadm(country = "ARG", level = 2, path = VAR_DIR)
ba_partidos <- gadm_arg[gadm_arg$NAME_1 == "Buenos Aires", ]
cat("Partidos bonaerenses:", nrow(ba_partidos), "\n")

tifs_prob <- list.files(OUT_DIR, pattern = "_idoneidad_prob.*\\.tif$", recursive = TRUE, full.names = TRUE)
cat("Rasters EMmean encontrados:", length(tifs_prob), "\n")

nombres_especies <- tifs_prob %>%
  basename() %>%
  gsub("_idoneidad_prob.*\\.tif", "", .) %>%
  gsub("_", " ", .)

resultados_todas_especies <- list()

for (i in seq_along(tifs_prob)) {
  sp_nombre <- nombres_especies[i]
  sp_name   <- gsub(" ", "_", sp_nombre)
  tif_path  <- tifs_prob[i]

  cat(sprintf("Procesando: '%s'...\n", sp_nombre))

  r <- rast(tif_path)
  ba_partidos_proj <- if (crs(r) != crs(ba_partidos)) project(ba_partidos, crs(r)) else ba_partidos

  r_crop <- tryCatch(crop(r, ba_partidos_proj), error = function(e) NULL)

  df_sp <- data.frame(
    partido = ba_partidos_proj$NAME_2,
    provincia = ba_partidos_proj$NAME_1,
    idoneidad_media = 0.0,
    idoneidad_maxima = 0.0,
    area_alta_idoneidad_km2 = 0.0,
    species = sp_nombre,
    stringsAsFactors = FALSE
  )

  if (!is.null(r_crop) && hasValues(r_crop)) {
    r_mask <- mask(r_crop, ba_partidos_proj)

    mean_vals <- terra::extract(r_mask, ba_partidos_proj, fun = mean, na.rm = TRUE)
    max_vals  <- terra::extract(r_mask, ba_partidos_proj, fun = max, na.rm = TRUE)

    max_val_raster <- max(values(r_mask), na.rm = TRUE)
    if (is.na(max_val_raster)) max_val_raster <- 0

    umbral <- if (max_val_raster > 10) 600 else 0.6

    r_high <- r_mask >= umbral
    area_raster <- cellSize(r_mask, unit = "km")
    r_high_area <- mask(area_raster, r_high, maskvalues = 0)

    area_vals <- terra::extract(r_high_area, ba_partidos_proj, fun = sum, na.rm = TRUE)

    df_sp$idoneidad_media  <- round(mean_vals[, 2], 4)
    df_sp$idoneidad_maxima <- round(max_vals[, 2], 4)
    df_sp$area_alta_idoneidad_km2 <- round(area_vals[, 2], 2)

    if (max_val_raster > 10) {
      df_sp$idoneidad_media  <- round(df_sp$idoneidad_media / 1000, 4)
      df_sp$idoneidad_maxima <- round(df_sp$idoneidad_maxima / 1000, 4)
    }

    df_sp$idoneidad_media[is.na(df_sp$idoneidad_media)] <- 0
    df_sp$idoneidad_maxima[is.na(df_sp$idoneidad_maxima)] <- 0
    df_sp$area_alta_idoneidad_km2[is.na(df_sp$area_alta_idoneidad_km2)] <- 0
  } else {
    cat("  [AVISO] Sin cobertura geográfica en Buenos Aires.\n")
  }

  resultados_todas_especies[[sp_name]] <- df_sp
}

df_consolidado <- bind_rows(resultados_todas_especies)

df_ranking_global <- df_consolidado %>%
  group_by(partido, provincia) %>%
  summarise(
    idoneidad_media_bgen   = round(mean(idoneidad_media), 4),
    idoneidad_max_bgen     = round(max(idoneidad_maxima), 4),
    area_alta_total_km2    = round(sum(area_alta_idoneidad_km2), 2),
    riqueza_especies_aptas = sum(idoneidad_media > 0.2),
    .groups = "drop"
  ) %>%
  mutate(
    indice_prioridad = round((idoneidad_media_bgen * 0.7) + ((riqueza_especies_aptas / length(tifs_prob)) * 0.3), 4)
  ) %>%
  arrange(desc(indice_prioridad), desc(area_alta_total_km2)) %>%
  mutate(rango_prioridad = row_number())

write_csv(df_ranking_global, file.path(RANK_DIR, "consolidado_ranking_partidos_EMmean.csv"))
cat("\nGuardado en:", file.path(RANK_DIR, "consolidado_ranking_partidos_EMmean.csv"), "\n")

# -----------------------------------------------------------------------------
# Comparación directa contra el ranking oficial (EMwmeanByTSS)
# -----------------------------------------------------------------------------
oficial <- read_csv(file.path(BASE_DIR, "outputs", "ranking_partidos", "consolidado_ranking_partidos.csv"),
                     show_col_types = FALSE) %>%
  select(partido, indice_prioridad_oficial = indice_prioridad, rango_oficial = rango_prioridad)

comparacion <- df_ranking_global %>%
  select(partido, indice_prioridad_emmean = indice_prioridad, rango_emmean = rango_prioridad) %>%
  left_join(oficial, by = "partido") %>%
  mutate(cambio_rango = rango_oficial - rango_emmean) %>%
  arrange(rango_oficial)

write_csv(comparacion, file.path(RANK_DIR, "comparacion_ranking_oficial_vs_emmean.csv"))
cat("Comparación guardada en:", file.path(RANK_DIR, "comparacion_ranking_oficial_vs_emmean.csv"), "\n\n")

cat("=== TOP 15 oficial vs EMmean ===\n")
print(comparacion %>% filter(rango_oficial <= 15) %>%
        select(partido, rango_oficial, rango_emmean, cambio_rango))

cat("\n=== Resumen de cambios de rango ===\n")
cat("Mediana |cambio_rango|:", median(abs(comparacion$cambio_rango), na.rm = TRUE), "\n")
cat("Máximo |cambio_rango|:", max(abs(comparacion$cambio_rango), na.rm = TRUE), "\n")
cat("Partidos con cambio de rango = 0:", sum(comparacion$cambio_rango == 0, na.rm = TRUE), "de", nrow(comparacion), "\n")
