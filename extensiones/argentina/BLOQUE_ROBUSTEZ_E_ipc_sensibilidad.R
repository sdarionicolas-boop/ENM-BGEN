# =============================================================================
#
#   BLOQUE ROBUSTEZ E (Argentina/tesis): sensibilidad del Ipc al modelo
#   corregido (thinning 10km + block CV), subconjunto de 4 especies.
#
#   NO es un Ipc oficial alternativo. Compara, para el MISMO subconjunto
#   de 4 especies, el ranking de partidos bonaerenses que resulta de usar
#   los rasters oficiales (CV aleatoria, 19 especies en el pipeline
#   original pero acá solo las 4 elegidas para que la comparación sea
#   apples-to-apples) vs. los rasters corregidos (BLOQUE_ROBUSTEZ_D).
#
#   El Ipc oficial que se reporta como principal en la tesis sigue siendo
#   el de BLOQUE8_partidos_bonaerenses.R (19 especies). Esto es evidencia
#   de sensibilidad, no un reemplazo.
#
# =============================================================================

library(terra)
library(geodata)
library(dplyr)
library(readr)
library(sf)

BASE_DIR       <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
VAR_DIR        <- file.path(BASE_DIR, "variables")
OUT_DIR_OFICIAL <- file.path(BASE_DIR, "outputs")
OUT_DIR_ROBUST  <- file.path(BASE_DIR, "outputs_robustez")

ESPECIES_SUBSET <- c("Celtis tala", "Schinus molle",
                     "Solanum pseudocapsicum", "Cortaderia selloana")
sp_names <- gsub(" ", "_", ESPECIES_SUBSET)

# -----------------------------------------------------------------------------
# 1. Límites de Buenos Aires (reutiliza caché de BLOQUE8 si ya existe)
# -----------------------------------------------------------------------------
gadm_arg <- gadm(country = "ARG", level = 2, path = VAR_DIR)
ba_partidos <- gadm_arg[gadm_arg$NAME_1 == "Buenos Aires", ]
cat("Partidos bonaerenses:", nrow(ba_partidos), "\n")

# -----------------------------------------------------------------------------
# 2. Función: calcular idoneidad/área alta por partido para un set de rasters
# -----------------------------------------------------------------------------
calcular_ranking <- function(dir_base, sp_names, etiqueta) {

  resultados <- list()

  for (sp_name in sp_names) {
    tif_path <- file.path(dir_base, sp_name, paste0(sp_name, "_idoneidad_prob.tif"))
    if (!file.exists(tif_path)) {
      cat("  [AVISO] No existe:", tif_path, "- se omite.\n")
      next
    }

    r <- rast(tif_path)
    ba_proj <- if (crs(r) != crs(ba_partidos)) project(ba_partidos, crs(r)) else ba_partidos

    r_crop <- crop(r, ba_proj)
    r_mask <- mask(r_crop, ba_proj)

    mean_vals <- terra::extract(r_mask, ba_proj, fun = mean, na.rm = TRUE)

    max_val_raster <- max(values(r_mask), na.rm = TRUE)
    if (is.na(max_val_raster)) max_val_raster <- 0
    umbral <- if (max_val_raster > 10) 600 else 0.6

    r_high <- r_mask >= umbral
    area_raster <- cellSize(r_mask, unit = "km")
    r_high_area <- mask(area_raster, r_high, maskvalues = 0)
    area_vals <- terra::extract(r_high_area, ba_proj, fun = sum, na.rm = TRUE)

    idoneidad_media <- round(mean_vals[, 2], 4)
    if (max_val_raster > 10) idoneidad_media <- round(idoneidad_media / 1000, 4)
    idoneidad_media[is.na(idoneidad_media)] <- 0

    area_alta <- round(area_vals[, 2], 2)
    area_alta[is.na(area_alta)] <- 0

    resultados[[sp_name]] <- data.frame(
      partido = ba_proj$NAME_2,
      idoneidad_media = idoneidad_media,
      area_alta_km2 = area_alta,
      species = sp_name
    )
  }

  df <- bind_rows(resultados)

  df %>%
    group_by(partido) %>%
    summarise(
      idoneidad_media_bgen   = round(mean(idoneidad_media), 4),
      area_alta_total_km2    = round(sum(area_alta_km2), 2),
      riqueza_especies_aptas = sum(idoneidad_media > 0.2),
      .groups = "drop"
    ) %>%
    mutate(
      indice_prioridad = round((idoneidad_media_bgen * 0.7) +
                                ((riqueza_especies_aptas / length(sp_names)) * 0.3), 4),
      modelo = etiqueta
    ) %>%
    arrange(desc(indice_prioridad)) %>%
    mutate(rango = row_number())
}

cat("\nCalculando ranking OFICIAL (subconjunto de 4 especies, CV aleatoria)...\n")
ranking_oficial <- calcular_ranking(OUT_DIR_OFICIAL, sp_names, "oficial_subset4")

cat("\nCalculando ranking CORREGIDO (thinning 10km + block CV)...\n")
ranking_robustez <- calcular_ranking(OUT_DIR_ROBUST, sp_names, "robustez_thin10km_blockCV")

# -----------------------------------------------------------------------------
# 3. Comparación
# -----------------------------------------------------------------------------
comparacion <- ranking_oficial %>%
  select(partido, indice_oficial = indice_prioridad, rango_oficial = rango) %>%
  inner_join(
    ranking_robustez %>% select(partido, indice_robustez = indice_prioridad, rango_robustez = rango),
    by = "partido"
  ) %>%
  mutate(cambio_rango = rango_oficial - rango_robustez) %>%
  arrange(rango_oficial)

cat("\n=== Top 15 oficial (subconjunto 4 especies) vs. rango bajo modelo corregido ===\n")
print(head(comparacion, 15))

top10_oficial   <- comparacion$partido[comparacion$rango_oficial <= 10]
top10_robustez  <- comparacion$partido[comparacion$rango_robustez <= 10]
overlap_top10   <- length(intersect(top10_oficial, top10_robustez))

correlacion <- cor(comparacion$rango_oficial, comparacion$rango_robustez, method = "spearman")

cat("\nSolapamiento Top 10 (oficial vs. corregido):", overlap_top10, "de 10\n")
cat("Correlación de Spearman entre ambos rankings:", round(correlacion, 3), "\n")

write_csv(comparacion, file.path("C:/Users/sdari/Desktop/BGEN/ENM_jacaranda/outputs_robustez",
                                  "ipc_sensibilidad_comparacion.csv"))

cat("\nBLOQUE ROBUSTEZ E completado.\n")
