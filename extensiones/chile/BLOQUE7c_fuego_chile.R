# =============================================================================
#
#   BLOQUE 7c (Chile): Cruce con MapBiomas Fuego Colección 1
#   ¿Cuánto del hábitat disponible para Jubaea chilensis ya se quemó?
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   Fuente: MapBiomas Fuego Chile Colección 1 — frecuencia de área
#   quemada acumulada 2013-2025 (30m). Descarga directa:
#   storage.googleapis.com/mapbiomas-public/initiatives/chile/fire/
#   collection1/frequency_burned_v1/mapbiomas_fire_chile_col1_
#   frequency_burned_2013_2025.tif
#
#   Conecta con una amenaza ya documentada en la literatura para
#   Jubaea chilensis: incendios recurrentes (sección 4 de la Memoria
#   Técnica). Usa el TERCER producto oficial de MapBiomas Chile
#   elegible según el reglamento (junto a land cover colección 2 y,
#   de forma exploratoria en BLOQUE7b, el cruce con áreas protegidas).
#
# =============================================================================

library(terra)
library(dplyr)
library(readr)

BASE_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR  <- file.path(BASE_DIR, "variables")
OUT_DIR  <- file.path(BASE_DIR, "outputs")
MB_DIR   <- file.path(OUT_DIR, "mapbiomas")

# -----------------------------------------------------------------------------
# 1. Cargar capas: cruce (idoneidad x uso suelo), protección, y fuego
# -----------------------------------------------------------------------------
cruce <- rast(file.path(MB_DIR, "Col2_2024", "Jubaea_chilensis_cruce.tif"))

pas <- vect(file.path(VAR_DIR, "wdpa_chile_aoi.geojson"))
pas_proj <- project(pas, crs(cruce))
pas_rast <- rasterize(pas_proj, cruce, field = 1, background = 0)
pas_rast <- mask(pas_rast, cruce)

fuego_raw <- rast(file.path(VAR_DIR, "fire_frequency_2013_2025.tif"))
fuego_a <- resample(crop(fuego_raw, cruce), cruce, method = "near")
fuego_a <- mask(fuego_a, cruce)
quemado_alguna_vez <- fuego_a > 0

cat("Resumen de frecuencia de incendios en la AOI:\n")
print(summary(values(fuego_a)))

# -----------------------------------------------------------------------------
# 2. Hábitat Compatible de alta idoneidad: cuánto se quemó, según protección
# -----------------------------------------------------------------------------
compatible_alta <- cruce == 41
area_km2_rast <- cellSize(cruce, unit = "km")
km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

df <- tibble(
  compatible = values(compatible_alta)[, 1],
  protegido  = values(pas_rast)[, 1],
  quemado    = values(quemado_alguna_vez)[, 1],
  frecuencia = values(fuego_a)[, 1]
) %>% filter(compatible)

resumen_fuego <- df %>%
  mutate(estado_proteccion = ifelse(protegido == 1, "Protegido", "Desprotegido")) %>%
  group_by(estado_proteccion) %>%
  summarise(
    km2_total     = round(n() * km2_por_pixel, 1),
    km2_quemado   = round(sum(quemado, na.rm = TRUE) * km2_por_pixel, 1),
    pct_quemado   = round(100 * mean(quemado, na.rm = TRUE), 1),
    frecuencia_media_quemas = round(mean(frecuencia[quemado], na.rm = TRUE), 2)
  )

cat("\n=== Hábitat Compatible de Alta Idoneidad: incendios 2013-2025, por estado de protección ===\n")
print(resumen_fuego)

total_km2 <- round(nrow(df) * km2_por_pixel, 1)
total_quemado_km2 <- round(sum(df$quemado, na.rm = TRUE) * km2_por_pixel, 1)
cat(sprintf("\nTotal hábitat compatible: %.1f km2. Quemado al menos una vez (2013-2025): %.1f km2 (%.1f%%)\n",
            total_km2, total_quemado_km2, 100 * total_quemado_km2 / total_km2))

write_csv(resumen_fuego, file.path(MB_DIR, "Col2_2024", "fuego_resumen.csv"))

# -----------------------------------------------------------------------------
# 3. Detalle por área protegida (palmares específicos)
# -----------------------------------------------------------------------------
pas_id_rast <- rasterize(pas_proj, cruce, field = "NAME", background = NA)
lv <- levels(pas_id_rast)[[1]]

df2 <- tibble(
  id         = values(pas_id_rast)[, 1],
  compatible = values(compatible_alta)[, 1],
  quemado    = values(quemado_alguna_vez)[, 1]
) %>%
  filter(!is.na(id), compatible) %>%
  left_join(lv, by = c("id" = "ID")) %>%
  group_by(nombre = NAME) %>%
  summarise(
    km2_compatible = round(n() * km2_por_pixel, 1),
    km2_quemado    = round(sum(quemado, na.rm = TRUE) * km2_por_pixel, 1),
    pct_quemado    = round(100 * mean(quemado, na.rm = TRUE), 1)
  ) %>%
  arrange(desc(km2_compatible))

cat("\n=== Palmares/áreas protegidas: hábitat compatible y % quemado ===\n")
print(df2, n = 15)
write_csv(df2, file.path(MB_DIR, "Col2_2024", "fuego_por_area_protegida.csv"))

cat("\nBLOQUE 7c completado.\n")
