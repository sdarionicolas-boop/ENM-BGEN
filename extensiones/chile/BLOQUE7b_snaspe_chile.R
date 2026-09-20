# =============================================================================
#
#   BLOQUE 7b (Chile): Cruce con Áreas Protegidas (SNASPE + otras figuras)
#   ¿Cuánto del hábitat disponible para Jubaea chilensis ya está protegido?
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   Fuente de áreas protegidas: WDPA (World Database on Protected Areas,
#   UNEP-WCMC/IUCN), filtrado a Chile y a la AOI de Jubaea chilensis,
#   obtenido vía Google Earth Engine (WCMC/WDPA/current/polygons).
#   Incluye Parques Nacionales, Reservas Nacionales y Monumentos
#   Naturales (SNASPE/CONAF), Santuarios de la Naturaleza (MMA),
#   sitios Ramsar y Reservas de la Biósfera UNESCO-MAB.
#
#   Pregunta:
#   Del hábitat "Compatible" (vegetación natural) de alta idoneidad
#   climática, ¿qué fracción está dentro de un área oficialmente
#   protegida y qué fracción está desprotegida y en riesgo?
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
# 1. Cargar el raster de cruce (idoneidad x uso de suelo) ya generado
#    para Colección 2 (resultado oficial)
# -----------------------------------------------------------------------------
cruce <- rast(file.path(MB_DIR, "Col2_2024", "Jubaea_chilensis_cruce.tif"))

# -----------------------------------------------------------------------------
# 2. Cargar y rasterizar las áreas protegidas (WDPA) a la misma grilla
# -----------------------------------------------------------------------------
pas <- vect(file.path(VAR_DIR, "wdpa_chile_aoi.geojson"))
cat("Áreas protegidas cargadas:", nrow(pas), "\n")
print(table(pas$DESIG_ENG))

pas_proj <- project(pas, crs(cruce))
pas_rast <- rasterize(pas_proj, cruce, field = 1, background = 0)
pas_rast <- mask(pas_rast, cruce)

# -----------------------------------------------------------------------------
# 3. Cruce: hábitat Compatible de alta idoneidad (código 41) x protegido/no
# -----------------------------------------------------------------------------
compatible_alta <- cruce == 41

area_km2_rast <- cellSize(cruce, unit = "km")
km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

tab <- table(
  compatible_alta = values(compatible_alta)[, 1],
  protegido       = values(pas_rast)[, 1]
)
print(tab)

km2_compatible_protegido    <- tab["TRUE", "1"] * km2_por_pixel
km2_compatible_desprotegido <- tab["TRUE", "0"] * km2_por_pixel
km2_compatible_total        <- km2_compatible_protegido + km2_compatible_desprotegido

resumen <- tibble(
  categoria = c("Compatible (alta idoneidad) - Protegido", "Compatible (alta idoneidad) - Desprotegido"),
  km2 = round(c(km2_compatible_protegido, km2_compatible_desprotegido), 1),
  porcentaje = round(c(km2_compatible_protegido, km2_compatible_desprotegido) / km2_compatible_total * 100, 1)
)

cat("\n=== Hábitat Compatible de Alta Idoneidad: protegido vs desprotegido ===\n")
print(resumen)
cat("\nTotal hábitat Compatible de alta idoneidad:", round(km2_compatible_total, 1), "km²\n")

write_csv(resumen, file.path(MB_DIR, "Col2_2024", "snaspe_proteccion_resumen.csv"))

# -----------------------------------------------------------------------------
# 4. Detalle: qué áreas protegidas específicas concentran más hábitat compatible
# -----------------------------------------------------------------------------
cat("\n=== Superficie de hábitat compatible por área protegida (top 15) ===\n")
pas_id_rast <- rasterize(pas_proj, cruce, field = "NAME", background = NA)
ids <- values(pas_id_rast)[, 1]
compat_vals <- values(compatible_alta)[, 1]

df_detalle <- tibble(nombre = ids, es_compatible = compat_vals) %>%
  filter(!is.na(nombre), es_compatible) %>%
  count(nombre, name = "n_pixeles") %>%
  mutate(km2 = round(n_pixeles * km2_por_pixel, 1)) %>%
  arrange(desc(km2))

print(head(df_detalle, 15))
write_csv(df_detalle, file.path(MB_DIR, "Col2_2024", "snaspe_detalle_por_area.csv"))

cat("\nBLOQUE 7b completado.\n")
