# =============================================================================
#
#   BLOQUE 8 (Chile): Ranking de palmares por prioridad de conservación
#   Jubaea chilensis — agrupamiento espacial de registros de presencia
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   La literatura (PMC9370131) identifica 6 poblaciones/palmares
#   principales: Culimo, Petorca, Ocoa, Viña del Mar/Valparaíso,
#   Cocalán y Candelaria, concentrando el 96% de los individuos.
#   Como no se dispone de polígonos oficiales para cada palmar, se
#   agrupan los registros de presencia (GBIF, ya filtrados a rango
#   natural) mediante clustering espacial (k-means, k=6) y se
#   etiqueta cada clúster con el nombre de la población más cercana
#   documentada en la literatura, cuando la distancia lo permite.
#
#   Para cada clúster se calcula un índice de prioridad de
#   conservación que combina disponibilidad de hábitat, protección
#   legal y exposición a incendios — mismo enfoque de índice
#   compuesto que BLOQUE8_partidos_bonaerenses.R (Argentina).
#
# =============================================================================

library(terra)
library(dplyr)
library(readr)

BASE_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
DATA_DIR <- file.path(BASE_DIR, "data")
VAR_DIR  <- file.path(BASE_DIR, "variables")
OUT_DIR  <- file.path(BASE_DIR, "outputs")
MB_DIR   <- file.path(OUT_DIR, "mapbiomas")

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Cargar presencias y agrupar en 6 clústeres espaciales (k-means)
# -----------------------------------------------------------------------------
pres <- read_csv(file.path(DATA_DIR, "registros_unificados_geocod.csv"), show_col_types = FALSE)
cat("Registros de presencia:", nrow(pres), "\n")

km <- kmeans(pres[, c("lon", "lat")], centers = 6, nstart = 50)
pres$cluster <- km$cluster

# Poblaciones documentadas en la literatura (PMC9370131), para referencia
# de etiquetado por proximidad al centroide de cada clúster.
poblaciones_lit <- tibble(
  nombre_lit = c("Ocoa", "Cocalán", "Viña del Mar/Valparaíso", "Candelaria", "Petorca", "Culimo"),
  individuos = c(70308, 35500, 7200, 1900, 1300, 204)
)

centroides <- pres %>%
  group_by(cluster) %>%
  summarise(lon_c = mean(lon), lat_c = mean(lat), n_registros = n())

cat("\nCentroides de clústeres:\n")
print(centroides)

# -----------------------------------------------------------------------------
# 2. Para cada clúster: polígono envolvente (buffer sobre convex hull) y
#    extracción de hábitat compatible, protección y fuego
# -----------------------------------------------------------------------------
cruce   <- rast(file.path(MB_DIR, "Col2_2024", "Jubaea_chilensis_cruce.tif"))
pas     <- vect(file.path(VAR_DIR, "wdpa_chile_aoi.geojson"))
pas_proj <- project(pas, crs(cruce))
pas_rast <- rasterize(pas_proj, cruce, field = 1, background = 0)
pas_rast <- mask(pas_rast, cruce)

fuego_raw <- rast(file.path(VAR_DIR, "fire_frequency_2013_2025.tif"))
fuego_a <- resample(crop(fuego_raw, cruce), cruce, method = "near")
fuego_a <- mask(fuego_a, cruce)
quemado <- fuego_a > 0

compatible_alta <- cruce == 41
area_km2_rast <- cellSize(cruce, unit = "km")
km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

resumen_clusters <- list()

for (cl in sort(unique(pres$cluster))) {
  pts <- pres %>% filter(cluster == cl)
  v <- vect(pts, geom = c("lon", "lat"), crs = "EPSG:4326")

  # Envolvente convexa + buffer de 3 km para capturar hábitat circundante
  hull <- convHull(v)
  hull_buf <- buffer(hull, width = 3000)

  cruce_c      <- crop(cruce, hull_buf, mask = TRUE)
  compat_c     <- cruce_c == 41
  pas_c        <- crop(pas_rast, hull_buf, mask = TRUE)
  quemado_c    <- crop(quemado, hull_buf, mask = TRUE)

  km2_compatible <- sum(values(compat_c), na.rm = TRUE) * km2_por_pixel
  km2_protegido  <- sum(values(compat_c)[values(compat_c) == 1] *
                           values(pas_c)[values(compat_c) == 1], na.rm = TRUE) * km2_por_pixel
  pct_protegido  <- ifelse(km2_compatible > 0, round(100 * km2_protegido / km2_compatible, 1), NA)

  vals_compat <- values(compat_c)[, 1]
  vals_quem   <- values(quemado_c)[, 1]
  idx_compat  <- which(vals_compat == 1)
  pct_quemado <- ifelse(length(idx_compat) > 0,
                         round(100 * mean(vals_quem[idx_compat], na.rm = TRUE), 1), NA)

  centroide <- centroides %>% filter(cluster == cl)

  # Etiquetar por población documentada más cercana (si está a <60 km)
  dist_km <- NA_real_
  nombre_asignado <- paste0("Clúster ", cl)
  # (las coordenadas exactas de cada población de la literatura no están
  #  publicadas como tabla; el etiquetado definitivo se deja pendiente de
  #  verificación manual -- ver nota en el resumen)

  resumen_clusters[[as.character(cl)]] <- tibble(
    cluster            = cl,
    n_registros        = centroide$n_registros,
    lon_centroide      = round(centroide$lon_c, 3),
    lat_centroide      = round(centroide$lat_c, 3),
    km2_compatible     = round(km2_compatible, 1),
    km2_protegido      = round(km2_protegido, 1),
    pct_protegido      = pct_protegido,
    pct_quemado        = pct_quemado
  )

  cat(sprintf("\nClúster %d: %d registros, centroide (%.3f, %.3f)\n",
              cl, centroide$n_registros, centroide$lon_c, centroide$lat_c))
  cat(sprintf("  Hábitat compatible: %.1f km2 | Protegido: %.1f%% | Quemado: %.1f%%\n",
              km2_compatible, pct_protegido, pct_quemado))
}

tabla <- bind_rows(resumen_clusters)

# -----------------------------------------------------------------------------
# 3. Índice de prioridad de conservación
#    Combina: mucho hábitat compatible (bien) + poco protegido (prioridad) +
#    mucho quemado (urgencia). Normalizado 0-1 por variable.
# -----------------------------------------------------------------------------
norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

tabla <- tabla %>%
  mutate(
    hab_norm      = norm01(km2_compatible),
    desprot_norm  = norm01(100 - pct_protegido),
    riesgo_fuego_norm = norm01(pct_quemado),
    indice_prioridad = round(0.5 * hab_norm + 0.3 * desprot_norm + 0.2 * riesgo_fuego_norm, 3)
  ) %>%
  arrange(desc(indice_prioridad)) %>%
  mutate(rango_prioridad = row_number())

cat("\n\n=== Ranking de clústeres por índice de prioridad de conservación ===\n")
print(tabla %>% select(rango_prioridad, cluster, n_registros, km2_compatible, pct_protegido, pct_quemado, indice_prioridad))

write_csv(tabla, file.path(MB_DIR, "Col2_2024", "ranking_palmares.csv"))

cat("\nNOTA: Los clústeres se etiquetan por número, no por nombre de población\n")
cat("(Ocoa/Cocalán/Petorca/etc.), porque las coordenadas exactas de cada\n")
cat("población en PMC9370131 no están tabuladas en el texto disponible.\n")
cat("El clúster con centroide más cercano a la Reserva Palmas de Cocalán\n")
cat("puede identificarse cruzando lon_centroide/lat_centroide con el WDPA.\n")

cat("\nBLOQUE 8 (ranking palmares) completado.\n")
