# =============================================================================
#
#   BLOQUE 10 (Chile): Validación cruzada con sitios conocidos
#   ¿El modelo predice alta idoneidad exactamente donde se sabe que
#   hay palmares reales de Jubaea chilensis?
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   Complementa la validación estadística interna (AUCroc/TSS, 80/20)
#   con una validación externa e independiente: se extrae la
#   distribución de clases de idoneidad dentro de las áreas
#   protegidas creadas específicamente para conservar Jubaea
#   chilensis (WDPA), sitios que el modelo nunca vio como parte del
#   ajuste de hiperparámetros o del umbral de clasificación.
#
# =============================================================================

library(terra)
library(dplyr)
library(readr)

BASE_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR  <- file.path(BASE_DIR, "variables")
OUT_DIR  <- file.path(BASE_DIR, "outputs")

idon_clase <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_clase.tif"))
idon_prob  <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_prob.tif"))

pas <- vect(file.path(VAR_DIR, "wdpa_chile_aoi.geojson"))
pas_proj <- project(pas, crs(idon_clase))

# Sitios de validación: áreas protegidas creadas específicamente para
# Jubaea chilensis, o que son palmares documentados en la literatura
# (independientes de los registros de presencia usados para entrenar).
sitios_validacion <- c(
  "Palmas de Cocalán",
  "Área de Palma Chilena de Monte Aranda",
  "Palmar El Salto",
  "La Campana - Peñuelas",
  "Bosque Fray Jorge (RB)"
)

etiquetas_clase <- c("1" = "Insustentable", "2" = "Bajo", "3" = "Moderado", "4" = "Alto")

resultados <- list()

for (sitio in sitios_validacion) {
  poly <- pas_proj[pas_proj$NAME == sitio, ]
  if (nrow(poly) == 0) {
    cat("No encontrado:", sitio, "\n")
    next
  }

  clase_vals <- extract(idon_clase, poly)[, 2]
  prob_vals  <- extract(idon_prob, poly)[, 2]

  tabla_clase <- table(factor(clase_vals, levels = 1:4))
  pct_clase <- round(100 * tabla_clase / sum(tabla_clase), 1)

  resultados[[sitio]] <- tibble(
    sitio = sitio,
    n_pixeles = sum(tabla_clase),
    prob_media = round(mean(prob_vals, na.rm = TRUE), 3),
    pct_insustentable = pct_clase["1"],
    pct_bajo = pct_clase["2"],
    pct_moderado = pct_clase["3"],
    pct_alto = pct_clase["4"]
  )

  cat(sprintf("\n%s (%d px):\n", sitio, sum(tabla_clase)))
  cat(sprintf("  Idoneidad media (prob): %.3f\n", mean(prob_vals, na.rm = TRUE)))
  cat(sprintf("  %% Alto: %.1f%% | Moderado: %.1f%% | Bajo: %.1f%% | Insustentable: %.1f%%\n",
              pct_clase["4"], pct_clase["3"], pct_clase["2"], pct_clase["1"]))
}

tabla_final <- bind_rows(resultados) %>% arrange(desc(prob_media))
cat("\n\n=== Resumen de validación cruzada ===\n")
print(tabla_final)

write_csv(tabla_final, file.path(OUT_DIR, "mapbiomas", "Col2_2024", "validacion_sitios_conocidos.csv"))

# Comparación de referencia: idoneidad media en la AOI completa (fuera de estos sitios)
prob_aoi_media <- mean(values(idon_prob), na.rm = TRUE)
cat(sprintf("\nIdoneidad media promedio en TODA la AOI (referencia): %.3f\n", prob_aoi_media))
cat("Si los sitios de validación superan claramente esta referencia, el modelo\n")
cat("identifica correctamente hábitat real más allá del ajuste estadístico interno.\n")

cat("\nBLOQUE 10 completado.\n")
