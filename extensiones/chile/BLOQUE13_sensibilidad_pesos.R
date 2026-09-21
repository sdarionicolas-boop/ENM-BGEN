# =============================================================================
#
#   BLOQUE 13 (Chile): Análisis de sensibilidad de los pesos del índice
#   de prioridad de conservación (ranking de palmares)
#
#   Responde a la observación de revisión: los pesos 50/30/20 (hábitat
#   compatible / desprotección / riesgo de incendio) no estaban
#   justificados. Se recalcula el ranking con esquemas alternativos
#   para verificar si el orden de prioridad es robusto al esquema de
#   ponderación elegido.
#
# =============================================================================

library(dplyr)
library(readr)

ranking <- read_csv(
  "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs/mapbiomas/Col2_2024/ranking_palmares.csv",
  show_col_types = FALSE
)

norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

ranking <- ranking %>%
  mutate(
    hab_norm          = norm01(km2_compatible),
    desprot_norm      = norm01(100 - pct_protegido),
    riesgo_fuego_norm = norm01(pct_quemado)
  )

esquemas <- list(
  "50/30/20 (original)" = c(0.50, 0.30, 0.20),
  "40/40/20"             = c(0.40, 0.40, 0.20),
  "60/20/20"             = c(0.60, 0.20, 0.20),
  "34/33/33 (igual peso)"= c(1/3, 1/3, 1/3),
  "70/15/15 (solo hábitat)" = c(0.70, 0.15, 0.15)
)

resultado <- ranking %>% select(nombre_tentativo)

for (nombre_esq in names(esquemas)) {
  w <- esquemas[[nombre_esq]]
  indice <- w[1]*ranking$hab_norm + w[2]*ranking$desprot_norm + w[3]*ranking$riesgo_fuego_norm
  rango <- rank(-indice, ties.method = "first")
  resultado[[nombre_esq]] <- rango
}

cat("=== Rango de cada palmar bajo distintos esquemas de ponderación ===\n\n")
print(as.data.frame(resultado), row.names = FALSE)

write_csv(resultado, "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs/mapbiomas/Col2_2024/sensibilidad_pesos_ranking.csv")

cat("\n=== Top-2 bajo cada esquema ===\n")
for (nombre_esq in names(esquemas)) {
  top2 <- ranking$nombre_tentativo[order(resultado[[nombre_esq]])][1:2]
  cat(nombre_esq, ":", paste(top2, collapse = " | "), "\n")
}
