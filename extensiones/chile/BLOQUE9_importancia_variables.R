# =============================================================================
#
#   BLOQUE 9 (Chile): Importancia de variables ambientales
#   Jubaea chilensis — extraída del modelo ya entrenado (biomod2)
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   No reentrena nada: el modelo ya fue corrido con var.import = 3
#   (BLOQUE 3 de ENM_BGEN_pipeline.R). Este bloque solo carga el
#   objeto guardado y extrae/grafica la importancia ya calculada.
#
# =============================================================================

library(biomod2)
library(dplyr)
library(readr)
library(ggplot2)

OUT_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs"
MODELS_DIR <- "C:/Users/sdari/Documents"

setwd(MODELS_DIR)
myBiomodModelOut <- get(load("Jubaea.chilensis/Jubaea.chilensis.Jubaea_chilensis_current.models.out"))

varimp <- get_variables_importance(myBiomodModelOut)
cat("Estructura de importancia de variables:\n")
print(head(varimp))

# Promedio por variable, a través de todos los algoritmos/runs/PA
resumen_varimp <- varimp %>%
  as_tibble() %>%
  group_by(expl.var) %>%
  summarise(
    importancia_media = round(mean(var.imp, na.rm = TRUE), 4),
    importancia_sd     = round(sd(var.imp, na.rm = TRUE), 4)
  ) %>%
  arrange(desc(importancia_media))

cat("\n=== Importancia media de variables (todos los algoritmos) ===\n")
print(resumen_varimp)

write_csv(resumen_varimp, file.path(OUT_DIR, "importancia_variables_jubaea.csv"))

# Gráfico
p <- ggplot(resumen_varimp, aes(x = reorder(expl.var, importancia_media), y = importancia_media)) +
  geom_col(fill = "#2E7D32") +
  geom_errorbar(aes(ymin = pmax(importancia_media - importancia_sd, 0),
                     ymax = importancia_media + importancia_sd), width = 0.2) +
  coord_flip() +
  labs(
    title = "Importancia de variables ambientales — Jubaea chilensis",
    subtitle = "Promedio ± DE entre 5 algoritmos (GLM, GBM, RF, MAXNET, XGBoost)",
    x = NULL, y = "Importancia (correlación con permutación)"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(OUT_DIR, "importancia_variables_jubaea.png"), p, width = 8, height = 5, dpi = 150)

cat("\nGráfico guardado en outputs/importancia_variables_jubaea.png\n")
cat("BLOQUE 9 completado.\n")
