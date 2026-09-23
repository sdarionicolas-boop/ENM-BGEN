# =============================================================================
#
#   BLOQUE ROBUSTEZ A (Argentina/tesis): prueba de distancias de thinning
#   espacial (spThin), subconjunto de 4 especies representativas del BGEN
#
#   Replica el ejercicio de extensiones/chile/BLOQUE14a_test_thinning.R
#   sobre un subconjunto de especies de Argentina, en vez de las 19, para
#   no rehacer todos los modelos oficiales (ver decisión documentada en
#   la tesis, sección "Prueba de robustez").
#
#   Subconjunto elegido (criterios objetivos, ver metricas_todas_especies.csv
#   y presencias_thin.csv del pipeline oficial):
#     - Celtis tala             (N=133, mínimo; TSS oficial máximo 0.855)
#     - Schinus molle           (N=369, máximo)
#     - Solanum pseudocapsicum  (N=302; mayor gap RF calibración-validación, 0.400)
#     - Cortaderia selloana     (N=281; mayor hábitat perdido absoluto del paper)
#
# =============================================================================

library(spThin)
library(readr)
library(dplyr)

BASE_DIR <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
DATA_DIR <- file.path(BASE_DIR, "data")

ESPECIES_SUBSET <- c("Celtis tala", "Schinus molle",
                     "Solanum pseudocapsicum", "Cortaderia selloana")

pres <- read_csv(file.path(BASE_DIR, "presencias_todas_especies.csv"), show_col_types = FALSE)
cat("Registros totales (19 especies):", nrow(pres), "\n")

pres_sub <- pres %>% filter(species %in% ESPECIES_SUBSET)
cat("Registros del subconjunto de robustez:", nrow(pres_sub), "\n")
print(table(pres_sub$species))

distancias <- c(1, 2, 5, 10, 20, 30, 50)

for (sp in ESPECIES_SUBSET) {

  sp_pres <- pres_sub %>% filter(species == sp)
  cat("\n=== ", sp, " (N original =", nrow(sp_pres), ") ===\n")

  for (d in distancias) {
    set.seed(42)
    thinned <- thin(
      loc.data = sp_pres,
      lat.col = "latitude", long.col = "longitude", spec.col = "species",
      thin.par = d, reps = 5,
      locs.thinned.list.return = TRUE,
      write.files = FALSE, write.log.file = FALSE, verbose = FALSE
    )
    n_reps <- sapply(thinned, nrow)
    cat(sprintf("  Distancia %3d km -> retenidos: min=%d, max=%d, media=%.0f (%.1f%% retenido)\n",
                d, min(n_reps), max(n_reps), mean(n_reps), 100 * mean(n_reps) / nrow(sp_pres)))
  }
}

cat("\nBLOQUE ROBUSTEZ A completado.\n")
