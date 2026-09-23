# =============================================================================
#
#   BLOQUE ROBUSTEZ B (Argentina/tesis): dataset final thinned a 10 km,
#   subconjunto de 4 especies representativas del BGEN
#
#   Replica extensiones/chile/BLOQUE14b_thinning_final.R. Genera un CSV
#   thinned por especie más un combinado, para usar en BLOQUE_ROBUSTEZ_C.
#
# =============================================================================

library(spThin)
library(readr)
library(dplyr)

BASE_DIR <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
DATA_DIR <- file.path(BASE_DIR, "data")
ROBUST_DIR <- file.path(DATA_DIR, "robustez")
dir.create(ROBUST_DIR, showWarnings = FALSE, recursive = TRUE)

ESPECIES_SUBSET <- c("Celtis tala", "Schinus molle",
                     "Solanum pseudocapsicum", "Cortaderia selloana")

pres <- read_csv(file.path(BASE_DIR, "presencias_todas_especies.csv"), show_col_types = FALSE) %>%
  filter(species %in% ESPECIES_SUBSET)

resultado_combinado <- list()

for (sp in ESPECIES_SUBSET) {

  sp_pres <- pres %>% filter(species == sp)
  sp_name <- gsub(" ", "_", sp)

  cat("\n===", sp, "===\n")

  set.seed(42)
  thinned <- thin(
    loc.data = sp_pres,
    lat.col = "latitude", long.col = "longitude", spec.col = "species",
    thin.par = 10, reps = 20,
    locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE, verbose = FALSE
  )

  n_reps <- sapply(thinned, nrow)
  cat("Registros retenidos por repetición (20 reps, thin.par=10km):\n")
  print(n_reps)

  best <- thinned[[which.max(n_reps)]]
  cat("Repetición elegida: máximo", nrow(best), "registros\n")

  pres_thin_sp <- best %>%
    transmute(species = sp, longitude = Longitude, latitude = Latitude)

  out_path <- file.path(ROBUST_DIR, paste0(sp_name, "_thinned10km.csv"))
  write_csv(pres_thin_sp, out_path)
  cat("Guardado en", out_path, ":", nrow(pres_thin_sp), "registros\n")

  resultado_combinado[[sp]] <- pres_thin_sp
}

combinado <- bind_rows(resultado_combinado)
write_csv(combinado, file.path(ROBUST_DIR, "presencias_thin_robustez_10km.csv"))

cat("\nResumen final (N por especie, 10km thinning):\n")
print(combinado %>% count(species))

cat("\nBLOQUE ROBUSTEZ B completado.\n")
