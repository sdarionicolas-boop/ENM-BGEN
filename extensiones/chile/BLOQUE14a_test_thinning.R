library(spThin)
library(readr)
library(dplyr)

pres <- read_csv("data/registros_unificados_geocod.csv", show_col_types = FALSE)
cat("Registros originales:", nrow(pres), "\n\n")

pres_sp <- pres %>% mutate(SPEC = especie)

distancias <- c(1, 2, 5, 10, 20, 30, 50)

for (d in distancias) {
  set.seed(42)
  thinned <- thin(
    loc.data = pres_sp,
    lat.col = "lat", long.col = "lon", spec.col = "SPEC",
    thin.par = d, reps = 5,
    locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE, verbose = FALSE
  )
  n_reps <- sapply(thinned, nrow)
  cat(sprintf("Distancia %3d km -> registros retenidos: min=%d, max=%d, media=%.0f (de %d originales, %.1f%% retenido)\n",
              d, min(n_reps), max(n_reps), mean(n_reps), nrow(pres), 100*mean(n_reps)/nrow(pres)))
}
