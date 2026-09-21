library(spThin)
library(readr)
library(dplyr)

pres <- read_csv("data/registros_unificados_geocod.csv", show_col_types = FALSE) %>%
  mutate(SPEC = especie)

set.seed(42)
thinned <- thin(
  loc.data = pres,
  lat.col = "lat", long.col = "lon", spec.col = "SPEC",
  thin.par = 10, reps = 20,
  locs.thinned.list.return = TRUE,
  write.files = FALSE, write.log.file = FALSE, verbose = FALSE
)

n_reps <- sapply(thinned, nrow)
cat("Registros retenidos por repetición (20 reps, thin.par=10km):\n")
print(n_reps)

# Elegir la repetición con más registros (maximiza poder estadístico
# dentro del criterio de 10 km de distancia mínima)
best <- thinned[[which.max(n_reps)]]
cat("\nRepetición elegida: máximo", nrow(best), "registros\n")

pres_thin <- best %>%
  transmute(especie = "Jubaea chilensis", lat = Latitude, lon = Longitude)

write_csv(pres_thin, "data/registros_thinned_10km.csv")
cat("Guardado en data/registros_thinned_10km.csv:", nrow(pres_thin), "registros\n")

cat("\nRango espacial:\n")
cat("Lat:", round(range(pres_thin$lat), 3), "\n")
cat("Lon:", round(range(pres_thin$lon), 3), "\n")
