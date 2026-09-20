# =============================================================================
#
#   BLOQUE 11 (Chile): Proyección a futuro — CMIP6 2050 (SSP2-4.5)
#   Jubaea chilensis — ¿cómo cambia la idoneidad climática al 2041-2060?
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  1.0.0  |  Septiembre 2026
#
#   No reentrena el modelo: usa el ensemble ya ajustado (biomod2) y
#   proyecta sobre un nuevo stack ambiental construido con variables
#   bioclimáticas CMIP6 (modelo MPI-ESM1-2-HR, SSP2-4.5, 2041-2060,
#   WorldClim v2.1 downscaling) para las mismas 7 variables retenidas
#   por VIF en el modelo actual (bio02, bio03, bio08, bio14, bio15,
#   slope, aspect). Topografía (slope/aspect) no cambia con el clima.
#
# =============================================================================

library(terra)
library(biomod2)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs")
MODELS_DIR <- "C:/Users/sdari/Documents"

# -----------------------------------------------------------------------------
# 1. Construir el stack ambiental futuro (mismas 7 variables, misma AOI)
# -----------------------------------------------------------------------------
env_current <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
aoi <- ext(env_current)

cmip6_raw <- rast(file.path(VAR_DIR, "climate", "wc2.1_30s", "wc2.1_30s_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif"))
names(cmip6_raw) <- sprintf("bio%02d", 1:19)

cmip6_aoi <- crop(cmip6_raw, aoi)
cmip6_aoi <- resample(cmip6_aoi, env_current[[c("bio02","bio03","bio08","bio14","bio15")]], method = "bilinear")

env_future <- c(
  cmip6_aoi[[c("bio02", "bio03", "bio08", "bio14", "bio15")]],
  env_current[[c("slope", "aspect")]]   # topografía: no cambia
)
names(env_future) <- c("bio02", "bio03", "bio08", "bio14", "bio15", "slope", "aspect")

writeRaster(env_future, file.path(VAR_DIR, "env_stack_future_2050.tif"), overwrite = TRUE)
cat("Stack futuro (2050, SSP2-4.5) armado:", nlyr(env_future), "variables\n")

# -----------------------------------------------------------------------------
# 2. Cargar el ensemble ya ajustado y proyectar sobre el clima futuro
# -----------------------------------------------------------------------------
setwd(MODELS_DIR)
myBiomodEM <- get(load("Jubaea.chilensis/Jubaea.chilensis.Jubaea_chilensis_current.ensemble.models.out"))

myBiomodEF_future <- BIOMOD_EnsembleForecasting(
  bm.em          = myBiomodEM,
  new.env        = env_future,
  proj.name      = "FUTURE_2050",
  models.chosen  = "all",
  metric.binary  = "TSS",
  metric.filter  = "TSS"
)

cat("Proyección futura completada.\n")

# -----------------------------------------------------------------------------
# 3. Clasificar (mismos umbrales que el mapa actual) y comparar
# -----------------------------------------------------------------------------
prob_future_path <- file.path(MODELS_DIR, "Jubaea.chilensis", "proj_FUTURE_2050",
                               "proj_FUTURE_2050_Jubaea.chilensis_ensemble.tif")
prob_future <- rast(prob_future_path)
# EMwmeanByTSS suele ser la primera capa; usar la de mean/wmean por TSS si hay varias
if (nlyr(prob_future) > 1) prob_future <- prob_future[[1]]
prob_future <- prob_future / 1000  # biomod2 escala 0-1000

prob_current <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_prob.tif"))
clase_current <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_clase.tif"))

prob_future <- resample(prob_future, prob_current, method = "bilinear")

breaks <- c(0, 0.2, 0.4, 0.6, 1.0)
clase_future <- classify(prob_future, rcl = cbind(breaks[-length(breaks)], breaks[-1], 1:4),
                          include.lowest = TRUE)

dir.create(file.path(OUT_DIR, "cmip6"), showWarnings = FALSE)
writeRaster(prob_future, file.path(OUT_DIR, "cmip6", "Jubaea_chilensis_idoneidad_prob_2050.tif"), overwrite = TRUE)
writeRaster(clase_future, file.path(OUT_DIR, "cmip6", "Jubaea_chilensis_idoneidad_clase_2050.tif"), overwrite = TRUE)

area_km2_rast <- cellSize(clase_current, unit = "km")
km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

tab_current <- table(factor(values(clase_current)[,1], levels = 1:4))
tab_future  <- table(factor(values(clase_future)[,1], levels = 1:4))

resumen <- tibble(
  clase = c("Insustentable", "Bajo", "Moderado", "Alto"),
  km2_actual = round(as.numeric(tab_current) * km2_por_pixel, 1),
  km2_2050   = round(as.numeric(tab_future) * km2_por_pixel, 1)
) %>%
  mutate(cambio_km2 = round(km2_2050 - km2_actual, 1),
         cambio_pct = round(100 * (km2_2050 - km2_actual) / km2_actual, 1))

cat("\n=== Cambio de superficie de idoneidad climática: actual vs. 2050 (SSP2-4.5) ===\n")
print(resumen)

write_csv(resumen, file.path(OUT_DIR, "cmip6", "cambio_idoneidad_2050.csv"))

prob_media_actual <- mean(values(prob_current), na.rm = TRUE)
prob_media_futura  <- mean(values(prob_future), na.rm = TRUE)
cat(sprintf("\nIdoneidad media actual: %.3f | Idoneidad media 2050: %.3f (cambio: %+.1f%%)\n",
            prob_media_actual, prob_media_futura,
            100 * (prob_media_futura - prob_media_actual) / prob_media_actual))

cat("\nBLOQUE 11 completado.\n")
