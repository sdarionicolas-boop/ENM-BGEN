# =============================================================================
#   EJEMPLO MÍNIMO REPRODUCIBLE — una especie, de punta a punta
#
#   Corre el núcleo del pipeline ENM-BGEN (variables → VIF → biomod2 →
#   ensemble EMmean → mapa de idoneidad) sobre Celtis tala, y compara las dos
#   estrategias de validación cruzada: aleatoria (la del pipeline oficial) vs.
#   por bloques espaciales (la que usa la prueba de robustez).
#
#   Qué es y qué no es:
#   - Es una versión REDUCIDA para entender y verificar el método en pocos
#     minutos. NO reproduce los resultados de las 19 especies.
#   - Variables a 10 arcmin (~18 km) en vez de 30 arcsec (~1 km); topografía
#     reducida a elevación; 1 set de 500 pseudo-ausencias en vez de 2 de 3.000.
#   - Los registros de presencia (61, ya rarificados a 10 km) vienen en
#     datos/. Las capas ambientales se descargan de WorldClim al correr
#     (no se redistribuyen en este repo).
#
#   Uso:   Rscript ejemplo_una_especie.R      (desde cualquier carpeta)
#          o en RStudio: abrir el archivo, "Set Working Directory > To Source
#          File Location" y correrlo completo.
#   Requisitos: R >= 4.2 y los paquetes listados abajo (+ conexión a internet
#          la primera vez, para WorldClim).
#   Salidas: ejemplo_minimo/salida/  (ignorada por git)
# =============================================================================

paquetes <- c("terra", "geodata", "usdm", "biomod2", "ENMeval", "dplyr", "readr", "tidyr")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Instalá los paquetes que faltan: install.packages(c(",
                         paste(sprintf('"%s"', faltan), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(terra); library(geodata); library(usdm)
  library(biomod2); library(dplyr); library(readr)
})

# --- Ubicación del script ----------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
DIR <- if (length(f)) dirname(normalizePath(f)) else getwd()
if (!file.exists(file.path(DIR, "datos", "celtis_tala_presencias_10km.csv")))
  stop("No encuentro datos/celtis_tala_presencias_10km.csv. Corré el script desde ",
       "la carpeta ejemplo_minimo/ o con Rscript.")

CACHE  <- file.path(DIR, "cache");  dir.create(CACHE,  showWarnings = FALSE)
SALIDA <- file.path(DIR, "salida"); dir.create(SALIDA, showWarnings = FALSE)

# --- Parámetros (los del pipeline, reducidos) ---------------------------------
options(width = 120)
SEED           <- 123
RES_ARCMIN     <- 10        # resolución de las capas (pipeline: 0.5 = 30 arcsec)
PA_NB_ABSENCES <- 500       # pseudo-ausencias (pipeline: 2 sets × 3.000)
PA_SRE_QUANT   <- 0.05
UMBRAL_VIF     <- 5
BM_MODELS      <- c("GLM", "GBM", "RF", "MAXNET", "XGBOOST")
CLASES_RCL <- matrix(c(0, 0.2, 1,  0.2, 0.4, 2,  0.4, 0.6, 3,  0.6, 1, 4),
                     ncol = 3, byrow = TRUE)   # 1 Insustentable … 4 Alto

# --- 1. Presencias y área de estudio -------------------------------------------
pres <- read_csv(file.path(DIR, "datos", "celtis_tala_presencias_10km.csv"),
                 show_col_types = FALSE)
cat("Presencias:", nrow(pres), "(rarificadas a 10 km)\n")

aoi <- ext(min(pres$longitude) - 1, max(pres$longitude) + 1,
           min(pres$latitude)  - 1, max(pres$latitude)  + 1)

# --- 2. Variables ambientales (WorldClim 2.1, 10 arcmin) ------------------------
bio  <- crop(worldclim_global(var = "bio",  res = RES_ARCMIN, path = CACHE), aoi)
elev <- crop(worldclim_global(var = "elev", res = RES_ARCMIN, path = CACHE), aoi)
names(bio)  <- sprintf("bio%02d", 1:19)
names(elev) <- "elev"
env_all <- c(bio, elev)

# Selección por colinealidad (VIF ≤ 5), igual que el pipeline
vs      <- vifstep(env_all, th = UMBRAL_VIF)
env     <- exclude(env_all, vs)
cat("Variables retenidas por VIF ≤", UMBRAL_VIF, ":", paste(names(env), collapse = ", "), "\n")

# --- 3. Modelado con dos estrategias de validación cruzada ----------------------
setwd(SALIDA)   # biomod2 guarda sus carpetas en el directorio de trabajo
xy <- as.data.frame(pres[, c("longitude", "latitude")])

set.seed(SEED)
bm_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(xy)), expl.var = env, resp.xy = xy,
  resp.name = "Celtis.tala",
  PA.nb.rep = 1, PA.nb.absences = PA_NB_ABSENCES,
  PA.strategy = "sre", PA.sre.quant = PA_SRE_QUANT
)

modelar <- function(estrategia) {
  if (estrategia == "random") {
    # como el pipeline oficial: 3 repeticiones, 80/20 al azar
    BIOMOD_Modeling(bm_data, modeling.id = "random", models = BM_MODELS,
                    CV.strategy = "random", CV.nb.rep = 3, CV.perc = 0.8,
                    metric.eval = c("TSS", "ROC"), var.import = 0, seed.val = SEED)
  } else {
    # bloques espaciales: validación geográficamente independiente
    tabla <- bm_CrossValidation(bm_data, strategy = "block")
    BIOMOD_Modeling(bm_data, modeling.id = "block", models = BM_MODELS,
                    CV.strategy = "user.defined", CV.user.table = tabla,
                    metric.eval = c("TSS", "ROC"), var.import = 0, seed.val = SEED)
  }
}

resultados <- list(); metricas <- list()
for (est in c("random", "block")) {
  cat("\n=== Validación cruzada:", est, "===\n")
  resultados[[est]] <- modelar(est)
  metricas[[est]] <- get_evaluations(resultados[[est]]) %>%
    filter(metric.eval %in% c("TSS", "ROC", "AUCroc")) %>%
    mutate(metric.eval = ifelse(metric.eval == "ROC", "AUCroc", metric.eval)) %>%
    group_by(algo, metric.eval) %>%
    summarise(calibracion = round(mean(calibration, na.rm = TRUE), 3),
              validacion  = round(mean(validation,  na.rm = TRUE), 3), .groups = "drop") %>%
    mutate(cv = est)
}
met <- bind_rows(metricas)
write_csv(met, file.path(SALIDA, "metricas_random_vs_block.csv"))

cat("\n=== TSS de validación por algoritmo: aleatoria vs. por bloques ===\n")
comp <- met %>% filter(metric.eval == "TSS") %>%
  select(algo, cv, calibracion, validacion) %>%
  tidyr::pivot_wider(names_from = cv, values_from = c(calibracion, validacion))
print(as.data.frame(comp), row.names = FALSE)
cat(sprintf("\nTSS medio de validación — aleatoria: %.3f | por bloques: %.3f\n",
            mean(met$validacion[met$metric.eval == "TSS" & met$cv == "random"]),
            mean(met$validacion[met$metric.eval == "TSS" & met$cv == "block"])))

# --- 4. Ensemble EMmean y mapa de idoneidad (modelos con validación por bloques) -
# Sin filtro por umbral de métrica: con validación espacial ningún algoritmo
# alcanza el TSS ≥ 0,7 que exige el pipeline oficial, y el ensemble quedaría vacío.
bm_ens <- BIOMOD_EnsembleModeling(
  bm.mod = resultados[["block"]], models.chosen = "all", em.by = "all",
  em.algo = "EMmean", metric.select = "all", var.import = 0
)
bm_proj <- BIOMOD_EnsembleForecasting(
  bm.em = bm_ens, proj.name = "actual", new.env = env, models.chosen = "all"
)
preds <- get_predictions(bm_proj)
prob  <- preds[[grep("EMmean", names(preds))[1]]]
if (global(prob, "max", na.rm = TRUE)[1, 1] > 1) prob <- prob / 1000
clase <- classify(prob, CLASES_RCL)

writeRaster(prob,  file.path(SALIDA, "celtis_tala_idoneidad_prob.tif"),  overwrite = TRUE)
writeRaster(clase, file.path(SALIDA, "celtis_tala_idoneidad_clase.tif"), overwrite = TRUE)

png(file.path(SALIDA, "celtis_tala_idoneidad.png"), width = 1400, height = 1100, res = 150)
plot(prob, main = "Celtis tala — idoneidad (EMmean, ejemplo reducido)", range = c(0, 1))
points(pres$longitude, pres$latitude, pch = 21, bg = "white", cex = 0.6)
dev.off()

cat("\nListo. Resultados en:", SALIDA, "\n",
    " - metricas_random_vs_block.csv\n",
    " - celtis_tala_idoneidad_prob.tif / _clase.tif\n",
    " - celtis_tala_idoneidad.png\n")
