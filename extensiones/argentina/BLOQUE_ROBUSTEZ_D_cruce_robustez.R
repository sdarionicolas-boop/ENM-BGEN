# =============================================================================
#
#   BLOQUE ROBUSTEZ D (Argentina/tesis): cruce oficial vs. robustez +
#   generación de rasters de idoneidad corregidos (thinning 10km + block CV)
#   para el subconjunto de 4 especies.
#
#   Parte 1: tabla comparativa de métricas (oficial vs. corregido), a
#   partir de metricas_todas_especies.csv (oficial) y
#   metricas_robustez_thinning_blockCV.csv (BLOQUE_ROBUSTEZ_C).
#
#   Parte 2: para cada especie, construye el ensemble (EMwmean, mismo
#   criterio que el oficial) sobre el modelo thinned+blockCV y proyecta
#   sobre el territorio -> outputs_robustez/{especie}/{especie}_idoneidad_prob.tif
#   Estos rasters corregidos alimentan BLOQUE_ROBUSTEZ_E (sensibilidad del Ipc).
#
# =============================================================================

library(biomod2)
library(terra)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR_OFICIAL  <- file.path(BASE_DIR, "outputs")
OUT_DIR    <- file.path(BASE_DIR, "outputs_robustez")
MODELS_DIR <- "C:/Users/sdari/Documents"

# Umbral 0: a diferencia del pipeline oficial (TSS>=0.7, AUCroc>=0.9), acá se
# incluyen todos los modelos sin filtrar por calidad, igual que en
# extensiones/chile/BLOQUE14d_cruce_robustez.R (linea 35) -- el objetivo es
# diagnostico, y bajo block CV ningun algoritmo alcanza TSS 0.7 en ninguna de
# las 4 especies (ver metricas_robustez_thinning_blockCV.csv), asi que un
# umbral de produccion dejaria el ensemble vacio.
EM_METRIC_TSS <- 0
EM_METRIC_AUC <- 0

CLASES_RCL <- matrix(c(
  0.0, 0.2, 1,
  0.2, 0.4, 2,
  0.4, 0.6, 3,
  0.6, 1.0, 4
), ncol = 3, byrow = TRUE)

ESPECIES_SUBSET <- c("Celtis tala", "Schinus molle",
                     "Solanum pseudocapsicum", "Cortaderia selloana")

# -----------------------------------------------------------------------------
# Parte 1: tabla comparativa de métricas oficial vs. robustez
# -----------------------------------------------------------------------------
oficial <- read_csv(file.path(OUT_DIR_OFICIAL, "metricas_todas_especies.csv"), show_col_types = FALSE) %>%
  filter(species %in% ESPECIES_SUBSET) %>%
  rename(validation_oficial = validation_media, calibration_oficial = calibration_media)

robustez <- read_csv(file.path(OUT_DIR, "metricas_robustez_thinning_blockCV.csv"), show_col_types = FALSE) %>%
  rename(validation_robustez = validation_media, calibration_robustez = calibration_media)

comparacion <- oficial %>%
  select(species, algo, metric.eval, calibration_oficial, validation_oficial) %>%
  inner_join(
    robustez %>% select(species, algo, metric.eval, calibration_robustez, validation_robustez),
    by = c("species", "algo", "metric.eval")
  ) %>%
  mutate(diferencia_validation = round(validation_robustez - validation_oficial, 3))

cat("=== Comparación oficial vs. robustez (por especie y algoritmo) ===\n")
print(comparacion)

resumen_por_especie <- comparacion %>%
  filter(metric.eval == "TSS") %>%
  group_by(species) %>%
  summarise(
    TSS_oficial_medio   = round(mean(validation_oficial), 3),
    TSS_robustez_medio  = round(mean(validation_robustez), 3),
    caida_TSS           = round(TSS_oficial_medio - TSS_robustez_medio, 3),
    .groups = "drop"
  )

cat("\n=== Resumen TSS por especie (oficial vs. robustez) ===\n")
print(resumen_por_especie)

write_csv(comparacion, file.path(OUT_DIR, "cruce_robustez_metricas.csv"))
write_csv(resumen_por_especie, file.path(OUT_DIR, "cruce_robustez_resumen_por_especie.csv"))

# -----------------------------------------------------------------------------
# Parte 2: ensemble + proyección con el modelo corregido, por especie
# -----------------------------------------------------------------------------
env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
setwd(MODELS_DIR)

for (sp in ESPECIES_SUBSET) {

  sp_name   <- gsub(" ", "_", sp)
  resp_name <- paste0(sp_name, "_thin10km")
  resp_punto <- gsub("_", ".", resp_name)
  sp_out    <- file.path(OUT_DIR, sp_name)
  dir.create(sp_out, showWarnings = FALSE, recursive = TRUE)

  cat("\n=== Ensemble + proyección (robustez):", sp, "===\n")

  tryCatch({

    ruta_modelo <- file.path(getwd(), resp_punto,
                              paste0(resp_punto, ".thin10km_blockCV.models.out"))
    if (!file.exists(ruta_modelo)) stop("Modelo no encontrado: ", ruta_modelo)

    bm_mod <- get(load(ruta_modelo))

    bm_ens <- BIOMOD_EnsembleModeling(
      bm.mod               = bm_mod,
      models.chosen        = "all",
      em.by                = "all",
      metric.select        = c("TSS", "AUCroc"),
      metric.select.thresh = c(EM_METRIC_TSS, EM_METRIC_AUC),
      em.algo              = c("EMmean", "EMwmean"),
      var.import            = 0
    )

    bm_proj <- BIOMOD_EnsembleForecasting(
      bm.em         = bm_ens,
      proj.name     = "CURRENT_robustez",
      new.env       = env_final,
      models.chosen = "all",
      metric.filter = c("TSS", "AUCroc")
    )

    ens_preds <- get_predictions(bm_proj)
    idx <- grep("EMwmeanByTSS", names(ens_preds))
    if (length(idx) == 0) idx <- grep("EMwmean", names(ens_preds))[1]
    ens_wmean <- ens_preds[[idx]]

    ens_prob <- if (max(values(ens_wmean), na.rm = TRUE) > 1) ens_wmean / 1000 else ens_wmean
    names(ens_prob) <- sp_name
    ens_class <- classify(ens_prob, CLASES_RCL)

    writeRaster(ens_prob, file.path(sp_out, paste0(sp_name, "_idoneidad_prob.tif")), overwrite = TRUE)
    writeRaster(ens_class, file.path(sp_out, paste0(sp_name, "_idoneidad_clase.tif")), overwrite = TRUE)

    cat("  ✓ Raster corregido guardado para", sp, "\n")

  }, error = function(e) {
    cat("  ✗ ERROR en ensemble/proyección de", sp, ":", conditionMessage(e), "\n")
  })
}

cat("\nBLOQUE ROBUSTEZ D completado.\n")
