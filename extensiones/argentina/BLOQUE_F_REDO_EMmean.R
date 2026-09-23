# =============================================================================
#
#   BLOQUE 4 REDO: re-extraer el ensemble como EMmean (no ponderado) en vez
#   de EMwmeanByTSS, para las 19 especies del pipeline oficial.
#
#   NO reentrena nada: carga los .models.out ya ajustados (BLOQUE 3 oficial,
#   validación aleatoria) y solo recalcula el objeto de ensemble + proyección,
#   que es barato comparado con reentrenar biomod2 desde cero.
#
#   Motivo del cambio: BLOQUE_ROBUSTEZ_A-D mostró que RF tiene calibración
#   perfecta (1,0) en las 4 especies testeadas bajo el esquema oficial, pero
#   es de los peores algoritmos bajo validación espacial rigurosa (block CV).
#   EMwmeanByTSS pondera por el TSS de validación aleatoria -- inflado
#   precisamente para el modelo que peor generaliza -- así que le da MÁS
#   peso al algoritmo menos confiable. EMmean no depende de esa métrica.
#
#   Escribe a OUT_DIR_EMMEAN (carpeta paralela), NO pisa los rasters
#   oficiales en outputs/ hasta que se revisen y se decida promoverlos.
#
# =============================================================================

library(dplyr)
library(readr)
library(terra)
library(biomod2)

BASE_DIR       <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
MODELS_DIR     <- "C:/Users/sdari/Documents"
VAR_DIR        <- file.path(BASE_DIR, "variables")
OUT_DIR_EMMEAN <- file.path(BASE_DIR, "outputs_emmean")
dir.create(OUT_DIR_EMMEAN, showWarnings = FALSE, recursive = TRUE)

EM_METRIC_TSS <- 0.7
EM_METRIC_AUC <- 0.9

CLASES_RCL <- matrix(c(
  0.0, 0.2, 1,
  0.2, 0.4, 2,
  0.4, 0.6, 3,
  0.6, 1.0, 4
), ncol = 3, byrow = TRUE)

env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
pres_thin <- read_csv(file.path(BASE_DIR, "presencias_thin.csv"), show_col_types = FALSE)
especies  <- unique(pres_thin$species)

cat("Especies a procesar:", length(especies), "\n")

setwd(MODELS_DIR)

resumen <- list()

for (sp in especies) {

  cat("\n", rep("=", 60), "\n")
  cat("BLOQUE 4 REDO (EMmean):", sp, "\n")

  sp_name  <- gsub(" ", "_", sp)
  sp_punto <- gsub(" ", ".", sp)
  sp_out   <- file.path(OUT_DIR_EMMEAN, sp_name)
  dir.create(sp_out, showWarnings = FALSE, recursive = TRUE)

  resultado <- tryCatch({

    ruta_modelo <- file.path(
      getwd(), sp_punto,
      paste0(sp_punto, ".", sp_name, "_current.models.out")
    )

    if (!file.exists(ruta_modelo)) stop("Modelo no encontrado: ", ruta_modelo)

    bm_mod <- load(ruta_modelo)
    bm_mod <- get(bm_mod)

    bm_ens <- BIOMOD_EnsembleModeling(
      bm.mod               = bm_mod,
      models.chosen        = "all",
      em.by                = "all",
      metric.select        = c("TSS", "AUCroc"),
      metric.select.thresh = c(EM_METRIC_TSS, EM_METRIC_AUC),
      em.algo               = c("EMmean", "EMwmean"),
      var.import            = 3
    )

    bm_proj <- BIOMOD_EnsembleForecasting(
      bm.em         = bm_ens,
      proj.name     = "CURRENT_EMMEAN_REDO",
      new.env       = env_final,
      models.chosen = "all",
      metric.filter = c("TSS", "AUCroc")
    )

    ens_preds <- get_predictions(bm_proj)

    idx_mean <- grep("EMmean", names(ens_preds))
    if (length(idx_mean) == 0) stop("No se encontró EMmean en las predicciones")
    idx_mean <- idx_mean[1]

    idx_wmean <- grep("EMwmean", names(ens_preds))
    idx_wmean <- if (length(idx_wmean) > 0) idx_wmean[1] else NA

    ens_mean <- ens_preds[[idx_mean]]
    ens_mean_prob <- if (max(values(ens_mean), na.rm = TRUE) > 1) ens_mean / 1000 else ens_mean
    names(ens_mean_prob) <- sp_name

    ens_class <- classify(ens_mean_prob, CLASES_RCL)

    writeRaster(ens_mean_prob, file.path(sp_out, paste0(sp_name, "_idoneidad_prob_EMmean.tif")), overwrite = TRUE)
    writeRaster(ens_class,     file.path(sp_out, paste0(sp_name, "_idoneidad_clase_EMmean.tif")), overwrite = TRUE)

    # Para comparación: cuánto cambia el promedio de idoneidad vs el EMwmean oficial
    diff_media <- NA
    if (!is.na(idx_wmean)) {
      ens_wmean <- ens_preds[[idx_wmean]]
      ens_wmean_prob <- if (max(values(ens_wmean), na.rm = TRUE) > 1) ens_wmean / 1000 else ens_wmean
      diff_media <- round(mean(values(ens_mean_prob), na.rm = TRUE) - mean(values(ens_wmean_prob), na.rm = TRUE), 4)
    }

    cat("  OK -- media EMmean:", round(mean(values(ens_mean_prob), na.rm = TRUE), 4),
        " | diferencia vs EMwmean oficial:", diff_media, "\n")

    list(species = sp, status = "OK",
         media_EMmean = round(mean(values(ens_mean_prob), na.rm = TRUE), 4),
         diff_vs_EMwmean = diff_media)

  }, error = function(e) {
    cat("  ERROR en", sp, ":", conditionMessage(e), "\n")
    list(species = sp, status = "ERROR", media_EMmean = NA, diff_vs_EMwmean = NA)
  })

  resumen[[sp]] <- resultado
}

resumen_df <- bind_rows(resumen)
write_csv(resumen_df, file.path(OUT_DIR_EMMEAN, "resumen_redo_emmean.csv"))
cat("\n\nBLOQUE 4 REDO completado. Resumen en:", file.path(OUT_DIR_EMMEAN, "resumen_redo_emmean.csv"), "\n")
print(resumen_df)
