# =============================================================================
#
#   BLOQUE ROBUSTEZ C (Argentina/tesis): reentrena el ensemble con
#   thinning 10km + validación cruzada por bloques espaciales, subconjunto
#   de 4 especies representativas del BGEN.
#
#   Replica extensiones/chile/BLOQUE14c_modelo_thinned_blockCV.R. Todo lo
#   demás (algoritmos, variables, parámetros de pseudo-ausencias) se
#   mantiene IDÉNTICO al pipeline oficial (ENM_BGEN_pipeline.R, BLOQUE 3)
#   para que la comparación sea limpia. Única diferencia deliberada:
#   metric.eval usa "ROC" (no "AUCroc") porque así lo define el pipeline
#   oficial de Argentina en BLOQUE 3.
#
# =============================================================================

library(biomod2)
library(terra)
library(readr)
library(dplyr)

BASE_DIR   <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
DATA_DIR   <- file.path(BASE_DIR, "data")
ROBUST_DIR <- file.path(DATA_DIR, "robustez")
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs_robustez")
MODELS_DIR <- "C:/Users/sdari/Documents"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 123
PA_NB_REP      <- 2
PA_NB_ABSENCES <- 3000
PA_SRE_QUANT   <- 0.05
BM_MODELS      <- c("GLM", "GBM", "RF", "MAXNET", "XGBOOST")

ESPECIES_SUBSET <- c("Celtis tala", "Schinus molle",
                     "Solanum pseudocapsicum", "Cortaderia selloana")

env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
cat("Variables ambientales:", paste(names(env_final), collapse = ", "), "\n")

setwd(MODELS_DIR)

resumen_todas <- list()

for (sp in ESPECIES_SUBSET) {

  sp_name  <- gsub(" ", "_", sp)
  resp_name <- paste0(sp_name, "_thin10km")

  cat("\n", rep("=", 60), "\n")
  cat("Modelando (robustez):", sp, "\n")
  cat(rep("=", 60), "\n")

  pres_sp <- read_csv(file.path(ROBUST_DIR, paste0(sp_name, "_thinned10km.csv")),
                       show_col_types = FALSE)
  cat("  Presencias (thinned 10km):", nrow(pres_sp), "\n")

  tryCatch({

    # -------------------------------------------------------------------
    # 1. Formatear datos biomod2 (idéntico al oficial salvo el input)
    # -------------------------------------------------------------------
    set.seed(SEED)
    bm_data <- BIOMOD_FormatingData(
      resp.var       = rep(1, nrow(pres_sp)),
      expl.var       = env_final,
      resp.xy        = as.data.frame(pres_sp[, c("longitude", "latitude")]),
      resp.name      = resp_name,
      PA.nb.rep      = PA_NB_REP,
      PA.nb.absences = PA_NB_ABSENCES,
      PA.strategy    = "sre",
      PA.sre.quant   = PA_SRE_QUANT
    )

    # -------------------------------------------------------------------
    # 2. Validación cruzada por BLOQUES ESPACIALES (en vez de random 80/20)
    # -------------------------------------------------------------------
    cv_block <- bm_CrossValidation(
      bm.format = bm_data,
      strategy  = "block"
    )

    # -------------------------------------------------------------------
    # 3. Modelado con los mismos 5 algoritmos, usando el CV por bloques
    # -------------------------------------------------------------------
    bm_mod_block <- BIOMOD_Modeling(
      bm.format     = bm_data,
      modeling.id   = "thin10km_blockCV",
      models        = BM_MODELS,
      CV.strategy   = "user.defined",
      CV.user.table = cv_block,
      var.import    = 3,
      metric.eval   = c("TSS", "ROC"),
      seed.val      = SEED
    )

    cat("  ✓ Modelado thinning+blockCV completado:", sp, "\n")

    # -------------------------------------------------------------------
    # 4. Métricas por algoritmo
    # -------------------------------------------------------------------
    evals_block <- get_evaluations(bm_mod_block)

    resumen_sp <- evals_block %>%
      filter(metric.eval %in% c("ROC", "TSS")) %>%
      mutate(metric.eval = ifelse(metric.eval == "ROC", "AUCroc", metric.eval)) %>%
      group_by(algo, metric.eval) %>%
      summarise(
        calibration_media = round(mean(calibration, na.rm = TRUE), 3),
        validation_media  = round(mean(validation,  na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      mutate(species = sp)

    print(resumen_sp)
    resumen_todas[[sp]] <- resumen_sp

  }, error = function(e) {
    cat("  ✗ ERROR en", sp, ":", conditionMessage(e), "\n")
  })

}

tabla_final <- bind_rows(resumen_todas)
write_csv(tabla_final, file.path(OUT_DIR, "metricas_robustez_thinning_blockCV.csv"))
cat("\nGuardado:", file.path(OUT_DIR, "metricas_robustez_thinning_blockCV.csv"), "\n")

cat("\nBLOQUE ROBUSTEZ C completado.\n")
