# =============================================================================
#
#   BLOQUE 14c (Chile): Modelo de robustez — thinning espacial (10 km) +
#   validación cruzada por bloques espaciales (spatial block CV)
#
#   Responde a la observación de revisión: el modelo oficial usa
#   thinning a resolución de píxel (~1 km, no aborda autocorrelación
#   espacial) y validación cruzada aleatoria (80/20), lo que puede
#   inflar AUCroc/TSS. Este bloque reentrena el mismo ensemble con:
#   - Presencias sometidas a *thinning* espacial real (spThin,
#     distancia mínima 10 km): 2.091 -> 71 registros.
#   - Validación cruzada por bloques espaciales (biomod2
#     CV.strategy = "block", partición geográfica en cuadrantes,
#     sin fuga de información entre calibración y validación).
#
#   Todo lo demás (variables, algoritmos, pseudo-ausencias) se
#   mantiene IDÉNTICO al modelo oficial para que la comparación sea
#   limpia.
#
# =============================================================================

library(biomod2)
library(terra)
library(readr)
library(dplyr)

BASE_DIR   <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
DATA_DIR   <- file.path(BASE_DIR, "data")
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs")
MODELS_DIR <- "C:/Users/sdari/Documents"

SEED <- 123
PA_NB_REP <- 2
PA_NB_ABSENCES <- 3000

# -----------------------------------------------------------------------------
# 1. Presencias thinned (10 km) + stack ambiental (idéntico al oficial)
# -----------------------------------------------------------------------------
pres <- read_csv(file.path(DATA_DIR, "registros_thinned_10km.csv"), show_col_types = FALSE)
cat("Presencias (thinned 10km):", nrow(pres), "\n")

env_stack <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
cat("Variables ambientales:", paste(names(env_stack), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2. Formatear datos biomod2
# -----------------------------------------------------------------------------
set.seed(SEED)
bm_data <- BIOMOD_FormatingData(
  resp.name      = "Jubaea_chilensis_thin10km",
  resp.var       = rep(1, nrow(pres)),
  resp.xy        = pres[, c("lon", "lat")],
  expl.var       = env_stack,
  PA.nb.rep      = PA_NB_REP,
  PA.nb.absences = PA_NB_ABSENCES,
  PA.strategy    = "sre"
)

# -----------------------------------------------------------------------------
# 3. Validación cruzada por BLOQUES ESPACIALES (en vez de random 80/20)
# -----------------------------------------------------------------------------
cv_block <- bm_CrossValidation(
  bm.format = bm_data,
  strategy  = "block"
)

cat("\nEstructura de la tabla de CV por bloques:\n")
print(dim(cv_block))
print(head(cv_block))

# -----------------------------------------------------------------------------
# 4. Modelado con los mismos 5 algoritmos, usando el CV por bloques
# -----------------------------------------------------------------------------
setwd(MODELS_DIR)

bm_mod_block <- BIOMOD_Modeling(
  bm.format     = bm_data,
  modeling.id   = "thin10km_blockCV",
  models        = c("GLM", "GBM", "RF", "MAXNET", "XGBOOST"),
  CV.strategy   = "user.defined",
  CV.user.table = cv_block,
  var.import    = 3,
  metric.eval   = c("TSS", "AUCroc"),
  seed.val      = SEED
)

cat("\nModelado con thinning + block CV completado.\n")

# -----------------------------------------------------------------------------
# 5. Métricas: comparar contra el modelo oficial (random CV, sin thinning)
# -----------------------------------------------------------------------------
evals_block <- get_evaluations(bm_mod_block)

resumen_block <- evals_block %>%
  filter(metric.eval %in% c("AUCroc", "TSS")) %>%
  group_by(algo, metric.eval) %>%
  summarise(
    calibration_media = round(mean(calibration, na.rm = TRUE), 3),
    validation_media  = round(mean(validation,  na.rm = TRUE), 3),
    .groups = "drop"
  )

cat("\n=== Métricas por algoritmo (thinning + block CV) ===\n")
print(resumen_block)

resumen_general <- resumen_block %>%
  group_by(metric.eval) %>%
  summarise(
    media_validacion = round(mean(validation_media, na.rm = TRUE), 3),
    min_validacion   = round(min(validation_media,  na.rm = TRUE), 3),
    max_validacion   = round(max(validation_media,  na.rm = TRUE), 3)
  )

cat("\n=== Resumen general (thinning 10km + block CV) ===\n")
print(resumen_general)

cat("\n=== Comparación con el modelo oficial (random CV, sin thinning) ===\n")
cat("Oficial   -> AUCroc: 0.948 (0.925-0.974) | TSS: 0.741 (0.698-0.769)\n")
cat("Robustez  -> ver tabla arriba\n")

write_csv(resumen_block, file.path(OUT_DIR, "metricas_thinning_blockCV_por_algoritmo.csv"))
write_csv(resumen_general, file.path(OUT_DIR, "metricas_thinning_blockCV_resumen.csv"))

cat("\nBLOQUE 14c completado.\n")
