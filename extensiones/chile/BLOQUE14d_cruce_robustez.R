# =============================================================================
#
#   BLOQUE 14d: Cruce de robustez -- ¿los porcentajes del 6.2 cambian
#   si usamos el modelo thinned + block CV en vez del oficial?
#
#   Convierte la afirmación "las limitaciones no invalidan los
#   hallazgos centrales" de retórica a evidencia empírica.
#
# =============================================================================

library(biomod2)
library(terra)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs")
MB_DIR     <- file.path(OUT_DIR, "mapbiomas")
MODELS_DIR <- "C:/Users/sdari/Documents"

# -----------------------------------------------------------------------------
# 1. Cargar el modelo thinned + block CV (ya entrenado en BLOQUE14c)
# -----------------------------------------------------------------------------
setwd(BASE_DIR)
thin_models_out <- "Jubaea.chilensis.thin10km/Jubaea.chilensis.thin10km.thin10km_blockCV.models.out"
bm_mod <- get(load(thin_models_out))

# Construir el ensemble ponderado por TSS (mismo criterio que el oficial)
bm_em <- BIOMOD_EnsembleModeling(
  bm.mod             = bm_mod,
  em.by              = "all",
  em.algo            = c("EMmean", "EMwmean"),
  metric.select      = c("TSS"),
  metric.select.thresh = c(0),
  metric.eval        = c("TSS", "AUCroc"),
  var.import         = 0
)

env <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))

bm_proj <- BIOMOD_EnsembleForecasting(
  bm.em          = bm_em,
  new.env        = env,
  proj.name      = "CURRENT_thin10km",
  models.chosen  = "all",
  metric.binary  = "TSS",
  metric.filter  = "TSS"
)

# -----------------------------------------------------------------------------
# 2. Cargar el raster de idoneidad thinned y clasificarlo con el mismo umbral
# -----------------------------------------------------------------------------
prob_path <- "Jubaea.chilensis.thin10km/proj_CURRENT_thin10km/proj_CURRENT_thin10km_Jubaea.chilensis.thin10km_ensemble.tif"
cat("Cargando proyección desde:", prob_path, "\n")
prob_thin <- rast(prob_path)
if (nlyr(prob_thin) > 1) prob_thin <- prob_thin[[1]]
prob_thin <- prob_thin / 1000  # biomod2 escala 0-1000

breaks <- c(0, 0.2, 0.4, 0.6, 1.0)
clase_thin <- classify(prob_thin, rcl = cbind(breaks[-length(breaks)], breaks[-1], 1:4),
                        include.lowest = TRUE)

# -----------------------------------------------------------------------------
# 3. Cruce con MapBiomas Chile Col2 (2024) -- misma lógica que el 6.2 oficial
# -----------------------------------------------------------------------------
mb_raw <- rast(file.path(VAR_DIR, "chile_coverage2_2024.tif"))

cod_compatible   <- c(3, 59, 60, 67, 11, 12, 63, 66)
cod_restaurable  <- c(9, 15, 21)
cod_incompatible <- c(18, 24, 25)
cod_excluido     <- c(0, 23, 29, 33, 34, 61)

rcl_mb <- rbind(
  cbind(cod_compatible,   1),
  cbind(cod_restaurable,  2),
  cbind(cod_incompatible, 3),
  cbind(cod_excluido,     NA)
)
mb_cat <- classify(mb_raw, rcl_mb, others = NA)

# Alinear con el grid de idoneidad
rcl_idon <- matrix(c(1, 1, 2, 1, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
idon_3cat <- classify(clase_thin, rcl_idon)

mb_cropped <- crop(mb_cat, idon_3cat)
mb_aligned <- resample(mb_cropped, idon_3cat, method = "near")
mb_aligned <- mask(mb_aligned, idon_3cat)

cruce_thin <- idon_3cat * 10 + mb_aligned

area_km2 <- cellSize(cruce_thin, unit = "km")
km2_por_px <- mean(values(area_km2), na.rm = TRUE)

freq_cruce <- freq(cruce_thin, usenames = FALSE)

etiquetas <- c(
  "41" = "Alta - Compatible", "42" = "Alta - Restaurable", "43" = "Alta - Incompatible",
  "31" = "Moderada - Compatible", "32" = "Moderada - Restaurable", "33" = "Moderada - Incompatible",
  "11" = "Baja - Compatible", "12" = "Baja - Restaurable", "13" = "Baja - Incompatible"
)

resumen <- freq_cruce %>%
  as_tibble() %>%
  rename(codigo = value, n_pixeles = count) %>%
  filter(!is.na(codigo)) %>%
  mutate(codigo_str = as.character(codigo),
         descripcion = etiquetas[codigo_str],
         km2 = round(n_pixeles * km2_por_px, 0)) %>%
  filter(!is.na(descripcion))

cat("=== Cruce oficial vs. robustez (thinning + block CV) ===\n")
cat("\nSuperficie por categoría en el modelo THINNED:\n")
print(resumen %>% select(descripcion, km2) %>% arrange(desc(km2)))

# Extraer los 3 valores clave de alta idoneidad
km2_alta <- resumen %>% filter(codigo %in% c(41, 42, 43))
total_alta <- sum(km2_alta$km2)
km2_alta <- km2_alta %>% mutate(pct = round(100 * km2 / total_alta, 1))

cat("\n=== SÓLO clase de ALTA idoneidad ===\n")
print(km2_alta %>% select(descripcion, km2, pct))

cat("\n=== Comparación con el modelo oficial ===\n")
cat("OFICIAL   (random CV): Compatible 71.9% | Restaurable 5.1% | Incompatible 22.9% | Total 20,760 km2\n")
cat(sprintf("THINNED   (block CV):  Compatible %.1f%% | Restaurable %.1f%% | Incompatible %.1f%% | Total %d km2\n",
    km2_alta$pct[km2_alta$codigo == 41],
    ifelse(any(km2_alta$codigo == 42), km2_alta$pct[km2_alta$codigo == 42], 0),
    ifelse(any(km2_alta$codigo == 43), km2_alta$pct[km2_alta$codigo == 43], 0),
    total_alta))

write_csv(km2_alta, file.path(MB_DIR, "Col2_2024", "cruce_robustez_thinned.csv"))
cat("\nGuardado en cruce_robustez_thinned.csv\n")
