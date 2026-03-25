# =============================================================================
#
#   PIPELINE DE MODELADO DE NICHO ECOLÓGICO – BGEN / UNAJ
#   Banco de Germoplasma de Especies Nativas
#   Universidad Nacional Arturo Jauretche
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – ARG/19/G24 (GEF/PNUD) – Red ARGENA
#   Versión:  1.0.0  |  Marzo 2026
#
#   Descripción:
#   Pipeline completo para modelar la distribución potencial de especies
#   nativas de interés para el BGEN usando ensemble modeling (biomod2).
#   Cubre desde la limpieza de registros de presencia hasta la generación
#   de mapas interactivos de idoneidad de hábitat.
#
#   Especies: 19 taxa nativos de Argentina / región pampeana
#   Fuentes de datos: GBIF + iNaturalist
#   Variables: WorldClim v2.1 (BIO01–BIO19) + SRTM (elevación, pendiente, aspecto)
#   Algoritmos: GLM, GBM, RF, MAXNET, XGBOOST (ensemble ponderado por TSS)
#
#   Estructura:
#     BLOQUE 0 – Setup: librerías, rutas y parámetros globales
#     BLOQUE 1 – Presencias: limpieza, normalización y thinning espacial
#     BLOQUE 2 – Variables: descarga, colinealidad (VIF) y stack final
#     BLOQUE 3 – Modelado: biomd2, pseudo-ausencias, validación cruzada
#     BLOQUE 4 – Ensemble y proyección: combinación y mapas de idoneidad
#     BLOQUE 5 – Métricas: tabla resumen de AUCroc y TSS por especie
#     BLOQUE 6 – Visualización: mapa interactivo multiespecies con leaflet
#
#   Requisitos:
#     R >= 4.2
#     Ver BLOQUE 0 para lista de paquetes
#
#   Uso:
#     1. Editá BASE_DIR en el BLOQUE 0 con tu ruta local
#     2. Colocá el CSV de presencias en BASE_DIR/data/
#     3. Corré el script bloque por bloque
#
# =============================================================================


# =============================================================================
# BLOQUE 0: Setup – librerías, rutas y parámetros globales
# =============================================================================

# Paquetes necesarios. Descomentá para instalar en una corrida nueva:
# install.packages(c(
#   "dplyr", "stringr", "readr",      # manipulación de datos
#   "terra", "geodata",                # datos espaciales y descarga de variables
#   "usdm", "corrplot",                # análisis de colinealidad
#   "biomod2",                         # modelado de nicho ecológico
#   "leaflet", "leaflet.extras",       # mapas interactivos
#   "htmlwidgets"                      # exportar HTML
# ))

library(dplyr)
library(stringr)
library(readr)
library(terra)
library(geodata)
library(usdm)
library(corrplot)
library(biomod2)
library(leaflet)
library(leaflet.extras)
library(htmlwidgets)

# -----------------------------------------------------------------------------
# Rutas  ←  EDITÁ SOLO ESTA SECCIÓN
# -----------------------------------------------------------------------------

BASE_DIR  <- "C:/Users/sdari/Desktop/ENM_jacaranda"   # raíz del proyecto
MODELS_DIR <- "C:/Users/sdari/Documents"               # donde biomod2 guarda modelos
                                                        # (el working dir al modelar)

# Subcarpetas (se crean automáticamente)
DATA_DIR  <- file.path(BASE_DIR, "data")
VAR_DIR   <- file.path(BASE_DIR, "variables")
OUT_DIR   <- file.path(BASE_DIR, "outputs")

dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(VAR_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Parámetros globales
# -----------------------------------------------------------------------------

# Presencias
MINIMO_REGISTROS  <- 50        # umbral mínimo de registros brutos por especie

# Variables ambientales
N_SAMPLE_VIF      <- 5000      # puntos para el análisis de correlación/VIF
UMBRAL_VIF        <- 5         # umbral máximo de VIF aceptable
SEED              <- 123       # semilla de aleatoriedad (reproducibilidad)

# Pseudo-ausencias
PA_NB_REP         <- 2         # número de sets de pseudo-ausencias
PA_NB_ABSENCES    <- 3000      # puntos de pseudo-ausencia por set
PA_SRE_QUANT      <- 0.05      # cuantil SRE (5% = conservador)

# Modelado
BM_MODELS         <- c("GLM", "GBM", "RF", "MAXNET", "XGBOOST")
CV_NB_REP         <- 3         # repeticiones de validación cruzada
CV_PERC           <- 0.8       # proporción de datos para calibración

# Ensemble
EM_METRIC_TSS     <- 0.7       # umbral mínimo TSS para entrar al ensemble
EM_METRIC_AUC     <- 0.9       # umbral mínimo AUCroc para entrar al ensemble

# Clasificación de idoneidad
CLASES_RCL <- matrix(c(
  0.0, 0.2, 1,    # 1 = Insustentable
  0.2, 0.4, 2,    # 2 = Bajo
  0.4, 0.6, 3,    # 3 = Moderado
  0.6, 1.0, 4     # 4 = Alto
), ncol = 3, byrow = TRUE)

cat("Setup completado.\n")
cat("Proyecto en:", BASE_DIR, "\n")


# =============================================================================
# BLOQUE 1: Presencias – limpieza, normalización y thinning espacial
# =============================================================================
#
# Entrada:  data/registros_unificados_geocod.csv
#           Columnas requeridas: especie, lat, lon
#
# Salida:   data/presencias_limpias.csv   (post-limpieza, pre-thinning)
#           data/presencias_thin.csv      (post-thinning, insumo del modelado)
#
# Nota sobre thinning:
#   Se conserva un único registro por píxel por especie, usando la resolución
#   del stack ambiental final (env_stack_vif5.tif). Esto reduce el sesgo
#   espacial de muestreo sin perder cobertura geográfica.
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Leer datos crudos
# -----------------------------------------------------------------------------

ruta_cruda <- file.path(DATA_DIR, "registros_unificados_geocod.csv")
df <- read_csv(ruta_cruda, show_col_types = FALSE)

cat("Registros crudos leídos:", nrow(df), "\n")

# -----------------------------------------------------------------------------
# 1.2 Normalizar nombres científicos
#     Problema: "Jacaranda mimosifolia" y "Jacaranda mimosifolia D.Don"
#     son la misma especie. Se conserva solo el binomio (Género + epíteto).
# -----------------------------------------------------------------------------

df <- df %>%
  mutate(especie_clean = word(especie, 1, 2))

cat("Especies únicas tras normalización:", length(unique(df$especie_clean)), "\n")

# -----------------------------------------------------------------------------
# 1.3 Limpiar coordenadas
# -----------------------------------------------------------------------------

df <- df %>%
  filter(
    !is.na(lat), !is.na(lon),
    lat != 0, lon != 0,
    lat >= -90,  lat <= 90,
    lon >= -180, lon <= 180
  )

cat("Registros tras limpieza de coordenadas:", nrow(df), "\n")

# -----------------------------------------------------------------------------
# 1.4 Filtrar por mínimo de registros brutos
#     Justificación: con 50 brutos se esperan ~20-30 únicos post-thinning,
#     suficiente para calibrar un modelo con pseudo-ausencias.
# -----------------------------------------------------------------------------

conteo <- df %>%
  count(especie_clean, name = "n_registros") %>%
  arrange(desc(n_registros))

especies_validas <- conteo %>%
  filter(n_registros >= MINIMO_REGISTROS) %>%
  pull(especie_clean)

cat("\nEspecies que superan el umbral mínimo (", MINIMO_REGISTROS, "registros):",
    length(especies_validas), "\n")

cat("Especies descartadas:",
    paste(conteo$especie_clean[!conteo$especie_clean %in% especies_validas],
          collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 1.5 Exportar presencias limpias
# -----------------------------------------------------------------------------

pres_limpias <- df %>%
  filter(especie_clean %in% especies_validas) %>%
  select(species = especie_clean, longitude = lon, latitude = lat)

write_csv(pres_limpias, file.path(DATA_DIR, "presencias_limpias.csv"))

cat("\nPresencias limpias guardadas:", nrow(pres_limpias), "registros,",
    length(unique(pres_limpias$species)), "especies\n")

# -----------------------------------------------------------------------------
# 1.6 Thinning espacial (1 registro por píxel por especie)
#     Requiere que el stack ambiental final ya exista (BLOQUE 2).
#     Si corrés el pipeline por primera vez, ejecutá BLOQUE 2 antes
#     y volvé acá para el thinning.
# -----------------------------------------------------------------------------

env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
pres_limpias <- read_csv(file.path(DATA_DIR, "presencias_limpias.csv"),
                         show_col_types = FALSE)

pres_vect <- terra::vect(
  pres_limpias,
  geom = c("longitude", "latitude"),
  crs  = "EPSG:4326"
)

pres_limpias$cell_id <- cellFromXY(env_final[[1]], crds(pres_vect))

pres_thin <- pres_limpias %>%
  group_by(species) %>%
  filter(!duplicated(cell_id)) %>%
  ungroup() %>%
  select(-cell_id)

cat("\nRegistros antes del thinning:", nrow(pres_limpias), "\n")
cat("Registros después del thinning:", nrow(pres_thin), "\n")
cat("Retención:", round(nrow(pres_thin) / nrow(pres_limpias) * 100, 1), "%\n\n")

pres_thin %>%
  count(species, sort = TRUE) %>%
  print(n = 25)

write_csv(pres_thin, file.path(DATA_DIR, "presencias_thin.csv"))
cat("\nThinning completado. Archivo guardado en data/presencias_thin.csv\n")


# =============================================================================
# BLOQUE 2: Variables ambientales – descarga, colinealidad y stack final
# =============================================================================
#
# Fuentes:
#   - WorldClim v2.1: 19 variables bioclimáticas (BIO01-BIO19), ~1 km
#   - SRTM: modelo digital de elevación + pendiente + aspecto, ~1 km
#
# Salida:
#   variables/bio_AOI.tif        (19 variables bioclimáticas recortadas)
#   variables/topo_AOI.tif       (elevación + pendiente + aspecto)
#   variables/env_stack.tif      (stack completo, 22 variables)
#   variables/env_stack_vif5.tif (stack reducido por VIF ≤ 5, insumo final)
#
# Nota sobre VIF:
#   Se aplica un doble filtro: correlación de Pearson |r| ≥ 0.8 + VIF stepwise
#   con umbral ≤ 5. Esto elimina la multicolinealidad que puede inflar la
#   importancia de variables correlacionadas y desestabilizar los modelos.
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Área de interés (AOI)
#     Argentina completa + buffer de 1 grado sobre el bounding box
#     de los registros de presencia.
# -----------------------------------------------------------------------------

pres_limpias <- read_csv(file.path(DATA_DIR, "presencias_limpias.csv"),
                         show_col_types = FALSE)

aoi <- ext(
  min(pres_limpias$longitude) - 1,
  max(pres_limpias$longitude) + 1,
  min(pres_limpias$latitude)  - 1,
  max(pres_limpias$latitude)  + 1
)

cat("AOI definido:\n")
print(aoi)

# -----------------------------------------------------------------------------
# 2.2 Variables bioclimáticas WorldClim v2.1
# -----------------------------------------------------------------------------

wc      <- worldclim_country("argentina", var = "bio", path = VAR_DIR, version = "2.1")
wc_aoi  <- crop(wc, aoi)
names(wc_aoi) <- sprintf("bio%02d", 1:19)

writeRaster(wc_aoi, file.path(VAR_DIR, "bio_AOI.tif"), overwrite = TRUE)
cat("WorldClim descargado y recortado:", nlyr(wc_aoi), "capas\n")

# -----------------------------------------------------------------------------
# 2.3 Topografía: elevación + pendiente + aspecto (SRTM)
# -----------------------------------------------------------------------------

elev      <- elevation_30s(country = "ARG", path = VAR_DIR)
elev_aoi  <- crop(elev, aoi)
slope     <- terrain(elev_aoi, v = "slope",  unit = "degrees")
aspect    <- terrain(elev_aoi, v = "aspect", unit = "degrees")

names(elev_aoi) <- "altitude"
names(slope)    <- "slope"
names(aspect)   <- "aspect"

topo_stack <- c(elev_aoi, slope, aspect)

writeRaster(topo_stack, file.path(VAR_DIR, "topo_AOI.tif"), overwrite = TRUE)
cat("Topografía generada:", nlyr(topo_stack), "capas\n")

# -----------------------------------------------------------------------------
# 2.4 Stack ambiental completo (22 variables)
#     Resamplear topografía al grid de WorldClim antes de apilar.
# -----------------------------------------------------------------------------

topo_stack <- resample(topo_stack, wc_aoi, method = "bilinear")
env_stack  <- c(wc_aoi, topo_stack)

writeRaster(env_stack, file.path(VAR_DIR, "env_stack.tif"), overwrite = TRUE)
cat("Stack completo armado:", nlyr(env_stack), "variables\n")

# -----------------------------------------------------------------------------
# 2.5 Análisis de colinealidad: correlación de Pearson + VIF stepwise
# -----------------------------------------------------------------------------

set.seed(SEED)

sample_pts <- spatSample(
  env_stack, size = N_SAMPLE_VIF,
  method = "random", na.rm = TRUE, as.points = FALSE
)
sample_df <- na.omit(as.data.frame(sample_pts))

# Matriz de correlación visual
corrplot(
  cor(sample_df, method = "pearson"),
  method = "color", type = "upper",
  tl.cex = 0.6, main = "Correlación de Pearson – 22 variables"
)

# VIF stepwise con umbral conservador
vif_result  <- vifstep(sample_df, th = UMBRAL_VIF)
vars_finales <- vif_result@results$Variables

cat("\nVariables seleccionadas por VIF ≤", UMBRAL_VIF, ":\n")
print(vars_finales)
cat("Total variables retenidas:", length(vars_finales), "\n")

# Barplot de VIF de variables seleccionadas
barplot(
  vif_result@results$VIF,
  names.arg = vif_result@results$Variables,
  las = 2, col = "steelblue",
  main = paste("VIF de variables seleccionadas (umbral ≤", UMBRAL_VIF, ")"),
  ylab = "Valor VIF"
)

# Stack final
env_final <- env_stack[[vars_finales]]
writeRaster(env_final, file.path(VAR_DIR, "env_stack_vif5.tif"), overwrite = TRUE)

cat("\nStack final guardado en variables/env_stack_vif5.tif\n")
cat("Bloque 2 completado. Siguiente: BLOQUE 1 sección thinning, luego BLOQUE 3.\n")


# =============================================================================
# BLOQUE 3: Modelado – biomod2, pseudo-ausencias y validación cruzada
# =============================================================================
#
# Este bloque entrena los modelos individuales para cada especie.
# biomod2 guarda los resultados en el working directory (MODELS_DIR),
# no en OUT_DIR. Esto es un comportamiento propio del paquete.
#
# Entrada:  data/presencias_thin.csv
#           variables/env_stack_vif5.tif
#
# Salida (en MODELS_DIR):
#   {Especie.nombre}/{Especie.nombre}.{Especie_nombre}_current.models.out
#
# Tiempo estimado: 15–60 min dependiendo del hardware (19 especies × 5 algoritmos)
# =============================================================================

# biomod2 guarda modelos en el working directory
setwd(MODELS_DIR)

env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
pres_thin <- read_csv(file.path(DATA_DIR, "presencias_thin.csv"),
                      show_col_types = FALSE)
especies  <- unique(pres_thin$species)

cat("Especies a modelar:", length(especies), "\n")

for (sp in especies) {

  cat("\n", rep("=", 60), "\n")
  cat("Modelando:", sp, "\n")
  cat(rep("=", 60), "\n")

  sp_name <- gsub(" ", "_", sp)

  sp_pres <- pres_thin %>%
    filter(species == sp) %>%
    select(longitude, latitude)

  cat("  Registros de presencia:", nrow(sp_pres), "\n")

  tryCatch({

    # -----------------------------------------------------------------------
    # 3a. Formatear datos para biomod2
    #     Estrategia SRE: las pseudo-ausencias se ubican fuera del rango
    #     ambiental observado de la especie (más realista que aleatorio puro).
    # -----------------------------------------------------------------------

    bm_data <- BIOMOD_FormatingData(
      resp.var       = rep(1, nrow(sp_pres)),
      expl.var       = env_final,
      resp.xy        = as.data.frame(sp_pres),
      resp.name      = sp_name,
      PA.nb.rep      = PA_NB_REP,
      PA.nb.absences = PA_NB_ABSENCES,
      PA.strategy    = "sre",
      PA.sre.quant   = PA_SRE_QUANT
    )

    # -----------------------------------------------------------------------
    # 3b. Entrenamiento de modelos individuales
    #     Validación cruzada aleatoria: CV_NB_REP repeticiones,
    #     CV_PERC de los datos para calibración.
    # -----------------------------------------------------------------------

    bm_mod <- BIOMOD_Modeling(
      bm.format    = bm_data,
      models       = BM_MODELS,
      CV.strategy  = "random",
      CV.nb.rep    = CV_NB_REP,
      CV.perc      = CV_PERC,
      metric.eval  = c("TSS", "ROC"),
      var.import   = 3,
      modeling.id  = paste0(sp_name, "_current")
    )

    cat("  ✓ Modelos entrenados:", sp, "\n")

  }, error = function(e) {
    cat("  ✗ ERROR en", sp, ":", conditionMessage(e), "\n")
  })

}

cat("\nBloque 3 completado. Modelos guardados en:", MODELS_DIR, "\n")
cat("Siguiente: BLOQUE 4 – Ensemble y proyección\n")


# =============================================================================
# BLOQUE 4: Ensemble y proyección – combinación de modelos y mapas
# =============================================================================
#
# Toma los modelos entrenados en BLOQUE 3, los combina en un ensemble
# ponderado por TSS y proyecta la idoneidad sobre Argentina.
#
# Ensemble elegido: EMwmeanByTSS
#   Pondera cada modelo según su TSS de validación. Los modelos con mejor
#   desempeño tienen más peso en el resultado final. Más robusto que
#   la media simple (EMmean) porque penaliza los modelos con sobreajuste.
#
# Entrada:  MODELS_DIR/{especie}/...models.out
#           variables/env_stack_vif5.tif
#
# Salida (en OUT_DIR/{especie}/):
#   {especie}_idoneidad_prob.tif    (probabilidad continua 0–1)
#   {especie}_idoneidad_clase.tif   (categorías 1=Insustentable ... 4=Alto)
# =============================================================================

setwd(MODELS_DIR)

env_final <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
pres_thin <- read_csv(file.path(DATA_DIR, "presencias_thin.csv"),
                      show_col_types = FALSE)
especies  <- unique(pres_thin$species)

cat("Especies a procesar:", length(especies), "\n")

for (sp in especies) {

  cat("\n", rep("=", 60), "\n")
  cat("Ensemble + Proyección:", sp, "\n")

  sp_name  <- gsub(" ", "_", sp)
  sp_punto <- gsub(" ", ".", sp)
  sp_out   <- file.path(OUT_DIR, sp_name)
  dir.create(sp_out, showWarnings = FALSE)

  tryCatch({

    # -----------------------------------------------------------------------
    # 4a. Cargar modelo entrenado desde disco
    # -----------------------------------------------------------------------

    ruta_modelo <- file.path(
      getwd(), sp_punto,
      paste0(sp_punto, ".", sp_name, "_current.models.out")
    )

    if (!file.exists(ruta_modelo)) stop("Modelo no encontrado: ", ruta_modelo)

    bm_mod <- load(ruta_modelo)
    bm_mod <- get(bm_mod)

    # -----------------------------------------------------------------------
    # 4b. Ensemble modeling
    #     Solo entran modelos que superan los umbrales de calidad:
    #     TSS ≥ EM_METRIC_TSS  y  AUCroc ≥ EM_METRIC_AUC
    # -----------------------------------------------------------------------

    bm_ens <- BIOMOD_EnsembleModeling(
      bm.mod               = bm_mod,
      models.chosen        = "all",
      em.by                = "all",
      metric.select        = c("TSS", "AUCroc"),
      metric.select.thresh = c(EM_METRIC_TSS, EM_METRIC_AUC),
      em.algo              = c("EMmean", "EMwmean"),
      var.import           = 3
    )

    # -----------------------------------------------------------------------
    # 4c. Proyección sobre el territorio
    # -----------------------------------------------------------------------

    bm_proj <- BIOMOD_EnsembleForecasting(
      bm.em         = bm_ens,
      proj.name     = "CURRENT",
      new.env       = env_final,
      models.chosen = "all",
      metric.filter = c("TSS", "AUCroc")
    )

    # -----------------------------------------------------------------------
    # 4d. Extraer predicción EMwmeanByTSS
    #     Se prefiere el ensemble ponderado por TSS sobre la media simple,
    #     por su mayor robustez ante modelos con sobreajuste (ej: RF).
    # -----------------------------------------------------------------------

    ens_preds <- get_predictions(bm_proj)

    idx <- grep("EMwmeanByTSS", names(ens_preds))
    if (length(idx) == 0) idx <- grep("EMwmean", names(ens_preds))[1]

    ens_wmean <- ens_preds[[idx]]

    # Escalar a 0–1 si biomod2 devolvió valores en 0–1000
    ens_prob <- if (max(values(ens_wmean), na.rm = TRUE) > 1) {
      ens_wmean / 1000
    } else {
      ens_wmean
    }

    names(ens_prob) <- sp_name

    # -----------------------------------------------------------------------
    # 4e. Clasificar y exportar rasters
    # -----------------------------------------------------------------------

    ens_class <- classify(ens_prob, CLASES_RCL)

    writeRaster(
      ens_prob,
      file.path(sp_out, paste0(sp_name, "_idoneidad_prob.tif")),
      overwrite = TRUE
    )

    writeRaster(
      ens_class,
      file.path(sp_out, paste0(sp_name, "_idoneidad_clase.tif")),
      overwrite = TRUE
    )

    cat("  Rango de idoneidad:",
        round(min(values(ens_prob), na.rm = TRUE), 3), "–",
        round(max(values(ens_prob), na.rm = TRUE), 3), "\n")
    cat("  ✓ Completado:", sp, "\n")

  }, error = function(e) {
    cat("  ✗ ERROR en", sp, ":", conditionMessage(e), "\n")
  })

}

cat("\nBloque 4 completado. Rasters en:", OUT_DIR, "\n")

# Verificar rasters generados
tifs <- list.files(OUT_DIR, pattern = "_idoneidad_prob.tif", recursive = TRUE)
cat("Total de rasters generados:", length(tifs), "de", length(especies), "\n")


# =============================================================================
# BLOQUE 5: Métricas – tabla resumen de AUCroc y TSS por especie
# =============================================================================
#
# Extrae las métricas de validación de cada modelo entrenado y genera
# una tabla consolidada con el desempeño promedio por algoritmo y especie.
#
# Columnas del output:
#   species           – nombre de la especie
#   algo              – algoritmo (GLM, GBM, RF, MAXNET, XGBOOST)
#   metric.eval       – métrica (AUCroc o TSS)
#   calibration_media – desempeño promedio en datos de calibración
#   validation_media  – desempeño promedio en datos de validación (el relevante)
#
# Interpretación:
#   AUCroc > 0.9 = discriminación excelente
#   TSS    > 0.7 = buen desempeño (umbral recomendado para ENM)
#   Si calibration >> validation: sobreajuste (típico en RF)
#
# Salida:  outputs/metricas_todas_especies.csv
# =============================================================================

setwd(MODELS_DIR)

pres_thin <- read_csv(file.path(DATA_DIR, "presencias_thin.csv"),
                      show_col_types = FALSE)
especies  <- unique(pres_thin$species)

resultados_resumen <- list()

for (sp in especies) {

  sp_name  <- gsub(" ", "_", sp)
  sp_punto <- gsub(" ", ".", sp)

  ruta_modelo <- file.path(
    getwd(), sp_punto,
    paste0(sp_punto, ".", sp_name, "_current.models.out")
  )

  tryCatch({

    bm_mod <- load(ruta_modelo)
    bm_mod <- get(bm_mod)
    evals  <- get_evaluations(bm_mod)

    resumen <- evals %>%
      filter(metric.eval %in% c("AUCroc", "TSS")) %>%
      group_by(algo, metric.eval) %>%
      summarise(
        calibration_media = round(mean(calibration, na.rm = TRUE), 3),
        validation_media  = round(mean(validation,  na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      mutate(species = sp)

    resultados_resumen[[sp_name]] <- resumen
    cat("✓", sp, "\n")

  }, error = function(e) {
    cat("✗", sp, ":", conditionMessage(e), "\n")
  })

}

tabla_final <- bind_rows(resultados_resumen)

write_csv(tabla_final, file.path(OUT_DIR, "metricas_todas_especies.csv"))

cat("\nTabla de métricas guardada en outputs/metricas_todas_especies.csv\n")
cat("Resumen general:\n")

tabla_final %>%
  group_by(metric.eval) %>%
  summarise(
    media_validacion = round(mean(validation_media, na.rm = TRUE), 3),
    min_validacion   = round(min(validation_media,  na.rm = TRUE), 3),
    max_validacion   = round(max(validation_media,  na.rm = TRUE), 3)
  ) %>%
  print()


# =============================================================================
# BLOQUE 6: Visualización – mapa interactivo multiespecies con leaflet
# =============================================================================
#
# Genera un visor web interactivo (HTML autocontenido) que permite explorar
# la idoneidad de hábitat de las 19 especies como capas independientes.
#
# Características del mapa:
#   - Tres fondos cartográficos: claro (CartoDB), topográfico (Esri), satelital
#   - Panel de capas: cada especie se puede activar/desactivar independientemente
#   - Por defecto muestra Jacaranda mimosifolia al abrir
#   - Leyenda con las 4 clases de idoneidad
#   - Minimapa y escala
#
# Paleta de colores:
#   1 Insustentable → gris   (#BDBDBD)
#   2 Bajo          → amarillo (#FEE08B)
#   3 Moderado      → naranja (#F46D43)
#   4 Alto          → verde   (#1A9850)
#
# Salida:  outputs/mapa_interactivo_ENM.html
#          Archivo único, portable, sin dependencias externas.
# =============================================================================

tifs_clase <- list.files(
  OUT_DIR,
  pattern   = "_idoneidad_clase.tif",
  recursive = TRUE,
  full.names = TRUE
)

nombres_especies <- tifs_clase %>%
  basename() %>%
  gsub("_idoneidad_clase.tif", "", .) %>%
  gsub("_", " ", .)

cat("Cargando", length(tifs_clase), "rasters para el mapa...\n")

# Paleta de colores por clase de idoneidad
pal <- colorFactor(
  palette  = c("#BDBDBD", "#FEE08B", "#F46D43", "#1A9850"),
  levels   = c(1, 2, 3, 4),
  na.color = "transparent"
)

# Mapa base
m <- leaflet() %>%
  addProviderTiles("CartoDB.Positron",   group = "Claro") %>%
  addProviderTiles("Esri.WorldTopoMap",  group = "Topográfico") %>%
  addProviderTiles("Esri.WorldImagery",  group = "Satelital")

# Agregar cada especie como capa independiente
for (i in seq_along(tifs_clase)) {

  sp_nombre <- nombres_especies[i]
  r_proj    <- terra::project(rast(tifs_clase[i]), "EPSG:4326")

  cat("  Agregando:", sp_nombre, "\n")

  m <- m %>%
    addRasterImage(
      r_proj,
      colors  = pal,
      opacity = 0.75,
      group   = sp_nombre
    )
}

# Controles, leyenda, escala y minimapa
m <- m %>%
  addLayersControl(
    baseGroups    = c("Claro", "Topográfico", "Satelital"),
    overlayGroups = nombres_especies,
    options       = layersControlOptions(collapsed = TRUE)
  ) %>%
  hideGroup(nombres_especies) %>%                    # ocultar todas por defecto
  showGroup("Jacaranda mimosifolia") %>%             # mostrar Jacaranda al abrir
  addLegend(
    position = "bottomright",
    colors   = c("#BDBDBD", "#FEE08B", "#F46D43", "#1A9850"),
    labels   = c("Insustentable", "Bajo", "Moderado", "Alto"),
    title    = "Idoneidad de hábitat",
    opacity  = 0.9
  ) %>%
  addScaleBar(position = "bottomleft") %>%
  addMiniMap(toggleDisplay = TRUE)

# Exportar como HTML autocontenido
ruta_html <- file.path(OUT_DIR, "mapa_interactivo_ENM.html")

htmlwidgets::saveWidget(m, file = ruta_html, selfcontained = TRUE)

cat("\nMapa interactivo guardado en:", ruta_html, "\n")
cat("Abrilo en cualquier navegador.\n")

# =============================================================================
# FIN DEL PIPELINE
# =============================================================================
#
# Archivos generados:
#
#   data/
#     presencias_limpias.csv          – registros post-limpieza
#     presencias_thin.csv             – registros post-thinning (insumo del modelado)
#
#   variables/
#     bio_AOI.tif                     – 19 variables bioclimáticas WorldClim
#     topo_AOI.tif                    – elevación, pendiente, aspecto (SRTM)
#     env_stack.tif                   – stack completo (22 variables)
#     env_stack_vif5.tif              – stack reducido por VIF ≤ 5 (9 variables)
#
#   outputs/
#     metricas_todas_especies.csv     – AUCroc y TSS por especie y algoritmo
#     mapa_interactivo_ENM.html       – visor web de las 19 especies
#     {especie}/
#       {especie}_idoneidad_prob.tif  – probabilidad continua 0-1
#       {especie}_idoneidad_clase.tif – categorías 1 a 4
#
#   MODELS_DIR (ej: C:/Users/sdari/Documents)
#     {Especie.nombre}/               – modelos biomod2 (uno por especie)
#
# Para reproducir el análisis:
#   1. Editá BASE_DIR y MODELS_DIR en BLOQUE 0
#   2. Colocá el CSV de presencias en data/
#   3. Corré BLOQUE 0, BLOQUE 2, BLOQUE 1 (sección thinning)
#   4. Corré BLOQUE 3, BLOQUE 4, BLOQUE 5, BLOQUE 6
# =============================================================================
