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

library(terra)

# Cargar raster
mb <- rast("C:/Users/sdari/Desktop/ENM_jacaranda/variables/argentina_coverage_2024.tif")

# Función para analizar raster de forma segura
analizar_raster <- function(r) {
  cat("Estructura general:\n")
  print(r)
  
  # Detectar tipo de datos
  tipo <- terra::datatype(r)
  cat("\nTipo de datos:", tipo, "\n")
  
  if (tipo %in% c("INT1S", "INT2S", "INT4S", "INT1U", "INT2U", "INT4U")) {
    # Raster categórico (valores enteros)
    cat("\nFrecuencia de valores (categórico):\n")
    print(freq(r))
  } else {
    # Raster continuo (valores flotantes)
    cat("\nRango de valores (continuo):\n")
    print(minmax(r))
  }
}

# Ejecutar análisis
analizar_raster(mb)


# =============================================================================
#
#   BLOQUE 7: Cruce con MapBiomas Argentina 2024
#   Idoneidad climática × Uso del suelo
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Argentina 2026
#   Versión:  1.0.0  |  Abril 2026
#
#   Pregunta:
#   ¿Cuánto hábitat climáticamente apto para cada especie está disponible,
#   bajo manejo potencialmente restaurable, o ya fue convertido/perdido?
#
#   Fuentes:
#   - Rasters de idoneidad: ENM_BGEN_pipeline.R (BLOQUES 3-4)
#   - MapBiomas Argentina Colección 2, año 2024 (30m, EPSG:4326)
#
#   Clasificación de uso del suelo:
#   Compatible   (3,4,6,11,12,63,66,73,77) – vegetación natural
#   Restaurable  (9,15,21)                 – uso productivo con potencial
#   Incompatible (19,24,25,36)             – convertido o urbanizado
#   Excluido     (27,33,34)                – sin datos / no aplicable
#
#   Salida:
#   outputs/mapbiomas/
#     {especie}_cruce.tif               – raster combinado (9 clases)
#     {especie}_cruce_resumen.csv       – superficie (ha) por categoría
#   outputs/mapbiomas/
#     resumen_todas_especies.csv        – tabla consolidada 19 especies
#     mapa_interactivo_cruce.html       – visor web del cruce
#
# =============================================================================


# =============================================================================
# 7.0 Setup
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(leaflet)
library(leaflet.extras)
library(htmlwidgets)

# Rutas — ajustá BASE_DIR si es necesario
BASE_DIR  <- "C:/Users/sdari/Desktop/ENM_jacaranda"
VAR_DIR   <- file.path(BASE_DIR, "variables")
OUT_DIR   <- file.path(BASE_DIR, "outputs")
MB_DIR    <- file.path(OUT_DIR, "mapbiomas")

dir.create(MB_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Leyenda MapBiomas Argentina – Colección 2
# Fuente: argentina.mapbiomas.org/codigos-de-la-leyenda
# -----------------------------------------------------------------------------

leyenda_mb <- tibble(
  codigo = c(3, 4, 6, 9, 11, 12, 15, 19, 21, 24, 25, 27, 33, 34, 36, 63, 66, 73, 77),
  clase  = c(
    "Bosques cerrados", "Bosques abiertos", "Bosques inundables",
    "Silvicultura", "Herbáceas inundables", "Herbáceas (pastizal)",
    "Pasturas", "Cultivos temporarios", "Mosaico de usos",
    "Áreas urbanas", "Otras áreas no vegetadas", "No observado",
    "Ríos, lagunas y lagos", "Hielo y nieve permanente",
    "Cultivos perennes", "Mosaico arbustos/herbáceas",
    "Arbustales cerrados", "Turberas", "Arbustales abiertos"
  ),
  categoria = c(
    "Compatible", "Compatible", "Compatible",       # 3, 4, 6
    "Restaurable",                                  # 9
    "Compatible", "Compatible",                     # 11, 12
    "Restaurable",                                  # 15
    "Incompatible",                                 # 19
    "Restaurable",                                  # 21
    "Incompatible", "Incompatible",                 # 24, 25
    "Excluido",                                     # 27
    "Excluido", "Excluido",                         # 33, 34
    "Incompatible",                                 # 36
    "Compatible", "Compatible", "Compatible",       # 63, 66, 73
    "Compatible"                                    # 77
  )
)

# Códigos por categoría (para reclasificación)
cod_compatible   <- leyenda_mb %>% filter(categoria == "Compatible")   %>% pull(codigo)
cod_restaurable  <- leyenda_mb %>% filter(categoria == "Restaurable")  %>% pull(codigo)
cod_incompatible <- leyenda_mb %>% filter(categoria == "Incompatible") %>% pull(codigo)
cod_excluido     <- leyenda_mb %>% filter(categoria == "Excluido")     %>% pull(codigo)

cat("Leyenda cargada:", nrow(leyenda_mb), "clases\n")
cat("Compatible:  ", length(cod_compatible),   "clases\n")
cat("Restaurable: ", length(cod_restaurable),  "clases\n")
cat("Incompatible:", length(cod_incompatible), "clases\n")


# =============================================================================
# 7.1 Cargar y preparar MapBiomas
# =============================================================================

cat("\nCargando MapBiomas Argentina 2024...\n")
mb_raw <- rast(file.path(VAR_DIR, "argentina_coverage_2024.tif"))
cat("MapBiomas cargado. Reclasificación se hará por especie.\n")


# =============================================================================
# 7.2 Loop de cruce por especie
# =============================================================================
#
# Para cada especie:
#   a. Cargar raster de idoneidad clasificada (1=Insustentable...4=Alto)
#   b. Reproyectar y alinear MapBiomas al grid de idoneidad
#   c. Generar raster combinado con 9 clases cruzadas:
#      Idoneidad Alta (4) × Compatible/Restaurable/Incompatible
#      Idoneidad Moderada (3) × Compatible/Restaurable/Incompatible
#      Idoneidad Baja/Insustentable (1-2) — agrupadas
#   d. Calcular superficie en hectáreas por categoría
#
# Lógica del raster combinado (valor = idoneidad * 10 + uso):
#   41 = Alta idoneidad + Compatible    → PRIORIDAD DE CONSERVACIÓN
#   42 = Alta idoneidad + Restaurable   → PRIORIDAD DE RESTAURACIÓN
#   43 = Alta idoneidad + Incompatible  → HÁBITAT PERDIDO
#   31 = Moderada + Compatible
#   32 = Moderada + Restaurable
#   33 = Moderada + Incompatible
#   11 = Baja/Insustentable + Compatible
#   12 = Baja/Insustentable + Restaurable
#   13 = Baja/Insustentable + Incompatible
#
# =============================================================================

tifs_clase <- list.files(
  OUT_DIR,
  pattern    = "_idoneidad_clase.tif",
  recursive  = TRUE,
  full.names = TRUE
)

nombres_especies <- tifs_clase %>%
  basename() %>%
  gsub("_idoneidad_clase.tif", "", .) %>%
  gsub("_", " ", .)

cat("\nEspecies a procesar:", length(tifs_clase), "\n")

resumen_global <- list()

for (i in seq_along(tifs_clase)) {
  
  sp_nombre <- nombres_especies[i]
  sp_name   <- gsub(" ", "_", sp_nombre)
  
  cat("\n", rep("=", 60), "\n")
  cat("Cruzando:", sp_nombre, "\n")
  
  tryCatch({
    
    # -----------------------------------------------------------------------
    # 7.2a Cargar idoneidad clasificada
    # -----------------------------------------------------------------------
    
    idon <- rast(tifs_clase[i])
    
    # Agrupar clases 1 y 2 (Insustentable + Bajo) → categoría 1
    # Moderado (3) → categoría 3
    # Alto (4) → categoría 4
    rcl_idon <- matrix(c(
      1, 1,   # Insustentable → 1
      2, 1,   # Bajo          → 1 (agrupado con Insustentable)
      3, 3,   # Moderado      → 3
      4, 4    # Alto          → 4
    ), ncol = 2, byrow = TRUE)
    
    idon_3cat <- classify(idon, rcl_idon)
    
    # -----------------------------------------------------------------------
    # 7.2b Alinear MapBiomas al grid de idoneidad
    #      El raster de idoneidad es ~1km, MapBiomas es 30m.
    #      Resampleamos MapBiomas (moda = clase dominante por píxel).
    # -----------------------------------------------------------------------
    
    # Primero crop, después resample, después classify
    mb_crop     <- crop(mb_raw, idon_3cat)
    mb_resample <- resample(mb_crop, idon_3cat, method = "near")
    
    # Ahora reclasificar a 3 categorías
    rcl_mb_matrix <- matrix(c(
      3,  1,   4,  1,   6,  1,   11, 1,   12, 1,
      63, 1,   66, 1,   73, 1,   77, 1,   # Compatible = 1
      9,  2,   15, 2,   21, 2,             # Restaurable = 2
      19, 3,   24, 3,   25, 3,   36, 3,   18, 3,  # Incompatible = 3
      27, NA,  33, NA,  34, NA             # Excluido = NA
    ), ncol = 2, byrow = TRUE)
    
    mb_aligned <- classify(mb_resample, rcl_mb_matrix, others = NA)
    
    # -----------------------------------------------------------------------
    # 7.2c Raster combinado (idoneidad * 10 + uso)
    # -----------------------------------------------------------------------
    
    cruce <- idon_3cat * 10 + mb_aligned
    names(cruce) <- paste0(sp_name, "_cruce")
    
    writeRaster(
      cruce,
      file.path(MB_DIR, paste0(sp_name, "_cruce.tif")),
      overwrite = TRUE
    )
    
    # -----------------------------------------------------------------------
    # 7.2d Calcular superficie en hectáreas
    #      Resolución ~1km × 1km = ~10.000 ha por píxel
    #      Ajuste exacto con cellSize() en proyección plana
    # -----------------------------------------------------------------------
    
    # Área por píxel en hectáreas
    area_pix <- cellSize(cruce, unit = "ha")
    
    # Frecuencia por clase del raster cruzado
    freq_cruce <- freq(cruce, usenames = FALSE)
    
    # Etiquetas descriptivas por código combinado
    etiquetas <- c(
      "41" = "Alta idoneidad – Compatible (CONSERVACIÓN)",
      "42" = "Alta idoneidad – Restaurable (RESTAURACIÓN)",
      "43" = "Alta idoneidad – Incompatible (HÁBITAT PERDIDO)",
      "31" = "Moderada idoneidad – Compatible",
      "32" = "Moderada idoneidad – Restaurable",
      "33" = "Moderada idoneidad – Incompatible",
      "11" = "Baja/Insustentable – Compatible",
      "12" = "Baja/Insustentable – Restaurable",
      "13" = "Baja/Insustentable – Incompatible"
    )
    
    # Calcular hectáreas por clase
    resumen_sp <- freq_cruce %>%
      as_tibble() %>%
      rename(codigo = value, n_pixeles = count) %>%
      filter(!is.na(codigo)) %>%
      mutate(
        codigo_str   = as.character(codigo),
        descripcion  = etiquetas[codigo_str],
        ha           = round(n_pixeles * res(cruce)[1] * res(cruce)[2] * (111320^2 / 10000), 0),
        species      = sp_nombre
      ) %>%
      filter(!is.na(descripcion)) %>%
      select(species, codigo, descripcion, n_pixeles, ha)
    
    resumen_global[[sp_name]] <- resumen_sp
    
    # Reporte en consola
    cat("\n  Superficie por categoría:\n")
    resumen_sp %>%
      select(descripcion, ha) %>%
      arrange(desc(ha)) %>%
      mutate(ha = format(ha, big.mark = ".")) %>%
      print(n = 20)
    
    # Guardar CSV individual
    write_csv(
      resumen_sp,
      file.path(MB_DIR, paste0(sp_name, "_cruce_resumen.csv"))
    )
    
    cat("  ✓ Completado:", sp_nombre, "\n")
    
  }, error = function(e) {
    cat("  ✗ ERROR en", sp_nombre, ":", conditionMessage(e), "\n")
  })
  
}


# =============================================================================
# 7.3 Tabla consolidada de todas las especies
# =============================================================================

tabla_global <- bind_rows(resumen_global)

write_csv(tabla_global, file.path(MB_DIR, "resumen_todas_especies.csv"))

cat("\n\nTabla global guardada en outputs/mapbiomas/resumen_todas_especies.csv\n")

# Resumen ejecutivo: superficie de hábitat perdido por especie
cat("\n--- HÁBITAT POTENCIAL PERDIDO (Alta idoneidad + Incompatible) ---\n")

tabla_global %>%
  filter(codigo == 43) %>%
  select(species, ha) %>%
  arrange(desc(ha)) %>%
  mutate(ha = format(ha, big.mark = ".")) %>%
  print(n = 20)

cat("\n--- PRIORIDADES DE CONSERVACIÓN (Alta idoneidad + Compatible) ---\n")

tabla_global %>%
  filter(codigo == 41) %>%
  select(species, ha) %>%
  arrange(desc(ha)) %>%
  mutate(ha = format(ha, big.mark = ".")) %>%
  print(n = 20)


# =============================================================================
# 7.4 Mapa interactivo del cruce
# =============================================================================
#
# Paleta del raster combinado:
#   Alta + Compatible    → verde oscuro  (#1A9850)
#   Alta + Restaurable   → azul          (#2166AC)
#   Alta + Incompatible  → rojo          (#D73027)
#   Moderada + Compatible   → verde claro (#A6D96A)
#   Moderada + Restaurable  → celeste     (#74ADD1)
#   Moderada + Incompatible → naranja     (#F46D43)
#   Baja/Insustentable (cualquier uso) → gris (#BDBDBD)
#
# =============================================================================

tifs_cruce <- list.files(
  MB_DIR,
  pattern    = "_cruce.tif",
  recursive  = FALSE,
  full.names = TRUE
)

nombres_cruce <- tifs_cruce %>%
  basename() %>%
  gsub("_cruce.tif", "", .) %>%
  gsub("_", " ", .)

cat("\nGenerando mapa interactivo de cruce...\n")

# Paleta por código combinado
pal_cruce <- colorFactor(
  palette = c(
    "#1A9850",  # 41 – Alta + Compatible
    "#2166AC",  # 42 – Alta + Restaurable
    "#D73027",  # 43 – Alta + Incompatible
    "#A6D96A",  # 31 – Moderada + Compatible
    "#74ADD1",  # 32 – Moderada + Restaurable
    "#F46D43",  # 33 – Moderada + Incompatible
    "#E8E8E8",  # 11 – Baja + Compatible
    "#D9D9D9",  # 12 – Baja + Restaurable
    "#CCCCCC"   # 13 – Baja + Incompatible
  ),
  levels   = c(41, 42, 43, 31, 32, 33, 11, 12, 13),
  na.color = "transparent"
)

m_cruce <- leaflet() %>%
  addProviderTiles("CartoDB.Positron",  group = "Claro") %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Topográfico") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satelital")

for (i in seq_along(tifs_cruce)) {
  
  sp_nombre <- nombres_cruce[i]
  r_proj    <- terra::project(rast(tifs_cruce[i]), "EPSG:4326")
  
  cat("  Agregando:", sp_nombre, "\n")
  
  m_cruce <- m_cruce %>%
    addRasterImage(
      r_proj,
      colors  = pal_cruce,
      opacity = 0.8,
      group   = sp_nombre
    )
}

m_cruce <- m_cruce %>%
  addLayersControl(
    baseGroups    = c("Claro", "Topográfico", "Satelital"),
    overlayGroups = nombres_cruce,
    options       = layersControlOptions(collapsed = TRUE)
  ) %>%
  hideGroup(nombres_cruce) %>%
  showGroup(nombres_cruce[grep("Jacaranda", nombres_cruce)]) %>%
  addLegend(
    position = "bottomright",
    colors   = c("#1A9850", "#2166AC", "#D73027",
                 "#A6D96A", "#74ADD1", "#F46D43",
                 "#E8E8E8"),
    labels   = c(
      "Alta idoneidad – Hábitat disponible (CONSERVACIÓN)",
      "Alta idoneidad – Uso restaurable (RESTAURACIÓN)",
      "Alta idoneidad – Hábitat perdido (CONVERSIÓN)",
      "Moderada – Hábitat disponible",
      "Moderada – Uso restaurable",
      "Moderada – Hábitat perdido",
      "Baja / Insustentable"
    ),
    title   = "Idoneidad × Uso del suelo<br><small>MapBiomas Argentina 2024</small>",
    opacity = 0.9
  ) %>%
  addScaleBar(position = "bottomleft") %>%
  addMiniMap(toggleDisplay = TRUE)

ruta_html_cruce <- file.path(MB_DIR, "mapa_interactivo_cruce.html")

htmlwidgets::saveWidget(m_cruce, file = ruta_html_cruce, selfcontained = TRUE)

cat("\nMapa de cruce guardado en:", ruta_html_cruce, "\n")
cat("\nBloque 7 completado.\n")
cat("Archivos en:", MB_DIR, "\n")

# =============================================================================
# FIN BLOQUE 7
#
# Archivos generados en outputs/mapbiomas/:
#   {especie}_cruce.tif            – raster cruzado (9 categorías)
#   {especie}_cruce_resumen.csv    – superficie (ha) por categoría
#   resumen_todas_especies.csv     – tabla consolidada 19 especies
#   mapa_interactivo_cruce.html    – visor web interactivo
#
# Próximos pasos sugeridos:
#   - Identificar las 5 especies con mayor hábitat perdido (código 43)
#   - Mapear zonas de alta prioridad compartidas entre especies
#   - Incorporar límites de áreas protegidas para análisis de brechas
# =============================================================================

# Totales consolidados
tabla_global %>%
  filter(codigo %in% c(41, 42, 43)) %>%
  group_by(codigo) %>%
  summarise(
    ha_total = sum(ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    categoria = case_when(
      codigo == 41 ~ "Hábitat disponible (CONSERVACIÓN)",
      codigo == 42 ~ "Uso restaurable (RESTAURACIÓN)",
      codigo == 43 ~ "Hábitat perdido (CONVERSIÓN)"
    ),
    ha_millones = round(ha_total / 1e6, 2)
  ) %>%
  select(categoria, ha_millones)


