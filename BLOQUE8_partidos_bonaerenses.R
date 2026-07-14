# =============================================================================
#
#   BLOQUE 8: Análisis por Partido Bonaerense y Ranking de Colectas
#   Idoneidad climática por división político-administrativa (Prov. Bs. As.)
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ
#   Versión:  1.1.0  |  Julio 2026
#
#   Descripción:
#   Cruza las capas de idoneidad climática con los límites políticos de los
#   partidos de la Provincia de Buenos Aires (Argentina) para generar un
#   ranking territorial que guíe las campañas de recolección de germoplasma
#   y conservación ex situ.
#
#   Entrada:
#   - Rasters de idoneidad: outputs/{especie}/{especie}_idoneidad_prob.tif
#   - Límites GADM: descargados de forma dinámica (nivel 2 de Argentina)
#
#   Salida:
#   outputs/ranking_partidos/
#     {especie}_ranking_partidos.csv    - tabla de prioridades por partido
#     consolidado_ranking_partidos.csv  - tabla unificada multi-especie
#     ranking_partidos.geojson          - capa espacial para abrir en QGIS
#     mapa_prioridades_partidos.html    - visor web interactivo de prioridades
#
# =============================================================================

# 8.0 Setup
# -----------------------------------------------------------------------------
library(terra)
library(geodata)
library(dplyr)
library(readr)
library(sf)
library(leaflet)
library(htmlwidgets)

# Rutas - Apuntando al directorio BGEN de Argentina
BASE_DIR  <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
VAR_DIR   <- file.path(BASE_DIR, "variables")
OUT_DIR   <- file.path(BASE_DIR, "outputs")
RANK_DIR  <- file.path(OUT_DIR, "ranking_partidos")

dir.create(RANK_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n=====================================================================\n")
cat("INICIANDO BLOQUE 8 - ANÁLISIS DE PARTIDOS BONAERENSES Y COLECTAS\n")
cat("=====================================================================\n")
cat("Directorio de trabajo:", BASE_DIR, "\n")

# 8.1 Descargar y preparar límites políticos (Buenos Aires)
# -----------------------------------------------------------------------------
cat("Cargando/Descargando límites políticos GADM de Argentina (nivel 2)...\n")
gadm_arg <- gadm(country = "ARG", level = 2, path = VAR_DIR)

# Filtrar para obtener solo la provincia de Buenos Aires
ba_partidos <- gadm_arg[gadm_arg$NAME_1 == "Buenos Aires", ]
cat("Límites de Buenos Aires listos. Total partidos:", nrow(ba_partidos), "\n")

# 8.2 Buscar rasters de idoneidad en outputs/
# -----------------------------------------------------------------------------
# Patrón flexible para coincidir con '_idoneidad_prob.tif' e '_idoneidad_prob_current.tif'
tifs_prob <- list.files(
  OUT_DIR,
  pattern    = "_idoneidad_prob.*\\.tif$",
  recursive  = TRUE,
  full.names = TRUE
)

# Fallback a idoneidad de clase si no se encuentra el raster continuo de probabilidad
if (length(tifs_prob) == 0) {
  tifs_prob <- list.files(
    OUT_DIR,
    pattern    = "_idoneidad_clase.*\\.tif$",
    recursive  = TRUE,
    full.names = TRUE
  )
}

if (length(tifs_prob) == 0) {
  stop("No se encontraron rasters de idoneidad en la carpeta outputs/. Por favor ejecute el pipeline de modelado primero.")
}

nombres_especies <- tifs_prob %>%
  basename() %>%
  gsub("_idoneidad_prob.*\\.tif", "", .) %>%
  gsub("_idoneidad_clase.*\\.tif", "", .) %>%
  gsub("_", " ", .)

cat("Rasters de idoneidad encontrados:", length(tifs_prob), "\n")
for(s in seq_along(nombres_especies)) {
  cat(sprintf("  [%d] %s (Ruta: %s)\n", s, nombres_especies[s], basename(tifs_prob[s])))
}

# 8.3 Loop de extracción y cálculo por especie
# -----------------------------------------------------------------------------
resultados_todas_especies <- list()

# Creamos una copia del vector de partidos para acumular las estadísticas de cada especie
partidos_stats <- ba_partidos

for (i in seq_along(tifs_prob)) {
  sp_nombre <- nombres_especies[i]
  sp_name   <- gsub(" ", "_", sp_nombre)
  tif_path  <- tifs_prob[i]
  
  cat(sprintf("\nProcesando extracción espacial para: '%s'...\n", sp_nombre))
  
  r <- rast(tif_path)
  
  # Alinear CRS entre el raster y el vector
  if (crs(r) != crs(ba_partidos)) {
    cat("  Alineando proyecciones (CRS)...\n")
    ba_partidos_proj <- project(ba_partidos, crs(r))
  } else {
    ba_partidos_proj <- ba_partidos
  }
  
  # Intentar recortar el raster al extent de Buenos Aires para agilizar
  r_crop <- tryCatch({
    crop(r, ba_partidos_proj)
  }, error = function(e) {
    NULL
  })
  
  # Inicializar data frame de resultados para todos los partidos bonaerenses
  df_sp <- data.frame(
    partido = ba_partidos_proj$NAME_2,
    provincia = ba_partidos_proj$NAME_1,
    idoneidad_media = 0.0,
    idoneidad_maxima = 0.0,
    area_alta_idoneidad_km2 = 0.0,
    species = sp_nombre,
    stringsAsFactors = FALSE
  )
  
  # Si el raster intersecta y tiene valores en la Provincia de Buenos Aires
  if (!is.null(r_crop) && hasValues(r_crop)) {
    r_mask <- mask(r_crop, ba_partidos_proj)
    
    # Extraer valores promedio y máximos
    cat("  Calculando estadísticas de idoneidad...\n")
    mean_vals <- terra::extract(r_mask, ba_partidos_proj, fun = mean, na.rm = TRUE)
    max_vals  <- terra::extract(r_mask, ba_partidos_proj, fun = max, na.rm = TRUE)
    
    # Obtener el valor máximo para detectar si está en escala 0-1 o 0-1000
    max_val_raster <- max(values(r_mask), na.rm = TRUE)
    if (is.na(max_val_raster)) max_val_raster <- 0
    
    # Definir umbral de alta idoneidad (probabilidad > 0.6 o > 600)
    umbral <- if (max_val_raster > 10) 600 else 0.6
    
    cat("  Calculando superficie de alta idoneidad (umbral =", umbral, ")...\n")
    r_high <- r_mask >= umbral
    area_raster <- cellSize(r_mask, unit = "km")
    r_high_area <- mask(area_raster, r_high, maskvalues = 0)
    
    area_vals <- terra::extract(r_high_area, ba_partidos_proj, fun = sum, na.rm = TRUE)
    
    # Asignar resultados
    df_sp$idoneidad_media  <- round(mean_vals[, 2], 4)
    df_sp$idoneidad_maxima <- round(max_vals[, 2], 4)
    df_sp$area_alta_idoneidad_km2 <- round(area_vals[, 2], 2)
    
    # Normalizar de 0 a 1 si el raster de entrada usaba escala 0-1000 (biomod2)
    if (max_val_raster > 10) {
      df_sp$idoneidad_media  <- round(df_sp$idoneidad_media / 1000, 4)
      df_sp$idoneidad_maxima <- round(df_sp$idoneidad_maxima / 1000, 4)
    }
    
    # Limpiar NAs
    df_sp$idoneidad_media[is.na(df_sp$idoneidad_media)] <- 0
    df_sp$idoneidad_maxima[is.na(df_sp$idoneidad_maxima)] <- 0
    df_sp$area_alta_idoneidad_km2[is.na(df_sp$area_alta_idoneidad_km2)] <- 0
    
  } else {
    cat("  [AVISO] La especie no tiene cobertura geográfica en la Provincia de Buenos Aires.\n")
    cat("          Todos los estadísticos se inicializan en 0 para esta región.\n")
  }
  
  # Guardar archivo individual de ranking ordenado por idoneidad media
  df_ranking <- df_sp %>%
    arrange(desc(idoneidad_media), desc(area_alta_idoneidad_km2)) %>%
    mutate(rango = row_number())
  
  write_csv(df_ranking, file.path(RANK_DIR, paste0(sp_name, "_ranking_partidos.csv")))
  cat("  ✓ Archivo guardado:", paste0(sp_name, "_ranking_partidos.csv"), "\n")
  
  resultados_todas_especies[[sp_name]] <- df_sp
  
  # Agregar atributos a la copia del vector espacial
  col_name_media <- paste0(sp_name, "_id_media")
  col_name_area  <- paste0(sp_name, "_area_km2")
  
  partidos_stats[[col_name_media]] <- df_sp$idoneidad_media
  partidos_stats[[col_name_area]]  <- df_sp$area_alta_idoneidad_km2
}

# 8.4 Consolidación de prioridades multiespecie (Ranking Global)
# -----------------------------------------------------------------------------
cat("\nGenerando ranking consolidado de prioridad de colecta...\n")

df_consolidado <- bind_rows(resultados_todas_especies)

# Calcular estadísticas de síntesis por partido bonaerense
df_ranking_global <- df_consolidado %>%
  group_by(partido, provincia) %>%
  summarise(
    idoneidad_media_bgen   = round(mean(idoneidad_media), 4),
    idoneidad_max_bgen     = round(max(idoneidad_maxima), 4),
    area_alta_total_km2    = round(sum(area_alta_idoneidad_km2), 2),
    riqueza_especies_aptas = sum(idoneidad_media > 0.2),
    .groups = "drop"
  ) %>%
  mutate(
    # Índice de prioridad de colecta: combina idoneidad media (70%) y riqueza de especies aptas (30%)
    indice_prioridad = round((idoneidad_media_bgen * 0.7) + ((riqueza_especies_aptas / length(tifs_prob)) * 0.3), 4)
  ) %>%
  arrange(desc(indice_prioridad), desc(area_alta_total_km2)) %>%
  mutate(rango_prioridad = row_number())

write_csv(df_ranking_global, file.path(RANK_DIR, "consolidado_ranking_partidos.csv"))
cat("✓ Guardado ranking consolidado en: consolidado_ranking_partidos.csv\n")

# Mostrar top 10 partidos prioritarios en consola
cat("\nTop 10 Partidos Bonaerenses de Mayor Prioridad para Colecta:\n")
print(head(df_ranking_global %>% select(rango_prioridad, partido, idoneidad_media_bgen, riqueza_especies_aptas, indice_prioridad), 10))

# 8.5 Exportar capa espacial para QGIS (GeoJSON)
# -----------------------------------------------------------------------------
cat("\nExportando datos espaciales a GeoJSON para QGIS...\n")

# Hacer match por nombre de partido para conservar la geometría e incorporar los campos consolidados
partidos_stats$idoneidad_bgen_media <- df_ranking_global$idoneidad_media_bgen[match(partidos_stats$NAME_2, df_ranking_global$partido)]
partidos_stats$idoneidad_bgen_max   <- df_ranking_global$idoneidad_max_bgen[match(partidos_stats$NAME_2, df_ranking_global$partido)]
partidos_stats$area_bgen_total_km2  <- df_ranking_global$area_alta_total_km2[match(partidos_stats$NAME_2, df_ranking_global$partido)]
partidos_stats$riqueza_aptas        <- df_ranking_global$riqueza_especies_aptas[match(partidos_stats$NAME_2, df_ranking_global$partido)]
partidos_stats$indice_prioridad     <- df_ranking_global$indice_prioridad[match(partidos_stats$NAME_2, df_ranking_global$partido)]
partidos_stats$rango_prioridad      <- df_ranking_global$rango_prioridad[match(partidos_stats$NAME_2, df_ranking_global$partido)]

# Convertir a sf para escritura compatible
ba_sf <- st_as_sf(partidos_stats)

# Guardar GeoJSON
geojson_path <- file.path(RANK_DIR, "ranking_partidos.geojson")
st_write(ba_sf, geojson_path, delete_dsn = TRUE, quiet = TRUE)
cat("✓ Capa GeoJSON guardada en:", geojson_path, "\n")

# 8.6 Generar Visor Web Interactivo (Leaflet)
# -----------------------------------------------------------------------------
cat("Generando visor interactivo HTML (Leaflet)...\n")

# Reproyectar a WGS84 para leaflet
ba_sf_4326 <- st_transform(ba_sf, 4326)

# Paleta de colores para el mapa
pal <- colorNumeric(
  palette = "YlOrRd",
  domain = ba_sf_4326$indice_prioridad
)

# Mapa Leaflet interactivo
mapa <- leaflet(ba_sf_4326) %>%
  addProviderTiles("CartoDB.Positron", group = "Mapa Base Claro") %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Topográfico") %>%
  addPolygons(
    fillColor = ~pal(indice_prioridad),
    weight = 1.5,
    opacity = 1,
    color = "#888",
    dashArray = "3",
    fillOpacity = 0.65,
    group = "Prioridades de Colecta BGEN",
    highlightOptions = highlightOptions(
      weight = 3,
      color = "#444",
      dashArray = "",
      fillOpacity = 0.85,
      bringToFront = TRUE
    ),
    label = ~paste0(
      "<strong>Partido:</strong> ", NAME_2, 
      "<br><strong>Rango de Prioridad:</strong> #", rango_prioridad,
      "<br><strong>Índice Prioridad:</strong> ", round(indice_prioridad, 4),
      "<br><strong>Idoneidad Promedio:</strong> ", round(idoneidad_bgen_media, 4),
      "<br><strong>Riqueza Aptas:</strong> ", riqueza_aptas, " de ", length(tifs_prob)
    ) %>% lapply(htmltools::HTML),
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "13px",
      direction = "auto"
    )
  ) %>%
  addLegend(
    pal = pal,
    values = ~indice_prioridad,
    opacity = 0.75,
    title = "Índice Prioridad BGEN",
    position = "bottomright"
  ) %>%
  addLayersControl(
    baseGroups = c("Mapa Base Claro", "Topográfico"),
    overlayGroups = c("Prioridades de Colecta BGEN"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addScaleBar(position = "bottomleft")

# Guardar visor HTML auto-contenido
html_path <- file.path(RANK_DIR, "mapa_prioridades_partidos.html")
saveWidget(mapa, file = html_path, selfcontained = FALSE)
cat("✓ Visor interactivo HTML guardado en:", html_path, "\n")

cat("\n=====================================================================\n")
cat("BLOQUE 8 COMPLETADO CON ÉXITO. RESULTADOS EN outputs/ranking_partidos/\n")
cat("=====================================================================\n\n")
