# =============================================================================
#
#   BLOQUE 7 (Chile): Cruce con MapBiomas Chile
#   Idoneidad climática × Uso del suelo — Jubaea chilensis (Palma chilena)
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Chile 2026
#   Versión:  2.0.0  |  Septiembre 2026
#
#   Procesa DOS colecciones de MapBiomas Chile, igual que se hizo para
#   Argentina/Perú cuando hay más de una colección relevante:
#
#   - Colección 2 (2024) es el PRODUCTO OFICIAL EXIGIDO por el reglamento
#     del Premio MapBiomas Chile 2026 (requisito de elegibilidad #1:
#     "Producto land cover colección 2"). Este es el resultado que se
#     postula.
#   - Colección 1 (2022) se conserva como ANEXO COMPARATIVO/TEMPORAL
#     (ya estaba calculado antes de leer el reglamento en detalle;
#     se mantiene porque permite mostrar robustez del resultado en
#     dos años/colecciones distintas, no porque sea válido para la
#     postulación).
#
#   Pregunta:
#   ¿Cuánto hábitat climáticamente apto para Jubaea chilensis está
#   disponible como vegetación natural, bajo uso restaurable, o ya
#   fue convertido/perdido por uso antrópico?
#
#   Fuentes:
#   - Raster de idoneidad: ENM_BGEN_pipeline.R (BLOQUES 3-4), este proyecto
#   - MapBiomas Chile Colección 1, año 2022 (30m, EPSG:4326)
#     Descarga directa: storage.googleapis.com/mapbiomas-public/initiatives/
#     chile/coverage/chile_coverage_2022.tif
#   - MapBiomas Chile Colección 2, año 2024 (30m, EPSG:4326)
#     Asset Google Earth Engine: projects/mapbiomas-chile/assets/LULC/
#     COLLECTION-02/CLASSIFICATIONS/classification-final/clasificacion-final-2
#     (banda "classification_2024"). No disponible como descarga plana
#     pública al momento de este análisis (sep-2026); se obtuvo vía
#     Earth Engine Python API, exportada en tiles y mosaicada con terra.
#   - Leyenda oficial: chile.mapbiomas.org/codigos-de-la-leyenda/
#
#   Clasificación de uso del suelo (códigos verificados contra la leyenda
#   oficial de Colección 2 y contrastados con los códigos realmente
#   presentes en la AOI de Jubaea chilensis, en ambas colecciones):
#
#   Compatible   (3,59,60,67,11,12,63,66) – bosque, humedal, pastizal
#                                            natural, estepa, matorral
#   Restaurable  (9,15,21)                – silvicultura, pastura,
#                                            mosaico agropecuario
#                                            (código 21 solo aparece en
#                                             Colección 1; Colección 2
#                                             no usa esta clase mosaico)
#   Incompatible (18,24,25)               – agricultura, infraestructura,
#                                            otra área no vegetada (degradada)
#   Excluido     (0,23,29,33,34,61)       – sin datos, arena/playa/duna,
#                                            afloramiento rocoso, agua,
#                                            hielo/nieve, salar
#                                            (naturalmente no vegetado,
#                                             fuera del marco de conversión
#                                             antrópica; no aplica el
#                                             esquema Compatible/Restaurable/
#                                             Incompatible)
#
#   Salida:
#   outputs/mapbiomas/{Col2_2024,Col1_2022}/
#     Jubaea_chilensis_cruce.tif           – raster combinado (9 clases)
#     Jubaea_chilensis_cruce_resumen.csv   – superficie (km²) por categoría
#     mapa_interactivo_cruce_chile.html    – visor web del cruce
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

BASE_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR  <- file.path(BASE_DIR, "variables")
OUT_DIR  <- file.path(BASE_DIR, "outputs")
MB_DIR   <- file.path(OUT_DIR, "mapbiomas")

dir.create(MB_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Leyenda MapBiomas Chile (común a Colección 1 y Colección 2)
# Fuente: chile.mapbiomas.org/codigos-de-la-leyenda/
# Códigos verificados por inspección directa (freq()) de ambos rasters
# recortados a la AOI de Jubaea chilensis (29,9°S-35,3°S).
# -----------------------------------------------------------------------------

leyenda_mb <- tibble(
  codigo = c(
    # Compatible – vegetación natural conservada
    3, 59, 60, 67, 11, 12, 63, 66,
    # Restaurable – uso productivo con potencial de recuperación
    9, 15, 21,
    # Incompatible – convertido / degradado por uso antrópico
    18, 24, 25,
    # Excluido – naturalmente no vegetado, sin datos
    0, 23, 29, 33, 34, 61
  ),
  clase = c(
    # Compatible
    "Bosque", "Bosque primario", "Bosque secundario", "Bosque achaparrado",
    "Humedal", "Pastizal (natural)", "Estepa", "Matorral",
    # Restaurable
    "Silvicultura", "Pastura", "Mosaico agropecuario",
    # Incompatible
    "Agricultura", "Infraestructura", "Otra área no vegetada (degradada)",
    # Excluido
    "Sin datos", "Arena, playa y duna", "Afloramiento rocoso",
    "Río, lago u océano", "Hielo y nieve", "Salar"
  ),
  categoria = c(
    rep("Compatible", 8),
    rep("Restaurable", 3),
    rep("Incompatible", 3),
    rep("Excluido", 6)
  )
)

cod_compatible   <- leyenda_mb %>% filter(categoria == "Compatible")   %>% pull(codigo)
cod_restaurable  <- leyenda_mb %>% filter(categoria == "Restaurable")  %>% pull(codigo)
cod_incompatible <- leyenda_mb %>% filter(categoria == "Incompatible") %>% pull(codigo)
cod_excluido     <- leyenda_mb %>% filter(categoria == "Excluido")     %>% pull(codigo)

cat("Leyenda cargada:", nrow(leyenda_mb), "clases\n")
cat("Compatible:  ", length(cod_compatible), "clases\n")
cat("Restaurable: ", length(cod_restaurable), "clases\n")
cat("Incompatible:", length(cod_incompatible), "clases\n")
cat("Excluido:    ", length(cod_excluido), "clases\n")

rcl_mb <- rbind(
  cbind(cod_compatible,   1),
  cbind(cod_restaurable,  2),
  cbind(cod_incompatible, 3),
  cbind(cod_excluido,     NA)
)


procesar_coleccion <- function(nombre_col, archivo_col, titulo_mapa_col) {
  cat("\n\n*************************************************************\n")
  cat("INICIANDO PROCESO PARA:", nombre_col, "\n")
  cat("*************************************************************\n")

  MB_COL_DIR <- file.path(MB_DIR, nombre_col)
  dir.create(MB_COL_DIR, showWarnings = FALSE, recursive = TRUE)

  # 7.1 Cargar y reclasificar MapBiomas Chile
  cat("\nCargando", titulo_mapa_col, "...\n")
  mb_raw <- rast(file.path(VAR_DIR, archivo_col))
  mb_cat <- classify(mb_raw, rcl_mb, others = NA)
  names(mb_cat) <- "uso_suelo"
  cat("MapBiomas reclasificado a 3 categorías.\n")

  # 7.2 Cruce con idoneidad climática
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
      idon <- rast(tifs_clase[i])
      rcl_idon <- matrix(c(1, 1, 2, 1, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
      idon_3cat <- classify(idon, rcl_idon)

      mb_cropped <- crop(mb_cat, idon_3cat)
      mb_aligned <- resample(mb_cropped, idon_3cat, method = "near")
      mb_aligned <- mask(mb_aligned, idon_3cat)

      cruce <- idon_3cat * 10 + mb_aligned
      names(cruce) <- paste0(sp_name, "_cruce")
      writeRaster(cruce, file.path(MB_COL_DIR, paste0(sp_name, "_cruce.tif")), overwrite = TRUE)

      area_km2_rast <- cellSize(cruce, unit = "km")
      freq_cruce <- freq(cruce, usenames = FALSE)
      km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

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

      resumen_sp <- freq_cruce %>%
        as_tibble() %>%
        rename(codigo = value, n_pixeles = count) %>%
        filter(!is.na(codigo)) %>%
        mutate(
          codigo_str  = as.character(codigo),
          descripcion = etiquetas[codigo_str],
          km2         = round(n_pixeles * km2_por_pixel, 0),
          species     = sp_nombre
        ) %>%
        filter(!is.na(descripcion)) %>%
        select(species, codigo, descripcion, n_pixeles, km2)

      resumen_global[[sp_name]] <- resumen_sp

      cat("\n  Superficie por categoría:\n")
      print(resumen_sp %>% select(descripcion, km2) %>% arrange(desc(km2)) %>% mutate(km2 = format(km2, big.mark = ".")), n = 20)
      write_csv(resumen_sp, file.path(MB_COL_DIR, paste0(sp_name, "_cruce_resumen.csv")))
      cat("  ✓ Completado:", sp_nombre, "\n")
    }, error = function(e) {
      cat("  ✗ ERROR en", sp_nombre, ":", conditionMessage(e), "\n")
    })
  }

  tabla_global <- bind_rows(resumen_global)
  write_csv(tabla_global, file.path(MB_COL_DIR, "resumen_todas_especies.csv"))

  # 7.3 Mapa interactivo del cruce
  tifs_cruce <- list.files(MB_COL_DIR, pattern = "_cruce.tif", full.names = TRUE)
  nombres_cruce <- tifs_cruce %>% basename() %>% gsub("_cruce.tif", "", .) %>% gsub("_", " ", .)

  pal_cruce <- colorFactor(
    palette = c("#1A9850", "#2166AC", "#D73027", "#A6D96A", "#74ADD1", "#F46D43", "#E8E8E8", "#D9D9D9", "#CCCCCC"),
    levels = c(41, 42, 43, 31, 32, 33, 11, 12, 13), na.color = "transparent"
  )
  m_cruce <- leaflet() %>%
    addProviderTiles("Esri.WorldGrayCanvas", group = "Claro") %>%
    addProviderTiles("Esri.WorldTopoMap", group = "Topográfico") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satelital")

  for (i in seq_along(tifs_cruce)) {
    sp_nombre <- nombres_cruce[i]
    r_proj <- terra::project(rast(tifs_cruce[i]), "EPSG:4326", method = "near")
    m_cruce <- m_cruce %>% addRasterImage(r_proj, colors = pal_cruce, opacity = 0.8, group = sp_nombre, method = "ngb", maxBytes = Inf)
  }
  m_cruce <- m_cruce %>%
    addLayersControl(baseGroups = c("Claro", "Topográfico", "Satelital"), overlayGroups = nombres_cruce, options = layersControlOptions(collapsed = TRUE)) %>%
    hideGroup(nombres_cruce) %>%
    showGroup(nombres_cruce[grep("Jubaea", nombres_cruce)]) %>%
    addLegend(
      position = "bottomright",
      colors = c("#1A9850", "#2166AC", "#D73027", "#A6D96A", "#74ADD1", "#F46D43", "#E8E8E8"),
      labels = c("Alta idoneidad – Hábitat disponible", "Alta idoneidad – Uso restaurable", "Alta idoneidad – Hábitat perdido", "Moderada – Hábitat disponible", "Moderada – Uso restaurable", "Moderada – Hábitat perdido", "Baja / Insustentable"),
      title = paste0("Idoneidad × Uso del suelo<br><small>", titulo_mapa_col, "</small>"), opacity = 0.9
    ) %>%
    addScaleBar(position = "bottomleft") %>%
    addMiniMap(toggleDisplay = TRUE)

  Sys.setenv(RSTUDIO_PANDOC = "C:/Users/sdari/AppData/Local/Programs/Quarto/bin/tools")
  ruta_html_cruce <- file.path(MB_COL_DIR, "mapa_interactivo_cruce_chile.html")
  htmlwidgets::saveWidget(m_cruce, file = ruta_html_cruce, selfcontained = TRUE)
  cat("\nMapa de cruce guardado en:", ruta_html_cruce, "\n")
}

procesar_coleccion(
  nombre_col      = "Col2_2024",
  archivo_col     = "chile_coverage2_2024.tif",
  titulo_mapa_col = "MapBiomas Chile Colección 2 (2024) — resultado oficial de la postulación"
)

procesar_coleccion(
  nombre_col      = "Col1_2022",
  archivo_col     = "chile_coverage_2022.tif",
  titulo_mapa_col = "MapBiomas Chile Colección 1 (2022) — anexo comparativo temporal"
)

# FIN BLOQUE 7 (Chile)
