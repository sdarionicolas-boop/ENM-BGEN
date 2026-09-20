# =============================================================================
#
#   BLOQUE 7 (Perú): Cruce con MapBiomas Perú 2024 (Colección 3)
#   Idoneidad climática × Uso del suelo — Mauritia flexuosa (Aguaje)
#
#   Autor:    Darío Nicolás Sánchez Leguizamón
#   Proyecto: Beca BIEI 2025 – BGEN/UNAJ | Premio MapBiomas Perú 2026
#   Versión:  2.0.0  |  Mayo 2026
#
#   Adaptación del BLOQUE7_mapbiomas.R (versión Argentina) a la leyenda
#   de MapBiomas Perú / Amazonia Colección 6.0.
#
#   CORRECCIÓN v2.0: Reemplazada leyenda de MapBiomas Argentina Col.2
#   por leyenda de MapBiomas Amazonia Col.6.0 (códigos distintos).
#   El script original fallaba porque los códigos 63,66,73,77 (Argentina)
#   no existen en Amazonia → todos los píxeles se volvían NA.
#
#   Pregunta:
#   ¿Cuánto hábitat climáticamente apto para cada especie está disponible,
#   bajo manejo potencialmente restaurable, o ya fue convertido/perdido?
#
#   Fuentes:
#   - Rasters de idoneidad: ENM_BGEN_pipeline.R (BLOQUES 3-4)
#   - MapBiomas Amazonia Colección 6.0, año 2023 (30m, EPSG:4326)
#     Fuente: amazon.mapbiomas.org
#
#   Clasificación de uso del suelo (Amazonia Col. 6.0):
#   Compatible   (3,4,6,11,12,13,68) – vegetación natural (bosque, humedal, pastizal)
#   Restaurable  (9,15,21,35)        – silvicultura, pasturas, mosaico agropecuario
#   Incompatible (18,24,25,30)       – cultivos, urbano, no vegetado, minería
#   Excluido     (23,27,29,33,34)    – agua, afloramientos, no observado
#
#   Salida:
#   outputs/mapbiomas/
#     {especie}_cruce.tif               – raster combinado (9 clases)
#     {especie}_cruce_resumen.csv       – superficie (km²) por categoría
#   outputs/mapbiomas/
#     resumen_todas_especies.csv        – tabla consolidada todas las especies
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
BASE_DIR <- "C:/Users/sdari/Desktop/ENM-BGEN-main"
VAR_DIR <- file.path(BASE_DIR, "variables")
OUT_DIR <- file.path(BASE_DIR, "outputs")
MB_DIR <- file.path(OUT_DIR, "mapbiomas")

dir.create(MB_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Leyenda MapBiomas Amazonia – Colección 6.0 (2023)
# Fuente: amazon.mapbiomas.org/codigos-de-la-leyenda
#
# NOTA: Esta leyenda es DIFERENTE a MapBiomas Argentina.
# Códigos como 63, 66, 73, 77 NO existen en Amazonia Col.6.0.
# Los códigos 13, 18, 30, 35, 68 SÍ existen en Amazonia pero
# no estaban en la leyenda Argentina → causaban todos los píxeles en NA.
# -----------------------------------------------------------------------------

leyenda_mb <- tibble(
  codigo = c(
    # Compatible – vegetación natural conservada
    3, 4, 6, 11, 12, 13, 68,
    # Restaurable – uso productivo con potencial de recuperación
    9, 15, 21, 35,
    # Incompatible – convertido / impermeabilizado
    18, 24, 25, 30,
    # Excluido – agua, afloramientos, sin datos
    23, 27, 29, 33, 34
  ),
  clase = c(
    # Compatible
    "Formación forestal (bosque)", "Formación savánica", "Formación inundable",
    "Campo inundado y pantano", "Formación campestre (pastizal natural)",
    "Otras formaciones no forestales", "Manglar / Apicum",
    # Restaurable
    "Silvicultura", "Pasturas", "Mosaico agropecuario", "Palma de aceite",
    # Incompatible
    "Agricultura/Cultivos temporarios", "Área urbanizada",
    "Otras áreas no vegetadas", "Minería",
    # Excluido
    "Playa, duna y arenal", "No observado",
    "Afloramiento rocoso", "Ríos, lagos y océano", "Glaciar"
  ),
  categoria = c(
    # Compatible
    "Compatible", "Compatible", "Compatible", "Compatible",
    "Compatible", "Compatible", "Compatible",
    # Restaurable
    "Restaurable", "Restaurable", "Restaurable", "Restaurable",
    # Incompatible
    "Incompatible", "Incompatible", "Incompatible", "Incompatible",
    # Excluido
    "Excluido", "Excluido", "Excluido", "Excluido", "Excluido"
  )
)

# Códigos por categoría (para reclasificación)
cod_compatible <- leyenda_mb %>%
  filter(categoria == "Compatible") %>%
  pull(codigo)
cod_restaurable <- leyenda_mb %>%
  filter(categoria == "Restaurable") %>%
  pull(codigo)
cod_incompatible <- leyenda_mb %>%
  filter(categoria == "Incompatible") %>%
  pull(codigo)
cod_excluido <- leyenda_mb %>%
  filter(categoria == "Excluido") %>%
  pull(codigo)

cat("Leyenda cargada:", nrow(leyenda_mb), "clases\n")
cat("Compatible:  ", length(cod_compatible), "clases\n")
cat("Restaurable: ", length(cod_restaurable), "clases\n")
cat("Incompatible:", length(cod_incompatible), "clases\n")



procesar_coleccion <- function(nombre_col, archivo_col, titulo_mapa_col) {
  cat("\n\n*************************************************************\n")
  cat("INICIANDO PROCESO PARA:", nombre_col, "\n")
  cat("*************************************************************\n")

  MB_COL_DIR <- file.path(MB_DIR, nombre_col)
  dir.create(MB_COL_DIR, showWarnings = FALSE, recursive = TRUE)

  # =============================================================================
  # 7.1 Cargar y preparar MapBiomas
  # =============================================================================
  cat("\nCargando", titulo_mapa_col, "...\n")
  mb_raw <- rast(file.path(VAR_DIR, archivo_col))

  rcl_mb <- rbind(
    cbind(cod_compatible,   1),
    cbind(cod_restaurable,  2),
    cbind(cod_incompatible, 3),
    cbind(cod_excluido,     NA)
  )
  mb_cat <- classify(mb_raw, rcl_mb, others = NA)
  names(mb_cat) <- "uso_suelo"
  cat("MapBiomas reclasificado a 3 categorías.\n")

  # =============================================================================
  # 7.2 Loop de cruce por especie
  # =============================================================================
  tifs_clase <- list.files(
    OUT_DIR,
    pattern    = "_idoneidad_clase_current.tif",
    recursive  = TRUE,
    full.names = TRUE
  )
  nombres_especies <- tifs_clase %>%
    basename() %>%
    gsub("_idoneidad_clase_current.tif", "", .) %>%
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

  # =============================================================================
  # 7.3 Tabla consolidada
  # =============================================================================
  tabla_global <- bind_rows(resumen_global)
  write_csv(tabla_global, file.path(MB_COL_DIR, "resumen_todas_especies.csv"))

  # =============================================================================
  # 7.4 Mapa interactivo
  # =============================================================================
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
    showGroup(nombres_cruce[grep("Mauritia", nombres_cruce)]) %>%
    addLegend(
      position = "bottomright",
      colors = c("#1A9850", "#2166AC", "#D73027", "#A6D96A", "#74ADD1", "#F46D43", "#E8E8E8"),
      labels = c("Alta idoneidad – Hábitat disponible", "Alta idoneidad – Uso restaurable", "Alta idoneidad – Hábitat perdido", "Moderada – Hábitat disponible", "Moderada – Uso restaurable", "Moderada – Hábitat perdido", "Baja / Insustentable"),
      title = paste0("Idoneidad × Uso del suelo<br><small>", titulo_mapa_col, "</small>"), opacity = 0.9
    ) %>%
    addScaleBar(position = "bottomleft") %>%
    addMiniMap(toggleDisplay = TRUE)
  
  ruta_html_cruce <- file.path(MB_COL_DIR, "mapa_interactivo_cruce.html")
  htmlwidgets::saveWidget(m_cruce, file = ruta_html_cruce, selfcontained = TRUE)
  cat("\nMapa de cruce guardado en:", ruta_html_cruce, "\n")
}

procesar_coleccion(
  nombre_col = "Amazonia_2023",
  archivo_col = "mapbiomas_collection60_integration_v1-classification_2023.tif",
  titulo_mapa_col = "MapBiomas Amazonia Col. 6.0 (2023)"
)

procesar_coleccion(
  nombre_col = "Peru_2024",
  archivo_col = "2024_coverage_lclu_16-1-1_4334c236-54e8-4f81-add2-6f8c31662401.tif",
  titulo_mapa_col = "MapBiomas Perú 2024"
)

# FIN BLOQUE 7
