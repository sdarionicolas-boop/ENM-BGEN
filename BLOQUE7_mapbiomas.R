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

# Reclasificar a 4 categorías numéricas:
# 1 = Compatible | 2 = Restaurable | 3 = Incompatible | 0 = Excluido (NA)
rcl_mb <- rbind(
  cbind(cod_compatible,   1),
  cbind(cod_restaurable,  2),
  cbind(cod_incompatible, 3),
  cbind(cod_excluido,     NA)
)

# Formato requerido por classify(): matriz de 3 columnas (from, to, becomes)
rcl_mb_matrix <- cbind(rcl_mb[, 1], rcl_mb[, 1], rcl_mb[, 2])

mb_cat <- classify(mb_raw, rcl_mb_matrix, others = NA)
names(mb_cat) <- "uso_suelo"

cat("MapBiomas reclasificado a 3 categorías.\n")


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

    mb_aligned <- resample(mb_cat, idon_3cat, method = "near")
    mb_aligned <- crop(mb_aligned, idon_3cat)
    mb_aligned <- mask(mb_aligned, idon_3cat)

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
