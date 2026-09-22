# =============================================================================
#
#   BLOQUE 15 (Chile): Proyección multi-GCM del ensemble a 2050
#   Cierra el flanco del "GCM único" (limitación 3 de la sección 7).
#
#   Ensamble de 6 GCMs cubriendo el rango likely IPCC + banda alta,
#   sin extremos IPCC-problemáticos (INM-CM5-0 muy bajo, hot models
#   UKESM/CanESM5). Todos disponibles a 30 arcsec en WorldClim,
#   verificado por HTTP HEAD.
#
#   GCMs y ECS aprox (equilibrium climate sensitivity, IPCC AR6):
#     MRI-ESM2-0    3.1 K   (media-baja)
#     MPI-ESM1-2-HR 3.0 K   (media)  -- YA descargado en BLOQUE 11
#     EC-Earth3-Veg 4.3 K   (alta)
#     IPSL-CM6A-LR  4.6 K   (alta)
#     ACCESS-CM2    4.7 K   (alta)
#     CNRM-CM6-1    4.9 K   (alta)
#
#   Sesgo residual reconocido: 4 de 6 GCMs tienen ECS > 4 K
#   (media likely IPCC ~3.2 K). Se declara explícitamente en 6.8.
#
#   Estrategia de espacio: descargar 1 GCM, recortar a AOI,
#   descartar el global (~9 GB) antes de bajar el siguiente.
#
# =============================================================================

library(terra)
library(biomod2)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/ENM_jubaea_chile"
VAR_DIR    <- file.path(BASE_DIR, "variables")
OUT_DIR    <- file.path(BASE_DIR, "outputs")
CMIP6_DIR  <- file.path(VAR_DIR, "cmip6_multi")
MODELS_DIR <- "C:/Users/sdari/Documents"

dir.create(CMIP6_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Config: GCMs a procesar (MPI ya está)
# -----------------------------------------------------------------------------
gcms <- tibble::tribble(
  ~modelo,           ~ECS,   ~url,
  "MRI-ESM2-0",      3.1,    "https://geodata.ucdavis.edu/cmip6/30s/MRI-ESM2-0/ssp245/wc2.1_30s_bioc_MRI-ESM2-0_ssp245_2041-2060.tif",
  "EC-Earth3-Veg",   4.3,    "https://geodata.ucdavis.edu/cmip6/30s/EC-Earth3-Veg/ssp245/wc2.1_30s_bioc_EC-Earth3-Veg_ssp245_2041-2060.tif",
  "IPSL-CM6A-LR",    4.6,    "https://geodata.ucdavis.edu/cmip6/30s/IPSL-CM6A-LR/ssp245/wc2.1_30s_bioc_IPSL-CM6A-LR_ssp245_2041-2060.tif",
  "ACCESS-CM2",      4.7,    "https://geodata.ucdavis.edu/cmip6/30s/ACCESS-CM2/ssp245/wc2.1_30s_bioc_ACCESS-CM2_ssp245_2041-2060.tif",
  "CNRM-CM6-1",      4.9,    "https://geodata.ucdavis.edu/cmip6/30s/CNRM-CM6-1/ssp245/wc2.1_30s_bioc_CNRM-CM6-1_ssp245_2041-2060.tif"
)

# El MPI ya está procesado -- lo cargamos directamente
mpi_aoi_path <- file.path(CMIP6_DIR, "MPI-ESM1-2-HR_AOI.tif")
if (!file.exists(mpi_aoi_path)) {
  cat("Recortando MPI-ESM1-2-HR (ya descargado en BLOQUE 11)...\n")
  mpi_raw <- rast(file.path(VAR_DIR, "climate/wc2.1_30s/wc2.1_30s_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif"))
  env_current <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
  aoi <- ext(env_current)
  mpi_aoi <- crop(mpi_raw, aoi)
  writeRaster(mpi_aoi, mpi_aoi_path, overwrite = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Descarga secuencial + recorte + borrado del global
# -----------------------------------------------------------------------------
for (i in seq_len(nrow(gcms))) {
  m <- gcms$modelo[i]
  aoi_path <- file.path(CMIP6_DIR, paste0(m, "_AOI.tif"))

  if (file.exists(aoi_path)) {
    cat(sprintf("[%s] recorte AOI ya existe, salto\n", m))
    next
  }

  tmp_global <- file.path(CMIP6_DIR, paste0("_TMP_", m, ".tif"))
  cat(sprintf("\n=== [%s] descargando... ===\n", m))
  t0 <- Sys.time()
  ret <- system2("curl", c("-s", "-o", shQuote(tmp_global), gcms$url[i]), timeout = 3600)
  if (ret != 0) stop("curl fallo en ", m)
  t1 <- Sys.time()
  cat(sprintf("  descarga: %.1f min\n", as.numeric(difftime(t1, t0, units = "mins"))))

  r <- rast(tmp_global)
  env_current <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))
  aoi <- ext(env_current)
  r_aoi <- crop(r, aoi)
  writeRaster(r_aoi, aoi_path, overwrite = TRUE)
  cat(sprintf("  recorte AOI guardado: %s\n", aoi_path))

  # Borrar el global para liberar espacio
  file.remove(tmp_global)
  cat(sprintf("  global borrado (liberado ~9GB)\n"))
}

# -----------------------------------------------------------------------------
# 3. Proyectar el ensemble ya entrenado sobre cada GCM
# -----------------------------------------------------------------------------
setwd(MODELS_DIR)
myBiomodEM <- get(load("Jubaea.chilensis/Jubaea.chilensis.Jubaea_chilensis_current.ensemble.models.out"))

env_current <- rast(file.path(VAR_DIR, "env_stack_vif5.tif"))

# Modelo actual y umbrales (para comparación)
prob_current <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_prob.tif"))
clase_current <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_clase.tif"))
area_km2_rast <- cellSize(clase_current, unit = "km")
km2_por_pixel <- mean(values(area_km2_rast), na.rm = TRUE)

breaks <- c(0, 0.2, 0.4, 0.6, 1.0)

# Función auxiliar: dado un GCM, arma el stack futuro, proyecta y devuelve
# cambio en superficie de "Alta" (clase 4).
proyectar_gcm <- function(modelo, aoi_path) {
  cat(sprintf("\n=== Proyectando ensemble sobre %s ===\n", modelo))
  cmip6 <- rast(aoi_path)
  names(cmip6) <- sprintf("bio%02d", 1:19)
  cmip6_al <- resample(cmip6, env_current[[c("bio02","bio03","bio08","bio14","bio15")]], method = "bilinear")

  env_future <- c(
    cmip6_al[[c("bio02","bio03","bio08","bio14","bio15")]],
    env_current[[c("slope","aspect")]]
  )
  names(env_future) <- c("bio02","bio03","bio08","bio14","bio15","slope","aspect")

  writeRaster(env_future,
              file.path(VAR_DIR, sprintf("env_stack_future_2050_%s.tif", modelo)),
              overwrite = TRUE)

  proj_name <- sprintf("FUTURE_2050_%s", gsub("-", "_", modelo))
  BIOMOD_EnsembleForecasting(
    bm.em         = myBiomodEM,
    new.env       = env_future,
    proj.name     = proj_name,
    models.chosen = "all",
    metric.binary = "TSS",
    metric.filter = "TSS"
  )

  prob_path <- file.path(MODELS_DIR, "Jubaea.chilensis",
                         sprintf("proj_%s", proj_name),
                         sprintf("proj_%s_Jubaea.chilensis_ensemble.tif", proj_name))
  prob_f <- rast(prob_path)
  if (nlyr(prob_f) > 1) prob_f <- prob_f[[1]]
  prob_f <- prob_f / 1000
  prob_f <- resample(prob_f, prob_current, method = "bilinear")

  clase_f <- classify(prob_f, rcl = cbind(breaks[-length(breaks)], breaks[-1], 1:4),
                       include.lowest = TRUE)

  tab_c <- table(factor(values(clase_current)[,1], levels = 1:4))
  tab_f <- table(factor(values(clase_f)[,1], levels = 1:4))

  tibble(
    GCM = modelo,
    km2_alta_actual = round(as.numeric(tab_c[4]) * km2_por_pixel, 0),
    km2_alta_2050   = round(as.numeric(tab_f[4]) * km2_por_pixel, 0),
    prob_media_actual = round(mean(values(prob_current), na.rm = TRUE), 3),
    prob_media_2050   = round(mean(values(prob_f), na.rm = TRUE), 3)
  ) %>%
    mutate(cambio_pct = round(100 * (km2_alta_2050 - km2_alta_actual) / km2_alta_actual, 1))
}

# Procesar los 6 GCMs (MPI + los 5 nuevos)
resultados <- list()

# MPI primero
resultados[["MPI-ESM1-2-HR"]] <- proyectar_gcm("MPI-ESM1-2-HR", mpi_aoi_path)

for (i in seq_len(nrow(gcms))) {
  m <- gcms$modelo[i]
  aoi_path <- file.path(CMIP6_DIR, paste0(m, "_AOI.tif"))
  if (file.exists(aoi_path)) {
    resultados[[m]] <- proyectar_gcm(m, aoi_path)
  }
}

# -----------------------------------------------------------------------------
# 4. Consolidación + estadísticas del ensamble
# -----------------------------------------------------------------------------
tabla <- bind_rows(resultados) %>%
  left_join(gcms %>% select(modelo, ECS) %>% add_row(modelo = "MPI-ESM1-2-HR", ECS = 3.0),
            by = c("GCM" = "modelo")) %>%
  arrange(ECS) %>%
  select(GCM, ECS, km2_alta_actual, km2_alta_2050, cambio_pct, prob_media_actual, prob_media_2050)

cat("\n\n=== TABLA CONSOLIDADA MULTI-GCM ===\n")
print(tabla, n = 20)

cambio <- tabla$cambio_pct
q <- quantile(cambio, c(0.25, 0.5, 0.75))
cat(sprintf("\n=== Estadísticas del ensamble (N=%d GCMs) ===\n", length(cambio)))
cat(sprintf("Mediana: %.1f%%\n", q[2]))
cat(sprintf("IQR:     %.1f%% a %.1f%%\n", q[1], q[3]))
cat(sprintf("Rango:   %.1f%% a %.1f%%\n", min(cambio), max(cambio)))
cat(sprintf("Media:   %.1f%%\n", mean(cambio)))

cat(sprintf("\nMPI-ESM1-2-HR (referencia individual): %.1f%%\n", tabla$cambio_pct[tabla$GCM == "MPI-ESM1-2-HR"]))

write_csv(tabla, file.path(OUT_DIR, "cmip6", "cambio_idoneidad_2050_multiGCM.csv"))

# Estadísticas resumen para usar en la memoria
resumen <- tibble(
  mediana = q[2],
  Q1 = q[1],
  Q3 = q[3],
  min = min(cambio),
  max = max(cambio),
  media = mean(cambio),
  N = length(cambio)
)
write_csv(resumen, file.path(OUT_DIR, "cmip6", "cambio_idoneidad_2050_estadisticas.csv"))

# -----------------------------------------------------------------------------
# 5. Gráfico
# -----------------------------------------------------------------------------
library(ggplot2)
tabla_plot <- tabla %>%
  arrange(cambio_pct) %>%
  mutate(GCM_lab = sprintf("%s (ECS %.1f)", GCM, ECS))

p <- ggplot(tabla_plot, aes(x = reorder(GCM_lab, cambio_pct), y = cambio_pct)) +
  geom_col(fill = "#B5502F") +
  geom_hline(yintercept = q[2], linetype = "dashed", color = "#2F5D3A", linewidth = 0.8) +
  annotate("text", x = 1, y = q[2] + 2,
           label = sprintf("Mediana %.1f%%", q[2]),
           color = "#2F5D3A", hjust = 0, size = 3.5) +
  coord_flip() +
  labs(
    title = sprintf("Cambio proyectado del hábitat de alta idoneidad a 2050 (SSP2-4.5)"),
    subtitle = sprintf("Ensamble multi-GCM (N=%d), rango ECS %.1f-%.1f K", nrow(tabla), min(tabla$ECS), max(tabla$ECS)),
    x = NULL, y = "Cambio en superficie de alta idoneidad (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(size = 10))

ggsave("C:/Users/sdari/Desktop/ENM-BGEN/extensiones/chile/resultados/cambio_idoneidad_2050_multiGCM.png",
       p, width = 9, height = 4.5, dpi = 150)
cat("\nGrafico guardado en resultados/cambio_idoneidad_2050_multiGCM.png\n")

cat("\nBLOQUE 15 completado.\n")
