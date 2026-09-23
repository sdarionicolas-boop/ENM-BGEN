library(terra)
library(dplyr)
library(readr)

BASE_DIR   <- "C:/Users/sdari/Desktop/BGEN/ENM_jacaranda"
OUT_OFICIAL <- file.path(BASE_DIR, "outputs")
OUT_EMMEAN  <- file.path(BASE_DIR, "outputs_emmean")

pres_thin <- read_csv(file.path(BASE_DIR, "presencias_thin.csv"), show_col_types = FALSE)
especies  <- unique(pres_thin$species)

resumen <- list()

for (sp in especies) {
  sp_name <- gsub(" ", "_", sp)

  tif_oficial <- file.path(OUT_OFICIAL, sp_name, paste0(sp_name, "_idoneidad_clase.tif"))
  tif_emmean  <- file.path(OUT_EMMEAN,  sp_name, paste0(sp_name, "_idoneidad_clase_EMmean.tif"))

  if (!file.exists(tif_oficial) || !file.exists(tif_emmean)) {
    cat("SKIP (falta archivo):", sp, "\n")
    next
  }

  r_of <- rast(tif_oficial)
  r_em <- rast(tif_emmean)

  px_area_km2 <- prod(res(r_of)) / 1e6  # asume CRS en metros; si es geográfico, ajustar

  area_of <- freq(r_of) %>% as.data.frame()
  area_em <- freq(r_em) %>% as.data.frame()

  alto_of <- sum(area_of$count[area_of$value == 4], na.rm = TRUE)
  alto_em <- sum(area_em$count[area_em$value == 4], na.rm = TRUE)

  cambio_pct <- if (alto_of > 0) round(100 * (alto_em - alto_of) / alto_of, 2) else NA

  resumen[[sp]] <- data.frame(
    species = sp,
    pixeles_alto_oficial = alto_of,
    pixeles_alto_emmean  = alto_em,
    cambio_pct_area_alta = cambio_pct
  )

  cat(sprintf("%-30s oficial=%d  emmean=%d  cambio=%s%%\n", sp, alto_of, alto_em, cambio_pct))
}

resumen_df <- bind_rows(resumen)
write_csv(resumen_df, file.path(BASE_DIR, "outputs_emmean", "comparacion_area_alta_emmean_vs_oficial.csv"))
cat("\nMediana cambio %:", median(resumen_df$cambio_pct_area_alta, na.rm = TRUE), "\n")
cat("Rango cambio %:", min(resumen_df$cambio_pct_area_alta, na.rm = TRUE), "a", max(resumen_df$cambio_pct_area_alta, na.rm = TRUE), "\n")
print(resumen_df)
