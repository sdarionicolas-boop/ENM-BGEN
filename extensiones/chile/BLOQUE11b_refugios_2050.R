library(terra)
library(dplyr)
library(readr)

OUT_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs"
MB_DIR  <- file.path(OUT_DIR, "mapbiomas")

clase_future <- rast(file.path(OUT_DIR, "cmip6", "Jubaea_chilensis_idoneidad_clase_2050.tif"))
cruce <- rast(file.path(MB_DIR, "Col2_2024", "Jubaea_chilensis_cruce.tif"))

# cruce = idoneidad_actual*10 + uso_suelo (1=Compatible,2=Restaurable,3=Incompatible)
uso_suelo_actual <- cruce %% 10
uso_suelo_actual <- resample(uso_suelo_actual, clase_future, method = "near")

alta_2050 <- clase_future == 4
alta_actual_tbl <- values(alta_2050)[, 1]
uso_tbl <- values(uso_suelo_actual)[, 1]

df <- tibble(alta_2050 = alta_actual_tbl, uso = uso_tbl) %>% filter(alta_2050 == 1, !is.na(uso))

tab <- df %>%
  mutate(categoria = case_when(
    uso == 1 ~ "Compatible (vegetación natural HOY)",
    uso == 2 ~ "Restaurable",
    uso == 3 ~ "Incompatible (ya convertido HOY)"
  )) %>%
  count(categoria) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("=== De los refugios climáticos de Alta idoneidad en 2050, uso de suelo ACTUAL ===\n")
print(tab)

write_csv(tab, file.path(OUT_DIR, "cmip6", "refugios_2050_uso_actual.csv"))
cat("\nBLOQUE 11b completado.\n")
