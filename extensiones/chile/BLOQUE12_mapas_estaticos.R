library(terra)
library(dplyr)
library(readr)
library(ggplot2)

OUT_DIR <- "C:/Users/sdari/Desktop/ENM_jubaea_chile/outputs"
MB_DIR  <- file.path(OUT_DIR, "mapbiomas")
DEST    <- "C:/Users/sdari/Desktop/ENM-BGEN/extensiones/chile/resultados"

# -----------------------------------------------------------------------------
# 1. Mapa estático: idoneidad climática
# -----------------------------------------------------------------------------
idon <- rast(file.path(OUT_DIR, "Jubaea_chilensis", "Jubaea_chilensis_idoneidad_clase.tif"))
idon_df <- as.data.frame(idon, xy = TRUE)
names(idon_df)[3] <- "clase"
idon_df$clase <- factor(idon_df$clase, levels = 1:4,
                         labels = c("Insustentable", "Bajo", "Moderado", "Alto"))

p1 <- ggplot(idon_df, aes(x = x, y = y, fill = clase)) +
  geom_raster() +
  scale_fill_manual(values = c("Insustentable" = "#D9D9D9", "Bajo" = "#FEE08B",
                                "Moderado" = "#F46D43", "Alto" = "#1A9850"),
                     name = "Idoneidad") +
  coord_equal() +
  labs(title = "Idoneidad climática — Jubaea chilensis",
       subtitle = "Ensemble modeling (biomod2), AUCroc 0,948 · TSS 0,741",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "right")

ggsave(file.path(DEST, "mapa_idoneidad_climatica.png"), p1, width = 6, height = 8, dpi = 150)
cat("Mapa 1 guardado.\n")

# -----------------------------------------------------------------------------
# 2. Mapa estático: cruce idoneidad x uso de suelo (Colección 2, oficial)
# -----------------------------------------------------------------------------
cruce <- rast(file.path(MB_DIR, "Col2_2024", "Jubaea_chilensis_cruce.tif"))
cruce_df <- as.data.frame(cruce, xy = TRUE)
names(cruce_df)[3] <- "codigo"

etiquetas <- c(
  "41" = "Alta — Compatible", "42" = "Alta — Restaurable", "43" = "Alta — Perdido",
  "31" = "Moderada — Compatible", "32" = "Moderada — Restaurable", "33" = "Moderada — Perdido",
  "11" = "Baja — Compatible", "12" = "Baja — Restaurable", "13" = "Baja — Perdido"
)
colores <- c(
  "Alta — Compatible" = "#1A9850", "Alta — Restaurable" = "#2166AC", "Alta — Perdido" = "#D73027",
  "Moderada — Compatible" = "#A6D96A", "Moderada — Restaurable" = "#74ADD1", "Moderada — Perdido" = "#F46D43",
  "Baja — Compatible" = "#E8E8E8", "Baja — Restaurable" = "#D9D9D9", "Baja — Perdido" = "#CCCCCC"
)
cruce_df$categoria <- etiquetas[as.character(cruce_df$codigo)]
cruce_df$categoria <- factor(cruce_df$categoria, levels = names(colores))
cruce_df <- cruce_df %>% filter(!is.na(categoria))

p2 <- ggplot(cruce_df, aes(x = x, y = y, fill = categoria)) +
  geom_raster() +
  scale_fill_manual(values = colores, name = NULL) +
  coord_equal() +
  labs(title = "Idoneidad × Uso del suelo — MapBiomas Chile Colección 2 (2024)",
       subtitle = "77% del hábitat de alta idoneidad disponible o restaurable · 22,9% perdido",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        legend.text = element_text(size = 8))

ggsave(file.path(DEST, "mapa_cruce_mapbiomas.png"), p2, width = 7, height = 8, dpi = 150)
cat("Mapa 2 guardado.\n")

# -----------------------------------------------------------------------------
# 3. Gráfico: ranking de palmares
# -----------------------------------------------------------------------------
ranking <- read_csv(file.path(MB_DIR, "Col2_2024", "ranking_palmares.csv"), show_col_types = FALSE)
ranking$etiqueta <- paste0(ranking$rango_prioridad, ". ", ranking$nombre_tentativo)

p3 <- ggplot(ranking, aes(x = reorder(etiqueta, indice_prioridad), y = indice_prioridad)) +
  geom_col(fill = "#B5502F") +
  coord_flip() +
  labs(title = "Ranking de prioridad de conservación — 6 palmares",
       subtitle = "Índice: 50% hábitat compatible + 30% desprotección + 20% riesgo de incendio",
       x = NULL, y = "Índice de prioridad") +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(size = 10))

ggsave(file.path(DEST, "ranking_palmares.png"), p3, width = 9.5, height = 4.5, dpi = 150)
cat("Mapa 3 (ranking) guardado.\n")

cat("\nBLOQUE 12 completado.\n")
