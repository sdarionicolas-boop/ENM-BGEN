# Extensión Perú — *Mauritia flexuosa* (Aguaje)

Adaptación del pipeline principal ([`ENM_BGEN_pipeline.R`](../../ENM_BGEN_pipeline.R)) para cruzar la idoneidad climática de *Mauritia flexuosa* con **MapBiomas Perú 2024 (Colección 3)**, postulada al Premio MapBiomas Perú 2026.

Mismo núcleo de modelado que la versión Argentina del repo; solo cambia la leyenda de reclasificación de uso del suelo, adaptada a las clases de MapBiomas Perú.

## Contenido

- [`BLOQUE7_mapbiomas_peru.R`](./BLOQUE7_mapbiomas_peru.R) — cruce de idoneidad climática × MapBiomas Perú 2024.
- [`Memoria_Tecnica_Premio_MapBiomas.md`](./Memoria_Tecnica_Premio_MapBiomas.md) — memoria técnica de la postulación.
- [`resultados/`](./resultados/) — visor interactivo, estadísticas y cartografía de salida:
  - [`mapa_interactivo_cruce_peru.html`](./resultados/mapa_interactivo_cruce_peru.html) — visor web interactivo del cruce.
  - `Mauritia_flexuosa_cruce_resumen.csv`, `resumen_todas_especies.csv` — estadísticas de superficie por categoría.
  - `screenshot_cruce.png`, `presentacion_inversor.png` — capturas de la cartografía automatizada (QGIS).

## Metodología

Ver la sección 5 de la [Memoria Técnica](./Memoria_Tecnica_Premio_MapBiomas.md) para el detalle completo. En resumen: reclasificación de MapBiomas Perú 2024 en tres categorías (Compatible / Restaurable / Incompatible), cruce algebraico con la capa de idoneidad climática (ensemble GLM+GBM+RF+MAXNET+XGBoost), y cuantificación de superficie (km²) por categoría.
