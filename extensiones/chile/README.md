# Extensión Chile — *Jubaea chilensis* (Palma chilena)

Adaptación del pipeline principal ([`ENM_BGEN_pipeline.R`](../../ENM_BGEN_pipeline.R)) para cruzar la idoneidad climática de *Jubaea chilensis* con **MapBiomas Chile Colección 1 (2022)**, postulada al Premio MapBiomas Chile 2026 (1.ª Edición), categoría Análisis Técnico o Científico.

Mismo núcleo de modelado que las versiones Argentina y Perú del repo; solo cambia la especie, el área de interés y la leyenda de reclasificación de uso del suelo, adaptada a las clases de MapBiomas Chile.

## Contenido

- [`BLOQUE7_mapbiomas_chile.R`](./BLOQUE7_mapbiomas_chile.R) — cruce de idoneidad climática × MapBiomas Chile Colección 1 (2022).
- [`Memoria_Tecnica_Premio_MapBiomas_Chile.md`](./Memoria_Tecnica_Premio_MapBiomas_Chile.md) — memoria técnica completa (introducción, objetivos, metodología, resultados, discusión).
- [`resultados/`](./resultados/) — visores interactivos, métricas y diagnósticos:
  - [`mapa_interactivo_jubaea_chile.html`](./resultados/mapa_interactivo_jubaea_chile.html) — visor de idoneidad climática.
  - [`mapa_interactivo_cruce_chile.html`](./resultados/mapa_interactivo_cruce_chile.html) — visor del cruce idoneidad × MapBiomas.
  - `Jubaea_chilensis_cruce_resumen.csv`, `cruce_resumen_todas_especies.csv` — superficie (km²) por categoría.
  - `metricas_modelo.csv` — AUCroc y TSS por algoritmo.
  - `correlacion_pearson.png`, `vif_barplot.png` — diagnósticos de colinealidad de variables ambientales.

## Resultado clave

Del hábitat de alta idoneidad climática para *Jubaea chilensis*: **92% persiste como vegetación natural o es restaurable, solo 4,3% fue convertido**. A diferencia de Argentina, la pérdida crítica de esta especie (2,5% de su población histórica) no se explica por conversión de uso de suelo, sino por presión directa sobre los individuos — ver la Memoria Técnica para el análisis completo.

## Metodología

Ver la [Memoria Técnica](./Memoria_Tecnica_Premio_MapBiomas_Chile.md) para el detalle completo: datos de presencia (GBIF, filtrados a rango natural), variables ambientales (WorldClim + SRTM), ensemble modeling (biomod2, 5 algoritmos), y reclasificación de MapBiomas Chile Colección 1 en categorías Compatible/Restaurable/Incompatible/Excluido.
