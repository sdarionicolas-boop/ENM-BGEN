# Extensión Chile — *Jubaea chilensis* (Palma chilena)

Adaptación del pipeline principal ([`ENM_BGEN_pipeline.R`](../../ENM_BGEN_pipeline.R)) para cruzar la idoneidad climática de *Jubaea chilensis* con **MapBiomas Chile Colección 2 (2024)**, postulada al Premio MapBiomas Chile 2026 (1.ª Edición), categoría Análisis Técnico o Científico.

Mismo núcleo de modelado que las versiones Argentina y Perú del repo; solo cambia la especie, el área de interés y la leyenda de reclasificación de uso del suelo, adaptada a las clases de MapBiomas Chile.

## Contenido

- [`BLOQUE7_mapbiomas_chile.R`](./BLOQUE7_mapbiomas_chile.R) — cruce de idoneidad climática × MapBiomas Chile, procesando **Colección 2 (2024, resultado oficial)** y **Colección 1 (2022, anexo comparativo)**.
- [`BLOQUE7b_snaspe_chile.R`](./BLOQUE7b_snaspe_chile.R) — cruce del hábitat compatible con áreas protegidas (WDPA/SNASPE): ¿cuánto del hábitat disponible ya está protegido?
- [`BLOQUE7c_fuego_chile.R`](./BLOQUE7c_fuego_chile.R) — cruce con el Producto Fuego de MapBiomas Chile: ¿cuánto del hábitat (protegido y desprotegido) ya se quemó?
- [`BLOQUE8_ranking_palmares.R`](./BLOQUE8_ranking_palmares.R) — clustering espacial de los registros de presencia en 6 palmares y ranking de prioridad de conservación (hábitat × protección × riesgo de incendio).
- [`BLOQUE9_importancia_variables.R`](./BLOQUE9_importancia_variables.R) — extrae la importancia de variables ya calculada por biomod2 (sin reentrenar).
- [`BLOQUE10_validacion_sitios_conocidos.R`](./BLOQUE10_validacion_sitios_conocidos.R) — validación externa: idoneidad predicha en palmares documentados en la literatura.
- [`BLOQUE11_proyeccion_cmip6.R`](./BLOQUE11_proyeccion_cmip6.R) — proyección del ensemble ya entrenado a 2050 (CMIP6, SSP2-4.5), sin reentrenar.
- [`BLOQUE11b_refugios_2050.R`](./BLOQUE11b_refugios_2050.R) — cruza los refugios climáticos de 2050 con el uso de suelo actual.
- [`Memoria_Tecnica_Premio_MapBiomas_Chile.md`](./Memoria_Tecnica_Premio_MapBiomas_Chile.md) — memoria técnica completa (introducción, objetivos, metodología, resultados, discusión, y análisis de sensibilidad metodológica entre colecciones).
- [`resultados/`](./resultados/) — visores interactivos, métricas y diagnósticos:
  - [`mapa_interactivo_jubaea_chile.html`](./resultados/mapa_interactivo_jubaea_chile.html) — visor de idoneidad climática.
  - [`Col2_2024/`](./resultados/Col2_2024/) — cruce oficial con MapBiomas Chile Colección 2 (2024).
  - [`Col1_2022/`](./resultados/Col1_2022/) — cruce con Colección 1 (2022), anexo comparativo.
  - `crosstab_col1_col2_alta_idoneidad.csv` — comparación píxel a píxel entre ambas colecciones.
  - `metricas_modelo.csv` — AUCroc y TSS por algoritmo.
  - `correlacion_pearson.png`, `vif_barplot.png` — diagnósticos de colinealidad de variables ambientales.

## Resultado clave (Colección 2, oficial)

Del hábitat de alta idoneidad climática para *Jubaea chilensis*: **77% persiste como vegetación natural o es restaurable, 22,9% fue convertido** (mayormente a agricultura). Un análisis de sensibilidad entre Colección 1 y Colección 2 muestra que el 95% de esa "pérdida" no es conversión real en 2 años, sino una reclasificación más granular de la antigua clase ambigua "mosaico agropecuario" — ver la Memoria Técnica, sección 6.3.

De ese hábitat disponible, **solo el 12,6% está dentro de un área protegida (WDPA/SNASPE) — el 87,4% no tiene ninguna figura de protección legal**. Las áreas protegidas con más hábitat compatible son exactamente los palmares documentados en la literatura (La Campana-Peñuelas, Fray Jorge, Palmas de Cocalán), lo que valida independientemente el modelo — ver sección 6.4 de la Memoria Técnica.

Un tercer cruce con el **Producto Fuego** de MapBiomas Chile agrega un matiz importante: las áreas protegidas se quemaron proporcionalmente **más** que las desprotegidas (12,4% vs 7,7%), y el palmar **"Palmar El Salto"** perdió el **80% de su hábitat compatible al fuego** entre 2013-2025 — ver sección 6.5.

Un **ranking de los 6 palmares** (por clustering espacial de presencias) identifica a **Cocalán y Petorca como máxima prioridad de conservación** (mucho hábitat, casi nada protegido), mientras que Ocoa (dentro del Parque Nacional La Campana) ya está mayoritariamente resguardada — ver sección 6.7.

La variable climática más importante es la **estacionalidad de la precipitación (bio15)**, coherente con el eje de megasequía que prioriza el comité evaluador. Y una **validación externa** con sitios documentados en la literatura confirma que el modelo predice correctamente alta idoneidad en los palmares reales (Palmar El Salto y Monte Aranda: 100% de su superficie en clase "Alta") — ver secciones 5.3 y 6.1.

Una **proyección a 2050** (CMIP6, SSP2-4.5) muestra que el hábitat de alta idoneidad climática se reduciría **46,2%** para 2041-2060, incluso en un escenario de emisiones moderado — y el 22,3% de los refugios climáticos que persistirían ya está convertido hoy. Ver sección 6.8.

## Nota sobre el acceso a datos

MapBiomas Chile Colección 2 no está disponible como descarga directa pública (a diferencia de Colección 1); se accede vía Google Earth Engine. El script incluye la lógica de reclasificación pero la descarga del raster fuente requiere autenticación propia con una cuenta de Earth Engine — ver la Memoria Técnica, sección 5.4, para el asset ID y el procedimiento.

## Metodología

Ver la [Memoria Técnica](./Memoria_Tecnica_Premio_MapBiomas_Chile.md) para el detalle completo: datos de presencia (GBIF, filtrados a rango natural), variables ambientales (WorldClim + SRTM), ensemble modeling (biomod2, 5 algoritmos), y reclasificación de MapBiomas Chile en categorías Compatible/Restaurable/Incompatible/Excluido.
