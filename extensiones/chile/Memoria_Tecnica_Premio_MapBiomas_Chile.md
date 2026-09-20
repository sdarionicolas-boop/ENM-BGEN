# Memoria Técnica - Premio MapBiomas Chile 2026 (1.ª Edición)

## 1. Título del Proyecto
**Idoneidad climática y disponibilidad real de hábitat para *Jubaea chilensis* (Palma chilena) mediante la integración de modelado de nicho ecológico con tres productos de MapBiomas Chile: cobertura de suelo, áreas protegidas y fuego.**

## 2. Autor(es)
Darío Nicolás Sánchez Leguizamón (Beca BIEI 2025 | Proyecto ARG/19/G24 | Red ARGENA)

## 3. Resumen Ejecutivo
*Jubaea chilensis* (palma chilena) es una palmera endémica de Chile central catalogada como Vulnerable por la UICN, con estudios recientes (2024) que argumentan su reclasificación a En Peligro Crítico. De la población preколombina estimada, hoy solo subsisten ~121.284 individuos (2,5%), concentrados en seis palmares relictos entre La Serena y el Maule.

Este trabajo aplica el mismo pipeline de ensemble modeling (biomod2) ya utilizado para 19 especies del BGEN-UNAJ en Argentina y para *Mauritia flexuosa* en Perú, ahora sobre *Jubaea chilensis*, integrando tres capas complementarias: **MapBiomas Chile Colección 2 (2024)** — el producto de land cover exigido por el requisito de elegibilidad N.º 1 del reglamento del Premio — para cuantificar disponibilidad real de hábitat; el **WDPA/SNASPE** para medir qué fracción de ese hábitat está formalmente protegida; y el **Producto Fuego Colección 1** de MapBiomas Chile para cuantificar exposición histórica a incendios. Como anexo comparativo se incluye además el mismo cruce de cobertura con la Colección 1 (2022), lo que permite evaluar la sensibilidad metodológica del resultado entre versiones sucesivas del producto MapBiomas.

Los cruces convergen en un mensaje único: el 77% del hábitat climáticamente apto sigue disponible, pero solo el 12,6% está protegido, y ni siquiera la protección formal garantiza resguardo frente a incendios — un palmar protegido específicamente para la especie ("Palmar El Salto") perdió el 80% de su hábitat compatible al fuego entre 2013 y 2025. Una proyección a futuro (CMIP6, SSP2-4.5, 2041-2060) agrega la dimensión más urgente: ese hábitat se reduciría un 46,2% adicional por cambio climático, incluso bajo un escenario de emisiones moderado.

## 4. Problema que Resuelve (Relevancia)
La literatura reciente atribuye el colapso poblacional de *Jubaea chilensis* a presión directa sobre los individuos (extracción de savia de palmeras decapitadas, sobrecosecha ilegal de semillas, depredación de semillas por *Rattus rattus* y de plántulas por conejos exóticos) más que a pérdida de hábitat por sí sola. Sin embargo, esta hipótesis no había sido contrastada cuantitativamente contra datos oficiales de cobertura de suelo.

Este estudio responde una pregunta concreta y accionable: **¿el hábitat climáticamente apto para la especie todavía existe como vegetación natural, o ya fue perdido por conversión de uso de suelo?** La respuesta determina si la estrategia de conservación debe priorizar la protección de individuos/semillas (si el hábitat está intacto) o además la restauración de hábitat (si está convertido).

## 5. Metodología

### 5.1. Datos de presencia
Se integraron 2.117 registros de ocurrencia de GBIF, filtrados al rango de distribución natural documentado en la literatura (29,9°S–35,3°S, franja costera mediterránea de Coquimbo a Maule), descartando 26 registros correspondientes a ejemplares cultivados fuera de rango (jardines y plazas, donde la especie se planta ampliamente como ornamental). Resultaron 2.091 registros válidos, sometidos a *thinning* espacial (1 registro por píxel).

### 5.2. Variables ambientales
19 variables bioclimáticas de WorldClim v2.1 y topografía SRTM (elevación, pendiente, aspecto) para Chile, reducidas mediante análisis VIF (≤5) para eliminar colinealidad.

### 5.3. Modelado de nicho ecológico (Ensemble Modeling)
Cinco algoritmos (GLM, GBM, RF, MAXNET, XGBoost) en R (`biomod2`), validación cruzada 80/20 con 3 repeticiones, ensamble ponderado por TSS. Métricas de validación:

| Métrica | Promedio | Rango |
|---------|----------|-------|
| AUCroc  | 0,948    | 0,925 – 0,974 |
| TSS     | 0,741    | 0,698 – 0,769 |

Ambas muy por encima de los umbrales de calidad del pipeline (TSS ≥ 0,7, AUCroc ≥ 0,9).

**Importancia de variables** (promedio ± DE entre los 5 algoritmos, ya calculada durante el ajuste del modelo, sin reentrenamiento adicional):

| Variable | Importancia media | Descripción |
|---|---:|---|
| bio15 | 0,506 ± 0,149 | Estacionalidad de la precipitación (coeficiente de variación) |
| bio08 | 0,269 ± 0,080 | Temperatura media del trimestre más lluvioso |
| bio14 | 0,090 ± 0,065 | Precipitación del mes más seco |
| bio03 | 0,075 ± 0,054 | Isotermalidad |
| bio02 | 0,053 ± 0,040 | Rango diurno medio de temperatura |
| slope | 0,049 ± 0,051 | Pendiente |
| aspect | 0,019 ± 0,032 | Orientación |

**La estacionalidad de la precipitación (bio15) domina ampliamente la idoneidad climática de la especie** — más del doble de importancia que la segunda variable. Esto es coherente con el eje de "megasequía y valles centrales" que el propio reglamento del Premio identifica como dinámica crítica para Chile: *Jubaea chilensis* no responde tanto a cuánta lluvia cae en total, sino a qué tan irregular es su distribución estacional, la métrica que la megasequía altera de forma más directa. Gráfico completo en [`resultados/importancia_variables_jubaea.png`](./resultados/importancia_variables_jubaea.png).

### 5.4. Integración con MapBiomas Chile Colección 2 (2024)
El reglamento del Premio MapBiomas Chile 2026 (requisito de elegibilidad N.º 1) exige el uso del **"Producto land cover colección 2"**. A diferencia de Argentina y Perú, la Colección 2 de Chile (publicada en 2025, serie 1999-2024) **no está disponible como descarga directa en el bucket público** de MapBiomas — su acceso oficial es exclusivamente vía **Google Earth Engine**, mediante el asset:

```
projects/mapbiomas-chile/assets/LULC/COLLECTION-02/CLASSIFICATIONS/classification-final/clasificacion-final-2
```

Se extrajo la banda `classification_2024` para el área de interés mediante la API de Python de Earth Engine (autenticación OAuth propia), exportada en 10 teselas (por límite de tamaño de descarga directa de GEE) y mosaicada con `terra::merge()` en R.

La capa se reclasificó según la leyenda oficial de MapBiomas Chile en tres dimensiones operativas, siguiendo el mismo criterio ya aplicado a Argentina y Perú:

| Categoría | Clases MapBiomas Chile | Códigos |
|-----------|------------------------|---------|
| **Compatible** | Bosque, humedal, pastizal natural, estepa, matorral | 3, 59, 60, 67, 11, 12, 63, 66 |
| **Restaurable** | Silvicultura, pastura, mosaico agropecuario | 9, 15, 21 |
| **Incompatible** | Agricultura, infraestructura, área degradada | 18, 24, 25 |
| *Excluido* | Agua, hielo/nieve, afloramiento rocoso, arena/duna, salar, sin datos | 0, 23, 29, 33, 34, 61 |

La categoría *Excluido* agrupa superficies naturalmente no vegetadas (alta cordillera, cuerpos de agua) que quedan fuera del esquema de conversión antrópica y no se computan como "pérdida" — una precisión metodológica necesaria en un país con una fracción tan grande de su territorio en la criósfera y zonas áridas, a diferencia de Argentina o Perú.

Luego se realizó el cruce algebraico con la capa de idoneidad climática reclasificada (Alta/Moderada/Baja), obteniendo 9 combinaciones posibles por celda. Como anexo, se repitió el mismo procedimiento con **Colección 1 (2022)** (descarga directa pública) para evaluar la estabilidad temporal/metodológica del resultado (ver sección 6.4).

## 6. Resultados

### 6.1. Validación cruzada con sitios conocidos

Las métricas AUCroc/TSS (sección 5.3) validan el modelo contra su propia muestra de entrenamiento/prueba (80/20). Como control independiente, se extrajo la idoneidad climática predicha dentro de sitios documentados como palmares reales de *Jubaea chilensis* — información externa al ajuste del modelo, ya que no se usó ningún límite administrativo en el entrenamiento:

| Sitio | Idoneidad media (prob.) | % superficie en clase "Alta" |
|---|---:|---:|
| Palmar El Salto | 0,788 | 100% |
| Área de Palma Chilena de Monte Aranda | 0,699 | 100% |
| La Campana – Peñuelas | 0,675 | 72,1% |
| Bosque Fray Jorge (RB) | 0,487 | 32,4% |
| Palmas de Cocalán | 0,453 | 33,3% |
| *Referencia: promedio en toda la AOI* | *0,249* | — |

**Los cinco sitios superan entre 1,8 y 3,2 veces la idoneidad media de referencia de toda el área de estudio**, y los dos palmares creados específicamente para la especie (Palmar El Salto, Monte Aranda) muestran 100% de su superficie en la clase de idoneidad más alta. Esto constituye una validación externa independiente de las métricas estadísticas internas: el modelo no solo ajusta bien sus propios datos de entrenamiento, sino que identifica correctamente hábitat real verificado en terreno por otros estudios. Tabla completa en [`resultados/Col2_2024/validacion_sitios_conocidos.csv`](./resultados/Col2_2024/validacion_sitios_conocidos.csv).

### 6.2. Superficie de hábitat de Alta Idoneidad Climática — Colección 2 (2024, resultado oficial)

| Categoría | Superficie | % |
|-----------|-----------:|---:|
| Compatible (conservación) | 14.937 km² | 71,9% |
| Restaurable (restauración) | 1.058 km² | 5,1% |
| Incompatible (hábitat perdido) | 4.765 km² | 22,9% |
| **Total idoneidad alta** | **20.760 km²** | 100% |

**El 77% del hábitat climáticamente óptimo para *Jubaea chilensis* persiste como vegetación natural o es restaurable; el 22,9% fue efectivamente convertido.**

### 6.3. Interpretación
Este resultado es relevante porque, aun siendo más conservador que una primera exploración con datos de 2022 (ver 6.3), sigue indicando que **la mayor parte del hábitat climáticamente apto está disponible**. La brecha entre 121.284 individuos actuales y el hábitat potencial disponible (miles de km²) señala que el factor limitante no es únicamente la disponibilidad de tierra, sino también la presión directa documentada en la literatura: extracción de savia, sobrecosecha de semillas y depredación por fauna exótica que colapsa la dispersión y el reclutamiento de nuevos individuos. La conservación de esta especie requiere, por lo tanto, una estrategia combinada: proteger tanto el 22,9% de hábitat ya convertido de mayor conversión adicional como los individuos reproductivos remanentes.

### 6.4. Sensibilidad metodológica: comparación Colección 1 (2022) vs. Colección 2 (2024)

Como control de robustez, se repitió el cruce con la Colección 1 (2022), obteniendo un resultado más optimista (92% disponible/restaurable, 4,3% perdido). Para entender esta diferencia se realizó un **cruce píxel a píxel** entre ambas colecciones dentro de la zona de alta idoneidad climática, en vez de asumir que se trata de conversión real de hábitat en dos años.

**Hallazgo:** de los 8.646 píxeles clasificados como "Mosaico Agropecuario" (código 21, categoría Restaurable) en Colección 1, la Colección 2 —que discontinuó esa clase ambigua y reclasifica con mayor granularidad— asigna esos mismos píxeles físicos a:

| Destino en Colección 2 | Píxeles | % | Categoría resultante |
|---|---:|---:|---|
| Agricultura (18) | 4.894 | 56,6% | Incompatible |
| Matorral (66) | 2.785 | 32,2% | Compatible |
| Pastizal (12) | 301 | 3,5% | Compatible |
| Infraestructura (24) | 217 | 2,5% | Incompatible |
| Pastura (15) | 172 | 2,0% | Restaurable |
| Silvicultura (9) | 109 | 1,3% | Restaurable |

El aumento neto de "hábitat perdido" entre colecciones es de 5.361 píxeles; **5.111 de ellos (95,3%) provienen directamente de esta reclasificación de la clase mosaico**, no de conversión real de uso de suelo ocurrida entre 2022 y 2024. Es decir: la Colección 2 no muestra que Chile perdió hábitat de palma en dos años — muestra que una fracción de lo que antes se veía como "uso mixto/restaurable" ambiguo es, con mayor resolución de clasificación, agricultura activa. Este es un resultado más preciso, no un empeoramiento del hábitat.

Esta comparación se documenta también como archivo de datos en [`resultados/crosstab_col1_col2_alta_idoneidad.csv`](./resultados/crosstab_col1_col2_alta_idoneidad.csv).

### 6.5. Brecha de protección: ¿cuánto del hábitat disponible está protegido?

El resultado de la sección 6.2 responde "cuánto hábitat hay disponible", pero no "cuánto de ese hábitat está legalmente resguardado de conversión futura". Para responder esto se cruzó la capa "Compatible × Alta idoneidad" (14.937 km²) con el **WDPA (World Database on Protected Areas, UNEP-WCMC/IUCN)**, filtrado a las 76 áreas protegidas de Chile que intersectan la AOI (Parques Nacionales, Reservas Nacionales y Monumentos Naturales del SNASPE/CONAF, Santuarios de la Naturaleza del MMA, sitios Ramsar y Reservas de la Biósfera UNESCO-MAB), obtenido vía Google Earth Engine.

| Categoría | Superficie | % |
|---|---:|---:|
| Compatible, alta idoneidad — **protegido** | 1.889 km² | 12,6% |
| Compatible, alta idoneidad — **desprotegido** | 13.048 km² | 87,4% |

**Solo el 12,6% del hábitat climáticamente óptimo y aún disponible está dentro de un área protegida. El 87,4% restante no tiene ninguna figura de protección legal** y queda expuesto a la misma presión agrícola que ya convirtió el 22,9% del hábitat históricamente apto (sección 6.2).

El detalle por área protegida (tabla completa en [`resultados/Col2_2024/snaspe_detalle_por_area.csv`](./resultados/Col2_2024/snaspe_detalle_por_area.csv)) confirma además la validez del modelo: las áreas que concentran más hábitat compatible son exactamente los sitios documentados en la literatura como palmares relictos —

| Área protegida | Hábitat compatible |
|---|---:|
| Parque Nacional La Campana – Peñuelas | 1.358 km² |
| Reserva de la Biósfera Bosque Fray Jorge | 397 km² |
| Reserva Nacional Lago Peñuelas | 76 km² |
| Palmas de Cocalán | 10 km² |
| Área de Palma Chilena de Monte Aranda | 5 km² |
| Palmar El Salto | 4 km² |

— es decir, el modelo de idoneidad climática predice correctamente alta idoneidad justo donde existen áreas protegidas creadas específicamente para esta especie ("Palmas de Cocalán", "Área de Palma Chilena de Monte Aranda", "Palmar El Salto"), lo que funciona como una validación adicional independiente de las métricas estadísticas (AUCroc/TSS) de la sección 5.3.

**Implicancia para gestión territorial:** la prioridad de conservación no es solo evitar más conversión de uso de suelo, sino ampliar la red de protección formal sobre el 87,4% de hábitat compatible que hoy no tiene ninguna figura legal — particularmente fuera de los núcleos ya protegidos de La Campana-Peñuelas y Fray Jorge.

### 6.6. Incendios: un tercer producto MapBiomas Chile aplicado al mismo hábitat

El reglamento del Premio también habilita el **Producto Fuego Colección 1** como fuente de datos elegible. Se utilizó la capa de **frecuencia de área quemada acumulada 2013-2025** (descarga directa pública) para cuantificar cuánto del hábitat compatible de alta idoneidad ya fue afectado por incendios — la misma amenaza que la literatura (sección 4) identifica como recurrente para la especie, ahora con datos espaciales concretos en vez de solo referencia bibliográfica.

| Estado de protección | Hábitat compatible | Quemado alguna vez (2013-2025) | % quemado |
|---|---:|---:|---:|
| Desprotegido | 13.048 km² | 1.010 km² | 7,7% |
| Protegido | 1.889 km² | 234 km² | 12,4% |
| **Total** | **14.937 km²** | **1.245 km²** | **8,3%** |

**Hallazgo relevante:** las áreas *protegidas* se quemaron proporcionalmente **más** que las desprotegidas (12,4% vs 7,7%). Esto no es contradictorio: la vegetación nativa continua que predomina dentro de las áreas protegidas (matorral, bosque esclerófilo) constituye mayor carga de combustible que el mosaico agrícola/urbano circundante. La implicancia práctica es directa: **la protección legal por sí sola no equivale a protección contra incendios** — se requiere manejo activo del combustible y cortafuegos dentro del propio sistema de áreas protegidas, no solo ampliar su superficie.

El detalle por área protegida (tabla completa en [`resultados/Col2_2024/fuego_por_area_protegida.csv`](./resultados/Col2_2024/fuego_por_area_protegida.csv)) identifica un caso crítico puntual:

| Área protegida | Hábitat compatible | % quemado 2013-2025 |
|---|---:|---:|
| **Palmar El Salto** | 3,6 km² | **80%** |
| Quebrada de La Plata | 1,4 km² | 50% |
| Lago Peñuelas | 75,8 km² | 44,8% |
| Cerro Santa Inés | 7,2 km² | 20% |
| La Campana – Peñuelas | 1.358 km² | 14,3% |
| Palmas de Cocalán | 10,1 km² | 14,3% |
| Fray Jorge, Monte Aranda, Roblería del Cobre, San Juan de Piche | — | 0% |

**"Palmar El Salto"** —una reserva creada específicamente para *Jubaea chilensis*— tuvo el **80% de su hábitat compatible quemado** en el período analizado. Es el hallazgo más urgente y accionable de este estudio: identifica un sitio puntual, pequeño y con nombre propio donde la intervención de manejo de fuego debería priorizarse de inmediato, en contraste con palmares como Fray Jorge o Monte Aranda que no registran quemas en el mismo período y podrían servir de referencia de buen manejo.

### 6.7. Ranking de palmares por prioridad de conservación

La literatura (PMC9370131) documenta seis poblaciones/palmares principales —Ocoa, Cocalán, Viña del Mar/Valparaíso, Candelaria, Petorca y Culimo— que concentran el 96% de los individuos, pero no publica sus coordenadas exactas en formato tabulado. Para construir un ranking accionable sin inventar límites administrativos, se agruparon los 2.091 registros de presencia mediante **clustering espacial (k-means, k=6)** y se calculó, para cada clúster, un índice de prioridad de conservación análogo al usado en `BLOQUE8_partidos_bonaerenses.R` (Argentina): combina hábitat compatible disponible (peso 50%), grado de desprotección (peso 30%) y exposición a incendios (peso 20%).

Cada clúster se etiquetó con la población documentada más cercana, indicando explícitamente el nivel de confianza de esa identificación (por coincidencia con polígonos WDPA conocidos, o por toponimia regional cuando no hay una referencia espacial exacta):

| Rango | Palmar (aprox.) | Confianza | Hábitat compatible | % Protegido | % Quemado | Índice prioridad |
|---|---|---|---:|---:|---:|---:|
| 1 | Cocalán | Media | 2.646 km² | 1,3% | 13,0% | 0,862 |
| 2 | Petorca | Media | 2.798 km² | 0,3% | 0% | 0,799 |
| 3 | Candelaria | Baja | 501 km² | 0% | 25,6% | 0,500 |
| 4 | Viña del Mar/Valparaíso | Alta | 1.120 km² | 38,9% | 18,0% | 0,391 |
| 5 | Culimo | Media | 581 km² | 0% | 0% | 0,317 |
| 6 | Ocoa (La Campana) | Alta | 1.074 km² | 63,2% | 7,6% | 0,184 |

**Ocoa queda en el último lugar de prioridad** precisamente porque ya está mayoritariamente protegida (63,2%, dentro del Parque Nacional La Campana-Peñuelas) — es la población que menos necesita intervención adicional. **Cocalán y Petorca encabezan el ranking**: concentran mucho hábitat compatible con protección casi nula, y deberían ser el foco de nuevas iniciativas de conservación formal. El clúster "Viña del Mar/Valparaíso" merece una salvedad: al coincidir con una zona urbana densamente muestreada por ciencia ciudadana (iNaturalist), su alto número de registros (1.422, el mayor de los seis) probablemente incluye ejemplares cultivados en plazas y jardines, no solo poblaciones silvestres — la misma limitación de sesgo de muestreo urbano discutida en la sección 5.1.

Tabla completa con centroides y metodología de etiquetado en [`resultados/Col2_2024/ranking_palmares.csv`](./resultados/Col2_2024/ranking_palmares.csv).

### 6.8. Proyección a futuro: cambio climático (CMIP6, 2050, SSP2-4.5)

Todos los resultados anteriores describen el presente. Como último análisis, se proyectó el ensemble ya entrenado (sin reentrenar) sobre un escenario climático futuro: **CMIP6, modelo MPI-ESM1-2-HR, SSP2-4.5 (escenario de emisiones moderado), horizonte 2041-2060**, manteniendo fijas las variables topográficas (pendiente, orientación) y actualizando las 5 variables bioclimáticas retenidas por VIF (bio02, bio03, bio08, bio14, bio15) a sus valores proyectados por WorldClim.

| Clase de idoneidad | Superficie actual | Superficie 2050 | Cambio |
|---|---:|---:|---:|
| Alto | 21.307 km² | 11.467 km² | **-46,2%** |
| Moderado | 13.945 km² | 14.388 km² | +3,2% |
| Bajo | 13.013 km² | 16.650 km² | +28,0% |
| Insustentable | 74.553 km² | 80.311 km² | +7,7% |

**La superficie de alta idoneidad climática caería un 46,2% para 2041-2060**, incluso bajo un escenario de emisiones moderado (SSP2-4.5, no el más pesimista). La idoneidad media de toda la AOI cae de 0,249 a 0,204 (-18,3%). Este es el hallazgo más severo del estudio: una especie que ya perdió el 97,5% de su población histórica, con solo 12,6% de su hábitat actual protegido, enfrenta además una contracción climática de casi la mitad de su hábitat óptimo en dos a tres décadas.

Un cruce adicional entre los refugios climáticos que persistirían en 2050 y el uso de suelo **actual** matiza el panorama: del área que seguiría siendo de alta idoneidad en 2050, el 75,1% es hoy vegetación natural disponible, pero el **22,3% ya está convertido** — es decir, incluso antes de llegar a 2050, una porción significativa de los futuros refugios climáticos ya no sería utilizable por la especie. Esto refuerza la urgencia de proteger hábitat compatible *ahora*: parte de lo que el clima futuro seguiría permitiendo, el uso de suelo presente ya lo está descartando.

Detalle completo en [`resultados/cambio_idoneidad_2050.csv`](./resultados/cambio_idoneidad_2050.csv) y [`resultados/refugios_2050_uso_actual.csv`](./resultados/refugios_2050_uso_actual.csv).

### 6.9. Visualización de resultados
- Visor interactivo de idoneidad climática: [`resultados/mapa_interactivo_jubaea_chile.html`](./resultados/mapa_interactivo_jubaea_chile.html)
- Visor interactivo del cruce, Colección 2 (oficial): [`resultados/Col2_2024/mapa_interactivo_cruce_chile.html`](./resultados/Col2_2024/mapa_interactivo_cruce_chile.html)
- Visor interactivo del cruce, Colección 1 (anexo): [`resultados/Col1_2022/mapa_interactivo_cruce_chile.html`](./resultados/Col1_2022/mapa_interactivo_cruce_chile.html)
- Brecha de protección (SNASPE/WDPA): [`resultados/Col2_2024/snaspe_proteccion_resumen.csv`](./resultados/Col2_2024/snaspe_proteccion_resumen.csv), [`resultados/Col2_2024/snaspe_detalle_por_area.csv`](./resultados/Col2_2024/snaspe_detalle_por_area.csv)
- Incendios 2013-2025: [`resultados/Col2_2024/fuego_resumen.csv`](./resultados/Col2_2024/fuego_resumen.csv), [`resultados/Col2_2024/fuego_por_area_protegida.csv`](./resultados/Col2_2024/fuego_por_area_protegida.csv)
- Ranking de palmares por prioridad: [`resultados/Col2_2024/ranking_palmares.csv`](./resultados/Col2_2024/ranking_palmares.csv)
- Validación con sitios conocidos: [`resultados/Col2_2024/validacion_sitios_conocidos.csv`](./resultados/Col2_2024/validacion_sitios_conocidos.csv)
- Importancia de variables: [`resultados/importancia_variables_jubaea.png`](./resultados/importancia_variables_jubaea.png), [`resultados/importancia_variables_jubaea.csv`](./resultados/importancia_variables_jubaea.csv)
- Diagnósticos de colinealidad: [`resultados/correlacion_pearson.png`](./resultados/correlacion_pearson.png), [`resultados/vif_barplot.png`](./resultados/vif_barplot.png)
- Proyección CMIP6 2050: [`resultados/cambio_idoneidad_2050.csv`](./resultados/cambio_idoneidad_2050.csv), [`resultados/refugios_2050_uso_actual.csv`](./resultados/refugios_2050_uso_actual.csv)

## 7. Conclusión
La integración de ensemble modeling con MapBiomas Chile Colección 2 permite aislar, por primera vez de forma cuantitativa, la disponibilidad real de hábitat climático para *Jubaea chilensis*. Los resultados indican que el 77% del hábitat climáticamente apto permanece disponible como vegetación natural o restaurable, mientras que el 22,9% ya fue convertido a usos incompatibles, mayormente agrícolas. El análisis de sensibilidad entre colecciones (sección 6.4) demuestra además que gran parte de la aparente diferencia frente a una primera exploración con Colección 1 no es conversión real sino una mejora en la resolución de clasificación de MapBiomas — un hallazgo metodológico relevante para cualquier estudio que combine series temporales de distintas colecciones de MapBiomas.

El cruce adicional con el WDPA (sección 6.5) agrega el dato más accionable del estudio: de ese hábitat disponible, **solo el 12,6% está dentro de un área protegida** — el 87,4% restante no tiene ninguna figura de protección legal. El cruce con el Producto Fuego (sección 6.6) matiza además qué significa "protegido": las áreas protegidas se quemaron proporcionalmente más que las desprotegidas (12,4% vs 7,7%), y un palmar específico —**Palmar El Salto**— tuvo el 80% de su hábitat compatible quemado en 2013-2025, mientras otros palmares (Fray Jorge, Monte Aranda) no registraron quemas.

Finalmente, la proyección a 2050 (sección 6.8) muestra que ninguna de estas medidas puede esperar: bajo un escenario de emisiones moderado (CMIP6, SSP2-4.5), el hábitat de alta idoneidad climática se reduciría **46,2%** para 2041-2060, y el 22,3% de los refugios climáticos que persistirían ya está convertido hoy. La ventana para actuar sobre el hábitat disponible actual se cierra en paralelo al propio cambio climático.

La conservación de esta especie crítica requiere, por lo tanto, una estrategia en cinco frentes: (1) ampliar la protección formal sobre el hábitat compatible hoy desprotegido, priorizando **Cocalán y Petorca** según el ranking de la sección 6.7; (2) frenar la conversión adicional a uso agrícola en las zonas de transición; (3) priorizar manejo activo de combustible y cortafuegos en palmares de alto riesgo de incendio como Palmar El Salto; (4) intervenir directamente sobre la dinámica poblacional — protección de individuos adultos, control de depredadores exóticos de semillas y plántulas, y facilitación activa de la regeneración en los seis palmares relictos identificados; y (5) incorporar la contracción climática proyectada a 2050 en la planificación de nuevas áreas protegidas, priorizando zonas que retengan alta idoneidad bajo cambio climático y no solo en el presente.

---
> **Repositorio de Código Abierto:** El código fuente (R), la lógica de reclasificación de la leyenda MapBiomas Chile, el script de extracción vía Google Earth Engine y el pipeline de modelado están disponibles públicamente en este mismo repositorio: [`BLOQUE7_mapbiomas_chile.R`](./BLOQUE7_mapbiomas_chile.R), [`BLOQUE7b_snaspe_chile.R`](./BLOQUE7b_snaspe_chile.R), [`BLOQUE7c_fuego_chile.R`](./BLOQUE7c_fuego_chile.R), [`BLOQUE8_ranking_palmares.R`](./BLOQUE8_ranking_palmares.R), [`BLOQUE9_importancia_variables.R`](./BLOQUE9_importancia_variables.R), [`BLOQUE10_validacion_sitios_conocidos.R`](./BLOQUE10_validacion_sitios_conocidos.R), [`BLOQUE11_proyeccion_cmip6.R`](./BLOQUE11_proyeccion_cmip6.R), [`BLOQUE11b_refugios_2050.R`](./BLOQUE11b_refugios_2050.R).

## Referencias
- The iconic *Jubaea chilensis* teeters on the edge of local extinction: a plea for enhanced conservation policies. *Biodiversity and Conservation* (2024). https://link.springer.com/article/10.1007/s10531-024-02929-3
- Genetic Diversity and Population Structure of *Jubaea chilensis*. PMC (2022). https://pmc.ncbi.nlm.nih.gov/articles/PMC9370131/
- Multiple Anthropogenic Pressures Lead to Seed Dispersal Collapse of *Jubaea chilensis*. *Frontiers in Ecology and Evolution* (2021). https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2021.719566/full
- Distribución, tamaño y estructura poblacional de *Jubaea chilensis* en "Las Palmas", Petorca. *Bosque* (SciELO). https://www.scielo.cl/pdf/bosque/v37n3/art07.pdf
- Ficha de clasificación de especies, Ministerio del Medio Ambiente de Chile. https://clasificacionespecies.mma.gob.cl/wp-content/uploads/2019/10/Jubaea_chilensis_14RCE_FINAL.pdf
- MapBiomas Chile, Colección 1 y 2 — Códigos de leyenda. https://chile.mapbiomas.org/codigos-de-la-leyenda/
- Premio MapBiomas Chile — 1.ª Edición, Bases del concurso (2026).
