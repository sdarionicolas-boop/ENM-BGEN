# Memoria Técnica - Premio MapBiomas Chile 2026 (1.ª Edición)

## 1. Título del Proyecto
**Idoneidad climática y disponibilidad real de hábitat para *Jubaea chilensis* (Palma chilena) mediante la integración de modelado de nicho ecológico y MapBiomas Chile Colección 1.**

## 2. Autor(es)
Darío Nicolás Sánchez Leguizamón (Beca BIEI 2025 | Proyecto ARG/19/G24 | Red ARGENA)

## 3. Resumen Ejecutivo
*Jubaea chilensis* (palma chilena) es una palmera endémica de Chile central catalogada como Vulnerable por la UICN, con estudios recientes (2024) que argumentan su reclasificación a En Peligro Crítico. De la población preколombina estimada, hoy solo subsisten ~121.284 individuos (2,5%), concentrados en seis palmares relictos entre La Serena y el Maule.

Este trabajo aplica el mismo pipeline de ensemble modeling (biomod2) ya utilizado para 19 especies del BGEN-UNAJ en Argentina y para *Mauritia flexuosa* en Perú, ahora sobre *Jubaea chilensis*, cruzando la idoneidad climática resultante con **MapBiomas Chile Colección 1 (2022)** para cuantificar cuánto del hábitat climáticamente apto persiste como vegetación natural, cuánto es restaurable y cuánto fue efectivamente convertido por uso antrópico.

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

### 5.4. Integración con MapBiomas Chile Colección 1 (2022)
La capa de cobertura y uso del suelo (`chile_coverage_2022.tif`, descarga directa desde `storage.googleapis.com/mapbiomas-public`) se reclasificó según la leyenda oficial de MapBiomas Chile en tres dimensiones operativas, siguiendo el mismo criterio ya aplicado a Argentina y Perú:

| Categoría | Clases MapBiomas Chile | Códigos |
|-----------|------------------------|---------|
| **Compatible** | Bosque, humedal, pastizal natural, estepa, matorral | 3, 59, 60, 67, 11, 12, 63, 66 |
| **Restaurable** | Silvicultura, pastura, mosaico agropecuario | 9, 15, 21 |
| **Incompatible** | Agricultura, infraestructura, área degradada | 18, 24, 25 |
| *Excluido* | Agua, hielo/nieve, afloramiento rocoso, arena/duna, salar, sin datos | 0, 23, 29, 33, 34, 61 |

La categoría *Excluido* agrupa superficies naturalmente no vegetadas (alta cordillera, cuerpos de agua) que quedan fuera del esquema de conversión antrópica y no se computan como "pérdida" — una precisión metodológica necesaria en un país con una fracción tan grande de su territorio en la criósfera y zonas áridas, a diferencia de Argentina o Perú.

Luego se realizó el cruce algebraico con la capa de idoneidad climática reclasificada (Alta/Moderada/Baja), obteniendo 9 combinaciones posibles por celda.

## 6. Resultados

### 6.1. Superficie de hábitat de Alta Idoneidad Climática

| Categoría | Superficie |
|-----------|-----------:|
| Compatible (conservación) | 12.931 km² |
| Restaurable (restauración) | 7.135 km² |
| Incompatible (hábitat perdido) | 892 km² |
| **Total idoneidad alta** | **20.958 km²** |

**El 92% del hábitat climáticamente óptimo para *Jubaea chilensis* persiste como vegetación natural o es restaurable; solo el 4,3% fue efectivamente convertido.**

### 6.2. Interpretación
Este resultado es metodológicamente relevante porque **contradice la hipótesis de pérdida de hábitat como causa principal del colapso poblacional**. El hábitat climáticamente apto está mayormente disponible — la brecha entre 121.284 individuos actuales y el hábitat potencial disponible (miles de km²) señala que el factor limitante no es la disponibilidad de tierra, sino la presión directa documentada en la literatura: extracción de savia, sobrecosecha de semillas y depredación por fauna exótica que colapsa la dispersión y el reclutamiento de nuevos individuos.

Esto reposiciona la prioridad de conservación: no se trata (únicamente) de proteger hábitat de la expansión agrícola, sino de proteger los individuos reproductivos y facilitar la regeneración natural en un hábitat que climática y territorialmente sigue siendo apto.

### 6.3. Visualización de resultados
- Visor interactivo de idoneidad climática: [`resultados/mapa_interactivo_jubaea_chile.html`](./resultados/mapa_interactivo_jubaea_chile.html)
- Visor interactivo del cruce con MapBiomas: [`resultados/mapa_interactivo_cruce_chile.html`](./resultados/mapa_interactivo_cruce_chile.html)
- Diagnósticos de colinealidad: [`resultados/correlacion_pearson.png`](./resultados/correlacion_pearson.png), [`resultados/vif_barplot.png`](./resultados/vif_barplot.png)

## 7. Conclusión
La integración de ensemble modeling con MapBiomas Chile Colección 1 permite aislar, por primera vez de forma cuantitativa, la disponibilidad real de hábitat climático para *Jubaea chilensis*. Los resultados indican que la conservación de esta especie crítica depende menos de la protección de hábitat frente a la conversión de uso de suelo — que ya está mayormente disponible — y más de intervenciones directas sobre la dinámica poblacional: protección de individuos adultos, control de depredadores exóticos de semillas y plántulas, y facilitación activa de la regeneración en los palmares relictos de Ocoa, Cocalán, Viña del Mar/Valparaíso, Candelaria, Petorca y Culimo.

---
> **Repositorio de Código Abierto:** El código fuente (R), la lógica de reclasificación de la leyenda MapBiomas Chile y el pipeline de modelado están disponibles públicamente en este mismo repositorio: [`BLOQUE7_mapbiomas_chile.R`](./BLOQUE7_mapbiomas_chile.R).

## Referencias
- The iconic *Jubaea chilensis* teeters on the edge of local extinction: a plea for enhanced conservation policies. *Biodiversity and Conservation* (2024). https://link.springer.com/article/10.1007/s10531-024-02929-3
- Genetic Diversity and Population Structure of *Jubaea chilensis*. PMC (2022). https://pmc.ncbi.nlm.nih.gov/articles/PMC9370131/
- Multiple Anthropogenic Pressures Lead to Seed Dispersal Collapse of *Jubaea chilensis*. *Frontiers in Ecology and Evolution* (2021). https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2021.719566/full
- Distribución, tamaño y estructura poblacional de *Jubaea chilensis* en "Las Palmas", Petorca. *Bosque* (SciELO). https://www.scielo.cl/pdf/bosque/v37n3/art07.pdf
- Ficha de clasificación de especies, Ministerio del Medio Ambiente de Chile. https://clasificacionespecies.mma.gob.cl/wp-content/uploads/2019/10/Jubaea_chilensis_14RCE_FINAL.pdf
- MapBiomas Chile, Colección 1 y 2 — Códigos de leyenda. https://chile.mapbiomas.org/codigos-de-la-leyenda/
