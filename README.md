# ENM-BGEN: Pipeline de Modelado de Nicho Ecológico

**Banco de Germoplasma de Especies Nativas (BGEN) – UNAJ**  
Beca BIEI 2025 | Proyecto ARG/19/G24 (GEF/PNUD) | Red ARGENA

---

## La historia de este repositorio

**Agosto de 2025.** Empieza una Beca de Iniciación a la Investigación (BIEI) de la Universidad Nacional Arturo Jauretche, dirigida por la Ing. Andrea Quinteros, en el Banco de Germoplasma de Especies Nativas (BGEN-Noreste; proyecto ARG/19/G24, GEF/PNUD; red ARGENA). El plan original no era este: era *"Especies nativas, a un click"*, una plataforma colaborativa para cargar datos fenológicos en campo. Esa primera meta se cumplió en la primera mitad de la beca.

**El giro.** Un banco de germoplasma necesita saber adónde ir a colectar semillas ahora, y una base fenológica propia necesita años de monitoreo. De común acuerdo con la dirección, la beca sumó datos históricos abiertos (GBIF, iNaturalist) y pasó al modelado de nicho ecológico: 19 especies nativas, un cruce con MapBiomas para saber cuánto hábitat sigue en pie y un índice para ordenar los 128 partidos bonaerenses según su prioridad de colecta. De esa pregunta práctica nació este pipeline (marzo de 2026). Se presentó en la UNAJ (Mes de la Ciencia), en un congreso internacional (UNASAM, mayo de 2026) y se postuló al Premio MapBiomas Argentina. El visor público está en [enm-bgen.lovable.app](https://enm-bgen.lovable.app).

**Agosto de 2026.** La beca terminó, con sus objetivos cumplidos. El trabajo no. Siguió por cuenta propia: el pipeline se aplicó a *Mauritia flexuosa* en Perú (en coautoría con Martín E. Cano, sumando una proyección climática y un módulo de potencial de cosecha para 382 comunidades nativas con visor web) y a *Jubaea chilensis* en Chile, y esa tercera aplicación obligó a preguntarse qué tan sólido era lo hecho. El informe de la beca daba por controlada la autocorrelación espacial con el thinning por píxel; las pruebas de Chile mostraron que no alcanzaba. Esas pruebas volvieron a Argentina, cambiaron el ensemble del pipeline y los resultados se recalcularon. Además, este pipeline es la base del Trabajo Final de Grado.

Lo que sigue es ese recorrido entero, con sus errores a la vista. El mapa del camino, las lecciones y el **estado de cada resultado** están en [`LABORATORIO.md`](./LABORATORIO.md).

---

## ¿Qué hace este pipeline?

Estima la distribución potencial de **19 especies nativas de Argentina** de interés para un banco de germoplasma, usando ensemble modeling con cinco algoritmos de aprendizaje automático. El resultado es un mapa de idoneidad de hábitat por especie, exportado como raster GeoTIFF y como visor web interactivo.

El pipeline cubre el flujo completo: desde registros de biodiversidad crudos hasta mapas listos para análisis de conservación.


---

## Especies modeladas

| | | |
|---|---|---|
| *Araujia sericifera* | *Ipomoea alba* | *Schinus molle* |
| *Austroeupatorium inulifolium* | *Jacaranda mimosifolia* | *Senna corymbosa* |
| *Ceiba speciosa* | *Passiflora caerulea* | *Solanum pseudocapsicum* |
| *Celtis tala* | *Phytolacca dioica* | *Syagrus romanzoffiana* |
| *Cortaderia selloana* | *Salpichroa origanifolia* | *Tecoma stans* |
| *Duranta erecta* | *Erythrina crista-galli* | *Tipuana tipu* |
| *Vachellia caven* | | |

---

## Estructura del repositorio
```
ENM-BGEN/
├── ENM_BGEN_pipeline.R             ← script maestro (BLOQUES 0-7: presencias → modelos → ensemble → mapas)
├── BLOQUE7_mapbiomas.R             ← cruce idoneidad × MapBiomas Argentina 2024
├── BLOQUE8_partidos_bonaerenses.R  ← ranking Ipc de los 128 partidos bonaerenses
├── Reporte_Tecnico_Ranking_Partidos.md / .pdf   ← reporte técnico del ranking Ipc
├── API_GBIF1.ipynb                 ← descarga de registros de presencia desde GBIF
├── convert_to_pdf.py               ← conversión del reporte técnico a PDF
├── .Renviron.example               ← variables de entorno con las rutas (ver "Configuración de rutas")
├── LABORATORIO.md                  ← recorrido del arco Argentina → Perú → Chile, lecciones y estado de cada resultado
├── ejemplo_minimo/                 ← una especie de punta a punta en ~1 min, con datos de muestra (ver "Probalo en un minuto")
├── README.md · LICENSE
│
├── extensiones/                    ← el mismo pipeline aplicado y puesto a prueba en otros contextos
│   ├── argentina/                  ← pruebas de robustez (BLOQUE_ROBUSTEZ_A-E) y corrección EMmean (BLOQUE_F-H)
│   ├── chile/                      ← Jubaea chilensis: MapBiomas Chile + WDPA + fuego + CMIP6 (BLOQUE7-15)
│   └── peru/                       ← Mauritia flexuosa: MapBiomas Perú
│
└── outputs/ranking_partidos/       ← único resultado versionado: ranking Ipc (CSV) con ensemble EMmean
```

Las carpetas de trabajo (no versionadas, se generan al correr el pipeline) son `data/` (presencias), `variables/` (capas ambientales) y `outputs/` (rasters y mapas), más los modelos de biomod2. Contienen archivos pesados (varios GB), por eso están en `.gitignore`: en el repo se versionan los scripts, los resultados tabulares y las pruebas de robustez.


---

## Metodología resumida

### Datos de presencia
- Fuentes: GBIF + iNaturalist
- Normalización de nombres científicos (binomio sin autoría)
- Limpieza de coordenadas
- Umbral mínimo: 50 registros por especie
- Thinning espacial: 1 registro por píxel por especie

### Variables ambientales
- WorldClim v2.1 (BIO01–BIO19, ~1 km, 1970–2000)
- SRTM: elevación + pendiente + aspecto (~1 km)
- Selección por VIF ≤ 5 → 9 variables retenidas: `bio02, bio03, bio08, bio09, bio13, bio14, bio15, slope, aspect`

### Modelado (biomod2)
- Pseudo-ausencias: estrategia SRE, 2 sets × 3.000 puntos
- Algoritmos: GLM, GBM, RF, MAXNET, XGBOOST
- Validación cruzada: 3 repeticiones, 80/20
- Ensemble: **EMmean** (media simple, no ponderada). *Se cambió desde `EMwmeanByTSS`: una prueba de robustez con thinning espacial + validación cruzada por bloques mostró que ponderar por TSS de validación aleatoria sobrepondera a Random Forest, el algoritmo que peor generaliza espacialmente (calibración perfecta = 1,0, firma clásica de sobreajuste). Detalle completo en `extensiones/argentina/BLOQUE_ROBUSTEZ_A-D` (diagnóstico sobre 4 especies) y `BLOQUE_F-H` (corrección aplicada a las 19 especies oficiales: mediana −0,71 % en área de hábitat, sin cambios en el Top 15 del ranking Ipc).*
- Umbrales de calidad: TSS ≥ 0.7 y AUCroc ≥ 0.9

### Desempeño (validación)

| Métrica | Promedio | Rango |
|---------|----------|-------|
| AUCroc | 0.96 | 0.93 – 0.98 |
| TSS | 0.78 | 0.70 – 0.91 |

### Clasificación de idoneidad

| Clase | Rango | Color |
|-------|-------|-------|
| 1 – Insustentable | 0.0 – 0.2 | ⬜ gris |
| 2 – Bajo | 0.2 – 0.4 | 🟨 amarillo |
| 3 – Moderado | 0.4 – 0.6 | 🟧 naranja |
| 4 – Alto | 0.6 – 1.0 | 🟩 verde |

---

## Requisitos
```
R >= 4.2
```

Paquetes (ver BLOQUE 0 del script para el comando de instalación):
```r
dplyr, stringr, readr      # manipulación de datos
terra, geodata             # datos espaciales
usdm, corrplot             # colinealidad
biomod2                    # modelado
leaflet, leaflet.extras, htmlwidgets  # visualización
```

---

## Uso

### 1. Clonar el repositorio
```bash
git clone https://github.com/sdarionicolas-boop/ENM-BGEN.git
cd ENM-BGEN
```

### 2. Configuración de rutas

Los scripts no tienen rutas fijas: leen las carpetas desde variables de entorno, así que no hace falta editar el código. Copiá el contenido de [`.Renviron.example`](./.Renviron.example) a tu `~/.Renviron` (en R: `usethis::edit_r_environ()`), reemplazá las rutas por las tuyas y reiniciá R:

| Variable | Qué es | La usan |
|---|---|---|
| `ENM_ARG_DIR` | carpeta del proyecto Argentina/BGEN (con `data/`, `variables/`, `outputs/`) | pipeline, BLOQUE7/8, `extensiones/argentina/` |
| `ENM_CHILE_DIR` | carpeta del proyecto Chile | `extensiones/chile/` |
| `ENM_PERU_DIR` | carpeta del proyecto Perú | `extensiones/peru/` |
| `ENM_MODELS_DIR` | carpeta donde biomod2 guarda los modelos (el *working directory* al modelar) | scripts que cargan modelos |
| `ENM_REPO_DIR` | raíz de tu clon de este repo | scripts que escriben dentro de `extensiones/` |
| `ENM_PANDOC_DIR` | (opcional) carpeta de pandoc, p. ej. la de Quarto | mapas HTML autocontenidos de Chile |

Si falta una variable, el script se detiene al inicio con un mensaje que dice cuál definir.

### 3. Colocar el archivo de presencias

Copiá tu CSV a `data/registros_unificados_geocod.csv`.  
Columnas requeridas: `especie`, `lat`, `lon`

### 4. Correr el pipeline

Ejecutá el script bloque por bloque en RStudio en este orden:
```
BLOQUE 0  →  BLOQUE 2  →  BLOQUE 1 (thinning)
          →  BLOQUE 3  →  BLOQUE 4
          →  BLOQUE 5  →  BLOQUE 6
```

---

## Probalo en un minuto

[`ejemplo_minimo/`](./ejemplo_minimo/) corre el núcleo del pipeline sobre una sola especie (*Celtis tala*) con datos de muestra incluidos y capas de WorldClim que se descargan solas. No requiere configurar rutas. Además de generar el mapa de idoneidad, compara la validación cruzada aleatoria con la de bloques espaciales, así que podés ver el efecto que documentan las pruebas de robustez (TSS medio por bloques ≈ 0,40, siempre por debajo del 0,45–0,65 de la validación aleatoria).

```bash
Rscript ejemplo_minimo/ejemplo_una_especie.R
```

---

## Adaptar el pipeline a otra región o a otras variables ambientales

El pipeline está escrito para WorldClim + SRTM, pero está pensado para poder sumar otras fuentes de variables ambientales (por ejemplo ENVIREM) o aplicarse a otra región. Las dos preguntas que más aparecen al adaptarlo:

### 1. "Al unir capas de distintas fuentes me quedan píxeles en NoData en los bordes (zona costera, por ejemplo)"

Esto pasa porque cada fuente de variables (WorldClim, ENVIREM, capas topográficas) suele tener su propia máscara de tierra/agua, y no siempre coinciden exactamente en el borde costero o en los límites del área de estudio. Al recortar cada stack por separado con su propia máscara y después unirlos (`c()` en terra), los píxeles donde una fuente tiene dato y la otra no quedan en NoData.

**Solución recomendada:**
1. Generá **una única máscara** de tierra para tu área de estudio (por ejemplo, a partir de una sola capa de referencia, como el DEM o una bio de WorldClim).
2. Recortá y alineá **todos** los stacks (WorldClim, ENVIREM, topografía, uso de suelo) a esa misma máscara y a la misma grilla (`resample()`) **antes** de unirlos con `c()`. Esto evita bordes desfasados entre capas de distinta fuente/resolución.
3. Si aun así quedan huecos puntuales (por ejemplo, un solo pixel aislado sin dato en alguna variable), podés rellenarlos con la media de la variable en el área de estudio (`focal()` con una ventana pequeña, o simplemente `values(r)[is.na(values(r))] <- mean(values(r), na.rm = TRUE)`).
4. **Si rellenás con la media, declaralo explícitamente en la metodología del trabajo final** — es una decisión metodológica, no un detalle técnico invisible: cuántos píxeles se rellenaron y con qué criterio.

Este es el mismo enfoque que usa `BLOQUE7_mapbiomas.R` para alinear MapBiomas (remuestreado a la resolución del stack de idoneidad) antes de cruzarlo — la lógica es la misma para cualquier par de capas de distinta fuente.

### 2. "Quiero correr el pipeline para varias especies, no solo una"

El pipeline ya está armado para eso: el loop de los BLOQUES 3 y 4 (entrenamiento y ensemble) recorre automáticamente todas las especies presentes en `data/presencias_thin.csv` (columna `species`) — no hace falta tocar el código para agregar especies, solo agregar sus registros de presencia al CSV de entrada con el mismo formato (`especie`, `lat`, `lon`). Cada especie se entrena y proyecta de forma independiente, y los resultados se guardan en `outputs/{especie}/`.

## BLOQUE 7: Cruce con MapBiomas Argentina 2024

Cruza los rasters de idoneidad de hábitat con la Colección 2 de MapBiomas Argentina (2024, 30m) para cuantificar cuánto hábitat climáticamente apto está disponible, bajo uso restaurable, o ya fue convertido.

### Clasificación de uso del suelo

| Categoría | Clases MapBiomas | Códigos |
|-----------|-----------------|---------|
| Compatible | Bosques, pastizales, arbustales, turberas | 3,4,6,11,12,63,66,73,77 |
| Restaurable | Silvicultura, pasturas, mosaico de usos | 9,15,21 |
| Incompatible | Cultivos, urbano, áreas degradadas | 19,24,25,36 |

### Resultados para 19 especies nativas

| Categoría | Superficie acumulada |
|-----------|---------------------|
| Hábitat disponible (conservación) | 171 millones de ha |
| Hábitat perdido (conversión) | 119 millones de ha |
| Uso restaurable (restauración) | 31 millones de ha |

### Archivos adicionales requeridos
- `argentina_coverage_2024.tif` – MapBiomas Argentina Colección 2, año 2024
  Descarga: [plataforma.argentina.mapbiomas.org](https://plataforma.argentina.mapbiomas.org)

> **¿Por qué BLOQUE 2 antes que BLOQUE 1 (thinning)?**  
> El thinning necesita el stack ambiental para asignar píxeles.  
> El stack se genera en BLOQUE 2.

---

## 🌎 Extensiones del proyecto

Este pipeline fue diseñado como un marco de trabajo reproducible y transfronterizo, no como un script aislado para una sola especie o país. La misma metodología (modelado de nicho ecológico + cruce algebraico con capas oficiales de uso de suelo MapBiomas) fue adaptada y validada en postulaciones adicionales:

- 🇵🇪 **Perú — *Mauritia flexuosa* (Aguaje):** cruce de idoneidad climática con **MapBiomas Perú 2024 (Colección 3)**, postulado al Premio MapBiomas Perú 2026. Ver [`extensiones/peru/`](./extensiones/peru/).
- 🇨🇱 **Chile — *Jubaea chilensis* (Palma chilena):** integra dos productos de MapBiomas Chile (cobertura Colección 2 y Producto Fuego) con áreas protegidas WDPA/SNASPE (fuente externa) y una proyección climática a 2050 con ensamble de 6 GCMs, postulado al Premio MapBiomas Chile 2026 (1.ª Edición). El 77% del hábitat de alta idoneidad climática persiste hoy, solo el 12,6% está protegido, un palmar específico ("Palmar El Salto") perdió el 80% de su hábitat al fuego, y la proyección a 2050 (CMIP6, SSP2-4.5) muestra alta incertidumbre inter-modelo (mediana -0,2%, rango -46,2% a +17,0%) — esa incertidumbre es el argumento para proteger ahora lo que ya es cierto. Ranking de los 6 palmares: Cocalán y Petorca, máxima prioridad. Ver [`extensiones/chile/`](./extensiones/chile/).

La arquitectura modular permite escalar el pipeline a nuevos países reutilizando el mismo núcleo de modelado (`ENM_BGEN_pipeline.R`) y adaptando solo la especie, el área de interés y la leyenda de reclasificación MapBiomas de cada país.

---

## Cómo citar

Si usás este pipeline en tu trabajo, podés citar el software:

> Sánchez Leguizamón, Darío Nicolás (2026). *ENM-BGEN: Pipeline de Modelado de Nicho Ecológico para Especies Nativas de Argentina, Perú y Chile* [Software]. Zenodo. https://doi.org/10.5281/zenodo.23149296

Ese DOI representa **todas las versiones** y siempre resuelve a la más reciente. Para citar una versión puntual (por ejemplo v1.1.0, `10.5281/zenodo.23598500`), usá el DOI que figura en la página de ese release en Zenodo. GitHub también arma la cita a partir de [`CITATION.cff`](./CITATION.cff) (botón "Cite this repository").

Trabajo de origen del pipeline:

> Sánchez Leguizamón, Darío Nicolás (2026). *Idoneidad de hábitat potencial para 19 especies de interés para un banco de germoplasma en Argentina: un enfoque de ensemble modeling.* Beca BIEI 2025 – BGEN/UNAJ. Proyecto ARG/19/G24 (GEF/PNUD).

---

## Contacto

Darío Nicolás Sánchez Leguizamón  
Banco de Germoplasma de Especies Nativas – UNAJ  

---

## Licencia y uso libre

Este proyecto se distribuye como código abierto bajo licencia MIT. La intención es que cualquier investigador, técnico o estudiante pueda adaptarlo a sus propias especies, regiones o bancos de germoplasma. No se incluyen archivos de datos pesados — todos los insumos se obtienen de fuentes abiertas (GBIF, iNaturalist, WorldClim, SRTM) o se generan al ejecutar el pipeline. ¡Sentite libre de explorar, modificar y compartir!