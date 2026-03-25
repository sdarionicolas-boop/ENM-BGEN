# ENM-BGEN: Pipeline de Modelado de Nicho Ecológico

**Banco de Germoplasma de Especies Nativas (BGEN) – UNAJ**  
Beca BIEI 2025 | Proyecto ARG/19/G24 (GEF/PNUD) | Red ARGENA

---

## ¿Qué hace este pipeline?

Estima la distribución potencial de **19 especies nativas de Argentina** de interés para un banco de germoplasma, usando ensemble modeling con cinco algoritmos de aprendizaje automático. El resultado es un mapa de idoneidad de hábitat por especie, exportado como raster GeoTIFF y como visor web interactivo.

El pipeline cubre el flujo completo: desde registros de biodiversidad crudos hasta mapas listos para análisis de conservación.

---

## Especies modeladas

| Especie | Especie |
|---------|---------|
| *Araujia sericifera* | *Passiflora caerulea* |
| *Austroeupatorium inulifolium* | *Phytolacca dioica* |
| *Ceiba speciosa* | *Salpichroa origanifolia* |
| *Celtis tala* | *Schinus molle* |
| *Cortaderia selloana* | *Senna corymbosa* |
| *Duranta erecta* | *Solanum pseudocapsicum* |
| *Erythrina crista-galli* | *Syagrus romanzoffiana* |
| *Ipomoea alba* | *Tecoma stans* |
| *Jacaranda mimosifolia* | *Tipuana tipu* |
| *Vachellia caven* | |

---

## Estructura del repositorio

```
ENM-BGEN/
├── ENM_BGEN_pipeline.R     ← script maestro (todo el pipeline)
├── README.md
│
├── data/                   ← datos de presencia (no incluidos en el repo)
│   ├── registros_unificados_geocod.csv   ← archivo de entrada
│   ├── presencias_limpias.csv            ← generado por BLOQUE 1
│   └── presencias_thin.csv               ← generado por BLOQUE 1
│
├── variables/              ← variables ambientales (generadas por BLOQUE 2)
│   ├── bio_AOI.tif
│   ├── topo_AOI.tif
│   ├── env_stack.tif
│   └── env_stack_vif5.tif
│
└── outputs/                ← resultados (generados por BLOQUES 4-6)
    ├── metricas_todas_especies.csv
    ├── mapa_interactivo_ENM.html
    └── {especie}/
        ├── {especie}_idoneidad_prob.tif
        └── {especie}_idoneidad_clase.tif
```

> **Nota:** Las carpetas `data/`, `variables/` y `outputs/` están en `.gitignore` porque contienen archivos pesados (rasters de varios GB). Solo se versiona el script.

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
- Ensemble: EMwmeanByTSS (ponderado por TSS de validación)
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
git clone https://github.com/tu-usuario/ENM-BGEN.git
cd ENM-BGEN
```

### 2. Configurar rutas

Editá las primeras líneas del BLOQUE 0 en `ENM_BGEN_pipeline.R`:

```r
BASE_DIR   <- "ruta/a/tu/carpeta/proyecto"
MODELS_DIR <- "ruta/donde/biomd2/guardara/modelos"
```

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

> **¿Por qué BLOQUE 2 antes que BLOQUE 1 (thinning)?**  
> El thinning necesita el stack ambiental para asignar píxeles.  
> El stack se genera en BLOQUE 2.

---

## Referencia

Si usás este pipeline en tu trabajo, podés citar:

> Sánchez Leguizamón, D.N. (2026). *Idoneidad de hábitat potencial para 19 especies de interés para un banco de germoplasma en Argentina: un enfoque de ensemble modeling.* Beca BIEI 2025 – BGEN/UNAJ. Proyecto ARG/19/G24 (GEF/PNUD).

---

## Contacto

Darío Nicolás Sánchez Leguizamón  
Banco de Germoplasma de Especies Nativas – UNAJ  
Proyecto ARG/19/G24 – GEF/PNUD – Red ARGENA
