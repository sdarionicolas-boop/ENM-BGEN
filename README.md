# ENM-BGEN: Pipeline de Modelado de Nicho Ecológico

**Banco de Germoplasma de Especies Nativas (BGEN) – UNAJ**  
Beca BIEI 2025 | Proyecto ARG/19/G24 (GEF/PNUD) | Red ARGENA

![R](https://img.shields.io/badge/R-4.2%2B-blue?logo=r&logoColor=white)
![biomod2](https://img.shields.io/badge/biomod2-ensemble%20modeling-success)

### ¿Qué hace este pipeline?

Estima la **distribución potencial** de 19 especies nativas de Argentina de interés para el BGEN mediante **ensemble modeling** (biomod2).  

Combina cinco algoritmos y genera:
- Mapas de idoneidad de hábitat (GeoTIFF continuo y clasificado)
- Visor web interactivo multiespecies (HTML con Leaflet)

**Flujo completo**: desde datos crudos de GBIF + iNaturalist hasta mapas listos para conservación.

### Especies modeladas

- Araujia sericifera
- Austroeupatorium inulifolium
- Ceiba speciosa
- Celtis tala
- Cortaderia selloana
- Duranta erecta
- Erythrina crista-galli
- Ipomoea alba
- Jacaranda mimosifolia
- Passiflora caerulea
- Phytolacca dioica
- Salpichroa origanifolia
- Schinus molle
- Senna corymbosa
- Solanum pseudocapsicum
- Syagrus romanzoffiana
- Tecoma stans
- Tipuana tipu
- Vachellia caven

### Estructura del repositorio

```bash
ENM-BGEN/
├── ENM_BGEN_pipeline.R          ← Script maestro
├── README.md
├── data/                        ← (en .gitignore)
│   ├── registros_unificados_geocod.csv
│   ├── presencias_limpias.csv
│   └── presencias_thin.csv
├── variables/                   ← (en .gitignore)
│   ├── bio_AOI.tif
│   ├── topo_AOI.tif
│   ├── env_stack.tif
│   └── env_stack_vif5.tif
└── outputs/                     ← (en .gitignore)
    ├── metricas_todas_especies.csv
    ├── mapa_interactivo_ENM.html
    └── {especie}/
        ├── {especie}_idoneidad_prob.tif
        └── {especie}_idoneidad_clase.tif

Nota: Las carpetas data/, variables/ y outputs/ están en .gitignore porque contienen archivos pesados (rasters de varios GB).
Metodología resumida
Datos de presencia

Fuentes: GBIF + iNaturalist
Normalización de nombres científicos
Thinning espacial: 1 registro por píxel por especie

Variables ambientales

WorldClim v2.1 (BIO01–BIO19)
SRTM: elevación + pendiente + aspecto
Selección por VIF ≤ 5 → 9 variables retenidas

Modelado (biomod2)

Pseudo-ausencias: SRE (2 sets × 3.000 puntos)
Algoritmos: GLM, GBM, RF, MAXNET, XGBOOST
Ensemble: EMwmeanByTSS (ponderado por TSS)

Desempeño (validación)




















MétricaPromedioRangoAUCroc0.960.93 – 0.98TSS0.780.70 – 0.91
Clasificación de idoneidad






























ClaseRangoColor1 – Insustentable0.0 – 0.2⬜ gris2 – Bajo0.2 – 0.4🟨 amarillo3 – Moderado0.4 – 0.6🟧 naranja4 – Alto0.6 – 1.0🟩 verde
Requisitos

R ≥ 4.2
Paquetes: dplyr, stringr, readr, terra, geodata, usdm, corrplot, biomod2, leaflet, leaflet.extras, htmlwidgets

Uso

Clonar el repositorioBash
git clone https://github.com/dnicolas_97/ENM-BGEN.git
cd ENM-BGEN

Editar BASE_DIR y MODELS_DIR en el BLOQUE 0 de ENM_BGEN_pipeline.R
Colocar tu CSV en data/registros_unificados_geocod.csv
Ejecutar en RStudio (orden obligatorio):text
BLOQUE 0 → BLOQUE 2 → BLOQUE 1 (thinning)
         → BLOQUE 3 → BLOQUE 4
         → BLOQUE 5 → BLOQUE 6

Licencia
Este proyecto se distribuye como código abierto bajo licencia MIT.
La intención es que cada usuario lo adapte a sus propias especies, regiones y necesidades de conservación.
No se incluyen archivos de datos pesados (están en .gitignore); todos los insumos se obtienen de fuentes abiertas (GBIF, iNaturalist, WorldClim) o se generan al ejecutar el pipeline.
¡Sentite libre de explorar, modificar y compartir!
Referencia
Sánchez Leguizamón, D.N. (2026). Idoneidad de hábitat potencial para 19 especies de interés para un banco de germoplasma en Argentina: un enfoque de ensemble modeling. Beca BIEI 2025 – BGEN/UNAJ. Proyecto ARG/19/G24 (GEF/PNUD).
Contacto
Darío Nicolás Sánchez Leguizamón
Becario BIEI 2025 – Banco de Germoplasma de Especies Nativas (BGEN) – UNAJ
