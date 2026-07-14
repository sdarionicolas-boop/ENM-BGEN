# Reporte Técnico: Análisis de Idoneidad Climática por Partido y Ranking de Priorización Territorial para la Colecta de Germoplasma

**Proyecto:** Beca de Iniciación a la Investigación (BIEI 2025) – UNAJ  
**Programa:** Banco de Germoplasma de Especies Nativas (BGEN-Noreste)  
**Autor:** Darío Nicolás Sánchez  
**Directora:** Ing. Andrea Quinteros  
**Fecha:** 13 de julio de 2026  

---

## 1. Introducción y Justificación
El planeamiento estratégico de las campañas de recolección de semillas para conservación *ex situ* requiere optimizar recursos logísticos y temporales en el campo. Los Modelos de Nicho Ecológico (ENM) estiman la distribución potencial teórica de las especies, pero sus resultados continuos de probabilidad a escala de píxel (~1 km²) deben integrarse con límites político-administrativos reales para ser de utilidad operativa para los colectores. 

Este documento presenta la metodología y los resultados del **Bloque 8** del pipeline de modelado, enfocado en cruzar los resultados de idoneidad climática de **19 especies nativas de interés para el BGEN** con la división por partidos de la Provincia de Buenos Aires (Argentina). El objetivo es jerarquizar y priorizar territorialmente las áreas de colecta, facilitando la toma de decisiones basada en datos científicos.

---

## 2. Metodología de Procesamiento Espacial

El análisis se codificó de manera automatizada en R mediante el script `BLOQUE8_partidos_bonaerenses.R` utilizando el siguiente flujo de trabajo:

```mermaid
graph TD
    A[Rasters de Idoneidad de 19 Especies] --> D[Cruce Espacial en R]
    B[Límites Políticos GADM Nivel 2 - ARG] --> E[Filtrado Prov. Buenos Aires]
    E --> D
    D --> F[Estadísticas de Idoneidad por Partido]
    D --> G[Cálculo de Área de Alta Idoneidad con cellSize]
    F --> H[Cálculo del Índice de Prioridad BGEN]
    G --> H
    H --> I[Ranking Consolidado CSV]
    H --> J[Capa Vectorial GeoJSON para QGIS]
    H --> K[Visor Leaflet Interactivo HTML]
end
```

### 2.1. Fuentes de Datos

1.  **Modelos de Idoneidad Climática:** 19 rasters de probabilidad continua (rango $0.0$ a $1.0$ o enteros de $0$ a $1000$ escalados por `biomod2`) correspondientes a las especies del BGEN (ej. *Jacaranda mimosifolia*, *Erythrina crista-galli*, *Vachellia caven*, etc.).
2.  **Límites Políticos:** Capas vectoriales obtenidas dinámicamente desde la base de datos global GADM (nivel 2 de divisiones administrativas de Argentina), filtrando los 128 partidos correspondientes a la Provincia de Buenos Aires.

### 2.2. Cálculo de Métricas Espaciales

Para cada especie y cada partido bonaerense se calcularon tres métricas clave:

*   **Idoneidad Climática Promedio ($\mu_{id}$):** Media de los valores de idoneidad de todos los píxeles contenidos dentro del polígono del partido. Representa el potencial general de establecimiento de la especie en el territorio.
*   **Idoneidad Climática Máxima ($Max_{id}$):** El valor del píxel con máxima idoneidad dentro del partido, útil para identificar refugios locales óptimos.
*   **Superficie con Alta Idoneidad ($km^2$):** Se definió como la superficie acumulada donde la idoneidad climática es $\ge 0.6$ (o $600/1000$). La superficie de las celdas se calculó mediante `cellSize(unit="km")` sobre el elipsoide de referencia WGS84 para evitar distorsiones cartográficas en altas latitudes.

### 2.3. Índice de Prioridad de Colecta BGEN ($I_{pc}$)
Para consolidar la información de las 19 especies nativas en una sola métrica de decisión territorial, se diseñó el **Índice de Prioridad de Colecta BGEN ($I_{pc}$)**. Este índice combina la idoneidad promedio de las especies y la riqueza de especies aptas dentro del partido mediante una fórmula de suma ponderada:

$$I_{pc} = (\bar{\mu}_{id} \times 0.7) + \left(\frac{R_{aptas}}{N_{total}} \times 0.3\right)$$

Donde:

*   $\bar{\mu}_{id}$: Promedio de las idoneidades medias de todas las especies presentes en el partido.
*   $R_{aptas}$ (Riqueza de Especies Aptas): Cantidad de especies que superan un umbral de idoneidad media de $0.2$ en el partido (lo que indica que el territorio es climáticamente propicio para su supervivencia).
*   $N_{total}$: Número total de especies analizadas ($19$).

---

## 3. Resultados y Ranking de Prioridad

La ejecución del script sobre los datos consolidados generó un ranking territorial completo. A continuación se detalla el **Top 15 de Partidos Bonaerenses Prioritarios para Colecta**:

| Rango | Partido | Idoneidad Media Global ($\bar{\mu}_{id}$) | Riqueza Especies Aptas | Área Alta Idoneidad Total ($km^2$) | Índice de Prioridad ($I_{pc}$) |
|:---:|:---|:---:|:---:|:---:|:---:|
| **1** | Quilmes | 0.8884 | 19 / 19 | 120.45 | **0.9219** |
| **2** | Vicente López | 0.8819 | 19 / 19 | 33.20 | **0.9173** |
| **3** | San Isidro | 0.8819 | 19 / 19 | 51.52 | **0.9173** |
| **4** | General San Martín | 0.8803 | 19 / 19 | 55.40 | **0.9162** |
| **5** | San Fernando | 0.8719 | 19 / 19 | 895.34 | **0.9103** |
| **6** | Tres de Febrero | 0.8703 | 19 / 19 | 45.60 | **0.9092** |
| **7** | Lomas de Zamora | 0.8690 | 19 / 19 | 87.21 | **0.9083** |
| **8** | Lanús | 0.8681 | 19 / 19 | 48.33 | **0.9077** |
| **9** | **Avellaneda** | **0.8653** | **19 / 19** | **52.12** | **0.9057** |
| **10** | Berazategui | 0.8619 | 19 / 19 | 185.30 | **0.9033** |
| **11** | Morón | 0.8603 | 19 / 19 | 55.90 | **0.9022** |
| **12** | San Miguel | 0.8584 | 19 / 19 | 82.10 | **0.9009** |
| **13** | Hurlingham | 0.8540 | 19 / 19 | 35.20 | **0.8978** |
| **14** | Ituzaingó | 0.8519 | 19 / 19 | 38.45 | **0.8963** |
| **15** | José C. Paz | 0.8490 | 19 / 19 | 50.12 | **0.8943** |

### 3.1. Discusión Espacial de los Resultados
Los resultados indican que los partidos costeros y del eje del Área Metropolitana de Buenos Aires (AMBA) norte y sur (Quilmes, Vicente López, San Isidro, Avellaneda) lideran el índice de prioridad. Esto se debe a dos factores biológicos y ambientales principales:
1.  **Influencia Costera y Efecto Térmico:** La ribera del Río de la Plata y el delta del Paraná actúan como corredores biológicos de clima templado-húmedo, mitigando las heladas intensas del interior y permitiendo la supervivencia de especies de afinidad subtropical.
2.  **Solapamiento de Nicho:** Las 19 especies seleccionadas comparten requerimientos de suelos con buen drenaje o formaciones ribereñas (marginales) típicas de la ecorregión de Espinal y Delta e Islas del Paraná, lo que concentra su idoneidad climática en el noreste bonaerense.

---

## 4. Estructura de Entregables Tecnológicos

Los productos cartográficos y estadísticos se organizaron en la carpeta `outputs/ranking_partidos/` con la siguiente estructura de archivos:

```text
outputs/ranking_partidos/
├── Araujia_sericifera_ranking_partidos.csv    ← Reporte individual de especie
├── ... (otros 18 archivos CSV por especie)
├── consolidado_ranking_partidos.csv           ← Matriz general de priorización (Top 15 arriba)
├── ranking_partidos.geojson                   ← Capa espacial para SIG (QGIS / ArcGIS)
├── mapa_prioridades_partidos.html             ← Visor Web interactivo Leaflet
└── mapa_prioridades_partidos_files/           ← Recursos JS/CSS asociados al mapa Leaflet
```

### 4.1. Características del Visor Web Leaflet (`mapa_prioridades_partidos.html`)
El mapa interactivo se diseñó para facilitar el acceso rápido desde dispositivos móviles en el campo sin necesidad de software SIG de escritorio:
*   **Simbología Coroplética:** Relleno de polígonos bajo una escala de color continuo (*YlOrRd* de ColorBrewer), donde los tonos rojos más oscuros representan los partidos de máxima prioridad.
*   **Información Dinámica al Pasar el Cursor (Hover labels):** Muestra de forma inmediata el nombre del partido, su rango de prioridad, el valor exacto de su índice $I_{pc}$, y la cantidad de especies aptas sobre las 19 totales.
*   **Doble Capa Base:** Permite alternar entre mapa claro simplificado (CartoDB) y mapa de topografía/relieve (Esri WorldTopoMap).

---

## 5. Conclusiones y Recomendaciones para Campo

1.  **Focalizar Viajes de Colecta:** Se recomienda priorizar Quilmes, el delta de San Fernando y la franja costera de Avellaneda y Berazategui para las salidas iniciales de campo, ya que concentran condiciones óptimas simultáneas para las 19 especies nativas bajo estudio.
2.  **Carga en Navegadores SIG Móviles:** El archivo `ranking_partidos.geojson` puede ser exportado de forma directa a aplicaciones móviles de mapeo (como QField o GPX Viewer) para que los colectores naveguen el territorio en tiempo real orientados por las prioridades.
3.  **Monitoreo Temporal:** La metodología propuesta sienta las bases para actualizar el índice incorporando proyecciones climáticas futuras (CMIP6), permitiendo identificar qué partidos mantendrán su idoneidad a largo plazo y cuáles se degradarán debido al calentamiento global.
