# Ejemplo mínimo reproducible

Una especie (*Celtis tala*), de punta a punta, en menos de un minuto: variables ambientales → selección por colinealidad (VIF) → modelos biomod2 → ensemble `EMmean` → mapa de idoneidad. Sirve para entender el método, verificar que tu instalación funciona y ver con tus propios ojos el hallazgo central de las pruebas de robustez del repo.

```bash
Rscript ejemplo_una_especie.R
```

(o abrí `ejemplo_una_especie.R` en RStudio, *Session > Set Working Directory > To Source File Location*, y corrélo completo). No hace falta configurar variables de entorno: este ejemplo es autocontenido.

## Qué necesitás

- R ≥ 4.2 y los paquetes `terra`, `geodata`, `usdm`, `biomod2`, `ENMeval`, `dplyr`, `readr`, `tidyr`.
- Conexión a internet la primera vez: las capas de WorldClim 2.1 se descargan a `cache/` (no se redistribuyen en este repositorio).
- Los 61 registros de presencia de *Celtis tala*, ya rarificados a 10 km, vienen en `datos/`.

## Qué hace

1. Define el área de estudio (±1° alrededor de las presencias) y descarga/recorta WorldClim a 10 arcmin.
2. Descarta variables colineales con VIF ≤ 5 (igual que el pipeline).
3. Entrena los 5 algoritmos del pipeline (GLM, GBM, RF, MAXNET, XGBoost) **dos veces**: con validación cruzada aleatoria 80/20 —la del pipeline oficial— y con validación por bloques espaciales.
4. Imprime y guarda la comparación del TSS de validación por algoritmo.
5. Construye el ensemble `EMmean` con los modelos de bloques y exporta el mapa.

## Qué deberías ver

La validación aleatoria **cambia bastante de una corrida a otra** (reparte los 61 registros al azar y biomod2 no respeta la semilla en ese paso); la de bloques es determinista. En 5 corridas seguidas el TSS medio fue:

| | Validación aleatoria | Validación por bloques |
|---|---:|---:|
| TSS medio de validación (5 algoritmos) | 0,45 – 0,65 (0,63 · 0,65 · 0,45 · 0,57 · 0,60) | 0,403 en las 5 |

La validación por bloques fue menor que la aleatoria en las 5 corridas, pero el tamaño de la brecha depende del azar de cada partición: con tan pocos registros, la propia estimación aleatoria es ruidosa. El AUCroc cae menos que el TSS (por bloques ≈ 0,78–0,82). Random Forest muestra calibración 1,000 en ambos esquemas: la firma de sobreajuste que motivó reemplazar el ensemble ponderado por TSS (`EMwmeanByTSS`) por `EMmean` en el pipeline.

Salidas, en `salida/` (ignorada por git): `metricas_random_vs_block.csv`, `celtis_tala_idoneidad_prob.tif`, `celtis_tala_idoneidad_clase.tif` y `celtis_tala_idoneidad.png`.

## Diferencias con el pipeline completo

| | Este ejemplo | Pipeline completo |
|---|---|---|
| Especies | 1 | 19 |
| Resolución de las capas | 10 arcmin (~18 km) | 30 arcsec (~1 km) |
| Topografía | solo elevación | elevación + pendiente + aspecto |
| Pseudo-ausencias | 1 set × 500 | 2 sets × 3.000 |
| Filtro de calidad del ensemble | ninguno | TSS ≥ 0,7 y AUCroc ≥ 0,9 |

El filtro de calidad se omite a propósito: bajo validación espacial ningún algoritmo llega a TSS ≥ 0,7 y el ensemble quedaría vacío. Es el mismo problema que documenta `extensiones/argentina/`.

## Límites

Es un ejemplo con una sola especie y capas gruesas: sirve para verificar el método y ver el efecto de la validación, **no** para sacar conclusiones de conservación sobre *Celtis tala*. Para eso, usá el pipeline completo (`ENM_BGEN_pipeline.R`).
