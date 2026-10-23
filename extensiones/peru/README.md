# Extensión Perú — *Mauritia flexuosa* (aguaje)

**Distribución potencial y estado de conservación del hábitat de *Mauritia flexuosa* en la Amazonía peruana: análisis integrado de nicho ecológico y datos MapBiomas Perú Colección 3**

Sánchez Leguizamón, D. N. y Cano, M. E. — postulado a la 3.ª edición del Premio MapBiomas Perú 2026.

Es la primera expansión del pipeline fuera de Argentina, y la más amplia en alcance: además del cruce con MapBiomas, suma una proyección climática a 2050 y un módulo socioeconómico que estima el potencial de cosecha sostenible de aguaje de comunidades nativas, con un visor web. El enlace que figura en el trabajo postulado (`plan-aguaje.netlify.app`) respondía 404 al 23 de octubre de 2026.

## Qué se hizo

| Componente | Detalle |
|---|---|
| Presencias | 22.551 registros de GBIF; tras limpieza y thinning a un registro por píxel (~1 km), 523 registros independientes. |
| Variables | WorldClim 2.1 + SRTM; 9 retenidas tras correlación y VIF ≤ 5: BIO02, BIO03, BIO08, BIO09, BIO13, BIO14, BIO15, pendiente y aspecto. |
| Modelado | biomod2, 5 algoritmos (GLM, GBM, RF, MAXNET, XGBoost), 2 sets de 3.000 pseudo-ausencias (SRE), validación cruzada 80/20 con 3 repeticiones; ensemble ponderado por TSS. |
| Idoneidad actual | ~196.400 km² de alta idoneidad (> 75 %) y ~120.300 km² de idoneidad moderada. |
| Cruce con MapBiomas Perú, Colección 3 (2024) | De la alta idoneidad, **94,7 % (178.949 km²) conserva vegetación natural**, 4,9 % (9.286 km²) está bajo uso restaurable y 0,4 % (811 km²) ya fue convertido. |
| Proyección climática | SSP2-4.5, MIROC6, 2041-2060: la alta idoneidad se contrae **81,5 %** (de ~196.400 a ~36.200 km²). |
| Comunidades nativas | Sobre 2.002 comunidades nativas (el total de registros del shapefile oficial de MIDAGRI/GEORURAL, verificado), **382 tienen potencial de cosecha sostenible positivo (~4.090 t/año)**; las de mayor potencial son Matsés (172 t/año), Puerto Elvira (135) y Santa María de Fátima (98). Usa zonas inundables, idoneidad y accesibilidad por la red vial nacional (5.365 segmentos). |
| Producto de transferencia | Visor web interactivo con idoneidad actual y futura, 401 comunidades con *popup*, red vial y panel de conclusiones, pensado para SERNANP, MINAM y organizaciones de conservación (hoy no disponible en línea; ver arriba). |

## Qué hay en este repositorio

- [`BLOQUE7_mapbiomas_peru.R`](./BLOQUE7_mapbiomas_peru.R): cruce de idoneidad × MapBiomas Perú 2024 (Compatible / Restaurable / Incompatible).
- [`resultados/`](./resultados/): visor HTML del cruce, resúmenes de superficie por categoría y capturas de la cartografía.
- [`Memoria_Tecnica_Premio_MapBiomas.md`](./Memoria_Tecnica_Premio_MapBiomas.md): una **versión preliminar** de una sola autoría. La versión postulada (v4b, con coautoría y los resultados de arriba) no está en el repositorio.

Los modelos y rasters de idoneidad (actual y futuro) existen en el proyecto local pero no se versionan por su peso. El código del módulo de cosecha por comunidad y el del visor web **no están en este repositorio**. El procedimiento está descrito en el trabajo: zonas inundables aptas (elevación < 300 m, pendiente < 5° y buffer de 1 km de ríos principales) intersectadas con el territorio de cada comunidad, y la fórmula `potencial (t/año) = área inundable (ha) × 0,13 t/ha × idoneidad media (0-1) × factor de accesibilidad`, con la accesibilidad por rangos de distancia a la red vial nacional (< 5, 5-10, 10-20, 20-50 y > 50 km; el trabajo no detalla los valores numéricos de cada rango). Los datos del visor (`mapa_data.js`, 13 de mayo de 2026) permiten verificar los resultados: **reproducen exactamente 382 comunidades con potencial positivo y 4.090 t/año (4.089,8)**, encabezadas por Matsés (172 t/año) y Puerto Elvira (135 t/año); el detalle por comunidad está en [`resultados/potencial_cosecha_comunidades.csv`](./resultados/potencial_cosecha_comunidades.csv). La verificación también mostró tres diferencias con los textos: (1) el tercer puesto del trabajo ("Santa María de Fátima", 98 t/año) no figura en esos datos, donde el tercero es Fray Pedro – Río Yavarí (123 t/año); (2) un borrador del trabajo menciona 94 comunidades con más de 10 t/año y 22 con más de 50 t/año, mientras que los datos dan 121 y 7; (3) 71 de los 2.002 registros de origen son duplicados exactos (23 con potencial positivo, unas 190 t/año contadas dos veces): sin ellos quedan 359 comunidades y ~3.900 t/año. El script que generó estos datos sigue sin estar disponible. El coeficiente de 0,13 t/ha se atribuye en el trabajo a "literatura y datos de campo de la región Loreto" sin citas; un informe técnico de mayo de 2026 para Pachamama / CCPP lo ubica en Pacaya Samiria.

## Limitaciones conocidas (a la luz de las pruebas del laboratorio)

Esta extensión se hizo antes de las pruebas de robustez de Chile y Argentina (ver [`LABORATORIO.md`](../../LABORATORIO.md)) y no se le aplicaron. En particular:

- **La proyección usa un solo modelo climático (MIROC6).** En Chile, la cifra de un único modelo (−46,2 %) resultó ser un caso atípico frente a un ensamble de seis (mediana −0,2 %). El −81,5 % de Perú no se puso a prueba de la misma manera: conviene leerlo como el resultado de un modelo, no como una proyección consensuada.
- **Validación aleatoria, sin thinning espacial real ni bloques.** Las métricas probablemente están infladas, como en las otras dos aplicaciones.
- **Ensemble ponderado por TSS** (`EMwmeanByTSS`), el esquema que luego se reemplazó por `EMmean` en el pipeline principal.
- **Sobre las métricas del ensemble.** El trabajo reporta AUCroc 0,96 y TSS 0,78 para el ensemble. Contra los modelos guardados, los valores por algoritmo coinciden con la tabla del trabajo (validación: AUCroc 0,86 a 0,97 y TSS 0,62 a 0,76; medias de los cinco: 0,935 y 0,70), pero el 0,96 / 0,78 del ensemble no se reproduce: el objeto de ensemble de biomod2 solo trae métricas de calibración (AUCroc 0,98, TSS 0,87) y no tiene validación propia.

## Metodología

El núcleo es el mismo del pipeline principal ([`ENM_BGEN_pipeline.R`](../../ENM_BGEN_pipeline.R)); lo que cambia es la leyenda de reclasificación de uso del suelo, adaptada a las clases de MapBiomas Perú. Mapeo que usa el script (y del que salen las superficies de arriba): Compatible 3, 4, 6, 11, 12, 13, 68; Restaurable 9, 15, 21, 35; Incompatible 18, 24, 25, 30; excluidas (agua, afloramientos, sin datos) 23, 27, 29, 33, 34. El texto del trabajo postulado lista 25 y 33 como compatibles, lo que no coincide con el script.
