# ENM-BGEN como laboratorio

Este repositorio no es un único trabajo. Es un **laboratorio**: un mismo pipeline de modelado de nicho ecológico (biomod2 + MapBiomas) aplicado en tres contextos, puesto a prueba contra sí mismo, y con el rastro completo de lo que se midió, lo que falló y lo que se corrigió. Lo que se ofrece es el proceso, además de los resultados.

Para recorrerlo en orden: [`ejemplo_minimo/`](./ejemplo_minimo/) (un minuto, una especie) → este documento → [`extensiones/chile/`](./extensiones/chile/) (la etapa más auditada) → [`extensiones/argentina/`](./extensiones/argentina/) (la auditoría aplicada al origen).

## El arco, en cuatro etapas

| Etapa | Qué se hizo | Dónde |
|---|---|---|
| 1. **Construcción** (Argentina) | Pipeline de 19 especies nativas para el Banco de Germoplasma (BGEN-UNAJ): presencias GBIF/iNaturalist, WorldClim + SRTM, 5 algoritmos, ensemble, cruce con MapBiomas Argentina 2024 y un índice de prioridad de colecta (Ipc) para los 128 partidos bonaerenses. | raíz del repo |
| 2. **Expansión** (Perú) | El mismo núcleo aplicado a *Mauritia flexuosa* (en coautoría con Martín E. Cano), con MapBiomas Perú Colección 3, una proyección climática a 2050 y un módulo socioeconómico: potencial de cosecha sostenible de aguaje para 382 comunidades nativas, con visor web público. La etapa más amplia en alcance. | [`extensiones/peru/`](./extensiones/peru/) |
| 3. **Auditoría** (Chile) | Aplicado a *Jubaea chilensis*, con la decisión de ponerlo a prueba: thinning espacial real, validación cruzada por bloques, comparación entre colecciones de MapBiomas, ensamble de 6 modelos climáticos, sensibilidad del ranking y revisión de qué cuenta como "validación". | [`extensiones/chile/`](./extensiones/chile/) |
| 4. **Vuelta al origen** (Argentina) | Las mismas pruebas sobre 4 de las 19 especies y una corrección aplicada a las 19 (ensemble `EMmean`). | [`extensiones/argentina/`](./extensiones/argentina/) |

## Lo que se aprendió

1. **La validación aleatoria infla las métricas cuando los registros están agregados.** Con thinning a 10 km y validación por bloques, el TSS medio cayó de 0,741 a 0,460 en Chile y de 0,766 a 0,340 en 4 especies de Argentina (una corrida por especie).
2. **Ponderar el ensemble por TSS aleatorio premia al sobreajuste.** Random Forest tiene calibración 1,0 y es de los peores algoritmos bajo validación espacial; el ensemble ponderado le daba más peso justamente a él. Cuál es el peor algoritmo varía por especie, así que la corrección fue estructural (media simple), no elegir un favorito.
3. **Una corrección puede importar poco, y eso también es un resultado.** Cambiar a `EMmean` movió el área de alta idoneidad una mediana de −0,71 % y no cambió ni un puesto del Top 15 del Ipc.
4. **Una validación "externa" puede ser circular.** En Chile, el 72,6 % de los registros de entrenamiento caía dentro de los polígonos usados para "validar", así que se reformuló como control de consistencia, no como validación independiente.
5. **Un solo modelo climático engaña.** Con uno, el hábitat de alta idoneidad de *Jubaea chilensis* caía 46,2 % a 2050; con seis, la mediana es −0,2 % (rango −46,2 % a +17,0 %). El hallazgo pasó a ser la incertidumbre.
6. **Cambiar de colección de MapBiomas no es conversión real.** El 95,3 % del aumento aparente de hábitat perdido entre las Colecciones 1 y 2 de Chile es reclasificación de una clase.
7. **Con pocos registros, hasta la estimación aleatoria es ruidosa.** En el ejemplo mínimo, el TSS aleatorio fue de 0,45 a 0,65 según la corrida; el de bloques, 0,403 siempre.

## Estado de cada resultado

Qué es oficial, qué es diagnóstico y qué está desactualizado, a 23 de octubre de 2026. Si vas a reutilizar una cifra, mirá primero esta tabla.

| Resultado | Estado | Nota |
|---|---|---|
| Métricas de validación, 19 especies (AUCroc 0,96, TSS 0,78) | **Oficial, con limitación** | Validación aleatoria; ver lección 1. |
| Rasters de idoneidad, 19 especies | **Oficial** (`EMmean` desde el 22/9/2026) | No se versionan (pesados); se regeneran con el pipeline. |
| Ranking Ipc de partidos (`outputs/ranking_partidos/`) | **Oficial, con limitación** | Con `EMmean`. Top 15 dominado por el conurbano (sesgo de muestreo ciudadano); el Top 10 es inestable bajo thinning + bloques (6 de 10 se mantienen, en 4 especies). |
| Cruce con MapBiomas Argentina (alta idoneidad: 170 / 117 / 31 Mha disponible / perdido / restaurable) | **Oficial** (`EMmean`, recalculado el 23/10/2026) | El paper presentado al Premio conserva sus cifras originales (171 / 119 / 31). Cambio por especie entre −7 % y +3 %; ver `extensiones/argentina/resultados/comparacion_cruce_mapbiomas_emmean_vs_original.csv`. |
| `mapa_interactivo_ENM.html` (BLOQUE 6) | **Oficial** (`EMmean`, regenerado el 23/10/2026) | No se versiona por su peso (7 MB); se regenera con el pipeline. |
| Pruebas de robustez de Argentina (`BLOQUE_ROBUSTEZ_A-E`) | **Diagnóstico** | 4 especies, una corrida cada una; la caída exacta tiene ruido. |
| Corrección `EMmean` (`BLOQUE_F-H`) | **Oficial** | Aplicada a las 19 especies. |
| Perú: cruce de *Mauritia flexuosa* con MapBiomas Perú Colección 3 (94,7 % de la alta idoneidad conserva vegetación natural) | **Histórico** | Versión enviada a la postulación, con el ensemble ponderado por TSS original; no se le aplicaron las pruebas de robustez. |
| Perú: proyección a 2050 (alta idoneidad −81,5 %) | **Histórico, sin contrastar** | Un solo modelo climático (MIROC6). En Chile, un solo modelo resultó atípico frente a seis; esta cifra no se puso a prueba así. |
| Perú: potencial de cosecha de 382 comunidades nativas (~4.090 t/año) y visor web | **Histórico, verificado con matices** | Los datos del visor reproducen exactamente 382 comunidades y 4.090 t/año. Pero 71 de los 2.002 registros están duplicados (sin ellos: 359 comunidades, ~3.900 t/año), el tercer puesto citado no figura en los datos y el código de cálculo no está disponible. El enlace del visor (`plan-aguaje.netlify.app`) hoy da 404. Ver [`extensiones/peru/`](./extensiones/peru/README.md). |
| Perú: métricas del ensemble (AUCroc 0,96 / TSS 0,78) | **A aclarar** | No se reproducen desde los modelos guardados; ver [`extensiones/peru/`](./extensiones/peru/README.md). |
| Chile: memoria técnica y cruce con MapBiomas Colección 2 | **Oficial, con limitaciones declaradas** | Enviada. Mantiene el ensemble ponderado original (declarado en su sección 7), por lo que difiere del pipeline actual. |
| Chile: proyección a 2050 con 6 modelos climáticos | **Oficial** | Mediana −0,2 %; la incertidumbre entre modelos es el hallazgo. |
| Chile: prueba de robustez (thinning 10 km + bloques) | **Diagnóstico** | Reentrenado con 71 registros; no reemplaza los resultados operativos. |
| `ejemplo_minimo/` | **Didáctico** | Capas gruesas y una especie: sirve para verificar el método, no para sacar conclusiones. |

## Qué no resuelve este laboratorio

- La autocorrelación espacial de fondo sigue sin corregirse en los resultados operativos de las 19 especies: se midió y se declara, pero no se reentrenó todo con validación espacial.
- Ninguno de estos resultados pasó por revisión de pares.
- Los datos de presencia y las capas ambientales no se redistribuyen; hay que reconstruirlos.
