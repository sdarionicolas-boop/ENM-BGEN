import os
import re
import markdown
import shutil
from playwright.sync_api import sync_playwright

MD_PATH = r"C:\Users\sdari\Desktop\BGEN\ENM_jacaranda\Reporte_Tecnico_Ranking_Partidos.md"
PDF_PATH = r"C:\Users\sdari\Desktop\BGEN\ENM_jacaranda\Reporte_Tecnico_Ranking_Partidos.pdf"
PDF_REPO_PATH = r"C:\Users\sdari\Desktop\ENM-BGEN-main\Reporte_Tecnico_Ranking_Partidos.pdf"

def replace_mermaid_with_flowchart(text):
    mermaid_pattern = r"```mermaid\s+graph TD.*?end\s+```"
    # Ensure no leading indentation inside each line of flow_html so Markdown doesn't treat it as code
    flow_html = """<div class="flowchart">
    <div class="flow-row">
        <div class="flow-col">
            <div class="flow-box input-box">
                <div class="box-title">Insumos de Distribución</div>
                <div class="box-text">19 rasters de probabilidad de idoneidad climática (biomod2)</div>
            </div>
        </div>
        <div class="flow-col">
            <div class="flow-box input-box">
                <div class="box-title">Límites Políticos</div>
                <div class="box-text">Capa vectorial GADM Nivel 2 (Argentina)</div>
            </div>
            <div class="flow-arrow">↓</div>
            <div class="flow-box proc-box">
                <div class="box-title">Filtrado Territorial</div>
                <div class="box-text">Selección de 128 partidos (Prov. Buenos Aires)</div>
            </div>
        </div>
    </div>
    <div class="flow-arrow-join">↓</div>
    <div class="flow-row">
        <div class="flow-box proc-box main-proc">
            <div class="box-title">Cruce Espacial en R (BLOQUE 8)</div>
            <div class="box-text">Alineación de CRS, recorte espacial de rasters y extracción de valores de celdas por partido</div>
        </div>
    </div>
    <div class="flow-arrow-split">↓</div>
    <div class="flow-row">
        <div class="flow-box metric-box">
            <div class="box-title">Estadísticas por Partido</div>
            <div class="box-text">Cálculo de Idoneidad Promedio (&mu;<sub>id</sub>) e Idoneidad Máxima (Max<sub>id</sub>)</div>
        </div>
        <div class="flow-box metric-box">
            <div class="box-title">Superficie de Alta Idoneidad</div>
            <div class="box-text">Área en km² con idoneidad &ge; 0.6 usando elipsoide (cellSize)</div>
        </div>
    </div>
    <div class="flow-arrow-join">↓</div>
    <div class="flow-row">
        <div class="flow-box index-box">
            <div class="box-title">Índice de Prioridad de Colecta BGEN (I<sub>pc</sub>)</div>
            <div class="box-text">Suma ponderada: 70% idoneidad media global + 30% riqueza de especies aptas</div>
        </div>
    </div>
    <div class="flow-arrow-join">↓</div>
    <div class="flow-row outputs-row">
        <div class="flow-box output-box">
            <div class="box-title">Ranking Consolidado</div>
            <div class="box-text">Matriz de priorización general (Top 15 y CSV unificado)</div>
        </div>
        <div class="flow-box output-box">
            <div class="box-title">Capa GeoJSON</div>
            <div class="box-text">Datos espaciales listos para SIG de escritorio (QGIS)</div>
        </div>
        <div class="flow-box output-box">
            <div class="box-title">Visor Web Leaflet</div>
            <div class="box-text">Mapa interactivo HTML móvil para uso en campo</div>
        </div>
    </div>
</div>"""
    text = re.sub(mermaid_pattern, flow_html, text, flags=re.DOTALL)
    return text

def parse_markdown_with_math(md_text):
    placeholders = {}
    counter = 0
    
    # 1. Match display math
    def display_repl(match):
        nonlocal counter
        placeholder = f"<!--MATH_DISPLAY_{counter}-->"
        placeholders[placeholder] = match.group(0)
        counter += 1
        return placeholder
        
    md_text = re.sub(r"\$\$.*?\$\$", display_repl, md_text, flags=re.DOTALL)
    
    # 2. Match inline math
    def inline_repl(match):
        nonlocal counter
        placeholder = f"<!--MATH_INLINE_{counter}-->"
        placeholders[placeholder] = match.group(0)
        counter += 1
        return placeholder
        
    md_text = re.sub(r"\$.*?\$", inline_repl, md_text)
    
    # 3. Parse markdown
    html = markdown.markdown(md_text, extensions=['tables'])
    
    # 4. Restore math
    for placeholder, math in placeholders.items():
        html = html.replace(placeholder, math)
        
    return html

def main():
    print("Leyendo archivo markdown...")
    if not os.path.exists(MD_PATH):
        print(f"Error: No existe el archivo Markdown en {MD_PATH}")
        return
        
    with open(MD_PATH, "r", encoding="utf-8") as f:
        md_content = f.read()
        
    # Normalizar retornos de carro
    md_content = md_content.replace('\r\n', '\n')
    
    print("Extrayendo metadatos y separando del cuerpo...")
    # Buscamos separar por el primer separador ---
    parts = md_content.split('\n---\n', 1)
    if len(parts) > 1:
        body_md = parts[1]
    else:
        body_md = md_content
        
    print("Reemplazando diagrama Mermaid con estructura HTML...")
    body_md = replace_mermaid_with_flowchart(body_md)
    
    print("Convirtiendo Markdown a HTML (preservando LaTeX)...")
    html_body = parse_markdown_with_math(body_md)
    
    # Creamos el HTML completo con estilos premium e integración de Google Fonts y MathJax
    html_document = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Reporte Técnico - BGEN</title>
    <!-- Carga de fuentes premium -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Outfit:wght@500;700;800&display=swap" rel="stylesheet">
    
    <!-- Configuración e inicio de MathJax -->
    <script>
    window.MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      }},
      options: {{
        ignoreHtmlClass: 'tex2jax_ignore',
        processHtmlClass: 'tex2jax_process'
      }}
    }};
    </script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            color: #2c3e35;
            line-height: 1.6;
            font-size: 10pt;
            background-color: #ffffff;
            margin: 0;
            padding: 0;
        }}
        
        /* Cabecera Premium */
        .report-header {{
            background: linear-gradient(135deg, #1e3d30 0%, #12241c 100%);
            color: #ffffff;
            padding: 35px 30px;
            border-radius: 12px;
            margin-bottom: 30px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        .header-badge {{
            display: inline-block;
            background-color: rgba(255,255,255,0.15);
            color: #94d2bd;
            font-size: 8pt;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 1.5px;
            padding: 4px 10px;
            border-radius: 20px;
            margin-bottom: 15px;
        }}
        .header-title {{
            font-family: 'Outfit', sans-serif;
            font-size: 18pt;
            font-weight: 800;
            line-height: 1.25;
            margin: 0 0 25px 0;
            color: #ffffff;
            border: none;
            padding: 0;
        }}
        .header-meta-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 15px;
            border-top: 1px solid rgba(255,255,255,0.15);
            padding-top: 20px;
        }}
        .meta-item {{
            display: flex;
            flex-direction: column;
        }}
        .meta-label {{
            font-size: 7.5pt;
            font-weight: 700;
            color: #94d2bd;
            letter-spacing: 1px;
            margin-bottom: 3px;
        }}
        .meta-val {{
            font-size: 9.5pt;
            color: #e5e5e5;
        }}
        
        /* Contenedor del documento */
        .content-container {{
            padding: 0 5px;
        }}
        
        /* Títulos de sección */
        h2 {{
            font-family: 'Outfit', sans-serif;
            font-size: 14pt;
            color: #1e3d30;
            margin-top: 30px;
            margin-bottom: 12px;
            border-bottom: 1px solid #e2e8f0;
            padding-bottom: 6px;
            page-break-after: avoid;
        }}
        h3 {{
            font-family: 'Outfit', sans-serif;
            font-size: 11pt;
            color: #2d5a27;
            margin-top: 20px;
            margin-bottom: 8px;
            page-break-after: avoid;
        }}
        
        p {{
            margin-top: 0;
            margin-bottom: 12px;
            text-align: justify;
        }}
        
        ul, ol {{
            margin-top: 0;
            margin-bottom: 15px;
            padding-left: 20px;
        }}
        li {{
            margin-bottom: 6px;
        }}
        
        /* Tablas */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            page-break-inside: avoid;
        }}
        th {{
            background-color: #1e3d30;
            color: #ffffff;
            font-family: 'Outfit', sans-serif;
            font-weight: 700;
            font-size: 8.5pt;
            padding: 8px 10px;
            border: 1px solid #d0d7de;
            text-align: center;
        }}
        td {{
            padding: 8px 10px;
            border: 1px solid #d0d7de;
            font-size: 8.5pt;
            color: #2c3e35;
        }}
        tr:nth-child(even) {{
            background-color: #f4f6f5;
        }}
        
        /* Alineación especial de columnas para la tabla de prioridades */
        table td:first-child, table th:first-child {{
            text-align: center;
            font-weight: bold;
        }}
        table td:nth-child(3), table th:nth-child(3),
        table td:nth-child(4), table th:nth-child(4),
        table td:nth-child(5), table th:nth-child(5),
        table td:nth-child(6), table th:nth-child(6) {{
            text-align: right;
        }}
        
        /* Ecuaciones matemáticas */
        .MathJax_Display, .formula-block {{
            background-color: #f4f6f5;
            border-left: 4px solid #2d5a27;
            padding: 12px 15px;
            margin: 18px 0;
            border-radius: 0 6px 6px 0;
            font-size: 10.5pt;
            overflow-x: auto;
        }}
        
        /* Diagrama de flujo */
        .flowchart {{
            display: flex;
            flex-direction: column;
            align-items: center;
            margin: 25px 0;
            gap: 0;
            page-break-inside: avoid;
        }}
        .flow-row {{
            display: flex;
            justify-content: center;
            gap: 15px;
            width: 100%;
        }}
        .flow-col {{
            display: flex;
            flex-direction: column;
            align-items: center;
        }}
        .flow-box {{
            width: 220px;
            padding: 10px 12px;
            border-radius: 8px;
            border: 1px solid #d0d7de;
            background-color: #ffffff;
            box-shadow: 0 2px 4px rgba(0,0,0,0.04);
            text-align: center;
        }}
        .box-title {{
            font-size: 9pt;
            font-weight: 700;
            color: #1e3d30;
            margin-bottom: 4px;
            font-family: 'Outfit', sans-serif;
        }}
        .box-text {{
            font-size: 7.5pt;
            color: #57606a;
            line-height: 1.3;
        }}
        .input-box {{
            border-top: 4.5px solid #005f73;
        }}
        .proc-box {{
            border-top: 4.5px solid #0a9396;
        }}
        .main-proc {{
            width: 340px;
            background-color: #f4f9f9;
        }}
        .metric-box {{
            border-top: 4.5px solid #94d2bd;
            width: 200px;
        }}
        .index-box {{
            border-top: 4.5px solid #ee9b00;
            width: 320px;
            background-color: #fefae0;
        }}
        .output-box {{
            border-top: 4.5px solid #ae2012;
            width: 150px;
        }}
        .flow-arrow {{
            font-size: 14pt;
            color: #0a9396;
            margin: 4px 0;
            font-weight: bold;
        }}
        .flow-arrow-join, .flow-arrow-split {{
            font-size: 14pt;
            color: #1e3d30;
            margin: 6px 0;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="report-header">
        <div class="header-badge">Programa BGEN &bull; UNAJ</div>
        <h1 class="header-title">Reporte Técnico: Análisis de Idoneidad Climática por Partido y Ranking de Priorización Territorial para la Colecta de Germoplasma</h1>
        <div class="header-meta-grid">
            <div class="meta-item">
                <span class="meta-label">PROYECTO</span>
                <span class="meta-val">Beca de Iniciación a la Investigación (BIEI 2025) – UNAJ</span>
            </div>
            <div class="meta-item">
                <span class="meta-label">FECHA DE EMISIÓN</span>
                <span class="meta-val">13 de julio de 2026</span>
            </div>
            <div class="meta-item">
                <span class="meta-label">PROGRAMA</span>
                <span class="meta-val">Banco de Germoplasma de Especies Nativas (BGEN-Noreste)</span>
            </div>
            <div class="meta-item">
                <span class="meta-label">INVESTIGADORES</span>
                <span class="meta-val">Darío Nicolás Sánchez (Autor) &bull; Ing. Andrea Quinteros (Directora)</span>
            </div>
        </div>
    </div>
    
    <div class="content-container tex2jax_process">
        {html_body}
    </div>
</body>
</html>"""

    print("Compilando PDF usando Playwright (Headless)...")
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            page = browser.new_page()
            
            # Cargar el HTML en el navegador
            page.set_content(html_document)
            
            # Esperar a que MathJax termine de renderizar todas las ecuaciones
            print("  Esperando a que MathJax cargue y renderice las formulas...")
            page.wait_for_selector(".MathJax", timeout=10000) # espera que aparezca al menos una ecuación renderizada
            
            # Evaluar y esperar el promise del ciclo de renderizado de MathJax
            page.evaluate("window.MathJax.startup.promise.then(() => window.MathJax.typesetPromise())")
            
            print("  Generando archivo PDF...")
            page.pdf(
                path=PDF_PATH,
                format="A4",
                print_background=True,
                margin={
                    "top": "2.5cm",
                    "bottom": "2.5cm",
                    "left": "2cm",
                    "right": "2cm"
                },
                display_header_footer=True,
                header_template='''
                    <div style="font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 7.5pt; width: 100%; margin-left: 2cm; margin-right: 2cm; display: flex; justify-content: space-between; border-bottom: 0.5px solid #e2e8f0; padding-bottom: 3px; color: #718096;">
                        <span>Reporte Tecnico: Distribucion Potencial y Priorizacion (BGEN)</span>
                        <span>Julio 2026</span>
                    </div>
                ''',
                footer_template='''
                    <div style="font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 7.5pt; width: 100%; margin-left: 2cm; margin-right: 2cm; display: flex; justify-content: space-between; border-top: 0.5px solid #e2e8f0; padding-top: 3px; color: #718096;">
                        <span>Banco de Germoplasma de Especies Nativas (BGEN-Noreste)</span>
                        <span>Pagina <span class="pageNumber"></span> de <span class="totalPages"></span></span>
                    </div>
                '''
            )
            browser.close()
            
        print(f"[OK] PDF generado con exito en: {PDF_PATH}")
        
        # Copiar al directorio local del repositorio si existe
        if os.path.exists(os.path.dirname(PDF_REPO_PATH)):
            shutil.copy2(PDF_PATH, PDF_REPO_PATH)
            print(f"[OK] Copia local actualizada en: {PDF_REPO_PATH}")
        else:
            print(f"[AVISO] No se pudo copiar a {PDF_REPO_PATH} porque el directorio no existe.")
            
    except Exception as e:
        print("[ERROR] durante la generacion del PDF:")
        print(e)

if __name__ == "__main__":
    main()
