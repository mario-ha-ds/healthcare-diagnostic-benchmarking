# Healthcare Diagnostic Benchmarking: Complexity Analysis in Medical Data 🏥

Este proyecto presenta un pipeline comparativo de minería de datos aplicado al diagnóstico médico. El objetivo central es evaluar cómo responden diferentes modelos (supervisados y no supervisados) ante dos naturalezas de datos opuestas: el **continuo biológico** frente a **marcadores bioquímicos discretos**.

## 🎯 Motivación del Proyecto
En el mundo real, los datasets médicos no siempre permiten una separación clara de clases. Este proyecto analiza el "gap" de rendimiento y complejidad entre:
1.  **Cardiovascular Disease (CVD):** Un escenario de alta entropía donde los factores de riesgo se solapan (continuum), dificultando el clustering y la clasificación.
2.  **Hepatitis C (HCV):** Un escenario con marcadores bioquímicos nítidos, utilizado aquí como benchmark de validación del pipeline.

## 📊 Los Datasets
* **CVD Risk Dataset ([Mendeley Data](https://data.mendeley.com/datasets/d9scg7j8fp/1)):** Datos de riesgo cardiovascular que requieren alta dimensionalidad (PCA) para ser explicados.
* **Hepatitis C Virus ([UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/503/hepatitis+c+virus+hcv+for+egyptian+patients)):** Datos de pacientes egipcios con estadios definidos de fibrosis y cirrosis.

## 🛠️ Stack Tecnológico
* **Lenguaje:** R 4.x
* **Gestión de Entorno:** `renv` (Reproducibilidad garantizada).
* **Librerías Clave:** `tidyverse`, `caret`, `VIM` (kNN Imputation), `C50`, `dbscan`, `factoextra`, `patchwork`.

## 📁 Estructura del Proyecto
```text
healthcare-diagnostic-benchmarking/
├── data/
│   ├── raw/                 # Datos originales (CVD y Hepatitis C)
│   └── clean/               # Datos procesados (.rds) para análisis
├── R/
│   └── 01_functions.R       # Motor de funciones modularizadas
├── notebooks/
│   ├── 01_eda_and_cleaning.Rmd  # Preprocesamiento y PCA Diagnóstico
│   ├── 02_clustering_study.Rmd  # Análisis no supervisado (K-means, DBSCAN)
│   └── 03_supervised_models.Rmd # Clasificación y Matrices de Coste
├── output/                  # Visualizaciones y gráficos exportados
└── renv.lock                # Bloqueo de versiones de librerías

## 🚀 Cómo Reproducir este Proyecto

Este proyecto utiliza `renv` para asegurar que las versiones de las librerías sean las mismas que las utilizadas en el desarrollo original.

1. Clona el repositorio:
```bash
git clone [https://github.com/mario-ha-ds/healthcare-diagnostic-benchmarking.git](https://github.com/mario-ha-ds/healthcare-diagnostic-benchmarking.git)

```


2. Abre el archivo `.Rproj` en RStudio.
3. Restaura el entorno (RStudio detectará `renv` automáticamente, si no, ejecuta en la consola):
```r
renv::restore()

```


4. Ejecuta los notebooks en orden (`01`, `02`, `03`).

## 🧠 Metodología Destacada

* **Tratamiento de Outliers:** Diferenciación entre errores de medida (CVD) y señales biológicas críticas (Hepatitis).
* **Diagnóstico de Complejidad:** Uso de PCA como métrica para cuantificar la dispersión de la información.
* **Safety-First Triage:** Implementación de **Matrices de Coste** en modelos C5.0, priorizando la Sensibilidad (Recall) sobre el Accuracy para minimizar Falsos Negativos en diagnósticos críticos.

---

**Autor:** [Tu Nombre]
**Contacto:** [Tu LinkedIn o Email]
