#' Load and Clean Column Names
#'
#' Reads a CSV file from a given path and standardizes column names 
#' using janitor::clean_names() for a consistent workflow.
#'
#' @param path String. The relative or absolute path to the .csv file.
#' @return A tibble with standardized column names.
#' @export
load_medical_data <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
}




#' Get the first n unique valid values from each column
#'
#' Iterates through a dataframe to extract unique values, ignoring
#' NAs and empty strings. This helps to understand the data nature 
#' (binary, categorical, or continuous).
#'
#' @param df The dataframe to inspect.
#' @param n The maximum number of unique values to return per column (default 10).
#' @return A list where each element contains unique values for a specific column.
#' @export
get_unique_valid_values <- function(df, n = 10) {
  purrr::map(df, function(col) {
    # 1. Convertimos a character para manejar tipos mixtos y espacios vacíos uniformemente
    col_char <- as.character(col)
    
    # 2. Filtro: Mantenemos valores que NO sean NA y que NO sean solo espacios en blanco
    valid_values <- col[!is.na(col_char) & trimws(col_char) != ""]
    
    # 3. Retornamos los primeros n valores únicos encontrados como vector
    result <- head(unique(valid_values), n)
    return(as.vector(result))
  })
}









#' Purga de Data Leakage y Variables Temporales
#'
#' Elimina atributos que contienen información del target o datos post-diagnóstico.
#' @param df_cvd Dataframe crudo de CVD.
#' @param df_hcv Dataframe crudo de HCV.
#' @return Una lista con ambos datasets purgados.
#' @export
purge_medical_leakage <- function(df_cvd, df_hcv) {
  
  # Purga CVD: Eliminamos scores calculados y redundancias de formato
  cvd_purged <- df_cvd %>% 
    dplyr::select(-cvd_risk_score, 
                  -blood_pressure_category, 
                  -blood_pressure_mm_hg)
  
  # Purga HCV: Somos puristas. Eliminamos TODO lo que no sea Baseline (Día 0)
  # Esto incluye cualquier medición en la semana 4, 12, etc.
  hcv_purged <- df_hcv %>% 
    dplyr::select(-alt4, -alt_12, -alt_24, -alt_36, -alt_48, -alt_after_24_w,
                  -rna_4, -rna_12, -rna_eot, -rna_ef)
  
  return(list(cvd = cvd_purged, hcv = hcv_purged))
}







#' Medical Data Typing and Factor Conversion
#'
#' Generic function to standardize data types. Converts specified columns to factors 
#' and ensures all other numerical columns are correctly typed as numeric.
#'
#' @param df The dataframe to process.
#' @param cat_cols Character vector. The names of the columns to be converted to factors.
#' @return A dataframe with optimized categorical and numerical types.
#' @export
type_medical_data <- function(df, cat_cols) {
  df %>%
    dplyr::mutate(
      # Convertimos a factor solo las columnas que le pedimos
      dplyr::across(dplyr::any_of(cat_cols), as.factor),
      # El resto de numéricas las aseguramos como double/numeric
      dplyr::across(where(is.numeric) & !dplyr::any_of(cat_cols), as.numeric)
    )
}








#' Detección y Visualización de Outliers (IQR)
#'
#' 1. Calcula outliers basándose en el rango intercuartílico (1.5 * IQR).
#' 2. Imprime una tabla con las variables afectadas y el conteo de outliers.
#' 3. Genera un gráfico de Boxplots marcando los puntos en rojo.
#'
#' @param df Dataframe con los datos.
#' @param title_suffix Texto opcional para el título del gráfico.
#' @return Un objeto ggplot con los boxplots.
#' @export
analyze_outliers <- function(df, title_suffix = "") {
  
  # 1. Preparar datos numéricos
  df_numeric <- df %>% dplyr::select(where(is.numeric))
  
  df_long <- df_numeric %>%
    tidyr::pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # 2. Detección Matemática
  outlier_summary <- df_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      q1 = quantile(value, 0.25, na.rm = TRUE),
      q3 = quantile(value, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - 1.5 * iqr,
      upper_bound = q3 + 1.5 * iqr,
      # Contamos cuántos valores están fuera de los límites
      n_outliers = sum(!is.na(value) & (value < lower_bound | value > upper_bound)),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_outliers > 0) %>%
    # Seleccionamos solo las columnas deseadas para el reporte
    dplyr::select(variable, n_outliers) %>%
    dplyr::arrange(desc(n_outliers))
  
  # 3. Imprimir Reporte en Consola
  cat("\n--- REPORTE DE OUTLIERS DETECTADOS:", title_suffix, "---\n")
  if (nrow(outlier_summary) > 0) {
    print(as.data.frame(outlier_summary))
  } else {
    cat(" No se han detectado outliers en ninguna variable numérica.\n")
  }
  
  # 4. Generar Gráfico
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = "", y = value)) +
    ggplot2::geom_boxplot(
      outlier.colour = "red", 
      outlier.shape = 16, 
      outlier.size = 2,
      fill = "lightblue",
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Análisis de Distribución y Outliers:", title_suffix),
      subtitle = "Puntos rojos indican valores fuera del rango 1.5 * IQR",
      x = NULL, y = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}











#' Imputación de Valores Ausentes mediante KNN
#'
#' @param df Dataframe con valores NA (normalmente cvd_clean).
#' @param k Número de vecinos para el cálculo de distancias (por defecto 5).
#' @return Dataframe completo con los valores ausentes imputados.
#' @export
impute_data_knn <- function(df, k = 5) {
  
  # Aplicamos kNN de la librería VIM
  # imp_var = FALSE evita la creación de columnas auxiliares de control
  df_imputed <- VIM::kNN(df, k = k, imp_var = FALSE)
  
  return(df_imputed)
}




#' Corrección de Errores Lógicos (Valores Negativos)
#'
#' Identifica valores negativos en columnas numéricas donde no deberían existir
#' (como medidas bioquímicas o antropométricas) y los sustituye por la mediana.
#' @param df Dataframe a procesar.
#' @return Dataframe con los errores corregidos.
#' @export
fix_logical_errors <- function(df) {
  df_fixed <- df
  
  # Identificamos columnas numéricas
  cols_num <- df_fixed %>% select(where(is.numeric)) %>% colnames()
  
  for (col in cols_num) {
    # Si hay algún valor negativo
    if (any(df_fixed[[col]] < 0, na.rm = TRUE)) {
      # Calculamos la mediana de los valores que SÍ son válidos (> 0)
      valid_median <- median(df_fixed[[col]][df_fixed[[col]] >= 0], na.rm = TRUE)
      
      # Aplicamos la corrección
      df_fixed[[col]] <- ifelse(df_fixed[[col]] < 0, valid_median, df_fixed[[col]])
      
      message(paste("Variable", col, ": Valores negativos corregidos con la mediana (", round(valid_median, 2), ")"))
    }
  }
  
  return(df_fixed)
}






#' Transformación Logarítmica de Carga Viral
#'
#' Aplica log10 a la variable rna_base para normalizar su distribución
#' y reducir el sesgo de magnitud en algoritmos de distancias.
#' @param df Dataset de HCV.
#' @return Dataset con la variable rna_base transformada a log10.
#' @export
transform_hcv_log <- function(df) {
  
  df_transformed <- df %>%
    # Usamos log10 por ser el estándar clínico para carga viral
    mutate(rna_base = log10(rna_base + 1)) %>%
    # Renombramos para que sea explícito en el resto del notebook
    rename(log_rna_base = rna_base)
  
  message("Variable rna_base transformada a log_rna_base (escala log10).")
  
  return(df_transformed)
}







#' Auditoría de Integridad y Corrección de Variables Derivadas
#'
#' Calcula el error porcentual en BMI y Waist-to-Height Ratio y los corrige.
#' @param df Dataframe de CVD (cvd_final).
#' @param tolerance Umbral de diferencia para considerar un valor como error (defecto 0.01).
#' @export
audit_and_fix_integrity <- function(df, tolerance = 0.01) {
  
  # 1. Cálculos de control
  df_audit <- df %>%
    dplyr::mutate(
      bmi_calc = weight_kg / (height_m^2),
      ratio_calc = abdominal_circumference_cm / height_cm,
      # Detectamos incongruencias fuera de la tolerancia (por redondeos)
      bmi_error = abs(bmi - bmi_calc) > tolerance,
      ratio_error = abs(waist_to_height_ratio - ratio_calc) > tolerance
    )
  
  # 2. Cálculo de estadísticas de error
  p_error_bmi <- mean(df_audit$bmi_error, na.rm = TRUE) * 100
  p_error_ratio <- mean(df_audit$ratio_error, na.rm = TRUE) * 100
  
  # 3. Reporte por consola
  cat("\n--- REPORTE DE INTEGRIDAD DE DATOS (CVD) ---\n")
  cat(sprintf("Incongruencias detectadas en BMI: %.2f%%\n", p_error_bmi))
  cat(sprintf("Incongruencias detectadas en Waist-to-Height Ratio: %.2f%%\n", p_error_ratio))
  
  # 4. Corrección: Sustituimos los valores originales por los calculados
  df_fixed <- df_audit %>%
    dplyr::mutate(
      bmi = bmi_calc,
      waist_to_height_ratio = ratio_calc
    ) %>%
    # Limpiamos las columnas auxiliares de la auditoría
    dplyr::select(-bmi_calc, -ratio_calc, -bmi_error, -ratio_error)
  
  cat("Variables recalculadas y corregidas con precisión matemática.\n")
  
  return(df_fixed)
}







#' Estandarización de Variables Numéricas (Z-score)
#'
#' @param df Dataframe con variables numéricas y categóricas.
#' @return Dataframe con las variables numéricas escaladas (media 0, sd 1).
#' @export
standardize_data <- function(df) {
  
  # Seleccionamos y escalamos solo las columnas numéricas
  df_scaled <- df %>%
    dplyr::mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
  
  return(df_scaled)
}








#' Despacho Final de Datasets por Modelo
#' 
#' Organiza los datasets ya procesados, eliminando categóricas donde no tocan.
#' @param cvd_sc Dataset CVD escalado.
#' @param hcv_sc Dataset HCV escalado (con log ya aplicado).
#' @param cvd_sn Dataset CVD limpio/imputado (escala original).
#' @param hcv_ty Dataset HCV tipado/limpio (escala original).
#' @return Lista con los 4 datasets finales.
#' @export
finalize_datasets <- function(cvd_sc, hcv_sc, cvd_sn, hcv_ty) {
  
  # 1. Datasets para CLUSTERING (Solo numéricas, ya escaladas)
  cvd_cluster <- cvd_sc %>% dplyr::select(where(is.numeric))
  hcv_cluster <- hcv_sc %>% dplyr::select(where(is.numeric))
  
  # 2. Datasets para SUPERVISADOS (Mantienen categóricas, escala original)
  # Solo los renombramos para seguir tu nomenclatura ideal
  cvd_super <- cvd_sn
  hcv_super <- hcv_ty
  
  return(list(
    cvd_cluster = cvd_cluster,
    hcv_cluster = hcv_cluster,
    cvd_super   = cvd_super,
    hcv_super   = hcv_super
  ))
}







#' Matriz de Correlación con Números (Clásica)
#'
#' @param df Dataframe a analizar.
#' @param title Título del gráfico.
#' @export
plot_correlation_matrix <- function(df, title = "Matriz de Correlación") {
  
  df_num <- df %>% dplyr::select(where(is.numeric))
  cor_matrix <- cor(df_num, use = "complete.obs")
  
  p <- ggcorrplot::ggcorrplot(cor_matrix, 
                              hc.order = TRUE, 
                              type = "lower", 
                              lab = TRUE, 
                              lab_size = 2.5, # Números pequeños para que quepan todos
                              method = "square", 
                              colors = c("#E46726", "white", "#6D9EC1"),
                              title = title,
                              outline.color = "white",
                              ggtheme = ggplot2::theme_minimal())
  print(p)
}