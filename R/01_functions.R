#' Load and Clean Medical Data
#'
#' Lee archivos médicos (CSV o DATA). Si detecta el dataset de Breast Cancer (wdbc),
#' inyecta automáticamente los nombres de las columnas clínicas siguiendo el estándar
#' de la UCI (Mean, SD, Worst).
#'
#' @param path String. Ruta al archivo .csv o .data.
#' @return A tibble con nombres estandarizados y limpios.
#' @export
load_medical_data <- function(path) {
  
  # Detectamos si es el archivo específico de Breast Cancer de la UCI
  is_wdbc <- grepl("wdbc", path)
  
  if (is_wdbc) {
    # Definimos las 10 características base
    base_feats <- c("radius", "texture", "perimeter", "area", "smoothness", 
                    "compactness", "concavity", "concave_points", "symmetry", "fractal_dimension")
    
    # Creamos los nombres con los sufijos clínicos solicitados
    wdbc_names <- c("id", "diagnosis", 
                    paste0(base_feats, "_mean"), 
                    paste0(base_feats, "_sd"), 
                    paste0(base_feats, "_worst"))
    
    # Cargamos el archivo sin cabecera (header = FALSE) e inyectamos el vector de nombres
    data <- readr::read_csv(path, col_names = wdbc_names, show_col_types = FALSE)
  } else {
    # Carga estándar para otros datasets con cabecera (como CVD)
    data <- readr::read_csv(path, show_col_types = FALSE)
  }
  
  # janitor::clean_names() asegura que todo sea snake_case consistente
  return(janitor::clean_names(data))
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









#' Purga Variables No Informativas y data leakage
#'
#' Esta función elimina variables que no aportan informacion o que compromenten al modelo
#' por ejemplo alguna variable repetida, alguna variable que aporte *data leakage* (atributos
#' derivados directa o indirectamente del target) y metadatos no predictivos
#' (por ejemplo, identificadores únicos), con el objetivo de garantizar
#' un pipeline de modelado válido y libre de información espuria.
#'
#' @param df_cvd Dataframe crudo del dataset de riesgo cardiovascular (CVD).
#' @param df_cancer Dataframe crudo del dataset de cáncer de mama.
#' @return Una lista con ambos datasets depurados y listos para el preprocesamiento.
#' @export
pre_filtering <- function(df_cvd, df_cancer) {
  
  # Purga CVD: Eliminamos scores, categorías redundantes
  # y representaciones complejas de variables clínicas
  cvd_purged <- df_cvd %>% 
    dplyr::select(
      -cvd_risk_score, 
      -blood_pressure_category, 
      -blood_pressure_mm_hg,
      -height_cm
    )
  
  # Purga Breast Cancer: Eliminamos el identificador único del paciente
  cancer_purged <- df_cancer %>% 
    dplyr::select(-id)
  
  return(list(cvd = cvd_purged, cancer = cancer_purged))
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
type_data <- function(df, cat_cols) {
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











#' Imputación de Valores Ausentes mediante KNN
#'
#' @param df Dataframe con valores NA.
#' @param k Número de vecinos para el cálculo de distancias (por defecto 5).
#' @return Dataframe completo con los valores ausentes imputados.
#' @export
impute_data_knn <- function(df, k = 5) {
  
  # Aplicamos kNN de la librería VIM
  # imp_var = FALSE evita la creación de columnas auxiliares de control
  df_imputed <- VIM::kNN(df, k = k, imp_var = FALSE)
  
  return(df_imputed)
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
      ratio_calc = abdominal_circumference_cm / (height_m * 100),
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














#' Prepare and Export Triple-Flavor Datasets
#'
#' Genera las tres variantes necesarias para el modelado (Clustering, Árboles, Logística).
#' Además, estandariza estéticamente todos los factores a mayúsculas.
#'
#' @param df Dataframe limpio y auditado.
#' @return Una lista con tres dataframes: $cluster, $trees, $logistic.
#' @export
prepare_model_variants <- function(df) {
  
  # 0. Paso Estético: Convertir niveles de factores a MAYÚSCULAS
  # Esto asegura consistencia visual en los gráficos futuros
  df_upper <- df %>%
    dplyr::mutate(dplyr::across(where(is.factor), ~ {
      # Obtenemos los niveles, los pasamos a mayúsculas y reasignamos
      levels(.x) <- toupper(levels(.x))
      .x
    }))
  
  # 1. Variante CLUSTERING: Solo numéricas + Estandarizadas (Z-Score)
  df_cluster <- df_upper %>%
    dplyr::select(where(is.numeric)) %>% # Eliminamos categóricas
    scale() %>%                          # Estandarizamos (media 0, sd 1)
    as.data.frame()                      # scale() devuelve matriz, lo devolvemos a df
  
  # 2. Variante TREES: Todo el dataset + Escala Original
  # Los árboles prefieren los datos "humanos" y manejan bien los factores
  df_trees <- df_upper
  
  # 3. Variante LOGÍSTICA: Factores + Numéricas Estandarizadas
  # La regresión necesita factores (para dummies) pero requiere números escalados
  df_logistic <- df_upper %>%
    dplyr::mutate(dplyr::across(where(is.numeric), scale))
  
  # Devolvemos una lista con los 3 sabores
  return(list(
    cluster  = df_cluster,
    trees    = df_trees,
    logistic = df_logistic
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












#' Limpieza Final de Redundancias para Logística
#'
#' Recibe los datasets de logística y elimina las colinealidades críticas.
#' @param df_cvd Dataset cvd_logistic.
#' @param df_cancer Dataset cancer_logistic.
#' @export
clean_redundancies_logistic <- function(df_cvd, df_cancer) {
  
  # 1. Poda CVD: Eliminamos precursores del BMI y Ratio + LDL estimado
  cvd_final <- df_cvd %>%
    dplyr::select(
      -dplyr::any_of(c("height_m", "height_cm", "weight_kg", 
                       "abdominal_circumference_cm", "estimated_ldl_mg_d_l"))
    )
  
  # 2. Poda Cáncer: Triaje basado en correlación > 0.85 y Wolberg
  # Eliminamos radios, perímetros, concavidades y medias redundantes
  vars_to_remove_cancer <- c(
    "radius_mean", "radius_sd", "radius_worst",
    "perimeter_mean", "perimeter_sd", "perimeter_worst",
    "concavity_mean", "concavity_sd", "concavity_worst",
    "area_mean", "texture_mean", "concave_points_mean", "compactness_mean"
  )
  
  cancer_final <- df_cancer %>% 
    dplyr::select(-dplyr::any_of(vars_to_remove_cancer))
  
  return(list(cvd = cvd_final, cancer = cancer_final))
}









#' Análisis de Componentes Principales "By Hand"
#' 
#' Ejecuta el flujo matemático completo del PCA sin librerías externas.
#' 
#' @param df Dataset numérico.
#' @param title Título para el Scree Plot.
#' @export
perform_manual_pca <- function(df, title = "PCA: Análisis de Varianza") {
  
  # 1. Aseguramos Centrado y Escalado (Matriz Z)
  # Aunque vengan escalados, por seguridad matemática lo aplicamos
  X <- as.matrix(df)
  X_c <- scale(X, center = TRUE, scale = TRUE)
  n <- nrow(X_c)
  
  # 2. Cálculo de la Matriz de Covarianzas (S)
  # S = (1/(n-1)) * Xt * X
  S <- (1 / (n - 1)) * (t(X_c) %*% X_c)
  
  # 3. Descomposición en Autovalores (lambda) y Autovectores (w)
  # Usamos la función base eigen() que resuelve Sw = lambda*w
  ev <- eigen(S)
  autovalores <- ev$values
  autovectores <- ev$vectors
  
  # Asignamos nombres para que sea interpretable
  colnames(autovectores) <- paste0("PC", 1:ncol(autovectores))
  rownames(autovectores) <- colnames(df)
  
  # 4. Proyección de los datos (Scores)
  # Z = Xc %*% W
  scores <- X_c %*% autovectores
  colnames(scores) <- paste0("PC", 1:ncol(scores))
  
  # 5. Cálculo de la Varianza Explicada (PVE)
  pve <- autovalores / sum(autovalores)
  pve_cum <- cumsum(pve)
  
  # --- Visualización (Scree Plot) ---
  plot_data <- data.frame(
    PC = 1:length(pve),
    Var = pve,
    CumVar = pve_cum
  )
  
  p <- ggplot(plot_data, aes(x = PC)) +
    geom_bar(aes(y = Var), stat = "identity", fill = "#2C3E50", alpha = 0.7) +
    geom_line(aes(y = CumVar), group = 1, color = "#E74C3C", size = 1) +
    geom_point(aes(y = CumVar), color = "#E74C3C", size = 2) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "#27AE60") +
    annotate("text", x = 2, y = 0.88, label = "UMBRAL 85%", color = "#27AE60", fontface = "bold") +
    labs(title = title,
         subtitle = "Implementación matemática manual (Eigen-decomposition)",
         x = "Componentes Principales", y = "% Varianza") +
    theme_minimal()
  
  print(p)
  
  # Devolvemos una lista con todo el rastro matemático
  return(list(
    autovalores = autovalores,
    loadings = autovectores,
    scores = scores,
    pve = pve,
    pve_cum = pve_cum
  ))
}