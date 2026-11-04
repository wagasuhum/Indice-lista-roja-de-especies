#### Load required packages - libraries to run the script----
packages <- c("rstudioapi","sf","dplyr","terra","red", "rredlist","ggplot2", "tibble", "reshape2")
lapply(packages,require,character.only=T)

# Set as working directory where the script is saved
dirfolder<- dirname(getSourceEditorContext()$path)
setwd(dirfolder)

#Load IUCN token----
token <- "QkqZt35of6w4oR6e11d2uME7ZPR9n61gdFsS"

# Cargar datos de especies endémicas
sp_endemicas <- read.csv2("C:/Users/walter.garcia/Downloads/RLI/RLI/RLI/data/table_2024-07-22-18-20-41.csv")
colnames(sp_endemicas)[2] <-"scientific_name"
sp_endemicas <- sp_endemicas[,2] 

# Obtener datos de especies de Colombia
IUCN_sp_col <- rredlist::rl_countries('CO', key = token)
IUCN_sp_df <- IUCN_sp_col$assessments

# Filtrar por especies endémicas
sp_endemicas_IUCN <- IUCN_sp_df %>%
  filter(taxon_scientific_name %in% sp_endemicas)

# Descargar historial de categorizaciones
IUCN_history <- lapply(seq_len(nrow(sp_endemicas_IUCN)), function(i) {
  name_split <- strsplit(sp_endemicas_IUCN$taxon_scientific_name[i], " ")[[1]]
  genus <- name_split[1]
  species <- paste(name_split[-1], collapse = " ")
  
  message("Descargando historial de ", genus, " ", species, " (", i, "/", nrow(sp_endemicas_IUCN), ")")
  
  res <- tryCatch({
    rredlist::rl_species(genus = genus, species = species, key = token)
  }, error = function(e) {
    message("⚠️ Error con ", genus, " ", species, ": ", e$message)
    return(NULL)
  })
  
  Sys.sleep(runif(1, 1, 2.5)) # pausa entre consultas
  return(res)
})

# Crear dataframe con el historial incluyendo clase taxonómica
IUCN_history_DF <- data.frame()

for (i in seq_along(IUCN_history)) {
  assessments <- IUCN_history[[i]]$assessments
  
  if (is.data.frame(assessments) && nrow(assessments) > 0) {
    # Agregar información taxonómica
    assessments$species <- IUCN_history[[i]]$taxon$scientific_name
    assessments$class <- IUCN_history[[i]]$taxon$class
    
    IUCN_history_DF <- rbind(IUCN_history_DF, assessments)
  }
}

# Limpiar y guardar datos
IUCN_history_DF_clean <- as.data.frame(lapply(IUCN_history_DF, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(y, collapse = ", "))
  } else {
    x
  }
}), stringsAsFactors = FALSE)

write.csv(IUCN_history_DF_clean, "IUCN_history_DF_with_class2.csv", row.names = FALSE)

# Crear matriz de cambios de categoría por especie
library(reshape2)

historyAssesment_data <- as.data.frame(IUCN_history_DF_clean)
historyAssesment_matrix <- dcast(
  historyAssesment_data,
  species + class ~ year_published,
  value.var = "red_list_category_code",
  fun.aggregate = function(x) unique(x)[1]
)

# Ajustar categorías IUCN
adjust_categories <- data.frame(
  Cat_IUCN = c("CR", "DD", "EN", "EN", "DD", "DD", "LC", "LC", "LC", "NT", "DD", "NT", "RE", "VU", "VU"),
  code = c("CR", "DD", "E", "EN", "I", "K", "LC", "LR/cd", "LR/lc", "LR/nt", "NA", "NT", "R", "V", "VU")
)

# Crear una lista para almacenar matrices por clase
class_matrices <- list()

# Obtener clases únicas
unique_classes <- unique(historyAssesment_matrix$class)
unique_classes <- unique_classes[!is.na(unique_classes)]

for (class_name in unique_classes) {
  # Filtrar especies de esta clase
  class_data <- historyAssesment_matrix %>%
    filter(class == class_name) %>%
    select(-class)
  
  # Convertir a matriz (especies en filas, años en columnas)
  class_matrix <- class_data %>%
    column_to_rownames("species") %>%
    as.matrix()
  
  
  # Ajustar categorías
  for(i in seq(nrow(adjust_categories))){
    class_matrix[which(class_matrix == adjust_categories[i,]$code, arr.ind = TRUE)] <- adjust_categories[i,]$Cat_IUCN
  }
  
  for(j in unique(adjust_categories$Cat_IUCN)){
    key <- c(tolower(j), toupper(j), j) %>% paste0(collapse = "|")
    class_matrix[which(grepl(key, class_matrix), arr.ind = TRUE)] <- j
  }
  
  class_matrix[which((!class_matrix %in% adjust_categories$Cat_IUCN) & !is.na(class_matrix), arr.ind = TRUE)] <- NA
  
  # Guardar en la lista
  class_matrices[[class_name]] <- class_matrix
}

process_class_matrix <- function(class_matrix, class_name) {
  
  # Crear dataframe
  df <- as.data.frame.matrix(class_matrix)
  
  # Función para reemplazar NAs con valor previo
  replace_na_with_previous <- function(df) {
    if (ncol(df) > 1) {
      for (col_idx in 2:ncol(df)) {
        current_col <- df[[col_idx]]
        prev_col <- df[[col_idx - 1]]
        
        # Reemplazar NA con valor de la columna anterior
        df[[col_idx]] <- ifelse(is.na(current_col), prev_col, current_col)
      }
    }
    return(df)
  }
  
  # Aplicar reemplazo de NAs
  df <- replace_na_with_previous(df)
  
  # Ordenar columnas por año (numéricamente)
  year_cols <- colnames(df)
  # Intentar convertir a numérico
  numeric_years <- suppressWarnings(as.numeric(year_cols))
  
  # Solo ordenar si podemos convertir a numérico
  if (all(!is.na(numeric_years))) {
    df <- df[, order(numeric_years)]
  }
  
  # Encontrar columna del año base (2000)
  target_year <- "2000"
  if (target_year %in% colnames(df)) {
    base_year_index <- which(names(df) == target_year)
    
    # Eliminar especies sin evaluación en el año base
    df_clean <- df[!is.na(df[, base_year_index]), , drop = FALSE]
    
    if (nrow(df_clean) > 0) {
      # Seleccionar columnas desde el año base en adelante
      if (base_year_index < ncol(df)) {
        df_filtered <- df_clean %>%
          select(all_of(names(df_clean)[base_year_index:ncol(df_clean)]))
      } else {
        df_filtered <- df_clean
      }
      
      # Calcular RLI para esta clase - MANEJO ESPECÍFICO PARA PAQUETE 'red'
      # Primero verificar qué devuelve rli()
      test_result <- tryCatch({
        red::rli(df_filtered, boot = FALSE)
      }, error = function(e) {
        message("  Error en red::rli(): ", e$message)
        return(NULL)
      })
      
      if (is.null(test_result)) {
        return(NULL)
      }
      
      # Determinar la estructura del resultado
      if (is.vector(test_result)) {
        # Si es un vector simple
        rli_values <- as.numeric(test_result)
        years <- colnames(df_filtered)
        
        # Verificar longitud
        if (length(rli_values) != length(years)) {
          message("  Advertencia: longitud de RLI (", length(rli_values), 
                  ") no coincide con años (", length(years), ")")
          
          # Si es más largo, tomar solo los primeros n años
          if (length(rli_values) > length(years)) {
            rli_values <- rli_values[1:length(years)]
          }
        }
        
        rli_df <- data.frame(
          Year = years,
          RLI = rli_values,
          Class = class_name,
          N_Especies = nrow(df_filtered),
          stringsAsFactors = FALSE
        )
        
      } else if (is.data.frame(test_result) || is.matrix(test_result)) {
        # Si es un dataframe o matriz
        test_result <- as.data.frame(test_result)
        
        # Buscar la columna con valores RLI
        rli_col <- NULL
        possible_names <- c("RLI", "rli", "value", "index", "V1")
        
        for (col_name in possible_names) {
          if (col_name %in% colnames(test_result)) {
            rli_col <- col_name
            break
          }
        }
        
        if (is.null(rli_col)) {
          # Tomar la primera columna
          rli_col <- colnames(test_result)[1]
        }
        
        years <- colnames(df_filtered)
        
        # Verificar longitud
        if (nrow(test_result) != length(years)) {
          message("  Advertencia: filas de resultado (", nrow(test_result), 
                  ") no coincide con años (", length(years), ")")
          
          # Si hay más filas, tomar solo las necesarias
          if (nrow(test_result) > length(years)) {
            test_result <- test_result[1:length(years), , drop = FALSE]
          }
        }
        
        rli_df <- data.frame(
          Year = years,
          RLI = as.numeric(test_result[[rli_col]]),
          Class = class_name,
          N_Especies = nrow(df_filtered),
          stringsAsFactors = FALSE
        )
        
      } else {
        message("  Formato de resultado no reconocido")
        return(NULL)
      }
      
      return(list(matrix = class_matrix, rli_data = rli_df, filtered_data = df_filtered))
    }
  } else {
    message("  Año 2000 no encontrado en los datos para ", class_name)
    
    # Usar todos los años disponibles
    if (ncol(df) > 0 && nrow(df) > 0) {
      # Filtrar especies que tienen al menos una evaluación
      df_clean <- df[apply(df, 1, function(x) any(!is.na(x))), , drop = FALSE]
      
      if (nrow(df_clean) > 0) {
        df_filtered <- df_clean
        
        # Calcular RLI
        test_result <- tryCatch({
          red::rli(df_filtered, boot = FALSE)
        }, error = function(e) {
          message("  Error en red::rli(): ", e$message)
          return(NULL)
        })
        
        if (is.null(test_result)) {
          return(NULL)
        }
        
        years <- colnames(df_filtered)
        
        if (is.vector(test_result)) {
          rli_values <- as.numeric(test_result)
          
          if (length(rli_values) != length(years)) {
            if (length(rli_values) > length(years)) {
              rli_values <- rli_values[1:length(years)]
            }
          }
          
          rli_df <- data.frame(
            Year = years,
            RLI = rli_values,
            Class = class_name,
            N_Especies = nrow(df_filtered),
            stringsAsFactors = FALSE
          )
          
          return(list(matrix = class_matrix, rli_data = rli_df, filtered_data = df_filtered))
        }
      }
    }
  }
  return(NULL)
}

# Función alternativa usando cálculo manual si red::rli falla
calculate_rli_manual <- function(df_filtered) {
  # Convertir categorías a pesos numéricos
  category_weights <- c(
    "CR" = 1.0,
    "EN" = 0.6667,
    "VU" = 0.3333,
    "NT" = 0.2,
    "LC" = 0,
    "DD" = NA,
    "NE" = NA,
    "RE" = NA
  )
  
  # Función para convertir categoría a peso
  category_to_weight <- function(cat) {
    ifelse(cat %in% names(category_weights), 
           category_weights[cat], 
           NA)
  }
  
  # Calcular RLI para cada año
  rli_values <- numeric(ncol(df_filtered))
  
  for (i in 1:ncol(df_filtered)) {
    year_data <- df_filtered[, i]
    
    # Convertir categorías a pesos
    weights <- sapply(year_data, function(x) {
      if (is.na(x)) return(NA)
      if (x %in% names(category_weights)) {
        return(category_weights[x])
      } else {
        return(NA)
      }
    })
    
    # Eliminar NAs
    valid_weights <- weights[!is.na(weights)]
    
    if (length(valid_weights) > 0) {
      # Calcular RLI: 1 - (suma de pesos / número de especies)
      rli_values[i] <- 1 - (sum(valid_weights) / length(valid_weights))
    } else {
      rli_values[i] <- NA
    }
  }
  
  return(rli_values)
}

# Versión simplificada de la función
process_class_matrix_simple <- function(class_matrix, class_name) {
  
  # Crear dataframe
  df <- as.data.frame.matrix(class_matrix)
  
  # Ordenar columnas por año
  year_cols <- colnames(df)
  numeric_years <- suppressWarnings(as.numeric(year_cols))
  
  if (all(!is.na(numeric_years))) {
    df <- df[, order(numeric_years)]
  }
  
  # Reemplazar NAs con valor anterior
  for (col_idx in 2:ncol(df)) {
    df[[col_idx]] <- ifelse(is.na(df[[col_idx]]), df[[col_idx-1]], df[[col_idx]])
  }
  
  # Buscar año 2000 o usar primer año disponible
  target_year <- "2000"
  if (target_year %in% colnames(df)) {
    base_year_index <- which(names(df) == target_year)
    df_clean <- df[!is.na(df[, base_year_index]), , drop = FALSE]
    
    if (base_year_index < ncol(df)) {
      df_filtered <- df_clean[, base_year_index:ncol(df_clean), drop = FALSE]
    } else {
      df_filtered <- df_clean
    }
  } else {
    # Usar primer año
    df_clean <- df[apply(df, 1, function(x) any(!is.na(x))), , drop = FALSE]
    df_filtered <- df_clean
  }
  
  if (nrow(df_filtered) == 0) {
    return(NULL)
  }
  
  # Calcular RLI (manualmente para mayor control)
  rli_values <- calculate_rli_manual(df_filtered)
  years <- colnames(df_filtered)
  
  rli_df <- data.frame(
    Year = years,
    RLI = rli_values,
    Class = class_name,
    N_Especies = nrow(df_filtered),
    stringsAsFactors = FALSE
  )
  
  return(list(matrix = class_matrix, rli_data = rli_df, filtered_data = df_filtered))
}

# Procesar cada clase con la versión simplificada
class_results <- list()
for (class_name in names(class_matrices)) {
  message("Procesando clase: ", class_name)
  
  # Primero intentar con la versión simplificada
  result <- process_class_matrix_simple(class_matrices[[class_name]], class_name)
  
  if (!is.null(result)) {
    class_results[[class_name]] <- result
    message("  ✓ Calculado RLI para ", nrow(result$rli_data), " años, ", 
            result$rli_data$N_Especies[1], " especies")
  } else {
    message("  ✗ No se pudo calcular RLI para ", class_name)
  }
}

# Si no hay resultados, intentar diagnóstico
if (length(class_results) == 0) {
  message("\nDiagnóstico del problema:")
  
  for (class_name in names(class_matrices)) {
    df <- as.data.frame.matrix(class_matrices[[class_name]])
    message("\nClase: ", class_name)
    message("  Dimensiones: ", nrow(df), " especies x ", ncol(df), " años")
    message("  Columnas (años): ", paste(colnames(df), collapse = ", "))
    
    # Verificar categorías presentes
    all_cats <- unique(as.vector(as.matrix(df)))
    message("  Categorías presentes: ", paste(all_cats[!is.na(all_cats)], collapse = ", "))
  }
}



# DEPURACIÓN: Ver qué está devolviendo red::rli para ACTINOPTERYGII
test_class <- "ACTINOPTERYGII"
test_matrix <- class_matrices[[test_class]]

if (!is.null(test_matrix)) {
  df <- as.data.frame.matrix(test_matrix)
  
  # Ordenar y limpiar como en la función
  year_cols <- colnames(df)
  numeric_years <- suppressWarnings(as.numeric(year_cols))
  if (all(!is.na(numeric_years))) {
    df <- df[, order(numeric_years)]
  }
  
  # Reemplazar NAs
  for (col_idx in 2:ncol(df)) {
    df[[col_idx]] <- ifelse(is.na(df[[col_idx]]), df[[col_idx-1]], df[[col_idx]])
  }
  
  # Filtrar año 2000
  if ("2000" %in% colnames(df)) {
    base_idx <- which(names(df) == "2000")
    df_clean <- df[!is.na(df[, base_idx]), , drop = FALSE]
    df_filtered <- df_clean[, base_idx:ncol(df_clean), drop = FALSE]
    
    message("\nDEPURACIÓN para ", test_class)
    message("df_filtered dimensiones: ", nrow(df_filtered), " x ", ncol(df_filtered))
    print(head(df_filtered[, 1:min(5, ncol(df_filtered))]))
    
    # Probar red::rli
    test_rli <- red::rli(df_filtered, boot = FALSE)
    message("\nResultado de red::rli:")
    print(str(test_rli))
    print(test_rli)
  }
}



# Probar con una sola clase
test_class <- "ACTINOPTERYGII"
if (test_class %in% names(class_matrices)) {
  test_result <- process_class_matrix_simple(class_matrices[[test_class]], test_class)
  
  if (!is.null(test_result)) {
    print(test_result$rli_data)
    
    # Gráfico rápido
    ggplot(test_result$rli_data, aes(x = Year, y = RLI, group = 1)) +
      geom_line(color = "blue") +
      geom_point() +
      labs(title = paste("RLI -", test_class)) +
      theme_minimal()
  }
}









#### Agrupar clases en categorías principales ####

# Primero, definimos un mapeo de clases a categorías principales
class_to_main_group <- list(
  # Anfibios
  "AMPHIBIA" = "Amphibians",
  
  # Aves
  "AVES" = "Birds",
  
  # Mamíferos
  "MAMMALIA" = "Mammals",
  
  # Reptiles (si tienes clase REPTILIA)
  "REPTILIA" = "Reptiles",
  
  # Peces
  "ACTINOPTERYGII" = "Fishes",
  "CHONDRICHTHYES" = "Fishes",
  
  # Plantas
  "MAGNOLIOPSIDA" = "Plants",
  "LILIOPSIDA" = "Plants", 
  "CYCADOPSIDA" = "Plants",
  "BRYOPSIDA" = "Plants",
  "JUNGERMANNIOPSIDA" = "Plants",
  "LYCOPODIOPSIDA" = "Plants",
  "SPHAGNOPSIDA" = "Plants",
  "ANTHOCEROTOPSIDA" = "Plants",
  
  # Invertebrados
  "ARTHONIOMYCETES" = "Invertebrates",
  "LECANOROMYCETES" = "Invertebrates",
  "AGARICOMYCETES" = "Invertebrates",
  "EUROTIOMYCETES" = "Invertebrates",
  
  # Otras clases se asignan a "Others" o se pueden incluir según necesidad
  "NOT ASSIGNED" = "Others"
)

# Función para agrupar matrices por categoría principal
group_matrices_by_main_category <- function(class_matrices, class_to_main_group) {
  main_group_matrices <- list()
  
  for (class_name in names(class_matrices)) {
    # Determinar la categoría principal
    if (class_name %in% names(class_to_main_group)) {
      main_group <- class_to_main_group[[class_name]]
    } else {
      main_group <- "Others"  # Categoría por defecto
    }
    
    # Agregar la matriz a la categoría principal
    if (main_group %in% names(main_group_matrices)) {
      # Si ya existe la categoría, combinar las matrices
      main_group_matrices[[main_group]] <- rbind(
        main_group_matrices[[main_group]],
        class_matrices[[class_name]]
      )
    } else {
      # Si no existe, crear nueva entrada
      main_group_matrices[[main_group]] <- class_matrices[[class_name]]
    }
  }
  
  return(main_group_matrices)
}

# Agrupar las matrices en categorías principales
main_group_matrices <- group_matrices_by_main_category(class_matrices, class_to_main_group)

# Mostrar las categorías creadas
message("Categorías principales creadas:")
print(names(main_group_matrices))
message("Número de especies por categoría:")
for (group_name in names(main_group_matrices)) {
  message("  ", group_name, ": ", nrow(main_group_matrices[[group_name]]), " especies")
}

#### Calcular RLI para cada categoría principal ####

main_group_results <- list()
for (group_name in names(main_group_matrices)) {
  message("\nProcesando categoría: ", group_name)
  
  result <- process_class_matrix_simple(main_group_matrices[[group_name]], group_name)
  
  if (!is.null(result)) {
    main_group_results[[group_name]] <- result
    message("  ✓ Calculado RLI para ", nrow(result$rli_data), " años, ", 
            result$rli_data$N_Especies[1], " especies")
  } else {
    message("  ✗ No se pudo calcular RLI para ", group_name)
  }
}

#### También calcular RLI para TODAS las especies juntas ####

# Combinar todas las matrices
all_species_matrix <- do.call(rbind, class_matrices)
all_result <- process_class_matrix_simple(all_species_matrix, "All species")

if (!is.null(all_result)) {
  main_group_results[["All species"]] <- all_result
  message("\n✓ Calculado RLI para todas las especies: ", 
          all_result$rli_data$N_Especies[1], " especies")
}

#### Crear una sola gráfica con todas las categorías ####

# Combinar todos los datos de RLI
all_rli_data_main <- data.frame()
for (group_name in names(main_group_results)) {
  group_rli <- main_group_results[[group_name]]$rli_data
  all_rli_data_main <- rbind(all_rli_data_main, group_rli)
}

# Verificar datos
message("\nDatos disponibles para gráfica:")
print(unique(all_rli_data_main$Class))

# Definir colores para cada categoría
group_colors <- c(
  "Amphibians" = "#FF6B6B",
  "Birds" = "#4ECDC4",
  "Mammals" = "#45B7D1",
  "Reptiles" = "#96CEB4",
  "Fishes" = "#FECA57",
  "Plants" = "#58B368",
  "Invertebrates" = "#9B59B6",
  "All species" = "#2C3E50",
  "Others" = "#95A5A6"
)

# Ordenar las categorías para la leyenda
ordered_groups <- c("All species", "Amphibians", "Birds", "Mammals", "Reptiles", 
                    "Fishes", "Plants", "Invertebrates", "Others")
all_rli_data_main$Class <- factor(all_rli_data_main$Class, levels = ordered_groups)

# Convertir Year a factor ordenado
all_rli_data_main$Year <- factor(
  all_rli_data_main$Year,
  levels = sort(unique(as.numeric(as.character(all_rli_data_main$Year))))
)

#### Crear la gráfica combinada ####
library(ggplot2)

combined_plot <- ggplot(all_rli_data_main, aes(x = Year, y = RLI, color = Class, group = Class)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  labs(
    title = "Red List Index by Major Taxonomic Groups",
    subtitle = "Colombian Endemic Species",
    x = "Year",
    y = "Red List Index (RLI)",
    color = "Taxonomic Group"
  ) +
  scale_color_manual(values = group_colors, breaks = ordered_groups) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# Mostrar la gráfica
print(combined_plot)

#### Guardar resultados ####
if (!dir.exists("result")) dir.create("result")

# Guardar la gráfica principal
ggsave("result/Combined_RLI_Plot.png", combined_plot, 
       width = 14, height = 8, dpi = 300)

# Guardar datos de RLI por categoría principal
write.csv(all_rli_data_main, "result/RLI_by_Main_Groups.csv", row.names = FALSE)

# Guardar resumen estadístico
summary_main <- all_rli_data_main %>%
  group_by(Class) %>%
  summarise(
    N_Species = first(N_Especies),
    First_Year = first(Year),
    First_RLI = round(first(RLI), 4),
    Last_Year = last(Year),
    Last_RLI = round(last(RLI), 4),
    Change = round(last(RLI) - first(RLI), 4),
    Percent_Change = round((last(RLI) - first(RLI)) / first(RLI) * 100, 2)
  ) %>%
  arrange(factor(Class, levels = ordered_groups))

print("\nResumen por categorías principales:")
print(summary_main)
write.csv(summary_main, "result/Summary_Main_Groups.csv", row.names = FALSE)

#### Crear gráfica alternativa: solo grupos principales (sin "Others") ####

# Filtrar para excluir "Others" si tiene pocos datos
if ("Others" %in% all_rli_data_main$Class) {
  main_groups_plot_data <- all_rli_data_main %>%
    filter(Class != "Others")
} else {
  main_groups_plot_data <- all_rli_data_main
}

# Crear gráfica sin "Others"
main_groups_only_plot <- ggplot(main_groups_plot_data, 
                                aes(x = Year, y = RLI, color = Class, group = Class)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  labs(
    title = "Red List Index by Major Taxonomic Groups",
    subtitle = "Colombian Endemic Species (Main Groups Only)",
    x = "Year",
    y = "Red List Index (RLI)",
    color = "Taxonomic Group"
  ) +
  scale_color_manual(values = group_colors[names(group_colors) != "Others"]) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(main_groups_only_plot)
ggsave("result/Main_Groups_Only_RLI_Plot.png", main_groups_only_plot, 
       width = 14, height = 8, dpi = 300)

#### Crear gráfica de área apilada (alternativa visual) ####

# Para gráfica de área, necesitamos pivotar los datos
area_plot_data <- all_rli_data_main %>%
  filter(Class != "All species") %>%  # Excluir "All species" para área apilada
  group_by(Year) %>%
  mutate(
    Proportion = N_Especies / sum(N_Especies),
    Weighted_RLI = RLI * Proportion
  ) %>%
  ungroup()

stacked_area_plot <- ggplot(area_plot_data, aes(x = Year, y = Weighted_RLI, fill = Class)) +
  geom_area(alpha = 0.7) +
  labs(
    title = "Weighted Red List Index by Taxonomic Group",
    subtitle = "Area proportional to number of species in each group",
    x = "Year",
    y = "Weighted RLI",
    fill = "Taxonomic Group"
  ) +
  scale_fill_manual(values = group_colors[!names(group_colors) %in% c("All species", "Others")]) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(stacked_area_plot)
ggsave("result/Stacked_Area_RLI_Plot.png", stacked_area_plot, 
       width = 12, height = 7, dpi = 300)

#### Informe final ####
cat("\n" , rep("=", 60), "\n", sep = "")
cat("RESUMEN EJECUTIVO\n")
cat(rep("=", 60), "\n")
cat("Total de especies analizadas:", all_result$rli_data$N_Especies[1], "\n")
cat("Período analizado:", min(as.numeric(as.character(all_rli_data_main$Year))), 
    "-", max(as.numeric(as.character(all_rli_data_main$Year))), "\n")
cat("Número de grupos taxonómicos:", length(unique(all_rli_data_main$Class)), "\n")
cat("\nCambio en RLI para todos las especies:\n")
cat("  Inicio (", all_result$rli_data$Year[1], "): ", 
    round(all_result$rli_data$RLI[1], 4), "\n", sep = "")
cat("  Final (", all_result$rli_data$Year[nrow(all_result$rli_data)], "): ", 
    round(all_result$rli_data$RLI[nrow(all_result$rli_data)], 4), "\n", sep = "")
cat("  Cambio total: ", 
    round(all_result$rli_data$RLI[nrow(all_result$rli_data)] - all_result$rli_data$RLI[1], 4), 
    "\n", sep = "")
cat(rep("=", 60), "\n")



# version -----------------------------------------------------------------


adjust_categories <- data.frame(
  Cat_IUCN = c("CR","EN","VU","NT","LC","DD","RE","EW","EX"),
  code = c("CR","EN","VU","NT","LC","DD","RE","EW","EX")
)


calculate_rli_iucn <- function(df) {
  
  weights <- c(
    "LC" = 0,
    "NT" = 1,
    "VU" = 2,
    "EN" = 3,
    "CR" = 4,
    "EW" = 5,
    "EX" = 5
  )
  
  # Mantener solo especies con al menos una evaluación válida
  valid_species <- apply(df, 1, function(x) any(x %in% names(weights)))
  df <- df[valid_species, , drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  N <- nrow(df)
  max_score <- 5 * N
  
  rli <- numeric(ncol(df))
  
  for (i in seq_len(ncol(df))) {
    cats <- df[, i]
    score <- sum(weights[cats], na.rm = TRUE)
    rli[i] <- 1 - (score / max_score)
  }
  
  return(rli)
}

process_class_matrix <- function(class_matrix, class_name) {
  
  df <- as.data.frame.matrix(class_matrix)
  
  # Ordenar columnas por año
  years <- suppressWarnings(as.numeric(colnames(df)))
  if (all(!is.na(years))) {
    df <- df[, order(years), drop = FALSE]
  }
  
  # Usar año base 2000 si existe
  if ("2000" %in% colnames(df)) {
    base <- which(colnames(df) == "2000")
    df <- df[, base:ncol(df), drop = FALSE]
  }
  
  # Calcular RLI
  rli_values <- calculate_rli_iucn(df)
  if (is.null(rli_values)) return(NULL)
  
  data.frame(
    Year = colnames(df),
    RLI = rli_values,
    Class = class_name,
    N_Especies = nrow(df),
    stringsAsFactors = FALSE
  )
}

class_rli_results <- list()

for (class_name in names(class_matrices)) {
  message("Procesando clase: ", class_name)
  
  rli_df <- process_class_matrix(class_matrices[[class_name]], class_name)
  
  if (!is.null(rli_df)) {
    class_rli_results[[class_name]] <- rli_df
    message("  ✓ ", nrow(rli_df), " años, ", rli_df$N_Especies[1], " especies")
  } else {
    message("  ✗ Sin datos válidos")
  }
}

all_rli_data_main <- do.call(rbind, class_rli_results)

summary(all_rli_data_main$RLI)



calculate_rli_iucn <- function(df) {
  
  weights <- c(
    "LC" = 0,
    "NT" = 1,
    "VU" = 2,
    "EN" = 3,
    "CR" = 4,
    "EW" = 5,
    "EX" = 5
  )
  
  # Mantener solo especies con AL MENOS UNA evaluación válida
  valid_species <- apply(df, 1, function(x) any(x %in% names(weights)))
  df <- df[valid_species, , drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  N <- nrow(df)
  max_score <- 5 * N
  
  rli <- numeric(ncol(df))
  
  for (i in seq_len(ncol(df))) {
    
    cats <- df[, i]
    
    # NO imputar categorías
    valid <- cats %in% names(weights)
    
    if (sum(valid) == 0) {
      rli[i] <- NA
    } else {
      score <- sum(weights[cats[valid]])
      rli[i] <- 1 - (score / max_score)
    }
  }
  
  return(rli)
}


process_class_matrix <- function(class_matrix, class_name) {
  
  df <- as.data.frame.matrix(class_matrix)
  
  # Ordenar años
  years <- suppressWarnings(as.numeric(colnames(df)))
  if (all(!is.na(years))) {
    df <- df[, order(years), drop = FALSE]
  }
  
  # Usar año base SOLO para seleccionar especies
  if ("2000" %in% colnames(df)) {
    base <- which(colnames(df) == "2000")
    keep <- !is.na(df[, base])
    df <- df[keep, , drop = FALSE]
    df <- df[, base:ncol(df), drop = FALSE]
  }
  
  rli_values <- calculate_rli_iucn(df)
  if (is.null(rli_values)) return(NULL)
  
  data.frame(
    Year = colnames(df),
    RLI = rli_values,
    Class = class_name,
    N_Especies = nrow(df),
    stringsAsFactors = FALSE
  )
}

table(class_matrices[["ACTINOPTERYGII"]], useNA = "ifany")




# version 3 ---------------------------------------------------------------

library(data.table)

# Cargar el archivo (más rápido que read.csv)
df <- fread("IUCN_history_DF_with_class2.csv", encoding = "UTF-8")

# Ver estructura
str(df)
head(df)
dim(df)

# Ver columnas disponibles
colnames(df)


calculate_rli_iucn <- function(df) {
  
  weights <- c(
    "LC" = 0,
    "NT" = 1,
    "VU" = 2,
    "EN" = 3,
    "CR" = 4,
    "EW" = 5,
    "EX" = 5
  )
  
  # Mantener solo especies con categorías válidas
  # Esto filtra especies que tienen AL MENOS UNA categoría válida en cualquier año
  valid_species <- apply(df, 1, function(x) any(x %in% names(weights)))
  df <- df[valid_species, , drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  N <- nrow(df)
  # El máximo teórico es cuando todas las especies están en EX/EW
  max_theoretical_score <- 5 * N
  
  rli <- numeric(ncol(df))
  
  for (i in seq_len(ncol(df))) {
    cats <- df[, i]
    
    # Identificar categorías válidas
    valid <- cats %in% names(weights)
    
    if (sum(valid) == 0) {
      rli[i] <- NA
    } else {
      # Solo sumar especies con categorías válidas en este año
      score <- sum(weights[cats[valid]])
      # Usar N_actual = número de especies con categoría válida en este año
      N_actual <- sum(valid)
      max_possible_score <- 5 * N_actual
      
      if (max_possible_score == 0) {
        rli[i] <- NA
      } else {
        rli[i] <- 1 - (score / max_possible_score)
      }
    }
  }
  
  return(rli)
}

# Versión alternativa si quieres usar el N constante (todas las especies)
calculate_rli_iucn_constant_N <- function(df) {
  
  weights <- c(
    "LC" = 0,
    "NT" = 1,
    "VU" = 2,
    "EN" = 3,
    "CR" = 4,
    "EW" = 5,
    "EX" = 5
  )
  
  # Filtrar especies: deben tener categoría válida en el PRIMER año (2000)
  if ("2000" %in% colnames(df)) {
    base_year <- "2000"
  } else {
    # Tomar el primer año disponible
    base_year <- colnames(df)[1]
  }
  
  keep <- df[, base_year] %in% names(weights)
  df <- df[keep, , drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  N <- nrow(df)
  max_score <- 5 * N
  
  rli <- numeric(ncol(df))
  
  for (i in seq_len(ncol(df))) {
    cats <- df[, i]
    
    # Para especies sin categoría en un año, usar la última disponible
    # Esto es importante para series temporales
    if (i > 1) {
      prev_cats <- df[, i-1]
      cats[!cats %in% names(weights)] <- prev_cats[!cats %in% names(weights)]
    }
    
    # Solo especies con categoría válida
    valid <- cats %in% names(weights)
    
    if (sum(valid) == 0) {
      rli[i] <- NA
    } else {
      score <- sum(weights[cats[valid]])
      # Usar N constante (todas las especies incluidas)
      rli[i] <- 1 - (score / max_score)
    }
  }
  
  return(rli)
}

# Para verificar la tabla de categorías
table_categories <- function(class_matrix, class_name) {
  df <- as.data.frame.matrix(class_matrices[[class_name]])
  
  cat("=== Clase:", class_name, "===\n")
  cat("Número total de especies:", nrow(df), "\n")
  cat("Número de años:", ncol(df), "\n\n")
  
  # Contar categorías por año
  for (year in colnames(df)) {
    cat("Año", year, ":\n")
    print(table(df[, year], useNA = "ifany"))
    cat("\n")
  }
  
  # Especies que no tienen categoría en 2000
  if ("2000" %in% colnames(df)) {
    missing_2000 <- sum(!df[, "2000"] %in% c("LC", "NT", "VU", "EN", "CR", "EW", "EX"))
    cat("Especies sin categoría válida en 2000:", missing_2000, "\n")
  }
}

# Ejecutar para verificar
table_categories(class_matrices[["ACTINOPTERYGII"]], "ACTINOPTERYGII")








# Verificar si tenemos datos en formato ancho o largo
has_year_cols <- any(grepl("^\\d{4}$", colnames(df)))

if (has_year_cols) {
  message("\n✅ Formato ANCHO detectado (una columna por año)")
  
  # Encontrar columnas que son años (solo dígitos)
  year_cols <- grep("^\\d{4}$", colnames(df), value = TRUE)
  message("Años encontrados: ", length(year_cols), " (de ", min(year_cols), " a ", max(year_cols), ")")
  
  # Identificar columnas de información (no años)
  info_cols <- setdiff(colnames(df), year_cols)
  message("Columnas de información: ", paste(info_cols, collapse = ", "))
  
  # Verificar si tenemos columna de clase
  if (any(c("class", "Class") %in% info_cols)) {
    class_col <- ifelse("class" %in% info_cols, "class", "Class")
    message("Columna de clase encontrada: '", class_col, "'")
    
    # Obtener clases únicas
    classes <- unique(df[[class_col]])
    message("Clases encontradas (", length(classes), "):")
    print(classes)
    
    # Crear lista de matrices por clase
    class_matrices <- list()
    
    for (cl in classes) {
      message("\nProcesando clase: ", cl)
      
      # Filtrar por clase
      class_data <- df[df[[class_col]] == cl, ]
      message("  Especies en esta clase: ", nrow(class_data))
      
      # Crear matriz (especies × años)
      # Usar primera columna de información como nombre de especie
      species_col <- info_cols[1]  # Ajusta según tu archivo
      rownames_vector <- make.names(class_data[[species_col]], unique = TRUE)
      
      # Crear matriz solo con columnas de años
      year_matrix <- as.matrix(class_data[, ..year_cols])
      rownames(year_matrix) <- rownames_vector
      
      # Almacenar en la lista
      class_matrices[[cl]] <- year_matrix
      
      message("  Matriz creada: ", nrow(year_matrix), " × ", ncol(year_matrix))
    }
    
  } else {
    message("⚠️ No se encontró columna de clase. Creando una sola matriz...")
    
    # Crear una sola matriz
    species_col <- info_cols[1]
    rownames_vector <- make.names(df[[species_col]], unique = TRUE)
    class_matrices <- list()
    class_matrices[["ALL"]] <- as.matrix(df[, ..year_cols])
    rownames(class_matrices[["ALL"]]) <- rownames_vector
  }
  
} else {
  # Formato largo (melted)
  message("\n📊 Formato LARGO detectado (species, year, category)")
  
  # Verificar columnas mínimas necesarias
  required_cols <- c("species", "year", "category", "class")
  if (all(required_cols %in% colnames(df))) {
    # Convertir a formato ancho
    message("Convirtiendo a formato ancho...")
    
    # Limpiar datos
    df_clean <- df %>%
      select(species, class, year, category) %>%
      distinct() %>%
      arrange(species, year)
    
    # Convertir a formato ancho
    df_wide <- df_clean %>%
      pivot_wider(
        id_cols = c(species, class),
        names_from = year,
        values_from = category,
        values_fill = NA
      )
    
    # Crear matrices por clase
    classes <- unique(df_wide$class)
    message("Clases encontradas (", length(classes), "):")
    print(classes)
    
    class_matrices <- list()
    year_cols <- setdiff(colnames(df_wide), c("species", "class"))
    
    for (cl in classes) {
      message("\nProcesando clase: ", cl)
      
      class_data <- df_wide %>%
        filter(class == cl)
      
      # Crear matriz
      species_names <- class_data$species
      year_data <- as.matrix(class_data[, year_cols])
      rownames(year_data) <- species_names
      
      class_matrices[[cl]] <- year_data
      message("  Especies: ", nrow(year_data))
      message("  Años: ", ncol(year_data))
    }
    
  } else {
    message("⚠️ Columnas requeridas no encontradas. Buscando alternativas...")
    
    # Intentar encontrar nombres alternativos
    species_col <- grep("species|name|scientific", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    year_col <- grep("year|date|assessment", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    category_col <- grep("category|status|red.list", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    class_col <- grep("class|group|taxon", colnames(df), value = TRUE, ignore.case = TRUE)[1]
    
    message("Columna de especies: ", species_col)
    message("Columna de año: ", year_col)
    message("Columna de categoría: ", category_col)
    message("Columna de clase: ", class_col)
  }
}

# Verificar que class_matrices existe
if (exists("class_matrices")) {
  message("\n✅ class_matrices creado exitosamente")
  message("Número de clases: ", length(class_matrices))
  message("Nombres de las clases: ", paste(names(class_matrices), collapse = ", "))
  
  # Ver detalles de cada matriz
  for (cl in names(class_matrices)) {
    mat <- class_matrices[[cl]]
    message("\nClase: ", cl)
    message("  Dimensiones: ", nrow(mat), " × ", ncol(mat))
    message("  Primeras columnas: ", paste(colnames(mat)[1:5], collapse = ", "))
    
    # Ver categorías únicas
    unique_cats <- unique(as.vector(mat))
    message("  Categorías únicas: ", paste(unique_cats[!is.na(unique_cats)][1:10], collapse = ", "))
  }
  
  # Ahora puedes ejecutar tu función
  source("table_categories.R")  # Si tienes la función en otro archivo
  # O ejecutar directamente:
  table_categories <- function(class_matrix, class_name) {
    df <- as.data.frame.matrix(class_matrix)
    
    cat("=== Clase:", class_name, "===\n")
    cat("Número total de especies:", nrow(df), "\n")
    cat("Número de años:", ncol(df), "\n\n")
    
    # Contar categorías por año (solo primeros 5 años para no saturar)
    for (year in colnames(df)[1:min(5, ncol(df))]) {
      cat("Año", year, ":\n")
      print(table(df[, year], useNA = "ifany"))
      cat("\n")
    }
    
    # Especies que no tienen categoría en 2000
    if ("2000" %in% colnames(df)) {
      valid_cats <- c("LC", "NT", "VU", "EN", "CR", "EW", "EX")
      missing_2000 <- sum(!df[, "2000"] %in% valid_cats & !is.na(df[, "2000"]))
      cat("Especies sin categoría válida en 2000:", missing_2000, "\n")
    }
  }
  
  # Ejecutar para la primera clase
  primera_clase <- names(class_matrices)[1]
  message("\n=== EJECUTANDO VERIFICACIÓN PARA: ", primera_clase, " ===")
  table_categories(class_matrices[[primera_clase]], primera_clase)
  
} else {
  message("\n❌ class_matrices no fue creado. Revisa la estructura de tus datos.")
  
  # Mostrar más detalles del dataframe
  message("\n=== RESUMEN DEL DATAFRAME ===")
  print(summary(df))
  
  message("\n=== MUESTRA DE DATOS ===")
  print(df[1:20, ])
}

# Ejecuta esto y comparte el resultado:
message("=== INFORMACIÓN CRÍTICA ===")
message("Nombres de columnas: ", paste(colnames(df), collapse = ", "))
message("Primeros valores de la primera fila:")
print(as.list(df[1, ]))
message("Estructura:")
str(df, list.len = 10)













# Cargar paquetes necesarios
library(data.table)
library(dplyr)
library(tidyr)

# Ver todas las columnas
message("=== COLUMNAS COMPLETAS ===")
print(colnames(df))

# Ver estructura completa
message("\n=== MUESTRA DE DATOS ===")
print(head(df, 20))

# Verificar si tenemos todas las columnas necesarias
message("\n=== BUSCANDO COLUMNAS CLAVE ===")

# Buscar columna de clase/taxonomía
class_col <- ifelse("class" %in% colnames(df), "class",
                    ifelse("Class" %in% colnames(df), "Class",
                           ifelse("taxon_class_name" %in% colnames(df), "taxon_class_name", NA)))

# Buscar columna de año
year_col <- ifelse("year" %in% colnames(df), "year",
                   ifelse("year_published" %in% colnames(df), "year_published",
                          ifelse("assessment_year" %in% colnames(df), "assessment_year", NA)))

# Buscar columna de categoría
cat_col <- ifelse("red_list_category_code" %in% colnames(df), "red_list_category_code",
                  ifelse("category" %in% colnames(df), "category",
                         ifelse("redlist_category" %in% colnames(df), "redlist_category", NA)))

# Buscar columna de nombre científico
species_col <- ifelse("taxon_scientific_name" %in% colnames(df), "taxon_scientific_name",
                      ifelse("species" %in% colnames(df), "species",
                             ifelse("scientific_name" %in% colnames(df), "scientific_name", NA)))

message("Columna de clase encontrada: ", class_col)
message("Columna de año encontrada: ", year_col)
message("Columna de categoría encontrada: ", cat_col)
message("Columna de especie encontrada: ", species_col)

# Transformar los datos al formato necesario
message("\n=== TRANSFORMANDO DATOS ===")

# 1. Seleccionar solo las columnas necesarias
df_clean <- df %>%
  select(
    species = !!species_col,
    class = !!class_col,
    year = !!year_col,
    category = !!cat_col
  ) %>%
  # Eliminar duplicados (si los hay)
  distinct() %>%
  # Ordenar
  arrange(species, year)

# Verificar los datos limpios
message("Datos limpios - Dimensiones: ", nrow(df_clean), " × ", ncol(df_clean))
print(head(df_clean, 20))

# 2. Verificar valores únicos
message("\n=== VALORES ÚNICOS ===")
message("Número de especies únicas: ", length(unique(df_clean$species)))
message("Clases únicas: ", paste(unique(df_clean$class), collapse = ", "))
message("Años únicos: ", paste(sort(unique(df_clean$year)), collapse = ", "))
message("Categorías únicas: ", paste(unique(df_clean$category), collapse = ", "))

# 3. Convertir a formato ancho (especies × años)
message("\n=== CREANDO FORMATO ANCHO ===")

df_wide <- df_clean %>%
  # Filtrar solo categorías válidas (opcional, puedes comentar esto)
  filter(category %in% c("LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD", "NE")) %>%
  # Pivotar a formato ancho
  pivot_wider(
    id_cols = c(species, class),
    names_from = year,
    values_from = category,
    values_fill = NA
  )

message("Formato ancho creado: ", nrow(df_wide), " × ", ncol(df_wide))

# 4. Crear matrices por clase
message("\n=== CREANDO MATRICES POR CLASE ===")

# Obtener años únicos (columnas numéricas)
year_cols <- setdiff(colnames(df_wide), c("species", "class"))
year_cols <- year_cols[grepl("^\\d{4}$", year_cols)]  # Solo columnas que son años
year_cols <- sort(as.numeric(year_cols))

# Mantener solo los años como caracteres
year_cols <- as.character(year_cols)

message("Años disponibles: ", length(year_cols), " (", min(year_cols), " - ", max(year_cols), ")")

# Obtener clases únicas
classes <- unique(df_wide$class)
message("Clases encontradas (", length(classes), "):")
print(classes)

# Crear la lista de matrices
class_matrices <- list()

for (cl in classes) {
  message("\nProcesando clase: ", cl)
  
  # Filtrar datos de esta clase
  class_data <- df_wide %>%
    filter(class == cl)
  
  message("  Número de especies: ", nrow(class_data))
  
  # Crear matriz
  species_names <- class_data$species
  year_data <- as.matrix(class_data[, year_cols])
  rownames(year_data) <- species_names
  
  # Almacenar en la lista
  class_matrices[[cl]] <- year_data
  
  message("  Matriz creada: ", nrow(year_data), " × ", ncol(year_data))
  
  # Verificar categorías únicas
  unique_cats <- unique(as.vector(year_data))
  unique_cats <- unique_cats[!is.na(unique_cats)]
  message("  Categorías únicas (", length(unique_cats), "): ", 
          paste(unique_cats[1:min(10, length(unique_cats))], collapse = ", "))
}

# 5. Verificar que class_matrices fue creado
message("\n=== VERIFICACIÓN FINAL ===")
message("Número de matrices creadas: ", length(class_matrices))
message("Clases en class_matrices: ", paste(names(class_matrices), collapse = ", "))

# 6. Función para verificar (adaptada a tu estructura)
table_categories <- function(class_matrix, class_name) {
  df_temp <- as.data.frame.matrix(class_matrix)
  
  cat("\n=== Clase:", class_name, "===\n")
  cat("Número total de especies:", nrow(df_temp), "\n")
  cat("Número de años:", ncol(df_temp), "\n")
  
  # Ver años disponibles
  cat("Años disponibles:", paste(colnames(df_temp), collapse = ", "), "\n\n")
  
  # Verificar si tenemos el año 2000
  if ("2000" %in% colnames(df_temp)) {
    cat("=== AÑO 2000 ===\n")
    tab_2000 <- table(df_temp[, "2000"], useNA = "ifany")
    print(tab_2000)
    cat("Total especies con datos en 2000:", sum(tab_2000), "\n\n")
  } else {
    cat("⚠️ Año 2000 no encontrado en los datos\n\n")
  }
  
  # Ver últimos años disponibles
  last_years <- tail(colnames(df_temp), 3)
  for (year in last_years) {
    cat("=== AÑO", year, "===\n")
    tab_year <- table(df_temp[, year], useNA = "ifany")
    print(tab_year)
    cat("Total especies con datos en", year, ":", sum(tab_year), "\n\n")
  }
  
  # Calcular porcentaje de NAs por año
  na_percent <- apply(df_temp, 2, function(x) round(mean(is.na(x)) * 100, 1))
  cat("=== PORCENTAJE DE DATOS FALTANTES POR AÑO ===\n")
  print(na_percent)
}

# 7. Probar con la primera clase
if (length(class_matrices) > 0) {
  primera_clase <- names(class_matrices)[1]
  message("\n=== PROBANDO CON LA CLASE: ", primera_clase, " ===")
  table_categories(class_matrices[[primera_clase]], primera_clase)
  
  # Guardar class_matrices para uso futuro
  saveRDS(class_matrices, "class_matrices.rds")
  message("\n✅ class_matrices guardado en 'class_matrices.rds'")
} else {
  message("❌ No se pudieron crear las matrices. Revisa los datos.")
}
