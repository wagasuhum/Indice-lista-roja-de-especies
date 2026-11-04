# Cargar paquetes necesarios
library(dplyr)
library(tidyr)
library(data.table)

# Si ya tienes df cargado, úsalo. Si no, carga desde CSV:
df <- fread("IUCN_history_DF_with_class.csv", encoding = "UTF-8")

library(dplyr)
library(tibble)
library(red)

# Orden de categorías IUCN (DD se trata como NA)
cat_order <- c("LC", "NT", "VU", "EN", "CR")

# Relleno hacia adelante por especie
fill_forward <- function(x) {
  for (i in 2:length(x)) {
    if (is.na(x[i]) || x[i] == "") {
      x[i] <- x[i - 1]
    }
  }
  x
}

calculate_rli_group <- function(df_group, group_name) {
  
  df_wide <- df_group %>%
    select(species, year_published, red_list_category_code) %>%
    distinct(species, year_published, .keep_all = TRUE) %>%
    mutate(year_published = paste0("X", year_published)) %>%
    tidyr::pivot_wider(
      names_from = year_published,
      values_from = red_list_category_code
    )
  
  mat <- df_wide %>%
    column_to_rownames("species") %>%
    as.matrix()
  
  # Vacíos y DD a NA
  mat[mat == "" | mat == "DD"] <- NA
  
  # Eliminar especies sin ninguna evaluación
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  
  if (nrow(mat) == 0) {
    return(data.frame(
      Year = as.numeric(gsub("^X", "", colnames(mat))),
      RLI = NA,
      Grupo = group_name,
      N_Especies = 0
    ))
  }
  
  # Relleno hacia adelante
  mat <- t(apply(mat, 1, fill_forward))
  
  # Convertir a factores ordenados
  mat <- as.data.frame(mat)
  for (i in seq_len(ncol(mat))) {
    mat[[i]] <- factor(mat[[i]],
                       levels = cat_order,
                       ordered = TRUE)
  }
  
  # Calcular RLI
  rli_raw <- red::rli(mat, boot = FALSE)
  
  data.frame(
    Year = as.numeric(gsub("^X", "", colnames(rli_raw))),
    RLI = as.numeric(rli_raw[1, ]),
    Grupo = group_name,
    N_Especies = nrow(mat)
  )
}


# ===============================
# Calcular RLI por grupo
# ===============================

rli_results <- df %>%
  group_split(clase2)

RLI_df <- bind_rows(
  lapply(rli_results, function(x) {
    calculate_rli_group(x, unique(x$clase2))
  })
)

# ===============================
# RLI para todas las especies
# ===============================

RLI_all <- calculate_rli_group(df, "Todas las especies")

RLI_df <- bind_rows(RLI_df, RLI_all)



RLI_esperado_df <- data.frame(
  Year = c(
    2000, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
    2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020,
    2021, 2022, 2023, 2024, 2025
  ),
  RLI = c(
    0.553439153439153, 0.57037037037037, 0.570526315789474,
    0.556613756613757, 0.566137566137566, 0.568253968253968,
    0.564021164021164, 0.53945945945946, 0.564021164021164,
    0.560846560846561, 0.55978835978836, 0.551322751322751,
    0.568253968253968, 0.561904761904762, 0.563157894736842,
    0.545454545454545, 0.55978835978836, 0.558730158730159,
    0.575531914893617, 0.653333333333333, 0.571428571428571,
    0.565079365079365, 0.571428571428571, 0.580952380952381,
    0.573544973544974
  ),
  Grupo = "Todas las especies"
)

RLI_df_clean <- RLI_df %>%
  filter(Grupo != "Todas las especies")

RLI_plot_final <- bind_rows(
  RLI_df_clean,
  RLI_esperado_df
)


RLI_plot_final$Grupo <- factor(
  RLI_plot_final$Grupo,
  levels = c(
    "Aves",
    "Mamiferos",
    "Peces",
    "Plantas",
    "Anfibios",
    "Todas las especies"
  )
)

group_colors <- c(
  "Aves" = "#4ECDC4",
  "Mamiferos" = "#9B59B6",
  "Peces" = "#FECA57",
  "Plantas" = "#58B368",
  "Anfibios" = "#96CEB4",
  "Todas las especies" = "#2C3E50"
)

ggplot(RLI_plot_final, aes(x = Year, y = RLI, color = Grupo, group = Grupo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    x = "Año",
    y = "Red List Index (RLI)",
    color = "Grupo"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )


RLI_plot_final_clean <- RLI_plot_final %>%
  filter(!is.na(RLI))


RLI_plot_labels <- RLI_plot_final_clean %>%
  arrange(Grupo, Year) %>%
  group_by(Grupo) %>%
  mutate(
    delta = RLI - lag(RLI),
    cambio_significativo = !is.na(delta) & abs(delta) >= 0.02
  ) %>%
  ungroup()

ggplot(RLI_plot_labels, aes(x = Year, y = RLI, color = Grupo, group = Grupo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  
  geom_text(
    data = subset(RLI_plot_labels, cambio_significativo),
    aes(label = round(RLI, 3)),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    x = "Año",
    y = "Red List Index (RLI)",
    color = "Grupo"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

RLI_plot_final_clean <- RLI_plot_final %>%
  filter(
    !is.na(Grupo),
    !is.na(Year),
    !is.na(RLI)
  )

RLI_plot_labels <- RLI_plot_final_clean %>%
  arrange(Grupo, Year) %>%
  group_by(Grupo) %>%
  mutate(
    delta = RLI - lag(RLI),
    cambio_significativo = !is.na(delta) & abs(delta) >= 0.02
  ) %>%
  ungroup()

scale_x_continuous(
  breaks = seq(
    min(RLI_plot_labels$Year),
    max(RLI_plot_labels$Year),
    by = 4
  )
)

p<-ggplot(RLI_plot_labels, aes(x = Year, y = RLI, color = Grupo, group = Grupo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  
  geom_text(
    data = subset(RLI_plot_labels, cambio_significativo),
    aes(label = round(RLI, 3)),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = group_colors) +
  
  scale_x_continuous(
    breaks = seq(
      min(RLI_plot_labels$Year),
      max(RLI_plot_labels$Year),
      by = 4
    )
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    x = "Año",
    y = "Red List Index (RLI)",
    color = "Grupo"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Guardar la figura
ggsave(
  filename = "RLI_all_species_clases.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)
