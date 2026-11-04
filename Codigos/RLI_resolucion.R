library(red)
library(dplyr)
library(tibble)
library(ggplot2)


df <- read.csv("C:/Users/walter.garcia/Documents/RLI/data/resolucion.csv", stringsAsFactors = FALSE, sep = ";")


# Species como nombres de fila
df_mat <- df %>%
  column_to_rownames("Species") %>%
  as.matrix()

fill_forward <- function(x) {
  for (i in 2:length(x)) {
    if (is.na(x[i]) || x[i] == "") {
      x[i] <- x[i - 1]
    }
  }
  x
}

df_filled <- t(apply(df_mat, 1, fill_forward))

df_clean <- df_filled[!is.na(df_filled[, 1]) & df_filled[, 1] != "", ]

cat_order <- c("VU", "EN", "CR")

df_factor <- as.data.frame(df_clean)

for (i in seq_len(ncol(df_factor))) {
  df_factor[[i]] <- factor(
    df_factor[[i]],
    levels = cat_order,
    ordered = TRUE
  )
}

rli_raw <- red::rli(df_factor, boot = FALSE)


RLI_df <- data.frame(
  Year = gsub("^X", "", colnames(rli_raw)),
  RLI  = as.numeric(rli_raw[1, ])
) %>%
  filter(Year != "Change/year") %>%
  mutate(Year = as.numeric(Year))


print(RLI_df)


library(ggplot2)

ggplot(RLI_df, aes(x = Year, y = RLI)) +
  geom_line(group = 1) +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()


RLI_plot_clean <- RLI_df %>%
  filter(!is.na(RLI)) %>%
  arrange(Year)

x_breaks <- sort(unique(RLI_plot_clean$Year))


p <- ggplot(RLI_plot_clean, aes(x = Year, y = RLI)) +
  geom_line(group = 1, linewidth = 1.1, color = "#2C3E50") +
  geom_point(size = 2.6, color = "#2C3E50") +
  geom_text(
    aes(label = round(RLI, 3)),
    vjust = -0.9,
    size = 3
  ) +
  scale_x_continuous(breaks = x_breaks) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    x = "Año",
    y = "Red List Index (RLI)"
  )

ggsave(
  filename = "RLI_all_species_2000_2025.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)


RLI_summary <- RLI_df %>%
  filter(!is.na(RLI)) %>%
  group_by(Year) %>%
  summarise(
    N_Species = first(N_Especies),
    First_Year = min(Year),
    First_RLI = RLI[Year == First_Year][1],
    Last_Year = max(Year),
    Last_RLI = RLI[Year == Last_Year][1],
    Change = Last_RLI - First_RLI,
    Percent_Change = (Change / First_RLI) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    First_RLI = round(First_RLI, 4),
    Last_RLI = round(Last_RLI, 4),
    Change = round(Change, 4),
    Percent_Change = round(Percent_Change, 2)
  )


# por clases --------------------------------------------------------------


library(red)
library(dplyr)
library(tibble)
library(ggplot2)

df <- read.csv(
  "C:/Users/walter.garcia/Documents/RLI/data/RLI_by_Group.csv",
  sep = ";",
  stringsAsFactors = FALSE
)

fill_forward <- function(x) {
  for (i in 2:length(x)) {
    if (is.na(x[i]) || x[i] == "") {
      x[i] <- x[i - 1]
    }
  }
  x
}


cat_order <- c("VU", "EN", "CR")


calculate_rli_group <- function(df_group, group_name) {
  
  mat <- df_group %>%
    select(Species, starts_with("X")) %>%
    column_to_rownames("Species") %>%
    as.matrix()
  
  # Rellenar hacia adelante
  mat <- t(apply(mat, 1, fill_forward))
  
  # Eliminar filas vacías
  mat <- mat[rowSums(is.na(mat) | mat == "") < ncol(mat), ]
  
  
  # Convertir a factores ordenados
  mat <- as.data.frame(mat)
  for (i in seq_len(ncol(mat))) {
    mat[[i]] <- factor(mat[[i]], levels = cat_order, ordered = TRUE)
  }
  
  # Calcular RLI
  rli_raw <- red::rli(mat, boot = FALSE)
  
  data.frame(
    Year = as.numeric(gsub("^X", "", colnames(rli_raw))),
    RLI  = as.numeric(rli_raw[1, ]),
    Grupo = group_name,
    N_Especies = nrow(mat)
  ) %>% filter(!is.na(Year))
}

rli_results <- df %>%
  group_split(clase2)

RLI_df <- bind_rows(
  lapply(rli_results, function(x) {
    calculate_rli_group(x, unique(x$clase2))
  })
)


print(RLI_df)

RLI_all <- calculate_rli_group(df, "Todas las especies")

RLI_df <- bind_rows(RLI_df, RLI_all)

RLI_plot <- RLI_df %>%
  filter(N_Especies > 0) %>%
  mutate(
    Year = as.numeric(Year),
    Grupo = factor(
      Grupo,
      levels = c(
        "Todas las especies",
        "Aves",
        "Invertebrados",
        "Peces",
        "Plantas",
        "Reptiles",
        "Anfibios",
        "Mamiferos"
      )
    )
  )


group_colors <- c(
  "Aves" = "#4ECDC4",
  "Invertebrados" = "#9B59B6",
  "Peces" = "#FECA57",
  "Plantas" = "#58B368",
  "Reptiles" = "#96CEB4",
  "Todas las especies" = "#2C3E50",
  "Anfibios"= "blue",
  "Mamiferos"="red"
)


ggplot(RLI_plot, aes(x = Year, y = RLI, color = Grupo, group = Grupo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Año",
    y = "Red List Index (RLI)",
    color = "Grupo"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )


RLI_plot_clean <- RLI_plot_clean %>%
  dplyr::filter(!is.na(Grupo))


years_real <- sort(unique(RLI_plot_clean$Year))


p <- ggplot(
  RLI_plot_clean,
  aes(
    x = factor(Year, levels = years_real),
    y = RLI,
    color = Grupo,
    group = Grupo
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  
  geom_text(
    aes(label = round(RLI, 3)),
    vjust = -0.9,
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


ggsave(
  filename = "RLI_por_grupo_2000_2025.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

library(dplyr)

RLI_summary_by_group <- RLI_plot %>%
  filter(!is.na(RLI)) %>%          # quitar NA
  arrange(Grupo, Year) %>%         # ordenar correctamente
  group_by(Grupo) %>%
  summarise(
    N_Species = first(N_Especies),
    First_Year = min(Year),
    First_RLI = RLI[Year == First_Year][1],
    Last_Year = max(Year),
    Last_RLI = RLI[Year == Last_Year][1],
    Change = Last_RLI - First_RLI,
    Percent_Change = (Change / First_RLI) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    First_RLI = round(First_RLI, 4),
    Last_RLI = round(Last_RLI, 4),
    Change = round(Change, 4),
    Percent_Change = round(Percent_Change, 2)
  )
print(RLI_summary_by_group)














