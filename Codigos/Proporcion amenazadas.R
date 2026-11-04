library(rredlist)
library(dplyr)
library(ggplot2)

token <- "QkqZt35of6w4oR6e11d2uME7ZPR9n61gdFsS"

# Especies evaluadas en Colombia
sp_col <- rl_countries("CO", key = token)

# Extraer dataframe
sp_col_df <- sp_col$assessments

str(sp_col_df)

sp_threatened <- sp_col_df %>%
  filter(red_list_category_code %in% c("CR","EN","VU"))

sp_count <- sp_threatened %>%
  group_by(year_published, red_list_category_code) %>%
  summarise(n = n(), .groups = "drop")

sp_prop <- sp_count %>%
  group_by(year_published) %>%
  mutate(
    total_year = sum(n),
    proportion = n / total_year
  )

sp_prop




ggplot(sp_prop,
       aes(x = year_published,
           y = proportion,
           fill = red_list_category_code)) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Año",
    y = "Proporción de especies",
    fill = "Categoría IUCN"
  )
