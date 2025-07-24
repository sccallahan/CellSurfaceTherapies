#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(usmap)
library(sf)
library(rnaturalearth)


#~~~~~~~~~~~~~~~~~~#
# read in studies ----
#~~~~~~~~~~~~~~~~~~#
studies <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

#~~~~~~~~~~~~~~~~~~#
# read in states/countries ----
#~~~~~~~~~~~~~~~~~~#
states <- readLines(here('data/raw/state_names.txt')) %>%
  .[-1]
countries <- read.fwf(here('data/raw/country-list.txt'), widths = c(5, 100), header = FALSE,
                      skip = 1, strip.white = TRUE, col.names = c('FIPS_ID', 'COUNTRY_NAME')) %>%
  filter(row_number() > 1)

#~~~~~~~~~~~~~~~~~~#
# extract states and countries ----
#~~~~~~~~~~~~~~~~~~#
## all curated trials ----
# extract states and countries
studies_trim <- studies %>% 
  filter(studies$NCT.Number %in% curated$ID) %>% 
  dplyr::select(c('NCT.Number', 'Locations')) %>% 
  mutate(Locations = tolower(Locations)) %>% 
  mutate(states_extract = sapply(str_extract_all(Locations, paste(tolower(states), collapse = "|")),
                               
                          function(x) if (length(x) > 0) paste(unique(x), collapse = "|") else NA)) %>%
  mutate(states_extract = na_if(states_extract, 'NA')) %>% 
  mutate(countries_extract = sapply(str_extract_all(Locations, paste(tolower(countries[[2]]), collapse = "|")),
                                    
                                 function(x) if (length(x) > 0) paste(unique(x), collapse = "|") else NA)) %>% 
  mutate(countries_extract = na_if(countries_extract, 'NA'))


# separate into different dataframes
states_data <- studies_trim %>% 
  dplyr::select(1,3) %>% 
  separate_rows(states_extract, sep = '\\|') %>% 
  filter(!is.na(states_extract))
length(unique(states_data$NCT.Number))

countries_data <- studies_trim %>% 
  dplyr::select(1,4) %>% 
  separate_rows(countries_extract, sep = '\\|') %>% 
  filter(!is.na(countries_extract))

### plotting, states ----
states_data_plot <- states_data %>% 
  count(states_extract)

states_data_plot$state <- statepop$abbr[match(states_data_plot$states_extract, tolower(statepop$full))]  # Convert to abbreviations

# Remove any NA (for unmatched states)
states_data_plot <- states_data_plot %>% filter(!is.na(state))

# Plot heatmap using usmap
plot_usmap(data = states_data_plot, values = "n", regions = "states", labels = TRUE) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Clinical Trials by State", fill = "Occurrences") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(here('output/figures/clinical_trials_map_states_curated.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_states_curated.svg'), width = 11)

### plotting, countries ----
countries_data_plot <- countries_data %>% 
  count(countries_extract)
# Load world country boundaries as an sf object
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  mutate(name = ifelse(name == 'United States of America', 'United States', name))

# Merge data with the world map
countries_data_plot <- world %>%
  mutate(name = tolower(name)) %>% 
  left_join(countries_data_plot, by = c("name" = "countries_extract"))

# Plot world map
ggplot(countries_data_plot) +
  geom_sf(aes(fill = n), color = "white", size = 0.1) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Clinical Trials by Country", fill = "Value") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
ggsave(here('output/figures/clinical_trials_map_countries_curated.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_countries_curated.svg'), width = 11)

## missed trials ----
misses <- read_delim(here('output/results/missed_trials_gpt4o_drug_ttddb.tsv'))
misses_clean <- misses %>% 
  dplyr::select(1) %>% 
  distinct()

### plotting, state ----
misses_state <- misses_clean %>% 
  left_join(states_data, by = c('ID' = 'NCT.Number')) %>% 
  count(states_extract)

misses_state$state <- statepop$abbr[match(misses_state$states_extract, tolower(statepop$full))]  # Convert to abbreviations

# Remove any NA (for unmatched states)
misses_state <- misses_state %>% filter(!is.na(state))

# Plot heatmap using usmap
plot_usmap(data = misses_state, values = "n", regions = "states", labels = TRUE) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Unmatched Clinical Trials by State", fill = "Occurrences") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(here('output/figures/clinical_trials_map_states_missed.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_states_missed.svg'), width = 11)


### plotting, countries ----
misses_country <- misses_clean %>% 
  left_join(countries_data, by = c('ID' = 'NCT.Number')) %>% 
  count(countries_extract)

world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  mutate(name = ifelse(name == 'United States of America', 'United States', name))

# Merge data with the world map
misses_country <- world %>%
  mutate(name = tolower(name)) %>% 
  left_join(misses_country, by = c("name" = "countries_extract"))

# Plot heatmap
ggplot(misses_country) +
  geom_sf(aes(fill = n), color = "white", size = 0.1) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Unmatched Clinical Trials by Country", fill = "Value") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
ggsave(here('output/figures/clinical_trials_map_countries_missed.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_countries_missed.svg'), width = 11)

## missed trial percent ----
percent_miss_states <- left_join(states_data_plot, misses_state, by = 'states_extract') %>% 
  rename_with(~ c('n_state', 'state', 'n_miss'), .cols = c(2:4)) %>% 
  select(-5) %>% 
  mutate(across(n_miss, ~replace_na(., 0))) %>% 
  mutate(percent = round((n_miss / n_state) * 100, digits = 2))

plot_usmap(data = percent_miss_states, values = "percent", regions = "states", labels = TRUE) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Unmatched Clinical Trials by State", fill = "Percent") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(here('output/figures/clinical_trials_map_states_missed_percent.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_states_missed_percent.svg'), width = 11)

countries_data_plot_edit <- countries_data_plot %>% 
  rename(n_country = n)
misses_country_edit <- misses_country %>% 
  rename(n_miss = n)
percent_miss_countries <- cbind(countries_data_plot_edit, n_miss = misses_country_edit$n_miss) %>% 
  mutate(n_miss_edit = ifelse(is.na(n_country) & is.na(n_miss), NA,
                              ifelse(!is.na(n_country) & is.na(n_miss), 0, n_miss))) %>% 
  mutate(percent = round((n_miss_edit / n_country)*100, digits = 2))

ggplot(percent_miss_countries) +
  geom_sf(aes(fill = percent), color = "white", size = 0.1) +
  scale_fill_gradient(low = "#008CFF", high = "#FF383F", na.value = "black") +
  theme_minimal() +
  labs(title = "Unmatched Clinical Trials by Country", fill = "Value") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
ggsave(here('output/figures/clinical_trials_map_countries_missed_percent.png'), bg = 'white', width = 11)
ggsave(here('output/figures/clinical_trials_map_countries_missed_percent.svg'), width = 11)
