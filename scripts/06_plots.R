#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(waffle)
library(colorR)
library(VennDiagram)
library(patchwork)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
grab_jsonl_results <- function(files) {
  all_data <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    con <- file(files[i], open = "r")
    df <- tryCatch(
      stream_in(con, verbose = FALSE),
      error = function(e) {
        warning(sprintf("Error reading file %s: %s", files[i], e$message))
        NULL
      },
      finally = close(con)
    )
    
    if (!is.null(df)) {
      # Rename the second column to response{n}
      colnames(df)[2] <- paste0("response", i)
      all_data[[i]] <- df
    }
  }
  
  # Remove NULL entries
  all_data <- Filter(Negate(is.null), all_data)
  
  if (length(all_data) == 0) {
    warning("No valid data frames to combine.")
    return(NULL)
  }
  
  # Perform left joins sequentially by 'id'
  combined_df <- all_data[[1]]
  
  for (j in 2:length(all_data)) {
    combined_df <- left_join(combined_df, all_data[[j]], by = "id")
  }
  
  return(combined_df)
}


get_totals <- function(df, name) {
  df %>%
    mutate(across(2:4, as.numeric)) %>%
    summarise(across(2:4, sum)) %>%
    mutate(Model = name, .before = 1)
}

#~~~~~~~~~~~~~~~~~~#
# read in trial matching and curated ----
#~~~~~~~~~~~~~~~~~~#
missed_trials_gpt_nofuzz <- read_tsv(here('output/results/missed_trials_gpt4o_nofuzz.tsv'))
matched_trials_gpt_nofuzz <- read_tsv(here('output/results/matched_trials_gpt4o_nofuzz.tsv'))
missed_trials_drug_ttddb <- read_tsv(here('output/results/missed_trials_drug_gpt4o_ttddb.tsv'))
matched_trials_drug_ttddb <- read_tsv(here('output/results/matched_trials_drug_gpt4o_ttddb.tsv'))
missed_trials_gpt_drug_ttddb <- read_tsv(here('output/results/missed_trials_gpt4o_drug_ttddb.tsv'))
matched_trials_gpt_drug_ttddb <- read_tsv(here('output/results/matched_trials_gpt4o_drug_ttddb.tsv'))

curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))
curated_split <- curated %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))

#~~~~~~~~~~~~~~~~~~#
# read in error list ----
#~~~~~~~~~~~~~~~~~~#
errors_annot <- read_tsv(file = here('output/results/missed_trials_gpt4o_drug_ttddb_annotated_update.tsv'))

#~~~~~~~~~~~~~~~~~~#
# cell surface data ----
#~~~~~~~~~~~~~~~~~~#
### curated ----
llama <- readRDS(here('output/results/curated_cell_surface_response_llama3-3_tempzero.rds'))
qwen <- readRDS(here('output/results/curated_cell_surface_response_qwen2-5_tempzero.rds'))

dat <- list.files(here('data/processed'), pattern = "^cell_surface_curated.*\\.jsonl$", full.names = TRUE)
gpt <- grab_jsonl_results(dat)
colnames(gpt) <- c('Target', 'surface1', 'surface2', 'surface3')

llama3_1 <- readRDS(here('output/results/curated_cell_surface_response_llama3-1_tempzero.rds'))
phi4 <- readRDS(here('output/results/curated_cell_surface_response_phi4_tempzero.rds'))
gemma <- readRDS(here('output/results/curated_cell_surface_response_gemma3-27b_tempzero.rds'))
granite <- readRDS(here('output/results/curated_cell_surface_response_granite_tempzero.rds'))
medllama <- readRDS(here('output/results/curated_cell_surface_response_medllama_tempzero.rds'))

### all genes ----
llama_allgenes <- readRDS(here('output/results/allgenes_cell_surface_response_llama3-3_tempzero.rds'))
gemma_allgenes <- readRDS(here('output/results/allgenes_cell_surface_response_gemma3-27b_tempzero.rds'))
qwen_allgenes <- readRDS(here('output/results/allgenes_cell_surface_response_qwen2-5_tempzero.rds'))
medllama_allgenes <- readRDS(here('output/results/allgenes_cell_surface_response_medllama_tempzero.rds'))

datAll <- list.files(here('data/processed'), pattern = "^cell_surface_allgenes.*\\.jsonl$", full.names = TRUE)
gpt_allgenes <- grab_jsonl_results(datAll)
colnames(gpt_allgenes) <- c('Target', 'surface1', 'surface2', 'surface3')

### literature calls ----
exp_calls <- read.csv(here('data/processed/cell_surface_experimental_all_genes.csv'))
compsets <- read.csv(here('data/processed/cell_surface_computational_all_genes.csv'))
result <- compsets %>%
  group_by(name) %>%
  summarize(count = n_distinct(set), .groups = "drop") %>%
  filter(count >= 4)

comp_calls <- result[,1,drop = TRUE]
combined_calls <- union(as.character(exp_calls$x), comp_calls)
combined_calls_df <- tibble(
  gene = combined_calls
)

#~~~~~~~~~~~~~~~~~~#
# figures ----
#~~~~~~~~~~~~~~~~~~#
## Figure 1 = workflow ----
## Figure 2 ----
### Therapy type guess figures ----
### Binary guess combined results pos/neg ----
binarytherapyguess <- data.frame(
  Group = rep(c('Llama 3.3', 'Llama 3.1:70b', 'Qwen2.5:72b', 'Gemma 3:27b', 'Phi-4', 'Granite 3.2:8b', 'GPT-4o', 'MedLlama3'), each =2),
  Category = rep(c("Known Positive", "Known Negative"), times = 8),
  Correct = c(811, 578,
              813, 572,
              797, 732,
              807, 635,
              800, 490,
              803, 618,
              808, 761,
              814, 162),
  Incorrect =c(814-811, 814-578,
               814-813, 814-572,
               814-797, 814-732,
               814-807, 814-635,
               814-800, 814-490,
               814-803, 814-618,
               814-808, 814-761,
               814-814, 814-162)
)


binarytherapyguess_long <- binarytherapyguess %>%
  mutate(total = Correct + Incorrect) %>% 
  mutate(percent = (Correct / total) * 100,
         percent_un = 100 - percent) %>% 
  pivot_longer(cols = c('percent', 'percent_un'))

binarytherapyguess_long <- binarytherapyguess_long %>%
  mutate(FillCategory = recode(name, "percent" = "Correct", "percent_un" = "Incorrect")) %>%
  mutate(FillCategory = factor(FillCategory, levels = rev(sort(unique(FillCategory))))) %>% 
  mutate(Category = factor(Category, levels = unique(Category)))

fill_colors <- c("Incorrect" = "lightgray", 
                 "Correct" = "#D56B79")

ggplot(binarytherapyguess_long, aes(x = Group, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Category, scales = 'free_y') +
  scale_fill_manual(name = "", values = fill_colors, breaks = c("Correct", "Incorrect")) +
  labs(title = "Cell Surface Therapy Matches by LLM (Binary)",
       x = "",
       y = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    panel.spacing = unit(1.2, "lines")
  )
ggsave(filename = here('output/figures/therapy_guesses_binary_by_LLM_percent.png'), bg = 'white', width = 12, height = 6)
ggsave(filename = here('output/figures/therapy_guesses_binary_by_LLM_percent.svg'), width = 12, height = 6)


### All LLM combined results PERCENTS ----
allLLMres <- data.frame(
  Group = rep(c('Llama 3.3', 'Llama 3.1:70b', 'Qwen2.5:72b', 'Gemma 3:27b', 'Phi-4', 'Granite 3.2:8b', 'GPT-4o', 'MedLlama3'), each =2),
  Category = rep(c("Targets", "Trials"), times = 8),
  Matched = c(740, 723, 739, 722, 692, 676, 714, 697, 649, 632, 634, 617, 744, 729, 584, 569),
  Unmatched =c(831-740, 814-723, 831-739, 814-722, 831-692, 814-676, 831-714, 814-697, 831-649, 814-632, 831-634, 814-617, 831-744, 814-729, 831-584, 814-569)
)

allLLMres_long <- allLLMres %>%
  mutate(total = Matched + Unmatched) %>% 
  mutate(percent = (Matched / total) * 100,
         percent_un = 100 - percent) %>% 
  dplyr::select(-Matched, -Unmatched, -total) %>% 
  pivot_longer(cols = c('percent', 'percent_un'))

allLLMres_long <- allLLMres_long %>%
  mutate(FillCategory = recode(name, "percent" = "Matched", "percent_un" = "Unmatched")) %>%
  mutate(FillCategory = factor(FillCategory, levels = rev(sort(unique(FillCategory)))))

fill_colors <- c("Unmatched" = "lightgray", 
                 "Matched" = "#80FFE8")

ggplot(allLLMres_long, aes(x = Group, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Category, scales = 'free_y') +
  scale_fill_manual(name = "", values = fill_colors) +
  labs(title = "Matches by LLM",
       x = "",
       y = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    panel.spacing = unit(1.2, "lines")
  )
ggsave(filename = here('output/figures/matched_trials_by_LLM_percent.png'), bg = 'white', width = 12, height = 6)
ggsave(filename = here('output/figures/matched_trials_by_LLM_percent.svg'), width = 12, height = 6)


### got-4o matches PERCENTS ----
# Sample data
trialmatches <- data.frame(
  Group = rep(c("NCI", "ClinicalTrials.gov", "Combined Approach"), each = 2),
  Category = rep(c("Targets", "Trials"), times = 3),
  Matches = c(628, 613, 624, 619, 744, 729),
  Unmatched = c(831-628, 814-613, 831-624, 814-619, 831-744, 814-729)
)

# Reshape data for ggplot
trialmatches_long <- trialmatches %>%
  mutate(total = Matches + Unmatched) %>% 
  mutate(percent = (Matches / total) * 100,
         percent_un = 100 - percent) %>% 
  dplyr::select(-Matches, -Unmatched, -total) %>% 
  pivot_longer(cols = c('percent', 'percent_un'))

trialmatches_long <- trialmatches_long %>%
  mutate(FillCategory = recode(name, "percent" = "Matched", "percent_un" = "Unmatched")) %>%
  mutate(FillCategory = ifelse(FillCategory == "Unmatched", "Unmatched", Group)) %>%
  mutate(FillCategory = factor(FillCategory, levels = c("Unmatched","NCI", "ClinicalTrials.gov", "Combined Approach")),
         Group = factor(Group, levels = c("ClinicalTrials.gov", "NCI", "Combined Approach")))

fill_colors <- c("Unmatched" = "lightgray", 
                 "NCI" = "#97D2FB", 
                 "ClinicalTrials.gov" = "#83BCFF", 
                 "Combined Approach" = "#80FFE8")

# Plot
ggplot(trialmatches_long, aes(x = Group, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Category) +
  scale_fill_manual(name = "", values = fill_colors) +
  labs(title = "GPT-4o Matches by Workflow",
       x = "",
       y = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    panel.spacing = unit(1.2, "lines")
  )
ggsave(filename = here('output/figures/matched_trials_by_workflow_percent.png'), bg = 'white', width = 10, height = 6)
ggsave(filename = here('output/figures/matched_trials_by_workflow_percent.svg'), width = 10, height = 6)

### Venns for method comparison ----
#### Targets ----
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = 628,   # Total elements in Set A, nci
  area2 = 624,   # Total elements in Set B, drugs
  cross.area = 508, # Overlapping elements between A and B
  category = c("NCI", "ClinicalTrials.gov"),
  lwd = 3,
  col = c("#97D2FB", "#83BCFF"),
  fill = c(alpha("#97D2FB",0.3), alpha("#83BCFF", 0.3)),
  cat.cex = 3,
  cex = 3,
  
  # Adjust category label positions
  cat.pos = c(315, 30),  # Angles to position labels outside circles
  cat.dist = c(0.05, 0.05)  # Distance of labels from the circles
)

#### png
png(here('output/figures/venn_diagram_targets.png'), width = 900, height = 800)
grid.newpage()
grid.draw(venn.plot)
dev.off()
#### svg
svg(here('output/figures/venn_diagram_targets.svg'), width = 12, height = 11)
grid.newpage()
grid.draw(venn.plot)
dev.off()

#### Trials ----
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = 613,   # Total elements in Set A, nci
  area2 = 619,   # Total elements in Set B, drugs
  cross.area = 503, # Overlapping elements between A and B
  category = c("NCI", "ClinicalTrials.gov"),
  lwd = 3,
  col = c("#97D2FB", "#83BCFF"),
  fill = c(alpha("#97D2FB",0.3), alpha("#83BCFF", 0.3)),
  cat.cex = 3,
  cex = 3,
  
  # Adjust category label positions
  cat.pos = c(315, 30),  # Angles to position labels outside circles
  cat.dist = c(0.05, 0.05)  # Distance of labels from the circles
)

#### png
png(here('output/figures/venn_diagram_trials.png'), width = 900, height = 800)
grid.newpage()
grid.draw(venn.plot)
dev.off()
#### svg
svg(here('output/figures/venn_diagram_trials.svg'), width = 12, height = 11) 
grid.newpage()
grid.draw(venn.plot)
dev.off()

### Matches by therapy ----
matched_trials_trim <- matched_trials_gpt_drug_ttddb %>% 
  dplyr::select(1, 5) %>% 
  distinct() %>% 
  mutate(method = rep('combined'))
table(curated$Class)
table(matched_trials_trim$Class)

nci_trials_trim <- matched_trials_gpt_nofuzz %>% 
  dplyr::select(1, 5) %>% 
  distinct() %>% 
  mutate(method = rep('nci'))
table(curated$Class)
table(nci_trials_trim$Class)

clin_trials_trim <- matched_trials_drug_ttddb %>% 
  dplyr::select(1, 5) %>% 
  distinct() %>% 
  mutate(method = rep('clinicaltrials'))
table(curated$Class)
table(clin_trials_trim$Class)

# need to make sure all counts have all possible therapies
# make sure every method has every therapy present
therapy_names <- c("ADC", "BiTE", "CAR-DC", "CAR-M", "CAR-NK", "CAR-T", "RPT")

ensure_all_therapies <- function(df, methodtype) {
  df <- as.data.frame(table(df$Class)) %>%
    mutate(method = methodtype)
  colnames(df)[1:2] <- c('Class', 'matched')
  
  # Create a full template with all therapies
  full_df <- data.frame(Class = therapy_names) %>%
    left_join(df, by = "Class") %>%
    mutate(matched = replace_na(matched, 0)) %>% 
    mutate(method = ifelse(is.na(method), methodtype, method))
  
  return(full_df)
}
# get counts
matched_counts <- ensure_all_therapies(matched_trials_trim, methodtype = 'combined')

nci_counts <- ensure_all_therapies(nci_trials_trim, methodtype = 'nci')

clin_counts <- ensure_all_therapies(clin_trials_trim, methodtype = 'clinicaltrials')

# total trial counts
total_counts <- as.data.frame(table(curated$Class))
colnames(total_counts)[1] <- c('Class')

# combine
match_data <- bind_rows(matched_counts, nci_counts, clin_counts)

match_data <- full_join(match_data, total_counts, by = "Class") %>%
  rename_with(~ c('total'), .cols = 4) %>% 
  mutate(percent = (matched / total) * 100,
         Unmatched = 100 - percent,
         Unmatched_num = total - matched)

# Create a palette function
color_palette <- colorRampPalette(c('#F72585', "#4CC9F0"))
colors <- color_palette(7)

# make long
match_data_long <- match_data %>%
  pivot_longer(cols = c(percent, Unmatched), names_to = "category", values_to = "value") %>% 
  mutate(FillCategory = ifelse(category == "Unmatched", "Unmatched", as.character(Class))) %>%
  mutate(FillCategory = factor(FillCategory, levels = c('Unmatched', 'ADC', 'BiTE', 'CAR-DC', 'CAR-M', 'CAR-NK', 'CAR-T', 'RPT')),
         category = factor(category, levels = c("Unmatched", "percent"))) %>% 
  mutate(method = factor(method, levels = c("combined", "nci", "clinicaltrials")))

fill_colors2 <- c(
  "ADC" = "#F72585",
  "BiTE" = "#DA4096",
  "CAR-DC" = "#BE5BA8",
  "CAR-M" = "#A177BA",
  "CAR-NK" = "#8492CC",
  "CAR-T" = "#68ADDE",
  "RPT" = "#4CC9F0",
  "Unmatched" = "gray80"
)

bar_labels <- match_data_long %>%
  group_by(Class, method) %>%
  summarise(num = unique(Unmatched_num), 
            value = sum(value), .groups = 'drop')

ggplot(match_data_long, aes(x = Class, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity") +
  geom_text(data = bar_labels,
            aes(x = Class, y = value, label = num),
            vjust = -0.5,
            size = 5,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "", values = fill_colors2) +
  labs(title = "Trial Matches by Therapy",
       x = "",
       y = "") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    strip.text = element_text(size = 16)
  ) +
  facet_wrap(~method, ncol = 1, labeller = as_labeller(c(
    'combined' = 'Combined Approach',
    'nci' = 'NCI',
    'clinicaltrials' = "ClinialTrials.gov")))

ggsave(filename = here('output/figures/matched_trials_by_therapy.png'), bg = 'white', height = 12, width = 10)
ggsave(filename = here('output/figures/matched_trials_by_therapy.svg'), height = 12, width = 10)

### Matches by therapy condensed ----
therapy_names_condensed <- c("ADC", "BiTE", "CAR-T", "Other CAR", "RPT")

ensure_all_therapies_condensed <- function(df, methodtype) {
  df <- as.data.frame(table(df$Class)) %>%
    mutate(method = methodtype)
  colnames(df)[1:2] <- c('Class', 'matched')
  
  # Create a full template with all therapies
  full_df <- data.frame(Class = therapy_names_condensed) %>%
    left_join(df, by = "Class") %>%
    mutate(matched = replace_na(matched, 0)) %>% 
    mutate(method = ifelse(is.na(method), methodtype, method))
  
  return(full_df)
}
# get counts
matched_counts_condensed <- matched_trials_trim %>% 
  mutate(Class = ifelse(Class == 'CAR-DC', 'Other CAR', 
                        ifelse(Class == 'CAR-M', 'Other CAR',
                               ifelse(Class == 'CAR-NK', 'Other CAR', Class))))
matched_counts_condensed <- ensure_all_therapies_condensed(matched_counts_condensed, methodtype = 'combined')

nci_counts_condensed <- nci_trials_trim %>% 
  mutate(Class = ifelse(Class == 'CAR-DC', 'Other CAR', 
                        ifelse(Class == 'CAR-M', 'Other CAR',
                               ifelse(Class == 'CAR-NK', 'Other CAR', Class))))
nci_counts_condensed <- ensure_all_therapies_condensed(nci_counts_condensed, methodtype = 'nci')

clin_counts_condensed <- clin_trials_trim %>% 
  mutate(Class = ifelse(Class == 'CAR-DC', 'Other CAR', 
                        ifelse(Class == 'CAR-M', 'Other CAR',
                               ifelse(Class == 'CAR-NK', 'Other CAR', Class))))
clin_counts_condensed <- ensure_all_therapies_condensed(clin_counts_condensed, methodtype = 'clinicaltrials')

# total trial counts
curated_condensed <- curated %>% 
  mutate(Class = ifelse(Class == 'CAR-DC', 'Other CAR', 
                        ifelse(Class == 'CAR-M', 'Other CAR',
                               ifelse(Class == 'CAR-NK', 'Other CAR', Class))))
total_counts_condensed <- as.data.frame(table(curated_condensed$Class))
colnames(total_counts_condensed)[1] <- c('Class')

# combine
match_data_condensed <- bind_rows(matched_counts_condensed, nci_counts_condensed, clin_counts_condensed)

match_data_condensed <- full_join(match_data_condensed, total_counts_condensed, by = "Class") %>%
  rename_with(~ c('total'), .cols = 4) %>% 
  mutate(percent = (matched / total) * 100,
         Unmatched = 100 - percent,
         Unmatched_num = total - matched)

# make long
match_data_long_condensed <- match_data_condensed %>%
  pivot_longer(cols = c(percent, Unmatched), names_to = "category", values_to = "value") %>% 
  mutate(FillCategory = ifelse(category == "Unmatched", "Unmatched", as.character(Class))) %>%
  mutate(FillCategory = factor(FillCategory, levels = c('Unmatched', 'ADC', 'BiTE', 'CAR-T', 'Other CAR', 'RPT')),
         category = factor(category, levels = c("Unmatched", "percent"))) %>% 
  mutate(method = factor(method, levels = c("combined", "nci", "clinicaltrials")))

fill_colors3 <- c(
  "ADC" = "#F72585",
  "BiTE" = "#DA4096",
  "Other CAR" = "#A177BA",
  "CAR-T" = "#68ADDE",
  "RPT" = "#4CC9F0",
  "Unmatched" = "gray80"
)

# total n
bar_labels_condensed <- match_data_long_condensed %>%
  group_by(Class, method) %>%
  summarise(num = unique(total),
            value = sum(value), .groups = 'drop')

ggplot(match_data_long_condensed, aes(x = Class, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity") +
  geom_text(data = bar_labels_condensed,
            aes(x = Class, y = value, label = num),
            vjust = -0.5,
            size = 5,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "", values = fill_colors3) +
  labs(title = "Trial Matches by Therapy",
       x = "",
       y = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),limits = c(0, 110), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    strip.text = element_text(size = 16)
  ) +
  facet_wrap(~method, ncol = 1, labeller = as_labeller(c(
    'combined' = 'Combined Approach',
    'nci' = 'NCI',
    'clinicaltrials' = "ClinialTrials.gov")))
ggsave(filename = here('output/figures/matched_trials_by_therapy_condensed_totals.png'), bg = 'white', height = 12, width = 6)
ggsave(filename = here('output/figures/matched_trials_by_therapy_condensed_totals.svg'), height = 12, width = 6)

### Matches by trial phase ----
matched_trials_trim_phase <- matched_trials_gpt_drug_ttddb %>% 
  dplyr::select(1, 3) %>% 
  distinct() %>% 
  mutate(method = rep('combined'))

nci_trials_trim_phase <- matched_trials_gpt_nofuzz %>% 
  dplyr::select(1, 3) %>% 
  distinct() %>% 
  mutate(method = rep('nci'))

clin_trials_trim_phase <- matched_trials_drug_ttddb %>% 
  dplyr::select(1, 3) %>% 
  distinct() %>% 
  mutate(method = rep('clinicaltrials'))

phase_counts <- as.data.frame(table(curated$Phase))
colnames(phase_counts)[1] <- c('Phase')

matched_phase_counts<- as.data.frame(table(matched_trials_trim_phase$Phase)) %>% 
  mutate(method = rep('combined'))
colnames(matched_phase_counts)[1] <- c('Phase')

nci_phase_counts <- as.data.frame(table(nci_trials_trim_phase$Phase)) %>% 
  mutate(method = rep('nci'))
colnames(nci_phase_counts)[1] <- c('Phase')

clin_phase_counts <- as.data.frame(table(clin_trials_trim_phase$Phase)) %>% 
  mutate(method = 'clinicaltrials')
colnames(clin_phase_counts)[1] <- c('Phase')

match_data_phase <- bind_rows(matched_phase_counts, nci_phase_counts, clin_phase_counts)

match_data_phase <- full_join(match_data_phase, phase_counts, by = "Phase") %>%
  rename_with(~ c('matched', 'total'), .cols = c(2,4)) %>% 
  mutate(percent = (matched / total) * 100,
         Unmatched = 100 - percent,
         Unmatched_num = total - matched)

match_data_phase_long <- match_data_phase %>%
  pivot_longer(cols = c(percent, Unmatched), names_to = "category", values_to = "value") %>% 
  mutate(FillCategory = ifelse(category == "Unmatched", "Unmatched", 'Matched')) %>%
  mutate(FillCategory = factor(FillCategory, levels = c("Unmatched", "Matched"))) %>% 
  mutate(method = factor(method, levels = c("combined", "nci", "clinicaltrials")))

# totals
bar_labels_phase <- match_data_phase_long %>%
  group_by(Phase, method) %>%
  summarise(num = unique(total),
            value = sum(value), .groups = 'drop')

fill_colors_phase <- c(
  "Matched" = "#FFCBA4",
  "Unmatched" = "gray80"
)

ggplot(match_data_phase_long, aes(x = Phase, y = value, fill = FillCategory)) +
  geom_bar(stat = "identity") +
  geom_text(data = bar_labels_phase,
            aes(x = Phase, y = value, label = num),
            vjust = -0.5,
            size = 5,
            inherit.aes = FALSE) +
  scale_fill_manual(name = "", values = fill_colors_phase) +
  labs(title = "Trial Matches by Phase",
       x = "",
       y = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 110), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
    strip.text = element_text(size = 16)
  ) +
  facet_wrap(~method, ncol = 1, labeller = as_labeller(c(
    'combined' = 'Combined Approach',
    'nci' = 'NCI',
    'clinicaltrials' = "ClinicalTrials.gov")))
ggsave(filename = here('output/figures/matched_trials_by_phase_totals.png'), bg = 'white', height = 12, width = 6)
ggsave(filename = here('output/figures/matched_trials_by_phase_totals.svg'), height = 12, width = 6)


## Figure 3 ----
### Therapy type errors ----
table(errors_annot$Target)
table(errors_annot$Class)

ggplot(data = errors_annot, aes(x = Class, fill = Class)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.5, size = 8) +
  labs(title = "Missed Trials by Therapy Type", y = "", x = "") +
  scale_fill_manual(values = c("ADC" = "#F72585",
                               "BiTE" = "#7209B7",
                               'CAR-NK' = "#3A0CA3",
                               'CAR-T' = "#4361EE",
                               'RPT' = "#4CC9F0"),
                    labels = c("ADC", "BiTE", 'CAR-NK', 'CAR-T', 'RPT'),
                    name = "") +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0),
    axis.title.x = element_blank()
  )
ggsave(filename = here('output/figures/missed_trials_by_therapy.png'), bg = 'white')
ggsave(filename = here('output/figures/missed_trials_by_therapy.svg'))

### Error type ----
capitalize_with_acronyms <- function(text, acronyms = c("LLM", "NA")) {
  # Convert to title case
  text_cap <- str_to_title(text)
  
  # Replace acronyms with their original uppercase form
  for (acronym in acronyms) {
    text_cap <- gsub(paste0("\\b", str_to_title(acronym), "\\b"), acronym, text_cap)
  }
  
  return(text_cap)
}

##### broader rename ----
errors_annot_split2 <- errors_annot %>% 
  separate_rows(Reason_category_mod, sep = ",") %>% 
  mutate(Reason_category_mod = str_trim(Reason_category_mod)) %>% 
  mutate(Reason_category_mod = str_replace_all(Reason_category_mod, "_", " ")) %>% 
  mutate(Reason_category_mod = capitalize_with_acronyms(Reason_category_mod)) %>% 
  mutate(Reason_category_mod = gsub('Nci', "NCI", Reason_category_mod)) %>% 
  mutate(error_class = ifelse(Reason_category_mod == 'Vague Therapy Name' | Reason_category_mod == 'Vague Therapy Description',
                        'Uncorrectable Error', 'Correctable Error')) %>% 
  mutate(Reason_category_mod = factor(Reason_category_mod, levels = c(
    'Vague Therapy Name', 'Vague Therapy Description', 'No NCI Match', 'String Match',
    'On-Topic', 'No Target Identified', 'Hallucination'
  )))

table(errors_annot_split2$Reason_category_mod)

ggplot(data = errors_annot_split2, aes(x = Reason_category_mod, fill = error_class)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.5, size = 8) +
  labs(title = "Error Types", y = "", x = "") +
  scale_fill_manual(values = c("Uncorrectable Error" = "#efbcd5",
                               "Correctable Error" = "#4b5267"),
                    labels = c("Correctable Error", "Uncorrectable Error"),
                    name = "") +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )
ggsave(filename = here('output/figures/missed_trials_by_error_gpt4o.png'), bg = 'white', width = 10, height = 6)
ggsave(filename = here('output/figures/missed_trials_by_error_gpt4o.svg'), width = 10, height = 6)


#### modified labels ----
tmp2 <- errors_annot_split2[, c(5,10)]
df_counts2 <- tmp2 %>%
  count(Class, Reason_category_mod)
ggplot(data = df_counts2, aes(fill=Class, values=n)) +
  geom_waffle(color = "white", size = 0.5, n_rows = 3) +
  facet_wrap(~Reason_category_mod, ncol=1) +
  scale_fill_manual(values = c("ADC" = "#F72585",
                               "BiTE" = "#7209B7",
                               'CAR-DC' = "#9D4EDD",
                               'CAR-NK' = "#3A0CA3",
                               'CAR-T' = "#4361EE",
                               'RPT' = "#4CC9F0"),
                    labels = c("ADC", "BiTE", 'CAR-DC', 'CAR-NK', 'CAR-T', 'RPT'),
                    name = "") +
  labs(title = "Therapy Type Composition of Each Error Type") +
  theme_void() +
  theme(
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(filename = here('output/figures/missed_trials_waffle_mod.png'), bg = 'white')
ggsave(filename = here('output/figures/missed_trials_waffle_mod.svg'))

## Figure 4 made in map_data ----
## Figure 5 ----
### Cell surface ----
#### curated list ----
modelNames <- c('gpt', 'llama', 'llama3_1', 'qwen', 'gemma', 'phi4', 'granite', 'medllama')

#clean up phi4 first
phi4 <- phi4 %>%
  mutate(across(2:4, ~ {
    n <- nchar(.)
    first <- substr(., 1, 1)
    last <- substr(., n, n)
    
    ifelse(
      n == 1 & grepl("^[0-9]$", .),
      ., 
      ifelse(
        grepl("^[0-9]$", first) & grepl("^[0-9]$", last),
        paste0(first, last),
        ifelse(grepl("^[0-9]$", first), first,
               ifelse(grepl("^[0-9]$", last), last, NA_character_))
      )
    )
  }))

modelTotals <- map2_dfr(mget(modelNames), modelNames, get_totals)
modelTotalsMelt <- pivot_longer(modelTotals, cols = 2:4)

fillReplicate <- c(
  "surface1" = "#00B0A2",
  "surface2" = "#85D3FF",
  "surface3" = "#005672"
)

modelNamesProper <- c(
  'Qwen2.5:72b',
  'Phi-4',
  'MedLlama3',
  'Llama 3.1:70b',
  'Llama 3.3',
  'Granite 3.2:8b',
  'GPT-4o',
  'Gemma 3:27b'
)

modelTotalsMelt$Model <- factor(modelTotalsMelt$Model,  levels = rev(sort(unique(modelTotalsMelt$Model))))

ggplot(modelTotalsMelt, aes(x = Model, y = value, color = name)) +
  geom_point(size = 4, position = position_dodge(0.2)) +
  scale_color_manual(values = fillReplicate,
                     labels = c('Replicate 1', 'Replicate 2', 'Replicate 3')) +
  labs(title = "Number of Correct Assignments",
       x = "",
       y = "",
       color = "Replicate") +
  scale_x_discrete(labels = modelNamesProper) +
  ylim(40, 80) +
  coord_flip() +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here('output/figures/cell_surface_assignment_curated_edited.png'), bg = 'white')
ggsave(here('output/figures/cell_surface_assignment_curated_edited.svg'))

#### all genes ----
# clean llama
llama_allgenes <- llama_allgenes %>%
  mutate(across(2:4, ~ {
    last <- substr(., nchar(.), nchar(.))
    ifelse(grepl("^[0-9]$", last), last, NA_character_)
  }))

# gpt no answer = 0
gpt_allgenes <- gpt_allgenes %>%
  mutate(across(2:4, ~ {
    last <- substr(., nchar(.), nchar(.))
    ifelse(grepl("^[0-9]$", last), last, 0)
  }))

modelNamesAll <- c('gpt_allgenes', 'llama_allgenes', 'qwen_allgenes', 'gemma_allgenes', 'medllama_allgenes')
modelTotalsAll <- map2_dfr(mget(modelNamesAll), modelNamesAll, get_totals)
modelTotalsMeltAll <- pivot_longer(modelTotalsAll, cols = 2:4)

fillReplicate <- c(
  "surface1" = "#00B0A2",
  "surface2" = "#85D3FF",
  "surface3" = "#005672"
)

modelNamesProperAll <- c(
  'Qwen2.5:72b',
  'MedLlama3',
  'Llama 3.3',
  'GPT-4o',
  'Gemma 3:27b'
)

modelTotalsMeltAll$Model <- factor(modelTotalsMeltAll$Model,  levels = rev(sort(unique(modelTotalsMeltAll$Model))))

ggplot(modelTotalsMeltAll, aes(x = Model, y = value, color = name)) +
  geom_point(size = 4, position = position_dodge(0.2)) +
  geom_hline(aes(yintercept = 4898, linetype = "Count"), color = "#D95F02") + # see below code (referenceGenes) for yintercept value
  scale_color_manual(values = fillReplicate,
                     labels = c('Replicate 1', 'Replicate 2', 'Replicate 3')) +
  scale_linetype_manual(name = "Literature Assignments", values = c("Count" = "dashed")) +
  labs(title = "Number of Cell Surface Assignments",
       x = "",
       y = "",
       color = "Replicate") +
  scale_x_discrete(labels = modelNamesProperAll) +
  ylim(3000, 20000) +
  coord_flip() +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(here('output/figures/cell_surface_assignment_allgenes.png'), bg = 'white')
ggsave(here('output/figures/cell_surface_assignment_allgenes.svg'))

# confusion matrix barplots
# for gpt-4o >=2 = 1, <2 = 0
referenceGenes <- combined_calls_df[combined_calls_df$gene %in% gemma_allgenes$symbol,] #4898 genes match
groundTruth <- ifelse(gemma_allgenes$symbol %in% referenceGenes$gene, 1, 0)
names(groundTruth) <- gemma_allgenes$symbol

gptFix <- gpt_allgenes %>%
  mutate(across(2:4, ~as.numeric(as.character(.)))) %>%
  transmute(result = ifelse(rowSums(across(2:4)) >= 2, 1, 0))
gptFix <- gptFix$result


groundTruthDF <- data.frame(
  gene = names(groundTruth),
  groundTruth = groundTruth,
  gemma = gemma_allgenes[,2],
  gpt = gptFix,
  llama = llama_allgenes[,2],
  qwen = qwen_allgenes[,2],
  medllama = medllama_allgenes[,2]
)

colnames(groundTruthDF) <- c('gene', 'groundTruth', 'gemma', 'gpt', 'llama', 'qwen', 'medllama')

get_confusion <- function(pred, truth) {
  pred <- as.integer(pred)
  truth <- as.integer(truth)
  tibble(
    TP = sum(pred == 1 & truth == 1),
    FP = sum(pred == 1 & truth == 0),
    TN = sum(pred == 0 & truth == 0),
    FN = sum(pred == 0 & truth == 1)
  )
}

modelLabels <- c('gemma', 'gpt', 'llama', 'qwen', 'medllama')

confMatrix <- map_dfr(modelLabels, ~
                            get_confusion(groundTruthDF[[.x]], groundTruthDF$groundTruth) %>% mutate(Group = .x)
) %>%
  pivot_longer(cols = -Group, names_to = "Metric", values_to = "Count")


confMatrix$Group <- factor(confMatrix$Group,
                               levels = rev(sort(unique(confMatrix$Group))))
customLabels <- c(
  gemma = "Gemma 3:27b",
  gpt = "GPT-4o",
  llama = "Llama 3.3",
  medllama = "MedLlama3",
  qwen = "Qwen2.5:72b"
)
ggplot(confMatrix, aes(x = Metric, y = Group, fill = Count)) +
  geom_tile() +
  geom_text(aes(label = Count), color = "white", size = 5) +
  scale_fill_gradient(low = "#32cbff", high = "#ef9cda") +
  # scale_fill_viridis_c(option = 'H') +
  scale_y_discrete(labels = customLabels) +
  labs(title = "Confusion Matrix Metrics",
       x = "",
       y = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
    )
ggsave(here('output/figures/cell_surface_confusionMatrix_allgenes.png'), bg = 'white')
ggsave(here('output/figures/cell_surface_confusionMatrix_allgenes.svg'))    
