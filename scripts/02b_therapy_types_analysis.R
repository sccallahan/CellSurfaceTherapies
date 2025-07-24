#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(purrr)
library(stringr)
library(rlang)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
cleanFinalColumn <- function(.data, col) {
  col <- ensym(col)
  
  .data %>%
    mutate(
      !!col := map_chr(
        !!col,
        ~ {
          if (is.na(.x) || trimws(.x) == "") return(NA_character_)
          
          parts <- str_split(.x, "\\|")[[1]]
          parts <- unique(parts)
          parts <- parts[!(parts == "NA" | nchar(parts) > 7)]
          
          if (length(parts) == 0) {
            NA_character_
          } else {
            paste(parts, collapse = "|")
          }
        }
      )
    )
}

processFile <- function(file, curated) {
  readRDS(file) %>%
    select(c(1, 25:27)) %>%
    mutate(
      across(
        c(response_zero, response_drug_zero, response_title_zero),
        ~ na_if(.x, "NA")
      ),
      final = pmap_chr(
        list(response_zero, response_drug_zero, response_title_zero),
        ~ paste(na.omit(c(...)), collapse = "|")
      ),
      final = str_replace_all(final, "\\s+", ""),
      final = na_if(final, "")
    ) %>%
    cleanFinalColumn(final) %>%
    mutate(
      binary = ifelse(is.na(final), 0, 1)
    ) %>%
    left_join(curated, by = c("NCT.Number" = "ID")) %>%
    select(1, 5, 6, 10)
}

processTrialList <- function(studies, trial_list){
  intervention_types <- c('DRUG', 'BIOLOGICAL', 'RADIATION', 'GENETIC', 'COMBINATION_PRODUCT')
  
  pancan <- studies %>% 
    dplyr::select(c(NCT.Number, Study.Title, Study.URL, Study.Status, Brief.Summary, Conditions, Interventions)) %>% 
    filter(str_detect(Interventions, paste0(intervention_types, collapse = "|")))
  
  possible_intervention_strings <- c('DRUG:', 'BIOLOGICAL:', 'RADIATION:', 'GENETIC:', 'COMBINATION_PRODUCT:')
  pattern <- paste0("(", paste(possible_intervention_strings, collapse = "|"), ") (.*)(?:\\||$)")
  
  pancan <- pancan %>% 
    filter(NCT.Number %in% trial_list$ID)
  intervention_list <- str_split(string = pancan$Interventions, pattern = "\\|")
  
  # this will maintain trial by row number in `pancan` object
  df <- tibble(data = intervention_list) %>%
    unnest_wider(data, names_sep = "_")
  filter_df <- df %>%
    mutate(across(everything(), ~ if_else(str_detect(.x,  paste0(intervention_types, collapse = "|")), .x, NA))) %>% # removes procedures, etc.
    mutate(across(everything(), ~ str_replace_all(., pattern = pattern, replacement = "\\2"))) %>% # keep just drug name
    mutate(id = row_number()) %>%  # keeps track of trial number by row
    pivot_longer(cols = -id, names_to = "variable", values_to = "drug") %>%
    filter(!is.na(drug)) %>%
    dplyr::select(id, drug)
  
  df_new_rows <- filter_df %>%
    rowwise() %>%
    mutate(
      drug = list(c(
        str_extract(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}"),  # Extract content
        str_remove(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}")  # Remove the extracted content from original string
      ))
    ) %>%
    unnest(cols = drug) %>%  # Unnest the list to create two rows
    ungroup() %>%   # Ensure ungrouping after rowwise operations
    filter(!is.na(drug))
  
  filter_df_mod <- df_new_rows %>% 
    mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
    mutate(drug = tolower(drug))
  
  filter_df_mod <- filter_df_mod %>% 
    mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
    separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
    mutate(drug = gsub("^antibody ", "", drug)) %>%
    mutate(drug = gsub("^anti ", "", drug)) %>% 
    mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
    mutate(drug = gsub(" ", "", drug)) %>% 
    mutate(drug = gsub("antibodydrugconjugate", "", drug))
  
  filter_df_mod <- filter_df_mod %>% distinct()
  
  return(list(
    pancan = pancan,
    final_result = filter_df_mod
  ))
}

processNCIResults <- function(nci_dd, model_res){
  nci <- readRDS(nci_dd)
  model_res <- readRDS(model_res)
  nciDD_targets <- left_join(nci, model_res)
  
  nciDD_targets <- nciDD_targets %>%
    mutate(across(everything(), ~ str_replace_all(., " ", ""))) %>%
    mutate(across(everything(), ~ ifelse(as.character(.) == "NA", NA, as.character(.))))
  
  df_new_rows <- nciDD_targets %>%
    rowwise() %>%
    mutate(
      drug = list(c(
        str_extract(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}"),  # Extract content
        str_remove(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}")  # Remove the extracted content from original string
      ))
    ) %>%
    unnest(cols = drug) %>%  # Unnest the list to create two rows
    ungroup() %>%   # Ensure ungrouping after rowwise operations
    filter(!is.na(drug))
  
  model_res_mod <- df_new_rows %>% 
    mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
    mutate(drug = tolower(drug))
  
  model_res_mod <- model_res_mod %>% 
    mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
    separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
    mutate(drug = gsub("^antibody ", "", drug)) %>% 
    mutate(drug = gsub("^anti ", "", drug)) %>% 
    mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
    mutate(drug = gsub(" ", "", drug)) %>% 
    mutate(drug = gsub("antibodydrugconjugate", "", drug))
  
  model_res_mod <- model_res_mod %>% distinct()
}

nciFinalMerge <- function(df) {
  df %>%
    select(c(1, 2, 26:27)) %>%
    mutate(
      across(
        c(response_zero, response_drug_zero),
        ~ na_if(.x, "NA")
      ),
      final = pmap_chr(
        list(response_zero, response_drug_zero),
        ~ paste(na.omit(c(...)), collapse = "|")
      ),
      final = str_replace_all(final, "\\s+", ""),
      final = na_if(final, "")
    ) %>%
    cleanFinalColumn(final) %>%
    mutate(
      binary = ifelse(is.na(final), 0, 1)
    ) %>% 
    distinct() %>% 
    group_by(across(1:2)) %>%
    summarise(
      across(c(response_zero, response_drug_zero), ~ paste(unique(na.omit(.x)), collapse = "|")),
      final = final %>%
        str_split("\\|") %>%
        unlist() %>%
        na.omit() %>%
        unique() %>%
        paste(collapse = "|"),
      binary = max(binary),
      .groups = "drop"
    )%>% 
    mutate(across(where(is.character), ~ na_if(.x, ""))) %>% 
    dplyr::select(-1) %>% 
    left_join(curated, by = c("NCT.Number" = "ID")) %>% 
    dplyr::select(1,4,5,9)
}

finalCalls <- function(nci_res, clin_res){
  final_calls <- nci_res %>%
    rename(
      final_nci = final,
      binary_nci = binary
    ) %>%
    left_join(clin_res) %>%
    mutate(
      collapsed = pmap_chr(
        list(final, final_nci),
        ~ {
          vals <- c(...)
          split_vals <- vals %>%
            na.omit() %>%
            str_split("\\|") %>%
            unlist() %>%
            unique()
          if (length(split_vals) == 0) NA_character_ else paste(sort(split_vals), collapse = "|")
        }
      ),
      binary_collapsed = pmax(binary, binary_nci, na.rm = TRUE)
    ) %>% 
    dplyr::select(1,4,7,8)
}

processFileNeg <- function(file) {
  readRDS(file) %>%
    select(c(1, 25:27)) %>%
    mutate(
      across(
        c(response_zero, response_drug_zero, response_title_zero),
        ~ na_if(.x, "NA")
      ),
      final = pmap_chr(
        list(response_zero, response_drug_zero, response_title_zero),
        ~ paste(na.omit(c(...)), collapse = "|")
      ),
      final = str_replace_all(final, "\\s+", ""),
      final = na_if(final, "")
    ) %>%
    cleanFinalColumn(final) %>%
    mutate(
      binary = ifelse(is.na(final), 0, 1)
    ) %>% 
    dplyr::select(1,5,6)
}

processTrialListNeg <- function(studies, trial_list){
  
  pancan <- studies %>% 
    dplyr::select(c(NCT.Number, Study.Title, Study.URL, Study.Status, Brief.Summary, Conditions, Interventions)) #%>% 
  pancan <- pancan %>% 
    filter(NCT.Number %in% trial_list$NCT.Number)
  intervention_list <- str_split(string = pancan$Interventions, pattern = "\\|")
  
  # this will maintain trial by row number in `pancan` object
  df <- tibble(data = intervention_list) %>%
    unnest_wider(data, names_sep = "_")
  filter_df <- df %>%
    mutate(across(everything(), ~ str_replace_all(., "\\b[A-Z]+:(.*?)\\b", "\\1"))) %>% # keep just drug name
    mutate(id = row_number()) %>%  # keeps track of trial number by row
    pivot_longer(cols = -id, names_to = "variable", values_to = "drug") %>%
    filter(!is.na(drug)) %>%
    dplyr::select(id, drug)
  
  df_new_rows <- filter_df %>%
    rowwise() %>%
    mutate(
      drug = list(c(
        str_extract(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}"),  # Extract content
        str_remove(drug, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}")  # Remove the extracted content from original string
      ))
    ) %>%
    unnest(cols = drug) %>%  # Unnest the list to create two rows
    ungroup() %>%   # Ensure ungrouping after rowwise operations
    filter(!is.na(drug))
  
  filter_df_mod <- df_new_rows %>% 
    mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
    mutate(drug = tolower(drug))
  
  filter_df_mod <- filter_df_mod %>% 
    mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
    separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
    mutate(drug = gsub("^antibody ", "", drug)) %>%
    mutate(drug = gsub("^anti ", "", drug)) %>% 
    mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
    mutate(drug = gsub(" ", "", drug)) %>% 
    mutate(drug = gsub("antibodydrugconjugate", "", drug))
  
  filter_df_mod <- filter_df_mod %>% distinct()
  
  return(list(
    pancan = pancan,
    final_result = filter_df_mod
  ))
}

#~~~~~~~~~~~~~~~~~~#
# read in data ----
#~~~~~~~~~~~~~~~~~~#
studies_ref <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ================= ACCURACY =================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~#
# read in/process clinical trial calls ----
#~~~~~~~~~~~~~~~~~~#
## llama3.3 ----
llama33_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-3_tempzero.rds'), curated = curated)
table(llama33_clin$binary)
llama33_clin %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## llama3.1 ----
llama31_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-1_tempzero.rds'), curated = curated)
table(llama31_clin$binary)
llama31_clin %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## qwen ----
qwen_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_qwen2-5_tempzero.rds'), curated = curated)
table(qwen_clin$binary)
qwen_clin %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## gemma ----
gemma_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_gemma3_tempzero.rds'), curated = curated)
table(gemma_clin$binary)
gemma_clin %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## granite ----
granite_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_granite_tempzero.rds'), curated = curated)
table(granite_clin$binary)
granite_clin %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## phi4 ----
phi_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_phi4_tempzero.rds'), curated = curated)
table(phi_clin$binary)
phi_clin %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

## medllama ----
medllama_clin <- processFile(file = here('data/processed/clinical_trial_therapy_guess_LLM_medllama_tempzero.rds'), curated = curated)
table(medllama_clin$binary)
medllama_clin %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~#
# read in/process NCI calls ----
#~~~~~~~~~~~~~~~~~~#
## This cleans all curated trials
trials_processed <- processTrialList(studies = studies_ref, trial_list = curated)
trials_processed_res <- trials_processed$final_result
pancan_res <- trials_processed$pancan

## Llama 3.3 ----
llama33_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                         model_res = here('data/processed/nci_therapy_guess_LLM_llama3-3_tempzero.rds'))

merge_first <- left_join(trials_processed_res, llama33_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

llama33_final <- nciFinalMerge(merge_second)

table(llama33_final$binary) 
llama33_final %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## Llama 3.1 ----
llama31_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                                 model_res = here('data/processed/nci_therapy_guess_LLM_llama3-1_tempzero.rds'))

merge_first <- left_join(trials_processed_res, llama31_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

llama31_final <- nciFinalMerge(merge_second)

table(llama31_final$binary)
llama31_final %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## qwen ----
qwen_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                                 model_res = here('data/processed/nci_therapy_guess_LLM_qwen2-5_tempzero.rds'))

merge_first <- left_join(trials_processed_res, qwen_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

qwen_final <- nciFinalMerge(merge_second)

table(qwen_final$binary)
qwen_final %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## gemma ----
gemma_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                              model_res = here('data/processed/nci_therapy_guess_LLM_gemma3_tempzero.rds'))

merge_first <- left_join(trials_processed_res, gemma_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

gemma_final <- nciFinalMerge(merge_second)

table(gemma_final$binary)
gemma_final %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## granite ----
granite_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                               model_res = here('data/processed/nci_therapy_guess_LLM_granite_tempzero.rds'))

merge_first <- left_join(trials_processed_res, granite_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

granite_final <- nciFinalMerge(merge_second)

table(granite_final$binary)
granite_final %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## phi4 ----
phi_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                               model_res = here('data/processed/nci_therapy_guess_LLM_phi4_tempzero.rds'))

merge_first <- left_join(trials_processed_res, phi_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

phi_final <- nciFinalMerge(merge_second)

table(phi_final$binary)
phi_final %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## medllama ----
medllama_nci <- processNCIResults(nci_dd = here('data/processed/nci_dd.rds'),
                               model_res = here('data/processed/nci_therapy_guess_LLM_medllama_tempzero.rds'))

merge_first <- left_join(trials_processed_res, medllama_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

medllama_final <- nciFinalMerge(merge_second)

table(medllama_final$binary)
medllama_final %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~#
# merge for final calls ----
#~~~~~~~~~~~~~~~~~~#
## llama 3.3 ----
llama33_combined <- finalCalls(nci_res = llama33_final,
                               clin_res = llama33_clin)

table(llama33_combined$binary_collapsed) 
llama33_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## llama 3.1 ----
llama31_combined <- finalCalls(nci_res = llama31_final,
                               clin_res = llama31_clin)

table(llama31_combined$binary_collapsed) 
llama31_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## qwen ----
qwen_combined <- finalCalls(nci_res = qwen_final,
                               clin_res = qwen_clin)

table(qwen_combined$binary_collapsed)
qwen_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## gemma ----
gemma_combined <- finalCalls(nci_res = gemma_final,
                               clin_res = gemma_clin)

table(gemma_combined$binary_collapsed)
gemma_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## granite ----
granite_combined <- finalCalls(nci_res = granite_final,
                             clin_res = granite_clin)

table(granite_combined$binary_collapsed)
granite_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## phi4 ----
phi_combined <- finalCalls(nci_res = phi_final,
                               clin_res = phi_clin)

table(phi_combined$binary_collapsed)
phi_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

## medllama3-v20 ----
medllama_combined <- finalCalls(nci_res = medllama_final,
                           clin_res = medllama_clin)

table(medllama_combined$binary_collapsed)
medllama_combined %>%
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NEGATIVE RESULTS/SPECIFICITY =================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~#
# read in/process negative clinical trial calls ----
#~~~~~~~~~~~~~~~~~~#
studies_guess_negative <- studies_ref[!studies_ref$NCT.Number %in% curated$ID, ]
set.seed(1337)
studies_guess_negative <- studies_guess_negative %>% 
  sample_n(814)

## llama 3.3 ----
llama33_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-3_tempzero_negative.rds'))
table(llama33_clin_neg$binary)
table(llama33_clin_neg$final)

## llama 3.1 ----
llama31_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-1_tempzero_negative.rds'))
table(llama31_clin_neg$binary)
table(llama31_clin_neg$final)

## qwen ----
qwen_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_qwen2-5_tempzero_negative.rds'))
table(qwen_clin_neg$binary)
table(qwen_clin_neg$final)

## gemma ----
gemma_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_gemma3_tempzero_negative.rds'))
table(gemma_clin_neg$binary)
table(gemma_clin_neg$final)

## granite ----
granite_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_granite_tempzero_negative.rds'))
table(granite_clin_neg$binary)
table(granite_clin_neg$final)

## phi4 ----
phi_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_phi4_tempzero_negative.rds'))
table(phi_clin_neg$binary)
table(phi_clin_neg$final)

## medllama ----
medllama_clin_neg <- processFileNeg(file = here('data/processed/clinical_trial_therapy_guess_LLM_medllama_tempzero_negative.rds'))
table(medllama_clin_neg$binary)
table(medllama_clin_neg$final)

#~~~~~~~~~~~~~~~~~~#
# read in/process NCI calls ----
#~~~~~~~~~~~~~~~~~~#
## This cleans all curated trials
trials_processed_neg <- processTrialListNeg(studies = studies_ref, trial_list = studies_guess_negative)
trials_processed_res_neg <- trials_processed_neg$final_result
pancan_res_neg <- trials_processed_neg$pancan

## Llama 3.3 ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, llama33_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

llama33_final_neg <- nciFinalMerge(merge_second)

table(llama33_final_neg$binary)

## Llama 3.1 ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, llama31_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

llama31_final_neg <- nciFinalMerge(merge_second)

table(llama31_final_neg$binary)

## qwen ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, qwen_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

qwen_final_neg <- nciFinalMerge(merge_second)

table(qwen_final_neg$binary)

## gemma ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, gemma_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

gemma_final_neg <- nciFinalMerge(merge_second)

table(gemma_final_neg$binary)

## granite ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, granite_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

granite_final_neg <- nciFinalMerge(merge_second)

table(granite_final_neg$binary)

## phi4 ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, phi_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

phi_final_neg <- nciFinalMerge(merge_second)

table(phi_final_neg$binary)

## medllama ----
# we use the same NCI results since that won't change
merge_first <- left_join(trials_processed_res_neg, medllama_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  # as_tibble %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

medllama_final_neg <- nciFinalMerge(merge_second)

table(medllama_final_neg$binary) 

#~~~~~~~~~~~~~~~~~~#
# merge for final calls ----
#~~~~~~~~~~~~~~~~~~#
## llama 3.3 ----
llama33_combined_neg <- finalCalls(nci_res = llama33_final_neg,
                               clin_res = llama33_clin_neg)

table(llama33_combined_neg$binary_collapsed)


## llama 3.1 ----
llama31_combined_neg <- finalCalls(nci_res = llama31_final_neg,
                                   clin_res = llama31_clin_neg)

table(llama31_combined_neg$binary_collapsed)

## qwen ----
qwen_combined_neg <- finalCalls(nci_res = qwen_final_neg,
                                   clin_res = qwen_clin_neg)

table(qwen_combined_neg$binary_collapsed)

## gemma ----
gemma_combined_neg <- finalCalls(nci_res = gemma_final_neg,
                                clin_res = gemma_clin_neg)

table(gemma_combined_neg$binary_collapsed)

## granite ----
granite_combined_neg <- finalCalls(nci_res = granite_final_neg,
                                 clin_res = granite_clin_neg)

table(granite_combined_neg$binary_collapsed)

## phi4 ----
phi_combined_neg <- finalCalls(nci_res = phi_final_neg,
                                   clin_res = phi_clin_neg)

table(phi_combined_neg$binary_collapsed)

## medllama ----
medllama_combined_neg <- finalCalls(nci_res = medllama_final_neg,
                               clin_res = medllama_clin_neg)

table(medllama_combined_neg$binary_collapsed)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MERGE AND SAVE =================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# names
keys <- c("llama3-3", "llama3-1", "qwen2-5", "gemma3-27b", "granite", "phi4", "medllama")

# df1 and df2 lists
df1_list <- list(llama33_combined, llama31_combined, qwen_combined,
                 gemma_combined, granite_combined, phi_combined, medllama_combined)
df2_list <- list(llama33_combined_neg, llama31_combined_neg, qwen_combined_neg,
                 gemma_combined_neg, granite_combined_neg, phi_combined_neg, medllama_combined_neg)

# Combine and save using Map()
Map(function(name, d1, d2) {
  merged <- bind_rows(
    d1 %>% mutate(groundTruth = 1),
    d2 %>% mutate(groundTruth = 0)
  )
  saveRDS(merged, here('output/results',
                       paste0("therapy_guess_merged_pos_neg_", name, "_tempzero.rds")))
}, keys, df1_list, df2_list)