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
library(data.table)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
cleanFinalColumn <- function(.data, col) {
  col <- ensym(col)  # tidy capture
  
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

grab_jsonl_results <- function(files) {
  all_data <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    con <- file(files[i], open = "r")
    all_data[[i]] <- tryCatch(
      stream_in(con, verbose = FALSE),
      error = function(e) {
        warning(sprintf("Error reading file %s: %s", files[i], e$message))
        NULL
      },
      finally = close(con)
    )
  }
  
  # Remove NULLs and bind all into one data frame
  combined_df <- rbindlist(all_data, fill = TRUE)
  
  return(combined_df)
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

processNCIResultsGPT <- function(nci_dd, model_res){
  nci <- readRDS(nci_dd)
  model_res <- model_res
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

nciFinalMergeGPT <- function(df) {
  df %>%
    dplyr::select(c(1, 2, 26:27)) %>%
    mutate(
      across(
        c(response_drug, response_def),
        ~ na_if(.x, "NA")
      ),
      final = pmap_chr(
        list(response_drug, response_def),
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
      across(c(response_drug, response_def), ~ paste(unique(na.omit(.x)), collapse = "|")),
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

finalCallsGPT <- function(nci_res, clin_res){
  final_calls <- nci_res %>%
    dplyr::rename(
      final_nci = final,
      binary_nci = binary
    ) %>%
    dplyr::select(-4) %>%  # remove Class
    left_join(gpt_clin, by = c('NCT.Number' = 'trial')) %>%
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
    dplyr::select(1,6,7,8)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# KNOWN POSITIVE ===================== ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~#
# read in clinical trial calls ----
#~~~~~~~~~~~~~~~~~~#
json_dir <- here('data/openai_chunks_therapies_curated_answers/')
json_files1 <- list.files(json_dir, pattern = ".*title.*.jsonl", full.names = TRUE)
json_files2 <- list.files(json_dir, pattern = ".*summary.*.jsonl", full.names = TRUE)
json_files3 <- list.files(json_dir, pattern = ".*intervention.*.jsonl", full.names = TRUE)

gpt1 <- grab_jsonl_results(json_files1) %>% 
  dplyr::rename(response_title = response)
gpt2 <- grab_jsonl_results(json_files2) %>% 
  dplyr::rename(response_summary = response)
gpt3 <- grab_jsonl_results(json_files3) %>% 
  dplyr::rename(response_intervention = response)

gpt_merge <- gpt1 %>% 
  left_join(gpt2) %>% 
  left_join(gpt3)

gpt_clin <- gpt_merge %>% 
    mutate(across(c(response_title, response_summary, response_intervention), ~ na_if(.x, 'NA')),
           final = pmap_chr(list(response_title, response_summary, response_intervention), ~ paste(na.omit(c(...)), collapse = "|"))) %>%
    mutate(final = str_replace_all(final, "\\s+", "")) %>%
    mutate(final = na_if(final, "")) %>%
    cleanFinalColumn(final) %>%
    mutate(binary = ifelse(is.na(final), 0, 1)) %>%
    left_join(curated, by = c('trial' = 'ID')) %>%
    dplyr::select(1,5,6,10)
table(gpt_clin$binary)
gpt_clin %>%
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~#
# read in nci ----
#~~~~~~~~~~~~~~~~~~#
## process trials first
trials_processed <- processTrialList(studies = studies_ref, trial_list = curated)
trials_processed_res <- trials_processed$final_result
pancan_res <- trials_processed$pancan

## gpt results
json_dir2 <- here('data/openai_chunks_therapies_nci_answers/')  # update this path
json_files1b <- list.files(json_dir2, pattern = ".*drug.*.jsonl", full.names = TRUE)
json_files2b <- list.files(json_dir2, pattern = ".*definition.*.jsonl", full.names = TRUE)

gpt1_nci <- grab_jsonl_results(json_files1b) %>% 
  dplyr::rename(response_drug = response) %>% 
  dplyr::rename(termId = trial) %>% 
  mutate(termId = as.integer(termId))
gpt2_nci <- grab_jsonl_results(json_files2b) %>% 
  dplyr::rename(response_def = response) %>% 
  dplyr::rename(termId = trial) %>% 
  mutate(termId = as.integer(termId))

gpt_merge <- left_join(final_dictionary_definitions, gpt1_nci) %>% 
  left_join(gpt2_nci)

# merge back to nci dd
gpt_nci <- processNCIResultsGPT(nci_dd = here('data/processed/nci_dd.rds'),
                                model_res = gpt_merge)

merge_first <- left_join(trials_processed_res, gpt_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

gpt_final <- nciFinalMergeGPT(merge_second)

table(gpt_final$binary)
gpt_final %>% 
  mutate(match = str_detect(final, fixed(Class)),
         exact = final == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~#
# merge for final calls ----
#~~~~~~~~~~~~~~~~~~#
gpt_combined <- finalCallsGPT(nci_res = gpt_final,
                               clin_res = gpt_clin)

table(gpt_combined$binary_collapsed)
gpt_combined %>% 
  mutate(match = str_detect(collapsed, fixed(Class)),
         exact = collapsed == Class ) %>%
  summarise(n_matches = sum(match, na.rm = TRUE),
            n_miss = 814 - sum(match, na.rm = TRUE),
            n_exact_matches = sum(exact, na.rm = TRUE),
            n_miss_exact = 814 - sum(exact, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# KNOWN NEGATIVE ===================== ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~#
# read in/process negative clinical trial calls ----
#~~~~~~~~~~~~~~~~~~#
studies_guess_negative <- studies_ref[!studies_ref$NCT.Number %in% curated$ID, ]
set.seed(1337)
studies_guess_negative <- studies_guess_negative %>% 
  sample_n(814)

json_dir_neg <- here('data/openai_chunks_therapies_curated_negative_answers/')
json_files1c <- list.files(json_dir_neg, pattern = ".*title.*.jsonl", full.names = TRUE)
json_files1d <- list.files(json_dir_neg, pattern = ".*summary.*.jsonl", full.names = TRUE)
json_files1e <- list.files(json_dir_neg, pattern = ".*intervention.*.jsonl", full.names = TRUE)

gpt1_neg <- grab_jsonl_results(json_files1c) %>% 
  dplyr::rename(response_title = response)
gpt2_neg <- grab_jsonl_results(json_files1d) %>% 
  dplyr::rename(response_summary = response)
gpt3_neg <- grab_jsonl_results(json_files1e) %>% 
  dplyr::rename(response_intervention = response)

gpt_merge_neg <- gpt1_neg %>% 
  left_join(gpt2_neg) %>% 
  left_join(gpt3_neg)

gpt_clin_neg <- gpt_merge_neg %>% 
  mutate(across(c(response_title, response_summary, response_intervention), ~ na_if(.x, 'NA')),
         final = pmap_chr(list(response_title, response_summary, response_intervention), ~ paste(na.omit(c(...)), collapse = "|"))) %>%
  mutate(final = str_replace_all(final, "\\s+", "")) %>%
  mutate(final = na_if(final, "")) %>%
  cleanFinalColumn(final) %>%
  mutate(binary = ifelse(is.na(final), 0, 1)) %>%
  left_join(curated, by = c('trial' = 'ID')) %>%
  dplyr::select(1,5,6,10)
table(gpt_clin_neg$binary)

#~~~~~~~~~~~~~~~~~~#
# read in/process negative NCI calls ----
#~~~~~~~~~~~~~~~~~~#
## This cleans all curated trials
trials_processed_neg <- processTrialListNeg(studies = studies_ref, trial_list = studies_guess_negative)
trials_processed_res_neg <- trials_processed_neg$final_result
pancan_res_neg <- trials_processed_neg$pancan

merge_first <- left_join(trials_processed_res_neg, gpt_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res_neg, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

gpt_final_neg <- nciFinalMergeGPT(merge_second)

table(gpt_final_neg$binary)

#~~~~~~~~~~~~~~~~~~#
# merge for final calls ----
#~~~~~~~~~~~~~~~~~~#
gpt_combined_neg <- finalCallsGPT(nci_res = gpt_final_neg,
                                   clin_res = gpt_clin_neg)

table(gpt_combined_neg$binary_collapsed)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MERGE AND SAVE ===================== ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
merged <- bind_rows(
  gpt_combined %>% mutate(groundTruth = 1),
  gpt_combined_neg %>% mutate(groundTruth = 0)
)
saveRDS(merged, here('output/results/therapy_guess_merged_pos_neg_gpt4o_tempzero_update.rds'))
