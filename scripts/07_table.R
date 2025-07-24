#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(data.table)
library(jsonlite)
library(org.Hs.eg.db)
library(limma)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
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

processTrialList <- function(studies){
  intervention_types <- c('DRUG', 'BIOLOGICAL', 'RADIATION', 'GENETIC', 'COMBINATION_PRODUCT')
  
  pancan <- studies %>% 
    dplyr::select(c(NCT.Number, Study.Title, Study.URL, Study.Status, Brief.Summary, Conditions, Interventions)) %>% 
    filter(str_detect(Interventions, paste0(intervention_types, collapse = "|")))
  
  possible_intervention_strings <- c('DRUG:', 'BIOLOGICAL:', 'RADIATION:', 'GENETIC:', 'COMBINATION_PRODUCT:')
  pattern <- paste0("(", paste(possible_intervention_strings, collapse = "|"), ") (.*)(?:\\||$)")
  
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
  nci <- nci_dd
  model_res <- model_res
  nciDD_targets <- left_join(nci, model_res)
  
  # will need to fix NAs
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
    dplyr::select(1,4,5) # keeps NCT, Therapy Type, Binary yes/no
}

finalCallsGPT <- function(nci_res, clin_res){
  final_calls <- nci_res %>%
    dplyr::rename(
      final_nci = final,
      binary_nci = binary
    ) %>%
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
    dplyr::select(1,9,10) %>% # keeps NCT, Therapy Type, Binary Yes/no
    dplyr::rename(therapy_type = collapsed,
                   therapy_type_binary = binary_collapsed)
}

cleanttdDB <- function(database){
  df_new_rows <- database %>%
    rowwise() %>%
    mutate(
      drug_name = list(c(
        str_extract(drug_name, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}"),  # Extract content
        str_remove(drug_name, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}")  # Remove the extracted content from original string
      ))
    ) %>%
    unnest(cols = drug_name) %>%  # Unnest the list to create two rows
    ungroup() %>%   # Ensure ungrouping after rowwise operations
    filter(!is.na(drug_name))
  
  
  df_new_rows2 <- df_new_rows %>%
    rowwise() %>%
    mutate(
      drug_alternate = list(c(
        str_extract(drug_alternate, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}"),  # Extract content
        str_remove(drug_alternate, "\\(.*?\\)|\\[.*?\\]|\\{.*?\\}")  # Remove the extracted content from original string
      ))
    ) %>%
    unnest(cols = drug_alternate) %>%  # Unnest the list to create two rows
    ungroup() %>%  # Ensure ungrouping after rowwise operations
    filter(!is.na(drug_alternate))
  
  ttdDB_lookup_mod <- df_new_rows2 %>% 
    mutate(drug_name = tolower((drug_name))) %>% 
    mutate(drug_alternate = tolower((drug_alternate))) %>% 
    distinct()
  
  ttdDB_lookup_mod <- ttdDB_lookup_mod %>% 
    mutate(drug_name = gsub("^anti ", "", drug_name)) %>% 
    mutate(drug_name =  gsub("[[:punct:]]", "", drug_name)) %>% 
    mutate(drug_name =  gsub(" ", "", drug_name)) %>% 
    mutate(drug_name =  gsub("antibodydrugconjugate", "", drug_name)) %>%
    mutate(drug_alternate = gsub("^anti ", "", drug_alternate)) %>% 
    mutate(drug_alternate = gsub("[[:punct:]]", "", drug_alternate)) %>% 
    mutate(drug_alternate =  gsub(" ", "", drug_name)) %>% 
    mutate(drug_alternate =  gsub("antibodydrugconjugate", "", drug_name))
  ttdDB_lookup_mod <- ttdDB_lookup_mod %>% 
    distinct()
}

processNCITargetCall <- function(nci_refined_calls){
  df_new_rows <- nci_refined_calls %>%
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
  
  gpt_refined_mod <- df_new_rows %>% 
    mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
    mutate(drug = tolower(drug))
  
  gpt_refined_mod <- gpt_refined_mod %>% 
    mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
    separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
    mutate(drug = gsub("^antibody ", "", drug)) %>% 
    mutate(drug = gsub("^anti ", "", drug)) %>% 
    mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
    mutate(drug = gsub(" ", "", drug)) %>% 
    mutate(drug = gsub("antibodydrugconjugate", "", drug))
}

processTrialDrugs <- function(trial_drugs_merged){
  df_new_rows <- trial_drugs_merged %>%
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
  
  trial_drugs_merged_mod <- df_new_rows %>% 
    mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
    mutate(drug = tolower(drug))
  
  trial_drugs_merged_mod <- trial_drugs_merged_mod %>% 
    mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
    separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
    mutate(drug = gsub("^antibody ", "", drug)) %>% 
    mutate(drug = gsub("^anti ", "", drug)) %>% 
    mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
    mutate(drug = gsub(" ", "", drug)) %>% 
    mutate(drug = gsub("antibodydrugconjugate", "", drug))
  
  trial_drugs_merged_mod <- trial_drugs_merged_mod %>%
    distinct() %>% 
    mutate(response_final = na_if(response_final, 'NA'))
  
  match_refine_trialdrug <- trial_drugs_merged_mod %>%
    dplyr::rename(nci_dd_target = response_final) %>% 
    left_join(., ttdDB_lookup_mod, by = c('drug' = 'drug_name'), relationship = 'many-to-many') %>%
    distinct() %>% 
    left_join(., ttdDB_lookup_mod, by = c('drug' = 'drug_alternate'), relationship = 'many-to-many') %>% 
    distinct() %>% 
    mutate(ttdDB_target = ifelse(is.na(target.x), target.y, target.x)) %>% 
    dplyr::select(-c(10,11, 13:16)) %>% 
    distinct() %>% 
    mutate(final_targets = ifelse(is.na(nci_dd_target) | is.na(ttdDB_target),
                                  coalesce(nci_dd_target, ttdDB_target),
                                  paste(nci_dd_target, ttdDB_target, sep = "|"))) %>%
    mutate(final_targets = str_replace_all(final_targets, "; ", "|")) %>%
    mutate(final_targets = gsub("Candi TMP1", 'TYMS', final_targets))
  
  known_trials_trialdrug_w_ttddb <- match_refine_trialdrug %>% 
    filter(!is.na(final_targets))
}

collapseTargetCalls <- function(df, NCTNumber, targets) {
  df %>%
    group_by(.data[[NCTNumber]]) %>%
    summarise(
      temp = {
        vals <- paste(.data[[targets]], collapse = "|")
        vals %>%
          str_split("\\|") %>%
          unlist() %>%
          unique() %>%
          sort() %>%
          paste(collapse = "|")
      },
      .groups = "drop"
    )
}

#~~~~~~~~~~~~~~~~~~#
# make cell surface literature calls ----
#~~~~~~~~~~~~~~~~~~#
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
# read in studies, nci dd, target calls ----
#~~~~~~~~~~~~~~~~~~#
# studies, drug dictionary parts
studies <- readRDS(here('data/processed/studies_20250525.rds'))
final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))
nci_dd <- readRDS(here('data/processed/nci_dd.rds'))

# target calls + db
trial_drugs <- readRDS(here('data/processed/clinical_trial_target_guess_LLM_gpt4o_drugname_refined.rds'))
gpt_refined <- readRDS(file = here('data/processed/nci_dd_merged_LLM_gpt4o_refined.rds'))
ttdDB_lookup <- readRDS(here('data/processed/ttdDB_lookup_final.rds'))

#~~~~~~~~~~~~~~~~~~#
# read/process NCI therapy assignment ----
#~~~~~~~~~~~~~~~~~~#
# grab files
json_dir_nci_type <- here('data/openai_chunks_therapies_nci_answers/')
json_files1b <- list.files(json_dir_nci_type, pattern = ".*drug.*.jsonl", full.names = TRUE)
json_files2b <- list.files(json_dir_nci_type, pattern = ".*definition.*.jsonl", full.names = TRUE)

gpt1_nci <- grab_jsonl_results(json_files1b) %>% 
  dplyr::rename(response_drug = response) %>% 
  dplyr::rename(termId = trial) %>% 
  mutate(termId = as.integer(termId))
gpt2_nci <- grab_jsonl_results(json_files2b) %>% 
  dplyr::rename(response_def = response) %>% 
  dplyr::rename(termId = trial) %>% 
  mutate(termId = as.integer(termId))

# process interventions
trials_processed <- processTrialList(studies = studies)
trials_processed_res <- trials_processed$final_result
pancan_res <- trials_processed$pancan

# merge back to dictionary definitions
gpt_merge <- left_join(final_dictionary_definitions, gpt1_nci) %>% 
  left_join(gpt2_nci)

# merge back to nci dd and process
gpt_nci <- processNCIResultsGPT(nci_dd = nci_dd,
                                model_res = gpt_merge)

# merge nci dd + clinical trials
merge_first <- left_join(trials_processed_res, gpt_nci, by = 'drug', relationship = 'many-to-many')
merge_first_clean <- merge_first %>% 
  distinct()

merge_second <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., merge_first_clean, by = 'id')

gpt_final_nci_therapy <- nciFinalMergeGPT(merge_second)

#~~~~~~~~~~~~~~~~~~#
# read/process clinical trial therapy assignment ----
#~~~~~~~~~~~~~~~~~~#
json_dir_clin_type <- here('data/openai_20250525_therapy_guess_answers/')
json_files1 <- list.files(json_dir_clin_type, pattern = ".*title.*.jsonl", full.names = TRUE)
json_files2 <- list.files(json_dir_clin_type, pattern = ".*summary.*.jsonl", full.names = TRUE)
json_files3 <- list.files(json_dir_clin_type, pattern = ".*intervention.*.jsonl", full.names = TRUE)

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
  dplyr::select(1,5,6) # keep NCT, Therapy type, binary yes/no

#~~~~~~~~~~~~~~~~~~#
# MERGE therapy assignment ----
#~~~~~~~~~~~~~~~~~~#
gpt_combined <- finalCallsGPT(nci_res = gpt_final_nci_therapy,
                              clin_res = gpt_clin)

#~~~~~~~~~~~~~~~~~~#
# process NCI gene target ----
#~~~~~~~~~~~~~~~~~~#
# process database + results from target calls
ttdDB_lookup_mod <- cleanttdDB(database = ttdDB_lookup)
gpt_refined_mod <- processNCITargetCall(gpt_refined)

# match back to clinical trial and integrate ttddb
match_a <- left_join(trials_processed_res, gpt_refined_mod, by = 'drug', relationship = 'many-to-many')
match_b <- match_a %>% 
  distinct()

match_refine <- match_b %>%
  dplyr::select(c(1:2, 5, 9, 11, 13, 17, 18, 21)) %>% 
  dplyr::select(c(1,2,4,8,9)) %>% # magic numbers that give id, drug, prettyURL, definition, response_final
  distinct() %>% 
  dplyr::rename(nci_dd_target = response_final) %>% 
  left_join(., ttdDB_lookup_mod, by = c('drug' = 'drug_name'), relationship = 'many-to-many') %>%
  distinct() %>% 
  dplyr::select(1:5, 7) %>% 
  distinct() %>% 
  mutate(ttdDB_target = target) %>% 
  mutate(final_targets = ifelse(is.na(nci_dd_target) | is.na(ttdDB_target),
                                coalesce(nci_dd_target, ttdDB_target),
                                paste(nci_dd_target, ttdDB_target, sep = "|"))) %>%
  mutate(final_targets = str_replace_all(final_targets, "; ", "|")) %>%
  mutate(final_targets = gsub("Candi TMP1", 'TYMS', final_targets)) %>% 
  dplyr::select(1,2,8) %>% 
  distinct() %>% 
  group_by(id, drug) %>%
  mutate(final_targets_merge = paste(final_targets, collapse = "|")) %>%
  ungroup() %>% 
  dplyr::select(1,2,4) %>% 
  distinct()

match_refine2 <- match_b %>%
  dplyr::select(c(1:2, 5, 9, 11, 13, 17, 18, 21)) %>% 
  dplyr::select(c(1,2,4,8,9)) %>% 
  distinct() %>% 
  dplyr::rename(nci_dd_target = response_final) %>% 
  left_join(., ttdDB_lookup_mod, by = c('drug' = 'drug_alternate'), relationship = 'many-to-many') %>% 
  distinct() %>% 
  dplyr::select(1:5, 7) %>% 
  distinct() %>% 
  mutate(ttdDB_target = target) %>% 
  mutate(final_targets = ifelse(is.na(nci_dd_target) | is.na(ttdDB_target),
                                coalesce(nci_dd_target, ttdDB_target),
                                paste(nci_dd_target, ttdDB_target, sep = "|"))) %>%
  mutate(final_targets = str_replace_all(final_targets, "; ", "|")) %>%
  mutate(final_targets = gsub("Candi TMP1", 'TYMS', final_targets)) %>% 
  dplyr::select(1,2,8) %>% 
  distinct() %>% 
  group_by(id, drug) %>%
  mutate(final_targets_merge = paste(final_targets, collapse = "|")) %>%
  ungroup() %>% 
  dplyr::select(1,2,4) %>% 
  distinct()

# recombine
match_refine3 <- left_join(match_refine, match_refine2) %>% 
  mutate(final_targets_merge = na_if(final_targets_merge, 'NA'))

# merge back to trials
pancan_merge <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., match_refine3, by = 'id')

# clean out NA
nci_targets_merged <- pancan_merge %>% 
  filter(!is.na(final_targets_merge))

#~~~~~~~~~~~~~~~~~~#
# process clinical trial gene target assignment ----
#~~~~~~~~~~~~~~~~~~#
trial_drugs_intermed <- left_join(trials_processed_res, trial_drugs, by = 'drug', relationship = 'many-to-one')
trial_drugs_merged <- rownames_to_column(pancan_res, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., trial_drugs_intermed, by = 'id')

trial_drugs_merged <-  trial_drugs_merged %>% 
  dplyr::rename(response_drugname = response) %>% 
  mutate(response_drugname = str_replace_all(response_drugname, " ", "")) %>%
  mutate(response_drugname = na_if(response_drugname, 'NA'))

clintrial_targets_merged <- processTrialDrugs(trial_drugs_merged)

#~~~~~~~~~~~~~~~~~~#
# MERGE gene target assignment ----
#~~~~~~~~~~~~~~~~~~#
nci_targets_collapsed <- collapseTargetCalls(nci_targets_merged, 'NCT.Number', 'final_targets_merge')
clintrial_targets_collapsed <- collapseTargetCalls(clintrial_targets_merged, 'NCT.Number', 'final_targets')
merged_targets <- full_join(nci_targets_collapsed, clintrial_targets_collapsed, by = "NCT.Number")

# merge columns keeping only unique strings
merged_targets_collapse <- merged_targets %>%
  mutate(targets = pmap_chr(list(temp.x, temp.y), function(x, y) {
    vals <- c(x, y)
    vals <- vals[!is.na(vals)]  # Remove NA
    vals <- unique(str_split(paste(vals, collapse = "|"), "\\|")[[1]])  # Split & dedup
    if (length(vals) == 0) NA_character_ else paste(sort(vals), collapse = "|")
  })) %>% 
  dplyr::select(1,4) # keep only NCT Number and final targets

#~~~~~~~~~~~~~~~~~~#
# MERGE therapy type and targets ----
#~~~~~~~~~~~~~~~~~~#
# initial join of therapies and targets
initial_table_res <- full_join(gpt_combined, merged_targets_collapse, by = 'NCT.Number') %>% # keep all trials
  dplyr::filter(therapy_type_binary == 1) %>%  # remove non cell surface therapy trials
  dplyr::select(-3) # remove binary column

# split genes into rows and find aliases
genes_corrected <- initial_table_res %>%
  separate_rows(targets, sep = "\\|") %>% 
  mutate(targets = na_if(targets, 'NA')) %>%  # if our search found multiple hit in nci
  mutate(targets = toupper(targets)) %>% 
  mutate(targets = str_replace(targets, 'CLDN18.2', 'CLDN18')) %>% 
  mutate(targets = str_replace(targets, 'PSMA', 'FOLH1')) %>% 
  mutate(limma = alias2SymbolTable(targets, species = 'Hs')) %>% 
  distinct()

get_multi_map <- function(alias) {
  # get each gene symbol
  symbols <- alias2Symbol(alias)
  # check for multi-mapping
  if (length(symbols) > 1) {
    return(paste(symbols, collapse = "|"))  # collapse
  } else {
    return(symbols)
  }
}
multi_map_results <- sapply(genes_corrected$targets, get_multi_map)
multi_map_df <- data.frame(
  alias = names(multi_map_results),
  multi_map = as.character(multi_map_results)
)

genes_corrected_final <- cbind(genes_corrected, multi_map_df) %>%
  separate_rows(multi_map, sep = "\\|")

table_stage2 <- genes_corrected_final %>%
  unique(.) %>%
  dplyr::rename(target_algo = multi_map) %>%
  mutate(target_algo = na_if(target_algo, 'character(0)')) %>%
  filter(!is.na(target_algo)) %>% 
  dplyr::select(1,2,6) # keeps NCT.Number, therapy type, multi-mapped targets

# filter for cell surface genes
table_cell_surface <- table_stage2 %>% 
  dplyr::filter(target_algo %in% combined_calls_df$gene)

# collapse genes back before re-join to clinical trial table
table_collapse_pre_join <- table_cell_surface %>% 
  group_by(NCT.Number) %>%
  summarize(
    therapy_type = dplyr::first(therapy_type),
    target = target_algo %>%
      str_split("\\|") %>%           
      unlist() %>%                  
      unique() %>%                  
      sort() %>%                     
      paste(collapse = "|"),       
    .groups = "drop"
  )

# rejoin to clinical trials table
table_final <- left_join(table_collapse_pre_join, studies, by = 'NCT.Number')

# save this table
write_delim(x = table_final, file = here('output/results/final_table.txt'))
