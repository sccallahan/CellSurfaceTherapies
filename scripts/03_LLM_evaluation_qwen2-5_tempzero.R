#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(parallel)
library(ollamar)
library(org.Hs.eg.db)
library(limma)


#~~~~~~~~~~~~~~~~~~#
# match qwen2.5:72b refined ----
#~~~~~~~~~~~~~~~~~~#
nciDD <- readRDS(here('data/processed/nci_dd_target_guess_LLM_qwen2-5_refined_tempzero.rds'))
nci <- readRDS(here('data/processed/nci_dd.rds'))
nciDD_targets <- left_join(nci, nciDD)

# will need to fix NAs
nciDD_targets <- nciDD_targets %>%
  mutate(across(everything(), ~ str_replace_all(., " ", ""))) %>%
  mutate(across(everything(), ~ ifelse(as.character(.) == "NA", NA, as.character(.))))

saveRDS(nciDD_targets, file = here('data/processed/nci_dd_merged_LLM_qwen2-5_refined_tempzero.rds'))

#~~~~~~~~~~~~~~~~~~#
# load target calls ----
#~~~~~~~~~~~~~~~~~~#
qwen_refined <- readRDS(file = here('data/processed/nci_dd_merged_LLM_qwen2-5_refined_tempzero.rds'))
ttdDB_lookup <- readRDS(here('data/processed/ttdDB_lookup_final.rds'))

#~~~~~~~~~~~~~~~~~~#
# load manual list ----
#~~~~~~~~~~~~~~~~~~#
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))
curated_split <- curated %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))

#~~~~~~~~~~~~~~~~~~#
# pull studies ----
#~~~~~~~~~~~~~~~~~~#
studies <- readRDS(here('data/processed/studies_818.rds'))

#~~~~~~~~~~~~~~~~~~#
# clean studies ----
#~~~~~~~~~~~~~~~~~~#
intervention_types <- c('DRUG', 'BIOLOGICAL', 'RADIATION', 'GENETIC', 'COMBINATION_PRODUCT')

pancan <- studies %>% 
  dplyr::select(c(NCT.Number, Study.Title, Study.URL, Study.Status, Brief.Summary, Conditions, Interventions)) %>% 
  filter(str_detect(Interventions, paste0(intervention_types, collapse = "|")))

possible_intervention_strings <- c('DRUG:', 'BIOLOGICAL:', 'RADIATION:', 'GENETIC:', 'COMBINATION_PRODUCT:')
pattern <- paste0("(", paste(possible_intervention_strings, collapse = "|"), ") (.*)(?:\\||$)")

# only select curated trials for analysis
pancan <- pancan %>% 
  filter(NCT.Number %in% curated$ID)
intervention_list <- str_split(string = pancan$Interventions, pattern = "\\|")

# this will maintain trial by row number in `pancan` object
df <- tibble(data = intervention_list) %>%
  unnest_wider(data, names_sep = "_")
filter_df <- df %>%
  mutate(across(everything(), ~ if_else(str_detect(.x,  paste0(intervention_types, collapse = "|")), .x, NA))) %>% # removes procedures, etc.
  mutate(across(everything(), ~ str_replace_all(., pattern = pattern, replacement = "\\2"))) %>% # keep just drug name
  mutate(id = row_number()) %>%  # keeps track of trial number by row
  pivot_longer(cols = -id, names_to = "variable", values_to = "drug") %>% # long format
  filter(!is.na(drug)) %>%
  dplyr::select(id, drug)

####
df_new_rows <- ttdDB_lookup %>%
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

####
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

####
df_new_rows <- qwen_refined %>%
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

qwen_refined_mod <- df_new_rows %>% 
  mutate(drug = str_replace(drug, "^ADC|-ADC", "")) %>% 
  mutate(drug = tolower(drug))

qwen_refined_mod <- qwen_refined_mod %>% 
  mutate(drug = str_replace_all(drug, "(injection of|for injection|cell injection|injection)", "")) %>% 
  separate_rows(drug, sep = "(,|;|/|\\+|(?<= )and(?= ))") %>% # Split by ';', '/', '+', or ' and '
  mutate(drug = gsub("^antibody ", "", drug)) %>% 
  mutate(drug = gsub("^anti ", "", drug)) %>% 
  mutate(drug = gsub("[[:punct:]]", "", drug)) %>%
  mutate(drug = gsub(" ", "", drug)) %>% 
  mutate(drug = gsub("antibodydrugconjugate", "", drug))

qwen_refined_mod <- qwen_refined_mod %>% distinct()


#~~~~~~~~~~~~~~~~~~#
# perform matches ----
#~~~~~~~~~~~~~~~~~~#
## qwen refined no fuzzy match ----
match5 <- left_join(filter_df_mod, qwen_refined_mod, by = 'drug', relationship = 'many-to-many')
match6 <- match5 %>% 
  distinct()

match_refine <- match6 %>%
  dplyr::select(c(1:2, 5, 9, 11, 13, 17, 18, 21)) %>% 
  dplyr::select(c(1,2,4,8,9)) %>% 
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

match_refine2 <- match6 %>%
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

match_refine3 <- left_join(match_refine, match_refine2) %>% 
  mutate(final_targets_merge = na_if(final_targets_merge, 'NA'))

### match back to trial list ----
pancan_mod_qwen_refined_nofuzz <- rownames_to_column(pancan, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., match_refine3, by = 'id')

known_trials_qwen_refined_nofuzz <- pancan_mod_qwen_refined_nofuzz %>% 
  filter(!is.na(final_targets_merge))

## qwen2.5:72b trial drug ----
trial_drugs <- readRDS(here('data/processed/clinical_trial_target_guess_LLM_qwen2-5_drugname_refined_tempzero.rds'))
trial_drugs_intermed <- left_join(filter_df, trial_drugs, by = 'drug', relationship = 'many-to-one')
trial_drugs_merged <- rownames_to_column(pancan, var = 'id') %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(., trial_drugs_intermed, by = 'id')
sum(curated$ID %in% trial_drugs_merged$NCT.Number)

trial_drugs_merged <-  trial_drugs_merged %>% 
  mutate(response_drugname = str_replace_all(response_drugname, " ", "")) %>%
  mutate(response_drugname = na_if(response_drugname, 'NA'))

sum(curated$ID %in% trial_drugs_merged$NCT.Number) 

## qwen2.5:72b trial drug w/ ttddb ----
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

sum(curated$ID %in% known_trials_trialdrug_w_ttddb$NCT.Number)

#~~~~~~~~~~~~~~~~~~#
# evaluate trial identification ----
#~~~~~~~~~~~~~~~~~~#

## Combo approach ----
# function to collect multi-maps
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

### qwen refined no fuzz ----
combo_compare_refined_nofuzz <- left_join(curated, known_trials_qwen_refined_nofuzz, by = c('ID' = 'NCT.Number'), relationship = 'one-to-many') %>%
  dplyr::select(1,4,5,14) %>%
  distinct()

# separate hits
genes_corrected3 <- combo_compare_refined_nofuzz %>%
  separate_rows(final_targets_merge, sep = "\\|") %>% 
  mutate(final_targets_merge = na_if(final_targets_merge, 'NA')) %>%  # if our search found multiple hit in nci
  mutate(final_targets_merge = toupper(final_targets_merge)) %>% 
  mutate(final_targets_merge = str_replace(final_targets_merge, 'CLDN18.2', 'CLDN18')) %>% 
  mutate(final_targets_merge = str_replace(final_targets_merge, 'PSMA', 'FOLH1')) %>% 
  mutate(limma = alias2SymbolTable(final_targets_merge, species = 'Hs')) %>% 
  distinct()

multi_map_results3 <- sapply(genes_corrected3$final_targets_merge, get_multi_map)
multi_map_df3 <- data.frame(
  alias = names(multi_map_results3),
  multi_map = as.character(multi_map_results3)
)

genes_corrected3 <- cbind(genes_corrected3, multi_map_df3) %>%
  separate_rows(multi_map, sep = "\\|")

# clean up data
final_compare3 <- genes_corrected3 %>%
  unique(.) %>%
  dplyr::rename(target_algo = multi_map) %>%
  mutate(target_algo = na_if(target_algo, 'character(0)')) %>%
  filter(!is.na(target_algo)) %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))


# filters for exact matches
combo_geneMatch3 <- final_compare3 %>%
  filter(str_detect(target_algo, Target))
length(unique(combo_geneMatch3$ID))

### llama trial drug w/ ttddb ----
combo_compare_trialdrug_ttddb <- left_join(curated, known_trials_trialdrug_w_ttddb, by = c('ID' = 'NCT.Number'), relationship = 'one-to-many') %>%
  dplyr::select(1,4,5,16) %>%
  distinct()

genes_corrected5 <- combo_compare_trialdrug_ttddb %>%
  separate_rows(final_targets, sep = "\\|") %>% 
  mutate(final_targets = na_if(final_targets, 'NA')) %>%  # if our search found multiple hit in nci
  mutate(final_targets = toupper(final_targets)) %>% 
  mutate(final_targets = str_replace(final_targets, 'CLDN18.2', 'CLDN18')) %>% 
  mutate(final_targets = str_replace(final_targets, 'PSMA', 'FOLH1')) %>% 
  mutate(limma = alias2SymbolTable(final_targets, species = 'Hs')) %>% 
  distinct()

multi_map_results5 <- sapply(genes_corrected5$final_targets, get_multi_map)
multi_map_df5 <- data.frame(
  alias = names(multi_map_results5),
  multi_map = as.character(multi_map_results5)
)

genes_corrected5 <- cbind(genes_corrected5, multi_map_df5) %>%
  separate_rows(multi_map, sep = "\\|")

# clean up data
final_compare5 <- genes_corrected5 %>%
  unique(.) %>%
  dplyr::rename(target_algo = multi_map) %>%
  mutate(target_algo = na_if(target_algo, 'character(0)')) %>%
  filter(!is.na(target_algo)) %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))

# filters for exact matches, also tells us how many trials we found successfully
combo_geneMatch5 <- final_compare5 %>%
  filter(str_detect(target_algo, Target))
length(unique(combo_geneMatch5$ID))

### naive combo llama refined nofuzz + trial drug w/ ttddb ----
combo_naive2 <- final_compare5 %>% 
  mutate(final_targets_merge = final_targets) %>% 
  dplyr::select(-final_targets) %>% 
  dplyr::select(1:3,7,4:6)
combo_naive2 <- rbind(combo_naive2, final_compare3)
combo_naive2 <- combo_naive2 %>% 
  distinct()

combo_geneMatch_naive2 <- combo_naive2 %>%
  filter(str_detect(target_algo, Target))
length(unique(combo_geneMatch_naive2$ID))

#~~~~~~~~~~~~~~~~~~#
# tally correct matches ----
#~~~~~~~~~~~~~~~~~~#
# qwen2.5 + ttddb
missed_trials_qwen_nofuzz <- curated_split[!curated_split$ID %in% combo_geneMatch3$ID,]
length(unique(missed_trials_qwen_nofuzz$ID))
matched_trials_qwen_nofuzz <- curated_split[curated_split$ID %in% combo_geneMatch3$ID,]
length(unique(matched_trials_qwen_nofuzz$ID))
# write_tsv(missed_trials_qwen_nofuzz, here('output/results/missed_trials_qwen2-5_nofuzz_tempzero.tsv'))
# write_tsv(matched_trials_qwen_nofuzz, here('output/results/matched_trials_qwen2-5_nofuzz_tempzero.tsv'))

# trial drug + ttddb
missed_trials_drug_ttddb <- curated_split[!curated_split$ID %in% combo_geneMatch5$ID,]
length(unique(missed_trials_drug_ttddb$ID))
matched_trials_drug_ttddb <- curated_split[curated_split$ID %in% combo_geneMatch5$ID,]
length(unique(matched_trials_drug_ttddb$ID))
# write_tsv(missed_trials_drug_ttddb, here('output/results/missed_trials_drug_qwen2-5_ttddb_tempzero.tsv'))
# write_tsv(matched_trials_drug_ttddb, here('output/results/matched_trials_drug_qwen2-5_ttddb_tempzero.tsv'))

# combined qwen/ttddb + trial drug/ttddb
missed_trials_qwen_drug_ttddb <- curated_split[!curated_split$ID %in% combo_geneMatch_naive2$ID,]
length(unique(missed_trials_qwen_drug_ttddb$ID))
matched_trials_qwen_drug_ttddb <- curated_split[curated_split$ID %in% combo_geneMatch_naive2$ID,]
length(unique(matched_trials_qwen_drug_ttddb$ID))
# write_tsv(missed_trials_qwen_drug_ttddb, here('output/results/missed_trials_qwen2-5_drug_ttddb_tempzero.tsv'))
# write_tsv(matched_trials_qwen_drug_ttddb, here('output/results/matched_trials_qwen2-5_drug_ttddb_tempzero.tsv'))
