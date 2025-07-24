#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(jsonlite)
library(tidyverse)
library(data.table)

#~~~~~~~~~~~~~~~~~~#
# chunk nci data ----
#~~~~~~~~~~~~~~~~~~#
nci <- readRDS(here('data/processed/nci_drug_definitions.rds'))

# have to chunk for API
chunk_size <- 40
num_chunks <- ceiling(nrow(nci) / chunk_size)

df_chunks <- split(nci, rep(1:num_chunks, each = chunk_size, length.out = nrow(nci)))

path <- here('data/openai_chunks_data/')
for (i in seq_along(df_chunks)) {
  chunk <- df_chunks[[i]]
  file_name <- sprintf("nci_drug_definitions_chunk_%d.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_custom_jsonl(chunk, file_name_full)
}

# non-batch ----
write_non_batch_jsonl <- function(df, file) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, 3])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a scientific researcher, tasked with determining the target gene of many thousands of drugs and compounds."
      ),
      list(
        role = "user",
        content = paste0(
          "Given a paragraph of text that describes the mechanism of action of the test drug or compound, return the target gene(s) for that drug or compound.",
          "The target gene(s) should be the specific gene(s) that the drug or compound directly interacts with or modulates to exert its effect. ",
          "Key phrases that indicate the target would be, for example, 'binds to' or 'targeting'. ",
          "If the drug or compound targets a receptor, return the gene symbol for the receptor, not the ligand. ",
          "We only are interested in the direct target(s) of the drug or compound. ",
          "Do not report other genes involved after the drug binds the primary target. ",
          "If there is no acceptable answer, return NA as the target name. ",
          "If possible, answer with only the gene name. ",
          "If there are multiple possible targets, you can also return the answer as each gene name separated by '|' without any spaces. ",
          "The gene name should be the correct HUGO gene symbol or the most common alias. ",
          "Respond with only the target gene and no additional explantion or text. ",
          "Do not explain your reasoning. ",
          "The paragraph to extract the target from is as follows: ", definition
        )
      )
    )
    
    record <- list(
      id = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

path <- here('data/openai_chunks_data/')
for (i in seq_along(df_chunks)) {
  chunk <- df_chunks[[i]]
  file_name <- sprintf("nci_drug_definitions_chunk_%d_non_batch.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_non_batch_jsonl(chunk, file_name_full)
}

#~~~~~~~~~~~~~~~~~~#
# target guess nci ----
#~~~~~~~~~~~~~~~~~~#
# RUN 01B_TARGET_GUESS_LLM_GPT4O_NONBATCH_LOOP_PT1.PY CHUNK 1 HERE ----

#~~~~~~~~~~~~~~~~~~#
# read in results and convert to dataframe ----
#~~~~~~~~~~~~~~~~~~#
json_dir <- here('data/openai_chunks_targets/')  # update this path
json_files <- list.files(json_dir, pattern = "*_chunk_[0-9]+\\.jsonl", full.names = TRUE)

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

tmp <- grab_jsonl_results(json_files)
tmp <- tmp %>% 
  mutate(id = as.integer(id))

tmp <- left_join(nci, tmp, by = c('termId' = 'id'))
# saveRDS(object = tmp, file = here('data/processed/nci_dd_target_guess_LLM_gpt4o.rds'))

tmp <- readRDS(here('data/processed/nci_dd_target_guess_LLM_gpt4o.rds'))
to_refine <- tmp %>% 
  mutate(response = str_replace_all(response, " ", "")) %>%
  mutate(response = ifelse(as.character(response) == "NA", NA, as.character(response)))

#~~~~~~~~~~~~~~~~~~#
# generate jsonl for refining answers ----
#~~~~~~~~~~~~~~~~~~#
df_chunks_refine <- split(tmp, rep(1:num_chunks, each = chunk_size, length.out = nrow(nci)))

write_non_batch_refine_jsonl <- function(df, file) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, 4])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a scientific researcher, tasked with determining gene names contained within of many thousands of provided strings."
      ),
      list(
        role = "user",
        content = paste0(
          "Given a string that can contain gene names in addition to other text, your job is to identify the gene names.",
          "If there is no gene name in the string, return NA. If there is no string provided, or if the string is empty, return NA.",
          "If possible, answer with only the gene names.",
          "If there are multiple gene names, you can also return the answer as each gene name separated by '|' without any spaces.",
          "The gene name should be the correct HUGO gene symbol or the most common alias.",
          "Respond with only the gene names and no additional explantion or text. Do not explain your reasoning.",
          "The string to evaluate for genes names is: ", definition
        )
      )
    )
    
    record <- list(
      id = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

path <- here('data/openai_chunks_data/')
for (i in seq_along(df_chunks_refine)) {
  chunk <- df_chunks_refine[[i]]
  file_name <- sprintf("nci_drug_definitions_refine_chunk_%d_non_batch.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_non_batch_refine_jsonl(chunk, file_name_full)
}

#~~~~~~~~~~~~~~~~~~#
# refine nci ----
#~~~~~~~~~~~~~~~~~~#
# RUN 01B_TARGET_GUESS_LLM_GPT4O_NONBATCH_LOOP_PT1.PY CHUNK 2 HERE ----

#~~~~~~~~~~~~~~~~~~#
# read in results and convert to dataframe ----
#~~~~~~~~~~~~~~~~~~#
json_dir <- here('data/openai_chunks_targets/')  # update this path
json_files <- list.files(json_dir, pattern = "*refine_chunk_[0-9]+\\.jsonl", full.names = TRUE)

tmp <- grab_jsonl_results(json_files)
tmp <- tmp %>% 
  mutate(id = as.integer(id))

tmp <- left_join(to_refine, tmp, by = c('termId' = 'id'))
tmp_final <- tmp %>% 
  mutate(response_final = ifelse(response == response_refine, response, paste(response, response_refine, sep = "|")))
# saveRDS(object = tmp_final, file = here('data/processed/nci_dd_target_guess_LLM_gpt4o_refined.rds'))

#~~~~~~~~~~~~~~~~~~#
# clinical trials drugs ----
#~~~~~~~~~~~~~~~~~~#
studies <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

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
  pivot_longer(cols = -id, names_to = "variable", values_to = "drug") %>%
  filter(!is.na(drug)) %>%
  dplyr::select(id, drug)

unique_drugs <- filter_df %>% 
  distinct(drug)

tmp <- unique_drugs

write_non_batch_drugs_jsonl <- function(df, file) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, 1])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a scientific researcher, tasked with determining the target gene of many thousands of drugs and compounds from clinical trials."
      ),
      list(
        role = "user",
        content = paste0(
          "Given a drug or compound name, return the target gene(s) for that drug or compound. The target is likely described in the clinical trial description or reporting of trial results.",
          "The target gene(s) should be the specific gene(s) that the drug or compound directly interacts with or modulates to exert its effect.",
          "Key phrases that indicate the target would be, for example, 'binds to' or 'targeting'.",
          "If the drug or compound targets a receptor, return the gene symbol for the receptor, not the ligand.",
          "We only are interested in the direct target(s) of the drug or compound. Do not report other genes involved after the drug binds the primary target.",
          "If there is no acceptable answer, return NA as the target name.",
          "If possible, answer with only the gene name.",
          "If there are multiple possible targets, you can also return the answer as each gene name separated by '|' without any spaces.",
          "The gene name should be the correct HUGO gene symbol or the most common alias.",
          "Respond with only the target gene and no additional explantion or text. Do not explain your reasoning.",
          "The drug to find a target for is: ", definition
        )
      )
    )
    
    record <- list(
      drug = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

chunk_size_drug <- 40
num_chunks_drug <- ceiling(nrow(tmp) / chunk_size_drug)

df_chunks_drug <- split(tmp, rep(1:num_chunks_drug, each = chunk_size_drug, length.out = nrow(tmp)))

path <- here('data/openai_chunks_data/')
for (i in seq_along(df_chunks_drug)) {
  chunk <- df_chunks_drug[[i]]
  file_name <- sprintf("clinical_trial_drugs_chunk_%d_non_batch.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_non_batch_drugs_jsonl(chunk, file_name_full)
}

#~~~~~~~~~~~~~~~~~~#
# target guess drugs ----
#~~~~~~~~~~~~~~~~~~#
# RUN 01B_TARGET_GUESS_LLM_GPT4O_NONBATCH_LOOP_PT1.PY CHUNK 3 HERE ----

#~~~~~~~~~~~~~~~~~~#
# read in results and convert to dataframe ----
#~~~~~~~~~~~~~~~~~~#
json_dir <- here('data/openai_chunks_targets/')  # update this path
json_files <- list.files(json_dir, pattern = "*drugs_chunk_[0-9]+\\.jsonl", full.names = TRUE)

tmp <- grab_jsonl_results(json_files)

# saveRDS(object = tmp, file = here('data/processed/clinical_trial_target_guess_LLM_gpt4o_drugname.rds'))

#~~~~~~~~~~~~~~~~~~#
# generate jsonl for refining answers ----
#~~~~~~~~~~~~~~~~~~#
tmp <- readRDS(here('data/processed/clinical_trial_target_guess_LLM_gpt4o_drugname.rds'))
to_refine <- tmp %>% 
  dplyr::rename(response = response_drugname) %>% 
  mutate(response = str_replace_all(response, " ", "")) %>%
  mutate(response = ifelse(as.character(response) == "NA", NA, as.character(response)))

write_non_batch_refine_drugs_jsonl <- function(df, file) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, 2])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a scientific researcher, tasked with determining gene names contained within of many thousands of provided strings."
      ),
      list(
        role = "user",
        content = paste0(
          "Given a string that can contain gene names in addition to other text, your job is to identify the gene names.",
          "If there is no gene name in the string, return NA. If there is no string provided, or if the string is empty, return NA.",
          "If possible, answer with only the gene names.",
          "If there are multiple gene names, you can also return the answer as each gene name separated by '|' without any spaces.",
          "The gene name should be the correct HUGO gene symbol or the most common alias.",
          "Respond with only the gene names and no additional explantion or text. Do not explain your reasoning.",
          "The string to evaluate for genes names is: ", definition
        )
      )
    )
    
    record <- list(
      drug = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

df_chunks_drug_refine <- split(to_refine, rep(1:num_chunks_drug, each = chunk_size_drug, length.out = nrow(tmp)))

path <- here('data/openai_chunks_data/')
for (i in seq_along(df_chunks_drug_refine)) {
  chunk <- df_chunks_drug_refine[[i]]
  file_name <- sprintf("clinical_trial_drugs_refine_chunk_%d_non_batch.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_non_batch_refine_drugs_jsonl(chunk, file_name_full)
}

#~~~~~~~~~~~~~~~~~~#
# refine drugs ----
#~~~~~~~~~~~~~~~~~~#
# RUN 01B_TARGET_GUESS_LLM_GPT4O_NONBATCH_LOOP_PT1.PY CHUNK 4 HERE ----

#~~~~~~~~~~~~~~~~~~#
# read in results and convert to dataframe ----
#~~~~~~~~~~~~~~~~~~#
json_dir <- here('data/openai_chunks_targets/') 
json_files <- list.files(json_dir, pattern = "*drugs_refine_chunk_[0-9]+\\.jsonl", full.names = TRUE)

tmp <- grab_jsonl_results(json_files)

tmp <- left_join(to_refine, tmp)
tmp_final <- tmp %>% 
  mutate(response_final = ifelse(response == response_refine, response, paste(response, response_refine, sep = "|")))
# saveRDS(object = tmp_final, file = here('data/processed/clinical_trial_target_guess_LLM_gpt4o_drugname_refined.rds'))
