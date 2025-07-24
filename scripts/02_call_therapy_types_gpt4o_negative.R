#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(data.table)

#~~~~~~~~~~~~~~~~~~#
# read in data ----
#~~~~~~~~~~~~~~~~~~#
studies_ref <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

#~~~~~~~~~~~~~~~~~~#
# guess therapy type ----
#~~~~~~~~~~~~~~~~~~#
studies_guess <- studies_ref[studies_ref$NCT.Number %in% curated$ID, ]

final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))

studies_guess_negative <- studies_ref[!studies_ref$NCT.Number %in% curated$ID, ]
set.seed(1337)
studies_guess_negative <- studies_guess_negative %>% 
  sample_n(814)

#~~~~~~~~~~~~~~~~~~#
# chunk data ----
#~~~~~~~~~~~~~~~~~~#
write_get_therapy_temp0_jsonl <- function(df, file, column) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, column])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a clinical trial researcher with expertise in oncology and cellular therapies."
      ),
      list(
        role = "user",
        content = paste0(
          "Given the following clinical trial summary, identify all applicable types of therapy being investigated.
          Choose only from the following options:
          Radiopharmaceutical or Radioligand (return as RPT); Antibody Drug Conjugate or ADC (return as ADC);
          BITE or bispecific T cell engager (return as BITE);
          CAR-T; CAR-NK; CAR-DC; CAR-M.
          Return CAR if it is a CAR-based therapy but the specific cell type is not clear.
          In cases where there is more than one therapy, return the answer separated by the pipe character (|).
          If the therapy is not one of the above, return NA. Do not include any explanation or commentary.
          The clinical trial summary is as follows: ", definition
        )
      )
    )
    
    record <- list(
      trial = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

write_get_therapy_drug_temp0_jsonl <- function(df, file, column) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    term_id <- as.character(df[i, 1])
    definition <- as.character(df[i, column])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a clinical trial researcher with expertise in oncology and cellular therapies."
      ),
      list(
        role = "user",
        content = paste0(
          " Given therapy or therapies used in the trial, identify all applicable types of therapy being investigated.
          Choose only from the following options:
          Radiopharmaceutical or Radioligand (return as RPT); Antibody Drug Conjugate or ADC (return as ADC);
          BITE or bispecific T cell engager (return as BITE);
          CAR-T; CAR-NK; CAR-DC; CAR-M.
          Return CAR if it is a CAR-based therapy but the specific cell type is not clear.
          In cases where there is more than one therapy, return the answer separated by the pipe character (|).
          If the therapy is not one of the above, return NA. Do not include any explanation or commentary.
          The therapy or therapies to label are as follows: ", definition
        )
      )
    )
    
    record <- list(
      trial = term_id,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

chunk_size <- 40
num_chunks <- ceiling(nrow(studies_guess_negative) / chunk_size)

df_chunks <- split(studies_guess_negative, rep(1:num_chunks, each = chunk_size, length.out = nrow(studies_guess_negative)))

path <- here('data/openai_chunks_therapies_curated_negative/')
for (i in seq_along(df_chunks)) {
  chunk <- df_chunks[[i]]
  file_name <- sprintf("clinical_trial_therapy_guess_summary_chunk_%d.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_get_therapy_temp0_jsonl(chunk, file_name_full, column = 6)
}

for (i in seq_along(df_chunks)) {
  chunk <- df_chunks[[i]]
  file_name <- sprintf("clinical_trial_therapy_guess_title_chunk_%d.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_get_therapy_temp0_jsonl(chunk, file_name_full, column = 2)
}

for (i in seq_along(df_chunks)) {
  chunk <- df_chunks[[i]]
  file_name <- sprintf("clinical_trial_therapy_guess_intervention_chunk_%d.jsonl", i)
  file_name_full <- paste0(path, file_name)
  write_get_therapy_drug_temp0_jsonl(chunk, file_name_full, column = 9)
}

#~~~~~~~~~~~~~~~~~~#
# answers ----
#~~~~~~~~~~~~~~~~~~#

# RUN 02_CALL_THERAPY_TYPES_GPT40_SUBMISSION_NEGATIVE.PY HERE