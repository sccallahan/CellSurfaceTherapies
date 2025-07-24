#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(ollamar)
library(lubridate)

#~~~~~~~~~~~~~~~~~~#
# read in data ----
#~~~~~~~~~~~~~~~~~~#
studies_ref <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))

#~~~~~~~~~~~~~~~~~~#
# use LLMs to guess therapy type ----
#~~~~~~~~~~~~~~~~~~#
studies_guess <- studies_ref[studies_ref$NCT.Number %in% curated$ID, ]

#~~~~~~~~~~~~~~~~~~#
## prompts ----
#~~~~~~~~~~~~~~~~~~#
get_therapy_count_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a clinical trial researcher with expertise in oncology and cellular therapies.
    Given the following clinical trial summary, identify all applicable types of therapy being investigated.
    Choose only from the following options:
    Radiopharmaceutical or Radioligand (return as RPT); Antibody Drug Conjugate or ADC (return as ADC);
    BITE or bispecific T cell engager (return as BITE);
    CAR-T; CAR-NK; CAR-DC; CAR-M.
    Return CAR if it is a CAR-based therapy but the specific cell type is not clear.
    In cases where there is more than one therapy, return the answer separated by the pipe character (|).
    If the therapy is not one of the above, return NA. Do not include any explanation or commentary.
    The clinical trial summary is as follows:"
    prompt <- paste0(context, text, sep = " ")
    
    # Generate a response using the specified model
    response <- generate(
      model = model,  # Replace with your chosen model
      prompt = prompt,
      output = "text",
      num_predict = -2,
      keep_alive = '10m',
      temperature = 0
    )
    
    # Print a message indicating completion for the sample
    message("Completed sample ", counter)
    
    # Increment the counter after processing
    counter <<- counter + 1
    
    return(response)
  }, error = function(e) {
    # Handle errors by returning a default response or a custom message
    message("Error in processing text: ", text)
    message("Error message: ", e$message)
    return("Error: Unable to process input")
  })
}

get_therapy_drugname_count_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a clinical trial researcher with expertise in oncology and cellular therapies.
    Given therapy or therapies used in the trial, identify all applicable types of therapy being investigated.
    Choose only from the following options:
    Radiopharmaceutical or Radioligand (return as RPT); Antibody Drug Conjugate or ADC (return as ADC);
    BITE or bispecific T cell engager (return as BITE);
    CAR-T; CAR-NK; CAR-DC; CAR-M.
    Return CAR if it is a CAR-based therapy but the specific cell type is not clear.
    In cases where there is more than one therapy, return the answer separated by the pipe character (|).
    If the therapy is not one of the above, return NA. Do not include any explanation or commentary.
    The therapy or therapies to label are as follows:"
    prompt <- paste0(context, text, sep = " ")
    
    # Generate a response using the specified model
    response <- generate(
      model = model,  # Replace with your chosen model
      prompt = prompt,
      output = "text",
      num_predict = -2,
      keep_alive = '10m',
      temperature = 0
    )
    
    # Print a message indicating completion for the sample
    message("Completed sample ", counter)
    
    # Increment the counter after processing
    counter <<- counter + 1
    
    return(response)
  }, error = function(e) {
    # Handle errors by returning a default response or a custom message
    message("Error in processing text: ", text)
    message("Error message: ", e$message)
    return("Error: Unable to process input")
  })
}

get_therapy_drugdef_count_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a clinical trial researcher with expertise in oncology and cellular therapies.
    Given the following drug description/drug definition, identify the type of therapy.
    Choose only from the following options:
    Radiopharmaceutical or Radioligand (return as RPT); Antibody Drug Conjugate or ADC (return as ADC);
    BITE or bispecific T cell engager (return as BITE);
    CAR-T; CAR-NK; CAR-DC; CAR-M.
    Return CAR if it is a CAR-based therapy but the specific cell type is not clear.
    In cases where there is more than one therapy, return the answer separated by the pipe character (|).
    If the therapy is not one of the above, return NA. Do not include any explanation or commentary.
    The drug description/definition is as follows:"
    prompt <- paste0(context, text, sep = " ")
    
    # Generate a response using the specified model
    response <- generate(
      model = model,  # Replace with your chosen model
      prompt = prompt,
      output = "text",
      num_predict = -2,
      keep_alive = '10m',
      temperature = 0
    )
    
    # Print a message indicating completion for the sample
    message("Completed sample ", counter)
    
    # Increment the counter after processing
    counter <<- counter + 1
    
    return(response)
  }, error = function(e) {
    # Handle errors by returning a default response or a custom message
    message("Error in processing text: ", text)
    message("Error message: ", e$message)
    return("Error: Unable to process input")
  })
}

#~~~~~~~~~~~~~~~~~~#
## labeling clinical trials ----
#~~~~~~~~~~~~~~~~~~#
### llama 3.3 ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-3_tempzero.rds'))

### llama 3.1 ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "llama3.1"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "llama3.1"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "llama3.1"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_llama3-1_tempzero.rds'))


### qwen ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_qwen2-5_tempzero.rds'))


### gemma ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_gemma3_tempzero.rds'))


### granite ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "granite3.2"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "granite3.2"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "granite3.2"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_granite_tempzero.rds'))


### phi ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "phi4"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "phi4"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "phi4"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_phi4_tempzero.rds'))

### medllama3-v20 ----
counter <- 1
start <- Sys.time()
studies_guess$response_zero <- sapply(studies_guess[[6]], function(desc) get_therapy_count_temp0(text = desc, model = "ahmgam/medllama3-v20"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_drug_zero <- sapply(studies_guess[[9]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "ahmgam/medllama3-v20"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
studies_guess$response_title_zero <- sapply(studies_guess[[2]], function(desc) get_therapy_count_temp0(text = desc, model = "ahmgam/medllama3-v20"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/clinical_trial_therapy_guess_LLM_medllama_tempzero.rds'))

#~~~~~~~~~~~~~~~~~~#
## labeling nci ----
#~~~~~~~~~~~~~~~~~~#
### llama 3.3 ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_llama3-3_tempzero.rds'))

### llama 3.1 ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "llama3.1:70b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "llama3.1:70b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_llama3-1_tempzero.rds'))

### qwen ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_qwen2-5_tempzero.rds'))

### gemma ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_gemma3_tempzero.rds'))

### granite ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "granite3.2"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "granite3.2"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_granite_tempzero.rds'))

### phi ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "phi4"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "phi4"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_phi4_tempzero.rds'))

### medllama3-v20 ----
counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_zero <- sapply(final_dictionary_definitions[[3]], function(desc) get_therapy_drugdef_count_temp0(text = desc, model = "ahmgam/medllama3-v20"))
stop <- Sys.time()
time <- stop - start
time

counter <- 1
start <- Sys.time()
final_dictionary_definitions$response_drug_zero <- sapply(final_dictionary_definitions[[2]], function(desc) get_therapy_drugname_count_temp0(text = desc, model = "ahmgam/medllama3-v20"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = studies_guess, file = here('data/processed/nci_therapy_guess_LLM_medllama_tempzero.rds'))