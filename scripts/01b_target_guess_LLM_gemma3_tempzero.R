#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(ollamar)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
# Get responses, print count to terminal to track
get_response2_count_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a scientific researcher, tasked with determining the target gene of many thousands of drugs and compounds.
    Given a paragraph of text that describes the mechanism of action of the test drug or compound, return the target gene(s) for that drug or compound. 
    The target gene(s) should be the specific gene(s) that the drug or compound directly interacts with or modulates to exert its effect.
    Key phrases that indicate the target would be, for example, 'binds to' or 'targeting'.
    If the drug or compound targets a receptor, return the gene symbol for the receptor, not the ligand.
    We only are interested in the direct target(s) of the drug or compound. Do not report other genes involved after the drug binds the primary target.
    If there is no acceptable answer, return NA as the target name. 
    If possible, answer with only the gene name. 
    If there are multiple possible targets, you can also return the answer as each gene name separated by '|' without any spaces.
    The gene name should be the correct HUGO gene symbol or the most common alias.
    Respond with only the target gene and no additional explantion or text. Do not explain your reasoning.
    The paragraph to extract the target from is as follows:"
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

refine_response2_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a scientific researcher, tasked with determining gene names contained within of many thousands of provided strings.
    Given a string that can contain gene names in addition to other text, your job is to identify the gene names. 
    If there is no gene name in the string, return NA. If there is no string provided, or if the string is empty, return NA. 
    If possible, answer with only the gene names. 
    If there are multiple gene names, you can also return the answer as each gene name separated by '|' without any spaces.
    The gene name should be the correct HUGO gene symbol or the most common alias.
    Respond with only the gene names and no additional explantion or text. Do not explain your reasoning.
    The string to evaluate for genes names is:"
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

get_response2_count_drugname_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a scientific researcher, tasked with determining the target gene of many thousands of drugs and compounds from clinical trials.
    Given a drug or compound name, return the target gene(s) for that drug or compound. The target is likely described in the clinical trial description or reporting of trial results. 
    The target gene(s) should be the specific gene(s) that the drug or compound directly interacts with or modulates to exert its effect.
    Key phrases that indicate the target would be, for example, 'binds to' or 'targeting'.
    If the drug or compound targets a receptor, return the gene symbol for the receptor, not the ligand.
    We only are interested in the direct target(s) of the drug or compound. Do not report other genes involved after the drug binds the primary target.
    If there is no acceptable answer, return NA as the target name. 
    If possible, answer with only the gene name. 
    If there are multiple possible targets, you can also return the answer as each gene name separated by '|' without any spaces.
    The gene name should be the correct HUGO gene symbol or the most common alias.
    Respond with only the target gene and no additional explantion or text. Do not explain your reasoning.
    The drug to find a target for is:"
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
# read in nci defintions ----
#~~~~~~~~~~~~~~~~~~#
final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))  # created in 01_scrape_drug_dictionary

tmp <- final_dictionary_definitions

#~~~~~~~~~~~~~~~~~~#
# read in trials ----
#~~~~~~~~~~~~~~~~~~#
studies <- readRDS(here('data/processed/studies_818.rds'))
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))

#~~~~~~~~~~~~~~~~~~#
# gemma3 NCI ----
#~~~~~~~~~~~~~~~~~~#
## First pass ----
counter <- 1
start <- Sys.time()
tmp$response_zero <- sapply(tmp[[3]], function(desc) get_response2_count_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(object = tmp, file = here('data/processed/nci_dd_target_guess_LLM_gemma3-27b_tempzero.rds'))

## refine answers ----
to_refine <- readRDS(here('data/processed/nci_dd_target_guess_LLM_gemma3-27b_tempzero.rds'))
to_refine <- to_refine %>% 
  mutate(response = str_replace_all(response_zero, " ", "")) %>%
  mutate(response = ifelse(as.character(response) == "NA", NA, as.character(response))) %>% 
  dplyr::select(-4)

counter <- 1
start <- Sys.time()
to_refine$response_refine <- sapply(to_refine[[4]], function(desc) refine_response2_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time
Sys.time()

to_refine_final <- to_refine %>% 
  mutate(response_final = ifelse(response == response_refine, response, paste(response, response_refine, sep = "|")))

saveRDS(object = to_refine_final, file = here('data/processed/nci_dd_target_guess_LLM_gemma3-27b_refined_tempzero.rds'))

#~~~~~~~~~~~~~~~~~~#
# gemma3 clin trials interventions ----
#~~~~~~~~~~~~~~~~~~#
## clean studies ----
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

## first pass ----
counter <- 1
start <- Sys.time()
tmp$response_drugname <- sapply(tmp[[1]], function(desc) get_response2_count_drugname_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time
Sys.time()

saveRDS(object = tmp, file = here('data/processed/clinical_trial_target_guess_LLM_gemma3-27b_drugname_tempzero.rds'))

## refine ----
counter <- 1
start <- Sys.time()
tmp$response_refine <- sapply(tmp[[2]], function(desc) refine_response2_temp0(text = desc, model = "gemma3:27b"))
stop <- Sys.time()
time <- stop - start
time
Sys.time()

tmp <- tmp %>% 
  mutate(response_final = ifelse(response_drugname == response_refine, response_drugname, paste(response_drugname, response_refine, sep = "|")))

saveRDS(object = tmp, file = here('data/processed/clinical_trial_target_guess_LLM_gemma3-27b_drugname_refined_tempzero.rds'))