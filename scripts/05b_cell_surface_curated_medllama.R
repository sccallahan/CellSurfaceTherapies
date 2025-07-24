#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(readr)
library(httr)
library(jsonlite)
library(ollamar)
library(lubridate)
library(tidyr)
library(stringr)

#~~~~~~~~~~~~~~~~~~#
# read in data ----
#~~~~~~~~~~~~~~~~~~#
curated <- read_tsv(file = '/mnt/scratch/group/szhao264/callahan/projects/trial_selector/data/raw/supplemental_trials_edit.tsv')
curated_split <- curated %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))

#~~~~~~~~~~~~~~~~~~#
# llm function ----
#~~~~~~~~~~~~~~~~~~#
call_cell_surface_engi_temp0 <- function(text, model = "llama3.3") {
  tryCatch({
    # provide context
    context <- "You are a scientific researcher tasked with determining whether a given gene’s protein product is present on the cell surface.
    A protein is considered a cell surface protein if it is located on the external side of the plasma membrane
    (also referred to as the cell membrane) at any point.
    This includes transmembrane proteins, GPI-anchored proteins, proteins with transient surface expression, and proteins
    that may be present on the cell surface under specific conditions or in particular cell types.
    Proteins associated only with internal cellular membranes (e.g., Golgi, ER) are not considered cell surface proteins
    unless they also appear on the plasma membrane. If the protein is present on the cell surface (even briefly or in specific contexts),
    respond with 1. If it is not present on the cell surface, respond with 0.
    Provide only the number (1 or 0) with no additional text or explanation. The gene in question is:"
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
    return(response)
  }, error = function(e) {
    # Handle errors by returning a default response or a custom message
    message("Error in processing text: ", text)
    message("Error message: ", e$message)
    return("Error: Unable to process input")
  })
}

#~~~~~~~~~~~~~~~~~~#
# curated trials ----
#~~~~~~~~~~~~~~~~~~#
# only need to call unique targets
tmp <- curated_split %>% 
  distinct(Target)

# gemma
tmp_medllama <- tmp
for (n in 1:3) {
  col_name <- paste0("surface", n)
  tmp_medllama[[col_name]] <- sapply(tmp_medllama[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "ahmgam/medllama3-v20")
  })
}
saveRDS(tmp_medllama,
        '/mnt/scratch/group/szhao264/callahan/projects/trial_selector/output/results/curated_cell_surface_response_medllama_tempzero.rds')
