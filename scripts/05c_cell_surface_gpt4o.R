#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(ollamar)
library(org.Hs.eg.db)
library(limma)
library(jsonlite)

#~~~~~~~~~~~~~~~~~~#
# read in curated data as test ----
#~~~~~~~~~~~~~~~~~~#
curated <- read_tsv(file = here('data/raw/supplemental_trials_edit.tsv'))
curated_split <- curated %>% 
  separate_rows(Target, sep = ",") %>% 
  mutate(Target = str_trim(Target))

#~~~~~~~~~~~~~~~~~~#
# llm function ----
#~~~~~~~~~~~~~~~~~~#
call_cell_surface_engi_temp0 <- function(text, model = "llama3.1:70b") {
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

write_curated_jsonl <- function(df, file) {
  con <- file(file, open = "w")
  on.exit(close(con))
  
  for (i in seq_len(nrow(df))) {
    gene <- as.character(df[i, 1])
    
    messages <- list(
      list(
        role = "system",
        content = "You are a scientific researcher tasked with determining whether a given gene’s protein product is present on the cell surface."
      ),
      list(
        role = "user",
        content = paste0(
          "A protein is considered a cell surface protein if it is located on the external side of the plasma membrane (also referred to as the cell membrane) at any point.",
          "This includes transmembrane proteins, GPI-anchored proteins, proteins with transient surface expression, and proteins that may be present on the cell surface under specific conditions or in particular cell types.",
          "Proteins associated only with internal cellular membranes (e.g., Golgi, ER) are not considered cell surface proteins unless they also appear on the plasma membrane.",
          "If the protein is present on the cell surface (even briefly or in specific contexts), respond with 1.",
          "If it is not present on the cell surface, respond with 0.",
          "Provide only the number (1 or 0) with no additional text or explanation.",
          "The gene in question is: ", gene
        )
      )
    )
    
    record <- list(
      id = gene,
      messages = messages
    )
    
    json_line <- toJSON(record, auto_unbox = TRUE)
    writeLines(json_line, con)
  }
}

write_curated_jsonl(df = tmp, file = here('data/raw/cell_surface_curated.jsonl'))

#~~~~~~~~~~~~~~~~~~#
# all genes ----
#~~~~~~~~~~~~~~~~~~#
## All genes -- 
all_genes <- suppressWarnings(suppressMessages(read_delim(file = here('data/raw/hgnc_complete_set.txt'), delim = '\t')))
gene_search <- all_genes %>% 
  filter(locus_group == 'protein-coding gene') %>% 
  dplyr::select(symbol) %>%
  distinct(.)

write_curated_jsonl(df = gene_search, file = here('data/raw/cell_surface_allgenes.jsonl'))
