#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(ollamar)
library(org.Hs.eg.db)
library(limma)

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

# llama 3.3
tmp_llama <- tmp
start <- Sys.time()
tmp_llama$surface1 <- sapply(tmp_llama[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

start <- Sys.time()
tmp_llama$surface2 <- sapply(tmp_llama[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

start <- Sys.time()
tmp_llama$surface3 <- sapply(tmp_llama[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "llama3.3"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(tmp_llama, here('output/results/curated_cell_surface_response_llama3-3_tempzero.rds'))

# qwen
tmp_qwen <- tmp
start <- Sys.time()
tmp_qwen$surface1 <- sapply(tmp_qwen[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

start <- Sys.time()
tmp_qwen$surface2 <- sapply(tmp_qwen[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

start <- Sys.time()
tmp_qwen$surface3 <- sapply(tmp_qwen[[1]], function(desc) call_cell_surface_engi_temp0(text = desc, model = "qwen2.5:72b"))
stop <- Sys.time()
time <- stop - start
time

saveRDS(tmp_qwen, here('output/results/curated_cell_surface_response_qwen2-5_tempzero.rds'))

# gemma
tmp_gemma <- tmp
for (n in 1:3) {
  col_name <- paste0("surface", n)
  tmp_gemma[[col_name]] <- sapply(tmp_gemma[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "gemma3:27b")
  })
}
saveRDS(tmp_gemma, here('output/results/curated_cell_surface_response_gemma3-27b_tempzero.rds'))

# phi4
tmp_phi4 <- tmp
for (n in 1:3) {
  col_name <- paste0("surface", n)
  tmp_phi4[[col_name]] <- sapply(tmp_phi4[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "phi4")
  })
}
saveRDS(tmp_phi4, here('output/results/curated_cell_surface_response_phi4_tempzero.rds'))

# granite
tmp_granite <- tmp
for (n in 1:3) {
  col_name <- paste0("surface", n)
  tmp_granite[[col_name]] <- sapply(tmp_granite[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "granite3.2")
  })
}
saveRDS(tmp_granite, here('output/results/curated_cell_surface_response_granite_tempzero.rds'))

# llama 3.1
tmp_llama3_1 <- tmp
for (n in 1:3) {
  col_name <- paste0("surface", n)
  tmp_llama3_1[[col_name]] <- sapply(tmp_llama3_1[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "llama3.1:70b")
  })
}
saveRDS(tmp_llama3_1, here('output/results/curated_cell_surface_response_llama3-1_tempzero.rds'))

#~~~~~~~~~~~~~~~~~~#
# all genes ----
#~~~~~~~~~~~~~~~~~~#
all_genes <- suppressWarnings(suppressMessages(read_delim(file = here('data/raw/hgnc_complete_set.txt'), delim = '\t')))
gene_search <- all_genes %>% 
  filter(locus_group == 'protein-coding gene') %>% 
  dplyr::select(symbol) %>% # masked with loading org.hs.eg.db
  distinct(.)
# saveRDS(gene_search, file = here('data/raw/allgenes.rds'))
# llama
all_genes_llama <- gene_search
for (n in 1:3) {
  col_name <- paste0("surface", n)
  all_genes_llama[[col_name]] <- sapply(all_genes_llama[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "llama3.3")
  })
}
saveRDS(all_genes_llama, here('output/results/allgenes_cell_surface_response_llama3-3_tempzero.rds'))

# gemma
all_genes_gemma <- gene_search
for (n in 1:3) {
  col_name <- paste0("surface", n)
  all_genes_gemma[[col_name]] <- sapply(all_genes_gemma[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "gemma3:27b")
  })
}
saveRDS(all_genes_gemma, here('output/results/allgenes_cell_surface_response_gemma3-27b_tempzero.rds'))

# qwen
all_genes_qwen <- gene_search
for (n in 1:3) {
  col_name <- paste0("surface", n)
  all_genes_qwen[[col_name]] <- sapply(all_genes_qwen[[1]], function(desc) {
    call_cell_surface_engi_temp0(text = desc, model = "qwen2.5:72b")
  })
}
saveRDS(all_genes_qwen, here('output/results/allgenes_cell_surface_response_qwen2-5_tempzero.rds'))