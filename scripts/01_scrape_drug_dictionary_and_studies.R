#~~~~~~~~~~~~~~~~~~#
# libraries ----
#~~~~~~~~~~~~~~~~~~#
library(here)
library(tidyverse)
library(httr)
library(jsonlite)
library(cli)

#~~~~~~~~~~~~~~~~~~#
# functions ----
#~~~~~~~~~~~~~~~~~~#
# These functions are taken from these repos:
# https://github.com/meerapatelmd/skyscraper
# https://github.com/meerapatelmd/secretary

ex_crawl_delay <- function(crawl_delay) {
  sp <- cli::make_spinner(which = "dots10",
                          template = "    {spin}")
  lapply(1:(crawl_delay*100), function(x) { sp$spin(); Sys.sleep(0.01) })
  sp$finish()
}

get_nci_dd <- function(crawl_delay = 5, size = 10000, letters, verbose = TRUE) {
  
  lttrs <- c(LETTERS, "%23")
  
  if (!missing(letters)) {
    lttrs <- toupper(letters)
  }
  
  output <- list()
  for (i in seq_along(lttrs)) {
    
    if (verbose) {
      secretary::typewrite_progress(iteration = i,
                                    total = length(lttrs))
      secretary::typewrite(secretary::blueTxt("Letter:"), lttrs[i])
    }
    
    ex_crawl_delay(crawl_delay = crawl_delay)
    response <- httr::GET(url = sprintf("https://webapis.cancer.gov/drugdictionary/v1/Drugs/expand/%s?", lttrs[i]),
                          query = list(size = size))
    parsed <- httr::content(x = response,
                            as = "text",
                            encoding = "UTF-8")
    results <- jsonlite::fromJSON(txt = parsed)
    
    
    output[[i]] <- results
    
  }
  
  names(output) <- lttrs
  
  # max_records <- results$meta$totalResults
  # batches <- ceiling(max_records/size)
  #
  # indices <-
  #         seq(from = 0,
  #             to = batches*size,
  #             by = size)
  # starting <<- indices[1:batches]
  #
  # output <- list()
  # for (i in 1:batches) {
  #         ex_crawl_delay(crawl_delay = crawl_delay)
  #         response <- httr::GET(url = "https://webapis.cancer.gov/drugdictionary/v1/Drugs",
  #                               query = list(size = size,
  #                                            from = starting[i]))
  #         parsed <- httr::content(x = response,
  #                                 as = "text",
  #                                 encoding = "UTF-8")
  #         output[[i]] <- jsonlite::fromJSON(txt = parsed)
  # }
  
  totalResults <-
    output %>%
    purrr::transpose() %>%
    purrr::pluck("meta") %>%
    purrr::transpose() %>%
    purrr::pluck("totalResults") %>%
    purrr::map(~ tibble::as_tibble_col(., column_name = "count")) %>%
    dplyr::bind_rows(.id = "letter")
  
  grandTotal <-
    output %>%
    purrr::transpose() %>%
    purrr::pluck("meta") %>%
    purrr::transpose() %>%
    purrr::pluck("totalResults") %>%
    unlist() %>%
    sum()
  
  results <-
    output %>%
    purrr::transpose() %>%
    purrr::pluck("results") %>%
    dplyr::bind_rows(.id = "letter")
  
  drugInfoSummaryLink <-
    results$drugInfoSummaryLink %>%
    dplyr::rename(uri_text = text)
  
  definition <- results$definition %>%
    dplyr::rename(html_text = text)
  
  
  aliases <-
    results$aliases %>%
    purrr::set_names(1:length(results$aliases)) %>%
    purrr::keep(~ !is.null(.)) %>%
    dplyr::bind_rows(.id = "rowid") %>%
    dplyr::mutate(rowid = as.integer(rowid)) %>%
    dplyr::select(rowid,
                  drug_type = type,
                  drug_name = name)
  
  results <-
    results %>%
    dplyr::select(-drugInfoSummaryLink,
                  -definition,
                  -aliases) %>%
    dplyr::bind_cols(drugInfoSummaryLink,
                     definition) %>%
    tibble::rowid_to_column(var = "rowid") %>%
    dplyr::left_join(aliases, by = "rowid") %>%
    dplyr::select(-rowid) %>%
    dplyr::distinct() %>%
    dplyr::rename(drug = name,
                  drug_name_type = type) %>%
    dplyr::transmute(ndd_datetime = Sys.time(),
                     letter,
                     preferredName,
                     termId,
                     drug,
                     firstLetter,
                     drug_name_type,
                     prettyUrlName,
                     nciConceptId,
                     nciConceptName,
                     uri,
                     uri_text,
                     html,
                     html_text,
                     drug_type,
                     drug_name)
  
  
  list(meta = list(grandTotal = grandTotal,
                   totalResults = totalResults),
       results = results)
  
}

drug_count <- function(size = 10000, crawl_delay = 5) {
  
  Sys.sleep(crawl_delay)
  
  response <- httr::GET(url = "https://webapis.cancer.gov/drugdictionary/v1/Drugs",
                        query = list(size = size,
                                     includeResourceTypes = "DrugTerm"))
  parsed <- httr::content(x = response,
                          as = "text",
                          encoding = "UTF-8")
  df <- jsonlite::fromJSON(txt = parsed)
  
  df$meta$totalResults
  
}

# other functions
pluck_and_bind <- function(lst) {
  tibble(
    termId = pluck(lst$termId),
    prettyUrlName = pluck(lst$prettyUrlName),
    definition = pluck(lst$definition$text)
  )
}

scrapeNCIDefinitions <- function(drug_names){
  if(is.null(drug_names)){
    stop("Must provide argument for 'drug_names' as a character vector")
  }
  base_url <- "https://webapis.cancer.gov/drugdictionary/v1/Drugs/"
  selected_results <- NULL
  for (i in seq_along(drug_names)){
    message(paste0('Working on drug: ', i))
    query_url <- paste0(
      base_url,
      drug_names[i])
    response <- httr::GET(url = query_url)
    if (status_code(response) != 200) {
      error_message <- httr::content(response, "text", encoding = "UTF-8")
      stop("Failed to fetch data: ", status_code(response), " - ", error_message)
    }
    parsed <- httr::content(x = response,
                            as = "text",
                            encoding = "UTF-8")
    results <- jsonlite::fromJSON(txt = parsed)
    tmp <- pluck_and_bind(results)
    selected_results <- bind_rows(selected_results, tmp)
  }
  return(selected_results)
}

scrapeStudieswithDates <- function(params){
  if(!is.list(params)){
    stop("Must provide argument for 'params' as a list")
  }
  # First construct a base url
  base_url <- "https://clinicaltrials.gov/api/v2/studies"
  studyList <- NULL    # Will hold the data
  i <- 1               # Will be used to iterate and edit pages 2+
  while(TRUE){
    # this will make the URL
    query_url <- paste0(
      base_url, "?",
      paste0(names(params), "=", gsub(params, pattern = " ", replacement = "+"), collapse = "&"))
    message(paste0('Conducting search using URL:', query_url))
    res <- GET(query_url)
    if (status_code(res) != 200) {
      error_message <- httr::content(res, "text", encoding = "UTF-8")
      stop("Failed to fetch data: ", status_code(res), " - ", error_message)
    }
    data <- httr::content(res, as = "parsed", type = "text/csv")
    data <- data[, -c(24,25,27,28)] # dates are causing a headaches and not very useful
    # print(head(data[, c(23, 24)]))
    
    # Process columns 23 and 24
    # Process columns 23 and 24
    for (col in c(23, 24)) {
      # Convert each column to character (even if it's already character, just to ensure)
      data[[col]] <- as.character(data[[col]])
      
      # Handle NA values separately (replace with placeholder string if needed)
      data[[col]][is.na(data[[col]])] <- "1900-01-01"  # Placeholder for missing values
      
      # Check and fix any "YYYY-MM" format (missing day)
      missing_day_format <- grepl("^\\d{4}-\\d{2}$", data[[col]])  # Detect "YYYY-MM"
      data[[col]][missing_day_format] <- paste0(data[[col]][missing_day_format], "-01")  # Add "-01" to those entries
      
      # If you just want to leave them as character strings, skip any date parsing step
      # No need for ymd() or any date parsing function here
      
      # Print to verify that it's read as a character string
      # message("First 10 values in column ", col, ": ")
      # print(head(data[[col]], 10))
    }
    
    # print(head(data[, c(23, 24)]))
    
    if(i > 1){        # pages 2+ don't have column names and instead make data the column names so we need to fix that
      tmp <- data.frame(data)
      tmp <- rbind(colnames(tmp), tmp)
      colnames(tmp) <- colnames(studyList)
    }
    else{             # The first page is fine as-is
      tmp <- data.frame(data)
    }
    studyList <- rbind(studyList, tmp)
    page_token <- res$headers$`x-next-page-token`  # grab the next page token
    # print(page_token)
    if(is.null(page_token)){  # if there's no token, break the loop
      break
    }
    else{
      params$pageToken <- page_token
    }
    i <- i + 1        # just tracks the page
  }
  return(studyList)
}

#~~~~~~~~~~~~~~~~~~#
# pull data ----
#~~~~~~~~~~~~~~~~~~#
# api entry is here: https://webapis.cancer.gov/drugdictionary/v1/Drugs/{drugname}
drug_dictionary <- get_nci_dd(verbose = FALSE)
drug_dataframe <- drug_dictionary$results
saveRDS(object = drug_dataframe, file = here('data/processed/nci_dd.rds'))

# there are duplicates based on naming, but IDs are tracked correctly
counts <- drug_count()

# these should match
length(unique(drug_dataframe$termId)) == counts

# isolate termIDs to pull definitions
to_search <- drug_dataframe %>% 
  distinct(termId, .keep_all = TRUE) %>% 
  select(termId, prettyUrlName)
id_search <- to_search$termId

#~~~~~~~~~~~~~~~~~~#
# match whole db ----
#~~~~~~~~~~~~~~~~~~#
drug_dataframe <- read_rds(here('data/processed/nci_dd.rds'))

to_search <- drug_dataframe %>% 
  distinct(termId, .keep_all = TRUE) %>% 
  select(termId, prettyUrlName)
id_search <- to_search$termId


final_dictionary_definitions <- scrapeNCIDefinitions(drug_names = id_search)
saveRDS(object = final_dictionary_definitions, file = here('data/processed/nci_drug_definitions.rds'))


final_dictionary_definitions <- readRDS(here('data/processed/nci_drug_definitions.rds'))

#~~~~~~~~~~~~~~~~~~#
# pull clinical trials ----
#~~~~~~~~~~~~~~~~~~#
params <- list(
  "format" = "csv",
  "query.cond" = "cancer OR solid tumor OR malignant neoplasm OR neoplasms OR neoplasm OR carcinoma OR neoplasms by site",
  "pageSize" = 1000
)

studies_raw <- scrapeStudieswithDates(params = params)
studies <- studies_raw %>%
  mutate(across(23:24, ~ str_remove(., "^(?i)x"))) %>%
  mutate(across(23:24, ~ str_remove(., ".\\.\\.\\.26$"))) %>%
  mutate(across(23:24, ~ ymd(.))) %>%
  filter(Start.Date != ymd('1900-01-01') | First.Posted != ymd('1900-01-1')) %>%
  filter(Start.Date <= ymd('2023-10-31') | First.Posted <= ymd('2023-10-31')) %>%
  dplyr::select(-c(23,24))

sum(curated$ID %in% studies$NCT.Number)
curated_out <- curated[!curated$ID %in% studies$NCT.Number,]
saveRDS(object = studies, file = here('data/processed/studies_818.rds'))

#~~~~~~~~~~~~~~~~~~#
# process ttdb ----
#~~~~~~~~~~~~~~~~~~#
ttdDB <- read.table(here('data/raw/P1-01-TTD_target_download.txt'),
                    col.names = c('target_id', 'descriptor', 'geneName_drugID', 'drug_name', 'phase'),
                    skip = 32,
                    sep = '\t',
                    fill = TRUE,
                    quote = "",
                    na.strings = "")
write.table(x = ttdDB, file = here('data/processed/ttDB_trimmed.tsv'), sep = '\t',
            col.names = TRUE, row.names = FALSE, quote = FALSE)

matchup <- ttdDB %>%
  filter(descriptor == "GENENAME") %>%
  select(target_id, geneName_drugID)
drug_matchup <- ttdDB %>%
  select(target_id, drug_name) %>%
  filter(!is.na(drug_name))
final_matchup <- full_join(drug_matchup, matchup, by = 'target_id') %>%
  relocate(target_id, .after = geneName_drugID)
saveRDS(final_matchup, file = here('output/results/ttdDB_lookup.rds'))
write.table(x = final_matchup, file = here('output/results/ttdDB_lookup.tsv'),
            row.names = FALSE, quote = FALSE, sep = '\t')
ttdDB_lookup <- readRDS(here('data/processed/ttdDB_lookup.rds'))

ttdDB_synonyms <- read.table(here('data/raw/P1-04-Drug_synonyms.txt'),
                              col.names = c('drug_id', 'synonyms', 'drug_name'),
                              skip = 32,
                              sep = '\t',
                              fill = TRUE,
                              quote = "",
                              na.strings = "")
 saveRDS(ttdDB_synonyms, file = here('data/processed/ttdDB_synonyms.rds'))
 write.table(x = ttdDB_synonyms, file = here('data/processed/ttdDB_synonyms.tsv'),
             row.names = FALSE, quote = FALSE, sep = '\t')
 ttdDB_synonyms <- readRDS(here('data/processed/ttdDB_synonyms.rds'))

# match drug id to drug names
linker <- ttdDB %>%
 filter(descriptor == "DRUGINFO") %>%
 select(3,4) %>%
 distinct()

# final table
# synonyms > linker > lookup
linker_final <- left_join(linker, ttdDB_synonyms, by = c('geneName_drugID' = 'drug_id')) %>%
 left_join(., ttdDB_lookup, by = c('drug_name.x' = 'drug_name')) %>%
 select(2,4,5) %>%
 setNames(c('drug_alternate', 'drug_name', 'target')) %>%
 filter(!is.na(target))
 write.table(linker_final, file = here('data/processed/ttdDB_lookup_final.tsv'),
             row.names = FALSE, quote = FALSE, sep = '\t')
 saveRDS(linker_final, file = here('data/processed/ttdDB_lookup_final.rds'))
 

#~~~~~~~~~~~~~~~~~~#
# cell surface ----
#~~~~~~~~~~~~~~~~~~#
## experimental data ----
cspa_entrez <- read.delim(here('data/raw/cst_lists/Bausch-Fluck2015_entrez.txt'),header=F)[,1]
cspa_sym <- read.delim(here('data/raw/cst_lists/Bausch-Fluck2015_symbol.txt'),header=F)[,1]
cspa <- mapIds(org.Hs.eg.db, as.character(cspa_entrez), "SYMBOL","ENTREZID")
cspa <- union(cspa,cspa_sym)
cspa <- union(cspa,alias2SymbolTable(cspa))
cspa <- cspa[!is.na(cspa)]

hpa_tab <- read.delim(here('data/raw/cst_lists/hpa_subcellular_location.tsv'),header=T,sep='\t')
rows_cst <- grepl('Plasma membrane',hpa_tab$Approved,fixed=T) | 
  grepl('Plasma membrane',hpa_tab$Enhanced,fixed=T) | 
  grepl('Plasma membrane',hpa_tab$Supported,fixed=T) | 
  grepl('Plasma membrane',hpa_tab$Main.location,fixed=T) | 
  grepl('Plasma membrane',hpa_tab$Additional.location,fixed=T)
hpa <- mapIds(org.Hs.eg.db, as.character(hpa_tab[rows_cst,'Gene']), "SYMBOL","ENSEMBL")
hpa <- union(hpa,hpa_tab[rows_cst,'Gene.name'])
hpa <- union(hpa,alias2SymbolTable(hpa))
hpa <- hpa[!is.na(hpa)]


## computational data ----
go_surface <- read.delim(here('data/raw/cst_lists/GO0009986.txt'),header=F,sep='\t')[,1]
go_surface_ext <- read.delim(here('data/raw/cst_lists/GO0009897.txt'),header=F,sep='\t')[,1]
go <- c(go_surface,go_surface_ext)
go <- union(go,alias2SymbolTable(go))
go <- go[!is.na(go)]

BauschFluck2018 <- read.delim(here('data/raw/cst_lists/Bausch-Fluck2018_symbol.txt'),header=F)[,1]
BauschFluck2018 <- union(BauschFluck2018,alias2SymbolTable(BauschFluck2018))
BauschFluck2018 <- BauschFluck2018[!is.na(BauschFluck2018)]

daCunha2009 <- read.delim(here('data/raw/cst_lists/daCunha2009_symbol.txt'),header=F)[,1]
daCunha2009 <- union(daCunha2009,alias2SymbolTable(daCunha2009))
daCunha2009 <- daCunha2009[!is.na(daCunha2009)]

Donnard2014 <- read.delim(here('data/raw/cst_lists/Donnard2014_symbol.txt'),header=F)[,1]
Donnard2014 <- union(Donnard2014,alias2SymbolTable(Donnard2014))
Donnard2014 <- Donnard2014[!is.na(Donnard2014)]

Lee2018 <- read.delim(here('data/raw/cst_lists/Lee2018_symbol.txt'),header=F)[,1]
Lee2018 <- union(Lee2018,alias2SymbolTable(Lee2018))
Lee2018 <- Lee2018[!is.na(Lee2018)]

Governa2022_ensembl <- read.delim(here('data/raw/cst_lists/Governa2022_SURFME_ensembl.txt'),header=F)[,1]
Governa2022_sym <- read.delim(here('data/raw/cst_lists/Governa2022_SURFME_symbol.txt'),header=F)[,1]
Governa2022 <- mapIds(org.Hs.eg.db, as.character(Governa2022_ensembl), "SYMBOL","ENSEMBL")
Governa2022 <- union(Governa2022,Governa2022_sym)
Governa2022 <- union(Governa2022,alias2SymbolTable(Governa2022))
Governa2022 <- Governa2022[!is.na(Governa2022)]

Hu2021_ensembl <- read.delim(here('data/raw/cst_lists/Hu2021_ensembl.txt'),header=F)[,1]
Hu2021_sym <- read.delim(here('data/raw/cst_lists/Hu2021_symbol.txt'),header=F)[,1]
Hu2021 <- mapIds(org.Hs.eg.db, as.character(Hu2021_ensembl), "SYMBOL","ENSEMBL")
Hu2021 <- union(Hu2021,Hu2021_sym)
Hu2021 <- union(Hu2021,alias2SymbolTable(Hu2021))
Hu2021 <- Hu2021[!is.na(Hu2021)]

# experimental
exp_calls <- union(cspa, hpa)
exp_calls <- unique(exp_calls)

write.csv(exp_calls, file = here('data/processed/cell_surface_experimental_all_genes.csv'),
          quote = FALSE, row.names = FALSE)

# computational
compsets <- list(go, BauschFluck2018, daCunha2009, Donnard2014, Lee2018, Governa2022, Hu2021)
names(compsets) <- c('go', 'BauschFluck2018', 'daCunha2009', 'Donnard2014', 'Lee2018', 'Governa2022', 'Hu2021')

names_df <- tibble::tibble(
  set = rep(names(compsets), lengths(compsets)),
  name = unlist(compsets)
)

write.csv(names_df, file = here('data/processed/cell_surface_computational_all_genes.csv'),
          quote = FALSE, row.names = FALSE)