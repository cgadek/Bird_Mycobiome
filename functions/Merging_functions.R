#Merging Functions####

require(phyloseq)
require(tidyverse)

# Concatenate unique values in a vector
concat_unique <- function(vec){
  uniq <- unique(as.character(vec))
  return(paste(uniq, collapse = "/"))
}

# Like psmelt, but only uses the otu_table and sample_data
ps_semi_melt <- function(ps){
  otu_table(ps) %>%
    data.frame(taxid = row.names(.)) %>%
    rename_with(function(x){gsub("X", "", x)}) %>%
    pivot_longer(!taxid, names_to = "sample_id", values_to = "abundance") %>%
    left_join(sample_data(ps) %>%
                data.frame(sample_id = row.names(.)),
              by = "sample_id")
}

# Function that summarizes a vector based on its class
summarise_vec <- function(vec){
  if(class(vec) %in% c("numeric", "integer", "logical")){
    return(mean(vec, na.rm = T))
  } else if (class(vec) %in% c("factor", "character")){
    return(concat_unique(vec))
  } else {
    stop("Error: unknown column type")
  }
}

# Converts a summary df to an otu_table
summ_to_otu_tbl <- function(summ){
  summ %>% 
    dplyr::select(taxid, sample_id, abundance) %>% 
    pivot_wider(names_from = "sample_id", values_from = "abundance") %>%
    column_to_rownames('taxid') %>%
    as.matrix() %>%
    otu_table(., taxa_are_rows = TRUE)
}

# Converts a summary df to sample_data
summ_to_sample_dat <- function(summ){
  summ %>% 
    dplyr::select(!c(taxid, abundance)) %>% 
    distinct() %>%
    column_to_rownames('sample_id') %>%
    sample_data()
}

# Function that merges phyloseq samples based on the names of one or more grouping factors
# present in sample_data(ps)
merge_ps_samples <- function(ps, grouping){
  
  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_table(ps) <- phyloseq::otu_table(t(otu_table(ps)), taxa_are_rows = T)
  }
  
  # Convert to long format
  ps_long <- ps_semi_melt(ps)
  
  # Summarise all columns
  summ <- ps_long %>%
    group_by(across(all_of(!!grouping))) %>%
    group_by(taxid, .add = T) %>%
    summarise(across(everything(), summarise_vec)) %>%
    ungroup()
  
  # Convert to otu_table and sample_data
  otu_tbl <- summ_to_otu_tbl(summ)
  sample_dat <- summ_to_sample_dat(summ)
  
  # Create new physeq object
  new_ps <- phyloseq(otu_tbl, sample_dat, tax_table(ps))
  return(new_ps)
}
