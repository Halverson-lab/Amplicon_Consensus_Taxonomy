#!/usr/bin/env Rscript

library(multidplyr) #parallel dplyr
library(readr)
library(tidyverse)
options("scipen"=999) # very important, keeps R from converting taxids to scientific notation

#### Environment setup ####
readRenviron("config.txt")

path_to_work_dir <- Sys.getenv("WORK_DIR")
setwd(path_to_work_dir)
if (!str_ends(path_to_work_dir, "/")){ path_to_work_dir <- paste0(path_to_work_dir, "/")}

path_to_data_dir <- Sys.getenv("DATABASE_DIR")
if (!str_ends(path_to_data_dir, "/")){ path_to_data_dir <- paste0(path_to_data_dir, "/")}

NA_thresh <- Sys.getenv("NA_THRESHOLD")

n_cores <- as.numeric(Sys.getenv("N_CORE"))

#Number of libraries
lib_num <- Sys.getenv("LIBRARY")
#list of names of all barcodes
barcode_list <- scan(paste0(path_to_work_dir,"barcode_names.txt"), what="", sep="\n")

contains_known_orgs <- Sys.getenv("KNOWN_ORG")

if(contains_known_orgs == TRUE){
  path_known_orgs <- Sys.getenv("KNOWN_ORG_SPECIES")
  # If sample includes added organisms the include a tsv of column labeled species, containing species name with taxid of each 
  # species names must match format if genus_species_taxid, for example "Pseudomonas_putida_303"
  known_orgs_df <- read_tsv(path_known_orgs)
  known_orgs_list <- known_orgs_df$species
  
  # additionally include a 2 colum csv
  path_known_orgs_samples <- Sys.getenv("KNOWN_ORG_SAMPLES")
  known_orgs_samples <- read_csv(path_known_orgs_samples)
  known_orgs_samples <- column_to_rownames(known_orgs_samples, var = "barcode")
} else {
  known_orgs_samples <- as.data.frame(barcode_list)
  known_orgs_samples$known_orgs_added <- F
  known_orgs_samples <- column_to_rownames(known_orgs_samples, var = "barcode_list")
}


# file path to the emu and sintax outputs
emu_path <- Sys.getenv("EMU_OUT")
sintax_path <- Sys.getenv("SINTAX_OUT")
laca_path <- Sys.getenv("LACA_OUT")

if (!str_ends(emu_path, "/")){ emu_path <- paste0(emu_path, "/")}
if (!str_ends(sintax_path, "/")){ sintax_path <- paste0(sintax_path, "/")}
if (!str_ends(laca_path, "/")){ laca_path <- paste0(laca_path, "/")}


# Read in the seqID to OTU table you generated
otu_df <- read_tsv(paste0(laca_path, "quant/seqID_to_otu.tsv"))

# Read in the taxonomy assignments for the otus
otu_taxonomy_df <- read_tsv(paste0(laca_path, "OTU_sintax.tsv"), col_names = FALSE)

# Read in EMU's taxonomy table for reference
raw_taxon_table <- read_tsv(paste0(path_to_data_dir,"taxonomy.tsv"), col_types = cols(.default = "c"))



#### Format dataframes, including the taxonomy dataframes. Do not edit, unless absolutely necessary ####

#Format emu taxonomy reference
# Paste the taxid id to the corresponding phylum level (Bactera_2, Pseudomonadota_1224, etc)
emu_taxid_taxonomy <- raw_taxon_table %>%
  mutate(domain = if_else(is.na(domain), NA, paste0(domain, "_", t_domain)),
         phylum = if_else(is.na(phylum), NA, paste0(phylum, "_", t_phylum)),
         class = if_else(is.na(class), NA, paste0(class, "_", t_class)),
         order = if_else(is.na(order), NA, paste0(order, "_", t_order)),
         family = if_else(is.na(family), NA, paste0(family, "_", t_family)),
         genus = if_else(is.na(genus), NA, paste0(genus, "_", t_genus)),
         species = if_else(is.na(species), NA, paste0(species, "_", t_species))
  ) %>%
  select(c(tax_id, domain, phylum, class, order, family, genus, species)) %>%   # keep only the selected columns
  mutate(across(!tax_id, ~ str_replace_all(., " ", "_"))) %>% # replace all spaces with underscores
  distinct() %>% # keep only distinct
  mutate(tax_id = paste0("taxid_", tax_id)) %>% # add taxid_ to the taxid numbers
  column_to_rownames(var = "tax_id") # make the taxid the rownames


# Get seqID to OTU from LACA
colnames(otu_df) <- c("count", "rep_cls", "cls", "otu", "read_id") #names the columns
otu_df <- otu_df[,c(4,5)] #just keep col otu and read_id
otu_df <- otu_df %>% mutate(otu = paste0("OTU_", otu)) #paste OTU_ in front of the OTU numbers


# Prep OTU to taxonomy table
colnames(otu_taxonomy_df) <- c("otu", "taxon_CI", "strand", "final_taxa") # name the columns
otu_taxonomy <- otu_taxonomy_df[,c(1,2)] # just keep the otu and the taxonomy assignment with confidence intervals
otu_taxonomy <- otu_taxonomy %>%
  mutate(domain = str_remove_all(str_extract(taxon_CI, "(?<=d:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","), # pull each level of taxonomy and remove the , 
         phylum = str_remove_all(str_extract(taxon_CI, "(?<=p:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","),
         class = str_remove_all(str_extract(taxon_CI, "(?<=c:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","),
         order = str_remove_all(str_extract(taxon_CI, "(?<=o:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","),
         family = str_remove_all(str_extract(taxon_CI, "(?<=f:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","),
         genus = str_remove_all(str_extract(taxon_CI, "(?<=g:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), ","),
         species = str_remove_all(str_extract(taxon_CI, "(?<=s:).*"), ",")) %>%
  select(-taxon_CI) %>%
  pivot_longer(cols = !otu, names_to = "taxon_level", values_to = "otu_taxa") %>% #pivot data to long form
  separate_wider_regex(otu_taxa, c(otu_taxa = ".*", "\\(", otu_taxa_CI = ".*")) %>% # split into the taxa and the CI for that taxa
  mutate(otu_taxa_CI = str_remove_all(otu_taxa_CI, pattern = "[^0-9.]"),  # remove anything not a number from the CI, like comma or parentheses
         otu_taxa_CI = na_if(otu_taxa_CI, ""), #na if empty or contains only _NA
         otu_taxa = na_if(otu_taxa, "_NA"))


#### Function to combine emu and sintax output, you can choose to edit thresholds if necessary ##
combine_taxonomy <- function(sample_id, sintax_df, emu_df, minimap_df){
  # Are there added orgs? 
  known_orgs <- known_orgs_samples[sample_id, 1]
  if(is.na(known_orgs)){known_orgs = FALSE} #if gnoto is empty then assign false
  
  #### Format sintax ####
  colnames(sintax_df) <- c("read_id", "taxon_CI", "strand", "final_taxa")
  sintax_taxa <- sintax_df[,c(1,2)] #keep only the read_id and the taxonomy levels with CI
  sintax_taxa <- sintax_taxa %>%
    mutate(domain = str_extract(taxon_CI, "(?<=d:)(.*?)\\([0-9]\\.[0-9]{2}\\)"), #extract each level of taxonomy
           phylum = str_extract(taxon_CI, "(?<=p:)(.*?)\\([0-9]\\.[0-9]{2}\\)"),
           class = str_extract(taxon_CI, "(?<=c:)(.*?)\\([0-9]\\.[0-9]{2}\\)"),
           order = str_extract(taxon_CI, "(?<=o:)(.*?)\\([0-9]\\.[0-9]{2}\\)"),
           family = str_extract(taxon_CI, "(?<=f:)(.*?)\\([0-9]\\.[0-9]{2}\\)"),
           genus = str_extract(taxon_CI, "(?<=g:)(.*?)\\([0-9]\\.[0-9]{2}\\)"),
           species = str_extract(taxon_CI, "(?<=s:).*"))
  
  ## if the read id contains info from sequencing then remove that info
  if(sum(str_detect(sintax_taxa$read_id, " ")) > 0){
    sintax_taxa <- sintax_taxa %>%
      mutate(read_id = gsub("[[:blank:]].*?$", "", read_id)) %>%
      select(-taxon_CI) # get rid of taxon_CI and run_info, they're not used
  } else {
    sintax_taxa <- sintax_taxa %>%
      select(-taxon_CI) # get rid of taxon_CI and run_info, they're not used
  }
  
  ## Add otu to sintax and format
  sintax_otu <- left_join(sintax_taxa, otu_df, by = "read_id") %>% #join keeping all values of the sintax df
    pivot_longer(cols = !c(read_id, otu), #pivot to long format
                 names_to = "taxon_level",
                 values_to = "sintax") %>%
    separate_wider_regex(sintax, c(sintax = ".*", "\\(", sintax_CI = ".*")) %>% #split the sintax taxonomy and CI into 2 columns
    mutate(sintax_CI = str_remove_all(sintax_CI, pattern = "[^0-9.]"), #clean up the CI again
           sintax_CI = na_if(sintax_CI, ""), #assign NA if empty or if _NA
           sintax = na_if(sintax, "_NA"))
  
  #### format minimap ####
  minimap_df_trim <- minimap_df %>% 
    group_by(QNAME) %>%
    mutate(max_AS = max(AS)) %>%
    filter(AS >= max_AS) %>%
    mutate(min_NM = min(NM)) %>%
    filter(NM <= min_NM) %>%
    mutate(percent_match = CMATCH/QLENGTH,
           percent_mismatch = N_MISMATCH/RLENGTH) %>%
    filter(percent_match >= 0.8 & percent_mismatch <= 0.2) %>%
    select(-c(max_AS, min_NM)) %>%
    separate_wider_delim(RNAME, delim = ":", names = c("taxid", "database", "seq_num")) %>%
    mutate(taxid = if_else(is.na(taxid), NA, paste0("taxid_", taxid)))


  #multipdlyr if nrows is greater than n_cores*1000
  if(nrow(minimap_df_trim) > (n_cores*1000)){
    #subset the taxonomy to only include the taxids needed, significantly speeds up working with large taxonomies
    minimap_tax_list <- unique(minimap_df_trim$taxid)
    minimap_emu_taxid_taxonomy <- emu_taxid_taxonomy[minimap_tax_list,]
    cluster_assign(cluster, minimap_emu_taxid_taxonomy = minimap_emu_taxid_taxonomy)
  
    #partition the dataframe between clusters
    minimap_part <- minimap_df_trim %>% group_by(taxid) %>% partition(cluster)
  
    minimap_with_tax <- minimap_part %>%
      mutate(
        domain = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "domain"]), #check if the taxid is NA, otherwise get the corresponding taxonomy from the taxonomy table
        phylum = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "phylum"]),
        class = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "class"]),
        order = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "order"]),
        family = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "family"]),
        genus = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "genus"]),
        species = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "species"])) %>%
      collect()
  } else {
    minimap_with_tax <- minimap_df_trim %>%
      mutate(
        domain = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "domain"]), 
        phylum = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "phylum"]),
        class = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "class"]),
        order = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "order"]),
        family = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "family"]),
        genus = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "genus"]),
        species = if_else(is.na(taxid), NA, minimap_emu_taxid_taxonomy[as.character(taxid), "species"])) 
  }

  
  ## if it scores equally well for 2+ diff species then species can't be assigned
  #keep top scoring assignment for each  
  minimap_with_tax_trim <- minimap_with_tax %>%
    group_by(QNAME, species) %>%
    mutate(max_AS = max(AS, na.rm = T)) %>%
    filter(AS >= max_AS) %>%
    mutate(min_NM = min(NM, na.rm = T)) %>%
    filter(NM <= min_NM) %>%
    mutate(min_N_MISMATCH = min(N_MISMATCH, na.rm = T)) %>%
    filter(N_MISMATCH <= min_N_MISMATCH) %>%
    select(-c(seq_num, database, taxid)) %>%
    distinct() %>%
    group_by(QNAME) %>%
    mutate(n=n())
  
  minimap_final <- filter(minimap_with_tax_trim, n == 1)
  
  #if two species still map equally well to a read then back it off to genus or family level
  if (max(minimap_with_tax_trim$n, na.rm = T) > 1) {
    #remove species to identify best matching by genus
    minimap_no_species <- minimap_with_tax_trim %>%
      filter(n>1) %>%
      mutate(species = as.character(NA)) %>%
      group_by(QNAME, genus) %>%
      mutate(max_AS = max(AS, na.rm = T)) %>%
      filter(AS >= max_AS) %>%
      mutate(min_NM = min(NM, na.rm = T)) %>%
      filter(NM <= min_NM) %>%
      mutate(min_N_MISMATCH = min(N_MISMATCH, na.rm = T)) %>%
      filter(N_MISMATCH <= min_N_MISMATCH) %>%
      distinct() %>%
      group_by(QNAME) %>%
      mutate(n=n())
    
    #bind to df
    minimap_final <- bind_rows(minimap_final, filter(minimap_no_species, n == 1))
    
    if (max(minimap_no_species$n, na.rm = T) > 1) {
      #remove genus to identify best matching by family
      minimap_no_genus <- minimap_no_species %>%
        filter(n>1) %>%
        mutate(genus = as.character(NA)) %>%
        group_by(QNAME, family) %>%
        mutate(max_AS = max(AS, na.rm = T)) %>%
        filter(AS >= max_AS) %>%
        mutate(min_NM = min(NM, na.rm = T)) %>%
        filter(NM <= min_NM) %>%
        mutate(min_N_MISMATCH = min(N_MISMATCH, na.rm = T)) %>%
        filter(N_MISMATCH <= min_N_MISMATCH) %>%
        distinct() %>%
        group_by(QNAME) %>%
        mutate(n=n())
      
      #bind to df
      minimap_final <- bind_rows(minimap_final, filter(minimap_no_genus, n == 1))
    }
    
  }

  
  minimap_final <- minimap_final %>%
    select(-c(min_N_MISMATCH, max_AS, min_NM, AS, NM, CINS, CDEL, n)) %>%
    pivot_longer(domain:species, names_to = "taxon_level", values_to = "minimap") %>%
    rename(read_id=QNAME)
  
  
  ## Format emu
  colnames(emu_df)[1] <- "read_id"
  # emu file passes back all possible assignments for each read. For each read keep the assignment with the highest confidence
  emu_long <- emu_df %>%
    pivot_longer(!read_id, names_to = "taxid", values_to = "confidence") %>%
    filter(!is.na(confidence)) %>%
    group_by(read_id) %>%
    summarize(n = n(),
              emu_max_CI = max(confidence),
              emu_tax_id = taxid[which.max(confidence)]) %>%
    mutate(emu_tax_id = if_else(is.na(emu_tax_id), NA, paste0("taxid_", emu_tax_id)))
  
  #subset the taxonomy to only include the taxids needed, significantly speeds up working with large taxonomies
  emu_tax_list <- unique(emu_long$emu_tax_id)
  subset_emu_taxid_taxonomy <- emu_taxid_taxonomy[emu_tax_list,]
  cluster_assign(cluster, subset_emu_taxid_taxonomy = subset_emu_taxid_taxonomy)
  
  #partition the dataframe between clusters
  emu_part <- emu_long %>% group_by(read_id) %>% partition(cluster)
  
  ## Add in emu taxonomy
  emu_with_taxon <- emu_part %>%
    mutate(
      domain = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "domain"]), #check if the taxid is NA, otherwise get the corresponding taxonomy from the taxonomy table
      phylum = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "phylum"]),
      class = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "class"]),
      order = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "order"]),
      family = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "family"]),
      genus = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "genus"]),
      species = if_else(is.na(emu_tax_id), NA, subset_emu_taxid_taxonomy[as.character(emu_tax_id), "species"])
    ) %>%
    select(-c(n, emu_max_CI, emu_tax_id)) %>% #remove these columns
    collect() 
  
  emu_with_taxon <- emu_with_taxon %>%
    left_join(otu_df, by = "read_id") %>% #join to the otu dataframe
    pivot_longer(!c(read_id, otu), names_to = "taxon_level", values_to = "emu") #pivot long
  
  
  ## Combine emu, sintax, and the otu dataframe
  all_assignments <- sintax_otu %>%
    full_join(emu_with_taxon, by = c("read_id", "taxon_level", "otu")) %>%
    left_join(otu_taxonomy, b = c("otu", "taxon_level")) %>%
    left_join(minimap_final) %>%
    ungroup() %>%
    mutate(comp_flags = 0) %>%
    mutate(comp_flags = if_else(!is.na(sintax) & !is.na(otu_taxa) & sintax == otu_taxa, (comp_flags + 1), comp_flags)) %>%
    mutate(comp_flags = if_else(!is.na(sintax) & !is.na(emu) & sintax == emu, comp_flags + 10, comp_flags)) %>%
    mutate(comp_flags = if_else(!is.na(sintax) & !is.na(minimap) & sintax == minimap, comp_flags + 100, comp_flags)) %>%
    mutate(comp_flags = if_else(!is.na(emu) & !is.na(otu_taxa) & emu == otu_taxa, comp_flags + 1000, comp_flags)) %>%
    mutate(comp_flags = if_else(!is.na(emu) & !is.na(minimap) & emu == minimap, comp_flags + 10000, comp_flags)) %>%
    mutate(comp_flags = if_else(!is.na(otu_taxa) & !is.na(minimap) & otu_taxa == minimap, comp_flags + 100000, comp_flags))
  
  # sintax == OTU <- 1
  # sintax == emu <- 10
  # sintax == minimap <- 100
  # emu == OTU <- 1000
  # emu == minimap <- 10000
  # OTU == minimap <- 100000
  # so 10110 means sintax = emu = minimap
  
  # decision tree
  # the different paths account for the fact that EMU is more accurate in systems containing known members, due to their expectation maximization algorithm
  
  if(contains_known_orgs){
    all_assignments <- all_assignments %>%
      mutate(final_taxon = case_when((comp_flags == 111111) ~ sintax, #all equal so pick one,
                                     (emu %in% known_orgs_list) ~ emu,
                                     (sintax_CI <= NA_thresh) ~ NA, # sintax has low CI then NA
                                     (comp_flags %in% c(10, 100, 1011, 10110, 100101)) ~ sintax, #sintax=emu, sintax=minimap, sintax=otu=emu, sintax=emu=minimap, sintax=otu=minimap
                                     (comp_flags == 111000) ~ emu, #emu=otu=minimap
                                     (comp_flags %in% c(0, 1100, 100000)) & (percent_match > 0.99) & (percent_mismatch < 0.01) ~ minimap, #minimap has a really good match,
                                     (comp_flags %in% c(1000, 1100)) & (otu_taxa_CI > NA_thresh) ~ emu, #emu = otu and minimap doesn't pass previous threshold
                                     is.na(emu) & is.na(minimap) & (sintax == otu_taxa) ~ sintax,
                                     is.na(emu) & is.na(minimap) & (sintax != otu_taxa) & (sintax_CI > otu_taxa_CI ) ~ sintax,
                                     is.na(emu) & is.na(minimap) & (sintax != otu_taxa) & (sintax_CI < otu_taxa_CI ) ~ otu_taxa,
                                     is.na(emu) & is.na(minimap) & is.na(otu_taxa) & (sintax_CI > 0.9) ~ sintax,
                                     .default = NA))
  }else{
    all_assignments <- all_assignments %>%
      mutate(final_taxon = case_when((comp_flags == 111111) ~ sintax, #all equal so pick one,
                                     #(emu %in% known_orgs_list) ~ emu,
                                     (sintax_CI <= NA_thresh) ~ NA, # sintax has low CI then NA
                                     (comp_flags %in% c(10, 100, 1011, 10110, 100101)) ~ sintax, #sintax=emu, sintax=minimap, sintax=otu=emu, sintax=emu=minimap, sintax=otu=minimap
                                     (comp_flags == 111000) ~ emu, #emu=otu=minimap
                                     (comp_flags %in% c(0, 1100, 100000)) & (percent_match > 0.99) & (percent_mismatch < 0.01) ~ minimap, #minimap has a really good match,
                                     (comp_flags %in% c(1000, 1100)) & (otu_taxa_CI > NA_thresh) ~ emu, #emu = otu and minimap doesn't pass previous threshold
                                     is.na(emu) & is.na(minimap) & (sintax == otu_taxa) ~ sintax,
                                     is.na(emu) & is.na(minimap) & (sintax != otu_taxa) & (sintax_CI > otu_taxa_CI ) ~ sintax,
                                     is.na(emu) & is.na(minimap) & (sintax != otu_taxa) & (sintax_CI < otu_taxa_CI ) ~ otu_taxa,
                                     is.na(emu) & is.na(minimap) & is.na(otu_taxa) & (sintax_CI > 0.9) ~ sintax,
                                     .default = NA))
  }
  
  #pivot out and keep the lowest level that's not NA then assign taxid to that level
  final_assignments <- all_assignments %>%
    select(c(read_id, taxon_level, final_taxon, otu)) %>% #keep these columns
    pivot_wider(names_from = taxon_level, values_from = final_taxon) %>% #pivot wide
    mutate(lowest_taxon = case_when(!is.na(species) ~ species, #identify the lowest taxid that isn't NA
                                    !is.na(genus) ~ genus,
                                    !is.na(family) ~ family,
                                    !is.na(order) ~ order,
                                    !is.na(class) ~ class,
                                    !is.na(phylum) ~ phylum,
                                    !is.na(domain) ~ domain,
                                    .default = NA),
           lowest_taxon_level = case_when(!is.na(species) ~ "species", #identify the level of the lowest taxid that isn't NA
                                          !is.na(genus) ~ "genus",
                                          !is.na(family) ~ "family",
                                          !is.na(order) ~ "order",
                                          !is.na(class) ~ "class",
                                          !is.na(phylum) ~ "phylum",
                                          !is.na(domain) ~ "domain",
                                          .default = NA)
    ) %>%
    group_by(domain, phylum, class, order, family, genus, species, lowest_taxon, lowest_taxon_level, otu) %>%
    summarize(read_counts = n()) %>% #summarize to reduce table size
    mutate(sample_id = sample_id) %>% #add column identifying the sample id
    separate_wider_regex(lowest_taxon, c(lowest_taxon = ".*", "_", lowest_taxid = ".*")) #split lowest taxon into name and taxid

  return(final_assignments)
}




#### Loop through all files ####

## nested for loops to loop through library (1-4) and barcode (1-72) to get sample_id for each
## Change the numbers to match yours, might need to modify the file names
#counter for number of samples

#create a cluster
cluster <- new_cluster(n_cores-1)
#pass tidyverse and the emu_taxid_taxonomy table to every cluster
cluster_library(cluster, "tidyverse")

#list all possible library and barcode combos
sample_list <- as.vector(outer(c(1:lib_num), barcode_list, paste, sep="_"))
#initialize empty dataframe
sample_id_df <- data.frame()

for(sample_id in sample_list){
  # check if emu file exists for sample
  if(file.exists(paste0(emu_path,"read_assignments/", sample_id, "_read-assignment-distributions.tsv")) & file.exists(paste0(sintax_path, sample_id, "_sintax.tsv"))) {
    # read in emu and sintax files
    emu_file <- read_tsv(paste0(emu_path,"read_assignments/", sample_id, "_read-assignment-distributions.tsv"), col_types = cols(.default = "d", `...1` = "c"))
    minimap_file <- read_csv(paste0(emu_path, "minimap2_aln_stats/", sample_id, "_aln_stats.csv"))
    sintax_file <- read_tsv(paste0(sintax_path, sample_id, "_sintax.tsv"), col_names = FALSE)
    
    # run the combine taxonomy function with the objects just generated and save the output
    output_df <- combine_taxonomy(sample_id, sintax_file, emu_file, minimap_file)
    output_df$sample_id <- sample_id
    
    sample_id_df <- rbind(sample_id_df, output_df)
  }
}


#### compare otu assignments and clean them up a little ####

# if >= 98% of an otu is assigned to one species then change the other <2% to match
otu_comp_df <- sample_id_df %>%
  group_by(otu, domain, phylum, class, order, family, genus, species, lowest_taxon, lowest_taxid, lowest_taxon_level) %>%
  summarize(n_samples=n(),
            total_reads = sum(read_counts)) %>%
  ungroup() %>%
  group_by(otu) %>%
  mutate(total_otu_reads = sum(total_reads),
         percent_otu = total_reads/total_otu_reads*100) %>%
  filter(percent_otu < 100 & percent_otu >= 98) %>%
  select(-c(n_samples:percent_otu)) %>%
  column_to_rownames(var = "otu")

otu_list <- rownames(otu_comp_df)

sample_id_df_otu_comp <- sample_id_df %>%
  mutate(
    domain = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "domain"], domain), 
    phylum = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "phylum"], phylum),
    class = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "class"], class),
    order = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "order"], order),
    family = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "family"], family),
    genus = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "genus"], genus),
    species = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "species"], species),
    lowest_taxon = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "lowest_taxon"], lowest_taxon),
    lowest_taxid = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "lowest_taxid"], lowest_taxid),
    lowest_taxon_level = if_else(otu %in% otu_list, otu_comp_df[as.character(otu), "lowest_taxon_level"], lowest_taxon_level)
  ) %>%
  group_by(lowest_taxon, lowest_taxid, lowest_taxon_level, otu, sample_id) %>%
  summarize(read_counts = sum(read_counts))


#### Convert to an otu table in microeco format, edit file names if you want ####

# table with taxid and otu, for use with alpha diversity measurements
final_taxid_OTU_table <- sample_id_df_otu_comp %>%
  mutate(tax_id = if_else(is.na(otu), paste0("taxid_", lowest_taxid), paste("taxid", lowest_taxid, otu, sep = "_"))) %>%
  ungroup() %>%
  select(c(tax_id, read_counts, sample_id)) %>% # keep only these columns
  pivot_wider(names_from = sample_id, values_from = read_counts, values_fill = 0) # pivot wide for otu matrix

# save the table
write_tsv(final_taxid_OTU_table, "abundance_table_with_OTU.tsv", col_names = T, eol = "\n") 


# table with taxid only, for non-alpha diversity uses
# just keep the lowest taxids and the read counts for each sample, classic otu table format
final_taxid_table <- sample_id_df_otu_comp %>%
  mutate(tax_id = paste("taxid", lowest_taxid, sep = "_")) %>%
  ungroup() %>%
  group_by(sample_id, tax_id) %>%
  summarize(counts = sum(read_counts)) %>%
  pivot_wider(names_from = sample_id, values_from = counts, values_fill = 0)

#save the table
write_tsv(final_taxid_table, "abundance_table.tsv", col_names = T, eol = "\n")
