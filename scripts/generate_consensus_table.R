#!/usr/bin/env Rscript


library(readr)
library(tidyverse)
options("scipen"=999) # very important, keeps R from converting taxids to scientific notation

#### Environment setup ####
readRenviron("config.txt")

path_to_work_dir <- Sys.getenv("WORK_DIR")

path_to_data_dir <- Sys.getenv("DATABASE_DIR")

path_to_gnoto <- Sys.getenv("GNOTO")

syncom <- Sys.getenv("SYNCOM")

if(syncom == TRUE){
  path_to_syncom <- Sys.getenv("SYNCOM_SPECIES")
  # If sample includes syncoms the include a tsv of column labeled species, containing species name with taxid of each member of the syncom 
  # species names must match format if genus_species_taxid, for example "Pseudomonas_putida_303"
  syncom_df <- read_tsv(path_to_syncom)
  syncom_list <- syncom_df$species
}

#Number of libraries
lib_num <- Sys.getenv("LIBRARY")
#list of names of all barcodes
barcode_list <- scan("barcode_names.txt", what="", sep="\n")

# specify working directory
setwd(path_to_work_dir)

# file path to the emu and sintax outputs, location is relative to workdir
emu_path <- paste0(path_to_work_dir, "/6_emu/read_assignments/")
sintax_path <- paste0(path_to_work_dir, "/7_sintax/")

# Read in the seqID to OTU table you generated
otu_df <- read_tsv(paste0(path_to_work_dir, "/5_laca/quant/seqID_to_otu.tsv"))

# Read in the taxonomy assignments for the otus
otu_taxonomy_df <- read_tsv(paste0(path_to_work_dir, "/5_laca/sintax_OTUs.tsv"), col_names = FALSE)

# Read in EMU's taxonomy table for reference
raw_taxon_table <- read_tsv(paste0(path_to_data_dir,"/taxonomy.tsv"))


#### Edit here - Load Sample info ####

# is sample gnotobiotic? 2 column csv with every barcode and TRUE/FALSE for gnotobiotic
gnotobiotic <- read_csv(path_to_gnoto)
gnotobiotic <- column_to_rownames(gnotobiotic, var = "barcode")



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
otu_df <- otu_df %>% mutate(otu = paste0("OTU_", otu)) #paster OTU_ in front of the OTU numbers


# Prep OTU to taxonomy table
colnames(otu_taxonomy_df) <- c("otu", "taxon_CI", "strand", "final_taxa") # name the columns
otu_taxonomy <- otu_taxonomy_df[,c(1,2)] # just keep the otu and the taxonomy assignment with confidence intervals
otu_taxonomy <- otu_taxonomy %>%
  mutate(domain = str_remove_all(str_extract(taxon_CI, "(?<=d:)(.*?)\\,"), ","), # pull each level of taxonomy and remove the , 
         phylum = str_remove_all(str_extract(taxon_CI, "(?<=p:)(.*?)\\,"), ","),
         class = str_remove_all(str_extract(taxon_CI, "(?<=c:)(.*?)\\,"), ","),
         order = str_remove_all(str_extract(taxon_CI, "(?<=o:)(.*?)\\,"), ","),
         family = str_remove_all(str_extract(taxon_CI, "(?<=f:)(.*?)\\,"), ","),
         genus = str_remove_all(str_extract(taxon_CI, "(?<=g:)(.*?)\\,"), ","),
         species = str_remove_all(str_extract(taxon_CI, "(?<=s:).*"), ",")) %>%
  select(-taxon_CI) %>%
  #pivot data to long form
  pivot_longer(cols = !otu,
               names_to = "taxon_level",
               values_to = "otu_taxa") %>%
  separate_wider_regex(otu_taxa, c(otu_taxa = ".*", "\\(", otu_taxa_CI = ".*")) %>% # split into the taxa and the CI for that taxa
  mutate(otu_taxa_CI = str_remove_all(otu_taxa_CI, pattern = "[^0-9.]"),  # remove anything not a number from the CI, like comma or parentheses
         otu_taxa_CI = na_if(otu_taxa_CI, ""), #na if empty or contains only _NA
         otu_taxa = na_if(otu_taxa, "_NA"))


#### Function to combine emu and sintax output, you can choose to edit thresholds if necessary ####

combine_taxonomy <- function(sample_id, sintax_df, emu_df){
  # is sample gnotobiotic 
  gnoto <- gnotobiotic[sample_id, 1]
  if(is.na(gnoto)){gnoto = FALSE} #if gnoto is empty then assign false
  
  ## Format sintax
  colnames(sintax_df) <- c("read_id", "taxon_CI", "strand", "final_taxa")
  sintax_taxa <- sintax_df[,c(1,2)] #keep only the read_id and the taxonomy levels with CI
  sintax_taxa <- sintax_taxa %>%
    mutate(domain = str_extract(taxon_CI, "(?<=d:)(.*?)\\,"), #extract each level of taxonomy
           phylum = str_extract(taxon_CI, "(?<=p:)(.*?)\\,"),
           class = str_extract(taxon_CI, "(?<=c:)(.*?)\\,"),
           order = str_extract(taxon_CI, "(?<=o:)(.*?)\\,"),
           family = str_extract(taxon_CI, "(?<=f:)(.*?)\\,"),
           genus = str_extract(taxon_CI, "(?<=g:)(.*?)\\,"),
           species = str_extract(taxon_CI, "(?<=s:).*"))
  
  ## if the read id contains info from sequencing then remove that info
  if(sum(str_detect(sintax_taxa$read_id, " ")) > 0){
    sintax_taxa <- sintax_taxa %>%
      separate_wider_regex(read_id, c(read_id = ".*?", "\\s+", run_info = ".*"))  %>% #pull the read id from the run info
      select(-c(taxon_CI, run_info)) # get rid of taxon_CI and run_info, they're not used
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
  
  ## Format emu
  colnames(emu_df)[1] <- "read_id"
  # emu file passes back all possible assignments for each read. For each read keep the assignment with the highest confidence
  emu_long <- emu_df %>%
    pivot_longer(!read_id, names_to = "taxid", values_to = "confidence") %>%
    filter(!is.na(confidence)) %>%
    group_by(read_id) %>%
    summarize(n = n(),
              emu_max_CI = max(confidence),
              emu_tax_id = taxid[which.max(confidence)])
  
  ## Add in emu taxonomy
  emu_with_taxon <- emu_long %>%
    mutate(emu_tax_id = if_else(is.na(emu_tax_id), NA, paste0("taxid_", emu_tax_id))) %>%
    mutate(
      domain = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "domain"]), #check if the taxid is NA, otherwise get the corresponding taxonomy from the taxonomy table
      phylum = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "phylum"]),
      class = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "class"]),
      order = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "order"]),
      family = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "family"]),
      genus = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "genus"]),
      species = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "species"])
    ) %>%
    select(-c(n, emu_max_CI, emu_tax_id)) %>% #remove these columns
    left_join(otu_df, by = "read_id") %>% #join to the otu dataframe
    pivot_longer(!c(read_id, otu), names_to = "taxon_level", values_to = "emu") #pivot long
  
  
  
  ## Combine emu, sintax, and the otu dataframe
  all_assignments <- sintax_otu %>%
    full_join(emu_with_taxon, by = c("read_id", "taxon_level", "otu")) %>%
    left_join(otu_taxonomy, b = c("otu", "taxon_level"))
  
  
  # decision tree
  # selecting the final assignment is based on whether the sample is gnotobiotic and then on if it is using a syncom
  # Due to the likelihood of cross contamination, syncom is marked by experiment and not by sample
  #   Path 1 - sample is gnotobiotic
  #   Path 2 - sample is not gnotobiotic but does contain a syncom
  #   Path 3 - sample is not gnotobiotic and does not contain a syncom
  # the different paths account for the fact that EMU is more accurate in systems containg known members, due to their expectation maximization algorithm
  # Confidence thresholds (0.3 and 0.7) were chosen by looking at zymo community samples, but cutoffs can be changed to preference
  if(gnoto == T){
    all_assignments <- all_assignments %>% 
      mutate(final_taxon = case_when((sintax_CI <= 0.3) & !(emu %in% syncom_list) ~ NA, # if sintax has low confidence and emu's assignment is not int the syncom list then NA
                                     (sintax == emu) ~ emu, #if emu and sintax agree then pick one
                                     (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax, #if emu is empty and sintax decently confident then sintax
                                     (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa, #if emu is empty, sintax is not cofident, and OTU is then use OTU
                                     (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA, # if emu is empty, and sintax and the otu are not confident then NA
                                     (sintax != emu) ~ emu, #if there is an emu assignment and it doesn't match sintax, choose emu
                                     .default = NA)) 
  }else{
    # experiment has a SynCom
    if(syncom == T){
      all_assignments <- all_assignments %>%
        mutate(final_taxon = case_when((sintax_CI <= 0.3) & !(emu %in% syncom_list) ~ NA, # if sintax has low confidence and emu's assignment is not int the syncom list then NA
                                       (sintax == emu) ~ sintax,  #if emu and sintax agree then pick one
                                       (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax, #if emu is empty and sintax decently confident then sintax
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa, #if emu is empty, sintax is not confident, and OTU is then use OTU
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA, # if emu is empty, and sintax and the otu are not confident then NA
                                       (sintax != emu) & (emu %in% syncom_list) ~ emu, # emu and sintax don't match, but emu's assignment is in the syncom, then choose emu
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI >= 0.7) ~ sintax, # emu and sintax don't agree, emu isn't in the syncom, and sintax is confident, then choose sintax
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax == otu_taxa) ~ sintax, # sintax and otu agree, but emu doesn't and its not a syncom, choose sintax
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu == otu_taxa) ~ emu, # emu and otu agree, choose them
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu != otu_taxa) ~ NA, # no one agrees and sintax is not confident, NA
                                       .default = NA)) 
    }else{
      #experiment does not have a SynCom
      all_assignments <- all_assignments %>%
        mutate(final_taxon = case_when((sintax_CI <= 0.3) ~ NA, # sintax has low CI then NA
                                       (sintax == emu) ~ sintax, #if emu and sintax agree then pick one
                                       (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax, #if emu is empty and sintax decently confident then sintax
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa, #if emu is empty, sintax is not confident, and OTU is then use OTU
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA, # if emu is empty, and sintax and the otu are not confident then NA
                                       (sintax != emu) & (sintax_CI >= 0.7) ~ sintax, # emu and sintax don't agree and sintax is confident, then choose sintax
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax == otu_taxa) ~ sintax, # sintax and otu agree, choose sintax
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu == otu_taxa) ~ emu, # emu and otu agree, choose emu
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu != otu_taxa) ~ NA, # no one agrees then NA
                                       .default = NA))
    }
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
    ## keep otu if lowest taxon is in family or genus, but not any others
    mutate(otu = if_else(lowest_taxon_level %in% c("family", "genus"), otu, NA)) %>%
    group_by(lowest_taxon, lowest_taxon_level, otu) %>%
    summarise(read_counts = n()) %>% #summarise to reduce table size
    mutate(sample_id = sample_id) %>% #add column identifying the sample id
    separate_wider_regex(lowest_taxon, c(lowest_taxon = ".*", "_", lowest_taxid = ".*")) #split lowest taxon into name and taxid
  
  return(final_assignments)
}




#### Edit here - Loop through all files ####

## nested for loops to loop through library (1-4) and barcode (1-72) to get sample_id for each
## Change the numbers to match yours, might need to modify the file names
#counter for number of samples
sample_num <- 0

for(plate_num in 1:lib_num){
  for(barcode in barcode_list){
    # check if emu file exists for sample, looks for a file in the following format 1-1_read-assignment-distributions.tsv
    if(file.exists(paste0(emu_path, plate_num, "_", barcode, "_read-assignment-distributions.tsv"))){
      
      # read in emu and sintax files
      emu_file <- read_tsv(paste0(emu_path, plate_num, "_", barcode, "_read-assignment-distributions.tsv"))
      sintax_file <- read_tsv(paste0(sintax_path, plate_num, "_", barcode, "_sintax.tsv"), col_names = FALSE)
      
      #save the sample id (1-1, etc.)
      sample_id <- paste0(plate_num, "_", barcode)
      
      # run the combine taxonomy function with the objects just generated and save the output
      output_df <- combine_taxonomy(sample_id, sintax_file, emu_file)
      
      # track number of samples
      sample_num <- sample_num + 1
      
      # if its the first sample then the output df becomes the dataframe
      if(sample_num == 1){
        sample_id_df <- output_df
      }else{
        # after the first, append the output df onto the previous
        sample_id_df <- rbind(sample_id_df, output_df)
      }
    }
  }
}
# after the loop there is a dataframe, sample_id_df, with the taxonomic assignment for every read

#### Convert to an otu table in microeco format, edit file names if you want ####

# table with taxid and otu, for use with alpha diversity measurements
final_taxid_OTU_table <- sample_id_df %>%
  mutate(tax_id = case_when(is.na(lowest_taxon) ~ NA, # NA if NA
                            is.na(otu) ~ paste("taxid", lowest_taxid, sep = "_"), # otu is NA then just use keep taxid
                            (lowest_taxon_level == "species") ~ paste("taxid", lowest_taxid, sep = "_"), # if lowest taxon is species level then don't include otu
                            !is.na(otu) ~ paste("taxid", lowest_taxid, otu, sep = "_"))) %>% # everything else paste lowest taxid plus otu
  ungroup() %>%
  select(c(tax_id, read_counts, sample_id)) %>% # keep only these columns
  pivot_wider(names_from = sample_id, values_from = read_counts) # pivot wide for otu matrix

# save the table
write_tsv(final_taxid_OTU_table, "microeco_taxid_OTU_table.tsv", col_names = T, eol = "\n") 



# table with taxid only, for non-alpha diversity uses
# just keep the lowest taxids and the read counts for each sample, classic otu table format
final_taxid_table <- sample_id_df %>%
  mutate(tax_id = paste("taxid", lowest_taxid, sep = "_")) %>%
  ungroup() %>%
  group_by(sample_id, tax_id) %>%
  summarize(counts = sum(read_counts)) %>%
  pivot_wider(names_from = sample_id, values_from = counts)

#save the table
write_tsv(final_taxid_table, "microeco_taxid_table.tsv", col_names = T, eol = "\n")
