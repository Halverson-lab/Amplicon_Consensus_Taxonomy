#--------Final pipeline to take EMU and sintax output and return normalized tables ready for analysis -----------------------

#### Environment setup - Edit here ####

library(readr)
library(tidyverse)
options("scipen"=999)

setwd("//las-dfs-01.las.iastate.edu/lss/research/larryh-lab/Ashley-P/Microbiome_sequencing")

#file path to the emu and sintax outputs
  emu_path <- "Pipeline_output_files/6_emu/read_assignments/"
  sintax_path <- "Pipeline_output_files/7_sintax/"

# Read in the seqID to OTU table you generated
  otu_df <- read_tsv("Pipeline_output_files/LACA_all/quant/seqID_to_otu.tsv")

# Read in the taxonomy assignments for the otus
  otu_taxonomy_df <- read_tsv("Pipeline_output_files/LACA_all/sintax_OTUs.tsv", col_names = FALSE)


# Read in taxonomy for emu to reference
  raw_taxon_table <- read_tsv("Pipeline_output_files/taxonomy.tsv")
  
  
#### Load Sample info - Edit here ####

# is sample gnotobiotic? 2 column csv with every barcode and TRUE/FALSE for gnotobiotic
gnotobiotic <- read_csv("Sample_info/gnotobiotic_sample_info.csv")
gnotobiotic <- column_to_rownames(gnotobiotic, var = "barcode")

# Does the experiment use SynComs?
syncom <- TRUE

# tsv of column labeled species, containing species name with taxid of each member of the syncom 
# for example "Pseudomonas_putida_303"
syncom_df <- read_tsv("marsc_species_list.tsv")
syncom_list <- syncom_df$species


#### prep emu taxonomy reference ####

emu_taxid_taxonomy <- raw_taxon_table %>%
  mutate(domain = if_else(is.na(kingdom), NA, paste0(kingdom, "_", t_kingdom)),
         phylum = if_else(is.na(phylum), NA, paste0(phylum, "_", t_phylum)),
         class = if_else(is.na(class), NA, paste0(class, "_", t_class)),
         order = if_else(is.na(order), NA, paste0(order, "_", t_order)),
         family = if_else(is.na(family), NA, paste0(family, "_", t_family)),
         genus = if_else(is.na(genus), NA, paste0(genus, "_", t_genus)),
         species = if_else(is.na(species), NA, paste0(species, "_", t_species))
  ) %>%
  select(c(tax_id, domain, phylum, class, order, family, genus, species)) %>%
  mutate(across(!tax_id, ~ str_replace_all(., " ", "_"))) %>%
  distinct()%>%
  mutate(tax_id = paste0("taxid_", tax_id)) %>%
  column_to_rownames(var = "tax_id")


#### Get seqID to OTU from LACA ####
colnames(otu_df) <- c("count", "rep_cls", "cls", "otu", "read_id")
otu_df <- otu_df[,c(4,5)]
otu_df <- otu_df %>%
  mutate(otu = paste0("OTU_", otu))


#### Prep OTU to taxonomy table ####

colnames(otu_taxonomy_df) <- c("otu", "taxon_CI", "strand", "final_taxa")
otu_taxonomy <- otu_taxonomy_df[,c(1,2)]
otu_taxonomy <- otu_taxonomy %>%
  mutate(domain = str_remove_all(str_extract(taxon_CI, "(?<=d:)(.*?)\\,"), ","),
         phylum = str_remove_all(str_extract(taxon_CI, "(?<=p:)(.*?)\\,"), ","),
         class = str_remove_all(str_extract(taxon_CI, "(?<=c:)(.*?)\\,"), ","),
         order = str_remove_all(str_extract(taxon_CI, "(?<=o:)(.*?)\\,"), ","),
         family = str_remove_all(str_extract(taxon_CI, "(?<=f:)(.*?)\\,"), ","),
         genus = str_remove_all(str_extract(taxon_CI, "(?<=g:)(.*?)\\,"), ","),
         species = str_remove_all(str_extract(taxon_CI, "(?<=s:).*"), ",")) %>%
  select(-taxon_CI) %>%
  pivot_longer(cols = !otu,
               names_to = "taxon_level",
               values_to = "otu_taxa") %>%
  separate_wider_regex(otu_taxa, c(otu_taxa = ".*", "\\(", otu_taxa_CI = ".*")) %>%
  mutate(otu_taxa_CI = str_remove_all(otu_taxa_CI, pattern = "[^0-9.]"), 
         otu_taxa_CI = na_if(otu_taxa_CI, ""),
         otu_taxa = na_if(otu_taxa, "_NA"))


#### Function to combine emu and sintax output ####

combine_taxonomy <- function(sample_id, sintax_df, emu_df){
  # is sample gnotobiotic
  gnoto <- gnotobiotic[sample_id, 1]
  if(is.na(gnoto)){gnoto = FALSE}
  
  ## Format sintax
  colnames(sintax_df) <- c("read_id", "taxon_CI", "strand", "final_taxa")
  sintax_taxa <- sintax_df[,c(1,2)]
  sintax_taxa <- sintax_taxa %>%
    mutate(domain = str_extract(taxon_CI, "(?<=d:)(.*?)\\,"),
           phylum = str_extract(taxon_CI, "(?<=p:)(.*?)\\,"),
           class = str_extract(taxon_CI, "(?<=c:)(.*?)\\,"),
           order = str_extract(taxon_CI, "(?<=o:)(.*?)\\,"),
           family = str_extract(taxon_CI, "(?<=f:)(.*?)\\,"),
           genus = str_extract(taxon_CI, "(?<=g:)(.*?)\\,"),
           species = str_extract(taxon_CI, "(?<=s:).*")) %>%
    separate_wider_regex(read_id, c(read_id = ".*?", "\\s+", run_info = ".*"))  %>%
    select(-c(taxon_CI, run_info)) 
  
  ## Add otu to sintax and format
  sintax_otu <- left_join(sintax_taxa, otu_df, by = "read_id") %>%
    pivot_longer(cols = !c(read_id, otu),
                 names_to = "taxon_level",
                 values_to = "sintax") %>%
    separate_wider_regex(sintax, c(sintax = ".*", "\\(", sintax_CI = ".*")) %>%
    mutate(sintax_CI = str_remove_all(sintax_CI, pattern = "[^0-9.]"), 
           sintax_CI = na_if(sintax_CI, ""),
           sintax = na_if(sintax, "_NA"))
  
  ## Format emu
  colnames(emu_df)[1] <- "read_id"
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
      domain = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "domain"]),
      phylum = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "phylum"]),
      class = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "class"]),
      order = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "order"]),
      family = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "family"]),
      genus = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "genus"]),
      species = if_else(is.na(emu_tax_id), NA, emu_taxid_taxonomy[as.character(emu_tax_id), "species"])
    ) %>%
    select(-c(n, emu_max_CI, emu_tax_id)) %>%
    left_join(otu_df, by = "read_id") %>%
    pivot_longer(!c(read_id, otu), names_to = "taxon_level", values_to = "emu")
  
  
  
  ## Combine emu and sintax
  all_assignments <- sintax_otu %>%
    full_join(emu_with_taxon, by = c("read_id", "taxon_level", "otu")) %>%
    left_join(otu_taxonomy, b = c("otu", "taxon_level"))
  
  
  # decision tree
  # selecting the final assignment is based on whether the sample is gnotobiotic and then on if it is using a syncom
  # Due to the likelihood of cross contamination, syncom is marked by experiment and not by sample
  #   Path 1 - sample is gnotobiotic
  #   Path 2 - sample is not gnotobiotic but does contain a syncom
  #   Path 3 - sample is not gnotobiotic and does not contain a syncom
  if(gnoto == T){
    all_assignments <- all_assignments %>%
      mutate(final_taxon = case_when((sintax_CI <= 0.3) & !(emu %in% syncom_list) ~ NA,
                                     (sintax == emu) ~ emu,
                                     (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax,
                                     (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa,
                                     (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA,
                                     (sintax != emu) ~ emu,
                                     .default = NA)) 
  }else{
    if(syncom == T){
      all_assignments <- all_assignments %>%
        mutate(final_taxon = case_when((sintax_CI <= 0.3) & !(emu %in% syncom_list) ~ NA,
                                       (sintax == emu) ~ sintax,
                                       (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax,
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa,
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA,
                                       (sintax != emu) & (emu %in% syncom_list) ~ emu ,
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI >= 0.7) ~ sintax,
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax == otu_taxa) ~ sintax,
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu == otu_taxa) ~ emu,
                                       (sintax != emu) & !(emu %in% syncom_list) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu != otu_taxa) ~ NA,
                                       .default = NA)) 
    }else{
      all_assignments <- all_assignments %>%
        mutate(final_taxon = case_when((sintax_CI <= 0.3) ~ NA,
                                       (sintax == emu) ~ sintax,
                                       (is.na(emu)) & (sintax_CI >= 0.7) ~ sintax,
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI >= 0.7) ~ otu_taxa,
                                       (is.na(emu)) & (sintax_CI < 0.7) & (otu_taxa_CI < 0.7) ~ NA,
                                       (sintax != emu) & (sintax_CI >= 0.7) ~ sintax,
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax == otu_taxa) ~ sintax,
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu == otu_taxa) ~ emu,
                                       (sintax != emu) & (sintax_CI < 0.7) & (sintax != otu_taxa) & (emu != otu_taxa) ~ NA,
                                       .default = NA))
    }
  }

  
  #pivot out and keep the lowest level that's not NA then assign taxid to that level
  final_assignments <- all_assignments %>%
    select(c(read_id, taxon_level, final_taxon, otu)) %>%
    pivot_wider(names_from = taxon_level, values_from = final_taxon) %>%
    mutate(lowest_taxon = case_when(!is.na(species) ~ species,
                                    !is.na(genus) ~ genus,
                                    !is.na(family) ~ family,
                                    !is.na(order) ~ order,
                                    !is.na(class) ~ class,
                                    !is.na(phylum) ~ phylum,
                                    !is.na(domain) ~ domain,
                                    .default = NA),
           lowest_taxon_level = case_when(!is.na(species) ~ "species",
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
    summarise(read_counts = n()) %>%
    mutate(sample_id = sample_id) %>%
    separate_wider_regex(lowest_taxon, c(lowest_taxon = ".*", "_", lowest_taxid = ".*"))
    
  return(final_assignments)
}




#### Loop through all files EDIT HERE####

## nested for loops, 1 through 4 and 1 through 72 to get sample_id for each
## Change the numbers to match yours, might need to modify the file names

for(plate_num in 1:4){
  for(barcode in 1:72){
    if(file.exists(paste0(emu_path, plate_num, "-HL0", barcode, "_read-assignment-distributions.tsv"))){
      emu_file <- read_tsv(paste0(emu_path, plate_num, "-HL0", barcode, "_read-assignment-distributions.tsv"))
      sintax_file <- read_tsv(paste0(sintax_path, plate_num, "-", barcode, "_sintax.tsv"), col_names = FALSE)
      sample_id <- paste0(plate_num, "-HL0", barcode)
      output_df <- combine_taxonomy(sample_id, sintax_file, emu_file)
      if(plate_num == 1 & barcode == 1){
        sample_id_df <- output_df
      }else{
        sample_id_df <- rbind(sample_id_df, output_df)
      }  
    }
  }
}


#### Convert to microeco format ####

#table with taxid + otu
final_taxid_OTU_table <- sample_id_df %>%
  mutate(tax_id = case_when(is.na(lowest_taxon) ~ NA,
                            is.na(otu) ~ paste("taxid", lowest_taxid, sep = "_"),
                            (lowest_taxon_level == "species") ~ paste("taxid", lowest_taxid, sep = "_"),
                            !is.na(otu) ~ paste("taxid", lowest_taxid, otu, sep = "_"))) %>%
  ungroup() %>%
  select(c(tax_id, read_counts, sample_id)) %>%
  pivot_wider(names_from = sample_id, values_from = read_counts)

write_tsv(final_taxid_OTU_table, "microeco_taxid_OTU_table.tsv", col_names = T, eol = "\n")



#table with taxid only
final_taxid_table <- sample_id_df %>%
  mutate(tax_id = paste("taxid", lowest_taxid, sep = "_")) %>%
  ungroup() %>%
  group_by(sample_id, tax_id) %>%
  summarize(counts = sum(read_counts)) %>%
  pivot_wider(names_from = sample_id, values_from = counts)

write_tsv(final_taxid_table, "microeco_taxid_table.tsv", col_names = T, eol = "\n")

