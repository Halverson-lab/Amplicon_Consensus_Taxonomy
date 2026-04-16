#!/usr/bin/env Rscript

library(readr)
library(tidyverse)
library(multidplyr)
options("scipen"=999)

#script to take in the minimap2 file and identify and label reads where the species can't be differentiated
args <- commandArgs(trailingOnly = TRUE)
readRenviron(args[1])
minimap_stats <- read_csv(args[2])

#### Environment setup ####
path_to_data_dir <- Sys.getenv("DATABASE_DIR")
setwd(path_to_data_dir)

n_cores <- Sys.getenv("N_CORE")

#read in taxonomy info since samtools can't handle full headers
taxonomy <- read_tsv("taxonomy.tsv", col_types = cols(.default = "c")) %>%
  distinct() %>%
  mutate(genus = if_else(genus == "Shigella", "Escherichia", genus),
         t_genus = if_else(t_genus == "620", "561", t_genus))

taxonomy_w_rownames <- taxonomy %>%
  column_to_rownames(var = "tax_id")

#if a read has more than 3 NA values in its taxonomy then toss it as problematic
no_tax_reads <- taxonomy %>%
  rowwise() %>%
  mutate(na_count = sum(is.na(c_across(phylum:species)))) %>%
  filter(na_count > 3)

no_tax_read_list <- append(no_tax_reads$tax_id, no_tax_reads$t_species)

#get a count of how many reads are in the database for each species
raw_seq_ids <- read_csv(paste0(path_to_data_dir, "seq_id_list.txt"), col_names = F)

all_seq_ids <- raw_seq_ids %>%
  separate_wider_delim(X1, delim = ":", names = c("tax_id", "db", "seqid")) %>%
  mutate(seq_name = paste(tax_id, db, seqid, sep = ":")) %>%
  left_join(taxonomy) %>%
  group_by(species) %>%
  mutate(total_species_seqs = n())

species_counts <- all_seq_ids %>%
  group_by(species, t_species) %>%
  summarize(total_species_seqs = n()) %>%
  filter(!is.na(t_species)) %>%
  column_to_rownames(var = "t_species")

#create a cluster
cluster <- new_cluster(as.integer(n_cores)-1)
#pass tidyverse and the emu_taxid_taxonomy table to every cluster
cluster_library(cluster, "tidyverse")
cluster_assign(cluster, taxonomy_w_rownames = taxonomy_w_rownames)
cluster_assign(cluster, species_counts = species_counts)

#partition the dataframe between clusters
minimap_part <- minimap_stats %>% group_by(QTAXID) %>% partition(cluster)

minimap_stats <- minimap_part %>%
  mutate(QTAXID = str_extract(QNAME, "[0-9].*?(?=:)"),
         RTAXID = str_extract(RNAME, "[0-9].*?(?=:)"),
         Q_T_SPECIES = taxonomy_w_rownames[as.character(QTAXID), "t_species"],
         R_T_SPECIES = taxonomy_w_rownames[as.character(RTAXID), "t_species"],
         QSPECIES = taxonomy_w_rownames[as.character(QTAXID), "species"],
         RSPECIES = taxonomy_w_rownames[as.character(RTAXID), "species"],
         QGENUS = taxonomy_w_rownames[as.character(QTAXID), "genus"],
         RGENUS = taxonomy_w_rownames[as.character(RTAXID), "genus"],
         QSPECIES_COUNT = species_counts[as.character(Q_T_SPECIES), "total_species_seqs"],
         RSPECIES_COUNT = species_counts[as.character(R_T_SPECIES), "total_species_seqs"]) %>%
  collect() %>%
  arrange(QSPECIES, RSPECIES)

minimap_failed_stats <- minimap_stats %>%
  filter(QSPECIES != RSPECIES) %>%
  filter(!(QTAXID %in% no_tax_read_list | RTAXID %in% no_tax_read_list))

problem_reads <- minimap_failed_stats %>%
  filter(QGENUS != RGENUS) %>%
  ungroup() %>%
  group_by(QNAME) %>%
  mutate(QNAME_freq_count = n()) %>%
  ungroup() %>%
  group_by(RNAME) %>%
  mutate(RNAME_freq_count = n()) %>%
  ungroup()

bad_reads_Q <- problem_reads %>%
  select(c(QNAME, QNAME_freq_count, QSPECIES_COUNT)) %>%
  distinct() %>%
  filter(QNAME_freq_count > 1) %>%
  rename(name=QNAME, freq_count=QNAME_freq_count, species_count=QSPECIES_COUNT)

bad_reads_R <- problem_reads %>%
  select(c(RNAME, RNAME_freq_count, RSPECIES_COUNT)) %>%
  distinct() %>%
  filter(RNAME_freq_count > 1) %>%
  rename(name=RNAME, freq_count=RNAME_freq_count, species_count=RSPECIES_COUNT)

bad_reads <- rbind(bad_reads_Q, bad_reads_R) %>%
  mutate(tax_id = str_extract(name, "[0-9].*?(?=:)"),
         t_species = taxonomy_w_rownames[as.character(tax_id), "t_species"],
         species = taxonomy_w_rownames[as.character(tax_id), "species"]) %>%
  group_by(name) %>%
  mutate(total_freq_count = sum(freq_count)) %>%
  select(-freq_count) %>%
  distinct() %>%
  group_by(species) %>%
  mutate(bad_count_in_species = length(unique(name)),
         percent_bad = bad_count_in_species/species_count) %>%
  filter(percent_bad < 0.05)


minimap_failed_stats_trimmed <- minimap_failed_stats %>%
  filter(!(QNAME %in% bad_reads$name | RNAME %in% bad_reads$name)) %>%
  ungroup()


taxid_pairs <- distinct(select(minimap_failed_stats_trimmed, c(Q_T_SPECIES, R_T_SPECIES))) %>%
  rename(QSPECIES=Q_T_SPECIES, RSPECIES=R_T_SPECIES)

group_list <- c()
index = 1

for (n in 1:nrow(taxid_pairs)) {
  # check to see if query is already in a group, add it if not
  q_id <- taxid_pairs[[n, "QSPECIES"]]
  r_id <- taxid_pairs[[n, "RSPECIES"]]
  
  if (n == 1){
    #initialize the list
    group_list[[paste0("group_", index)]] <- c(q_id, r_id)
  } else {
    existing_group = FALSE
    group_index = 0
    #see if the query or reference are in the lists already and record the group if it is
    for (i in seq_along(group_list)){
      if (q_id %in% group_list[[i]] | r_id %in% group_list[[i]]) {
        existing_group = TRUE
        group_index = i
      }
    }
    
    # if it was found, then add to existing group. otherwise add a new group
    if (existing_group) {
      group_list[[group_index]] <- c(group_list[[group_index]], q_id, r_id)
      group_list[[group_index]] <- unique(group_list[[group_index]])
    } else {
      index = index + 1
      group_list[[paste0("group_", index)]] <- c(q_id, r_id)
    }
  }
}


#extra check to see if there are any matches between groups and dereplicate if there are
derep_group_list <- c()
n <- 1
while ( n <= length(group_list) ){
  merge_group = FALSE
  i <- n + 1
  group_name <- paste0("group_", n)
  while ( i <= length(group_list) ) {
    #check for match
    elements_in_both <- intersect(group_list[[n]], group_list[[i]])
    # if there is a match
    if(length(elements_in_both) > 0){
      #merge the two groups
      derep_group_list[[n]] <- c(group_list[[n]], group_list[[i]])
      derep_group_list[[n]] <- unique(derep_group_list[[n]])
      merge_group <- TRUE
      
      #remove the group that got merged
      group_list[[n]] <- derep_group_list[[n]]
      group_list[[i]] <- NULL
      
      #adjust the count for the missing group
      i <- i - 1
    }
    i <- i + 1
  }
  if (merge_group == FALSE) {
    derep_group_list[[n]] <- group_list[[n]]
  }
  
  n <- n + 1
}
#name the new dereplicated groups
names(derep_group_list) <- paste0("group_", seq_along(derep_group_list))


#retrieve a list of all read ids and their taxid 
read_list <- minimap_stats[,c("QNAME", "QTAXID", "Q_T_SPECIES")] %>% 
  rename(seq_id=QNAME, taxid=QTAXID, t_species=Q_T_SPECIES) %>%
  rbind(rename(minimap_stats[,c("RNAME", "RTAXID", "R_T_SPECIES")], c(seq_id=RNAME, taxid=RTAXID, t_species=R_T_SPECIES))) %>%
  distinct() 


#add the group identifiers to this list
for (n in 1:nrow(read_list)) {
  taxid_id <- read_list[[n, "t_species"]]
  for (i in seq_along(derep_group_list)){
    if (taxid_id %in% derep_group_list[[i]]) {
      read_list[n, "group"] <- names(derep_group_list[i])
    }
  }
}


read_list_with_tax <- read_list %>%
  filter(!is.na(group)) %>%
  mutate(taxid = as.character(taxid)) %>%
  left_join(rename(taxonomy, taxid=tax_id))


#bad seqs in group
group_outliers <- read_list_with_tax %>%
  group_by(group) %>%
  mutate(total_seq_in_group = n()) %>%
  group_by(group, species) %>%
  mutate(n_reads_species = n(),
         percent_read_species = n_reads_species/total_seq_in_group) %>%
  select(-c(seq_id, taxid)) %>%
  distinct() %>%
  group_by(group, order) %>%
  mutate(percent_reads_order = sum(n_reads_species)/total_seq_in_group) %>%
  ungroup() %>%
  group_by(group, class) %>%
  mutate(percent_reads_class = sum(n_reads_species)/total_seq_in_group) %>%
  ungroup() %>%
  group_by(group, phylum) %>%
  mutate(percent_reads_phylum = sum(n_reads_species)/total_seq_in_group) %>%
  filter(percent_reads_order <= 0.01 | percent_reads_class <= 0.01 | percent_reads_phylum <= 0.01) %>%
  ungroup() %>%
  select(-c(total_seq_in_group, n_reads_species, percent_read_species, percent_reads_order, percent_reads_class , percent_reads_phylum))

read_list_with_tax_outliers_removed <- read_list_with_tax %>%
  anti_join(group_outliers)

group_outliers_removed <- anti_join(read_list_with_tax, read_list_with_tax_outliers_removed) %>%
  ungroup() %>%
  select(seq_id) %>%
  rename(name = seq_id)

#make group names based on taxonomy
rename_groups <- read_list_with_tax_outliers_removed %>%
  group_by(group) %>%
  mutate(across(where(is.character), ~ na_if(. , ""))) %>%
  mutate(n_uniq_domain=n_distinct(domain),
         n_uniq_phylum=n_distinct(phylum),
         n_uniq_class=n_distinct(class),
         n_uniq_order=n_distinct(order),
         n_uniq_family=n_distinct(family),
         n_uniq_genus=n_distinct(genus),
         n_uniq_species=n_distinct(species)) %>%
  mutate(group_taxa = case_when((n_uniq_genus < 4) & (!is.na(genus)) ~ genus,
                                (n_uniq_family < 4) & (!is.na(family)) ~ family,
                                (n_uniq_order < 4) & (!is.na(order)) ~ order,
                                (n_uniq_class < 4) & (!is.na(class)) ~ class,
                                (n_uniq_phylum < 4) & (!is.na(phylum)) ~ phylum),
         group_level = case_when(n_uniq_genus == 1 ~ "s",
                                 n_uniq_family == 1 ~ "g",
                                 n_uniq_order == 1 ~ "f",
                                 n_uniq_class == 1 ~ "o",
                                 n_uniq_phylum == 1 ~ "c",
                                 .default = "p")) %>%
  select(c(group, group_taxa, group_level)) %>%
  distinct() %>%
  group_by(group, group_level) %>%
  arrange(group_taxa) %>%
  summarize(group_name = paste(na.omit(group_taxa), collapse = "-")) %>%
  mutate(info_group_name = paste0(group,"_",group_level,"_", group_name))


read_list_with_new_group_name <- read_list_with_tax_outliers_removed %>%
  left_join(rename_groups) %>%
  ungroup()

# make one file for changing the headers of the relevant reads
# a file of which reads to remove 
# another file with the modified taxonomy, whatever that ends up being
# third is a file with all of the groups and who is in which (needs to be human and machine readable, so prob a csv)

#list of seqs to remove
seqs2remove <- bad_reads %>%
  ungroup() %>%
  select(name) %>%
  rbind(group_outliers_removed)

write.table(seqs2remove, paste0(path_to_data_dir, "sequences_to_be_removed.csv"), append = FALSE, sep = ",", dec = ".",
            row.names = F, col.names = F, quote = FALSE)

# function to convert to large number
trailing_zeros <- function(old_tax){
  taxid <- as.integer(old_tax)
  n_digits <- floor(log10(taxid)) + 1
  new_num <- (taxid * 10^(14 - n_digits))
  return(new_num)
}


#tsv for changing taxids
old2newIds <- read_list_with_new_group_name %>%
  mutate(group_num = as.integer(str_extract(group, "[0-9]*?$")),
         new_taxid = as.character(trailing_zeros(taxid) + group_num),
         new_id = seq_id) %>%
  separate_wider_delim(new_id, ":", names = c("old_tax", "new_db", "new_seq")) %>%
  unite("new_id", c(new_taxid, new_db, new_seq), sep = ":" ) %>%
  select(c(seq_id, new_id))

write_tsv(old2newIds, paste0(path_to_data_dir, "grouped-seq-headers.tsv"), col_name = FALSE, escape = "none")

#### human readable table of who is in which groups ####

#seq_id to group
seq_id_group <- read_list_with_new_group_name %>%
  select(c(seq_id, species, group, group_name)) %>%
  arrange(group, species)

write_csv(seq_id_group, paste0(path_to_data_dir, "seq_id_to_group_list.csv"))

#groups with species
group_to_species <- read_list_with_new_group_name %>%
  select(c(group, group_name, species)) %>%
  distinct() %>%
  arrange(species, group)

write_csv(group_to_species, paste0(path_to_data_dir, "group_to_species_list.csv"))


#### new taxonomies ####
#taxonomies to append to file
new_taxonomies <- read_list_with_new_group_name %>%
  mutate(group_num = as.integer(str_extract(group, "[0-9]*?$")),
         new_taxid = as.character(trailing_zeros(taxid) + group_num)) %>%
  select(c(taxid, domain, phylum, class, order, family, genus, species,
           t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species,
           new_taxid, group, group_level, group_name, group_num, info_group_name)) %>%
  mutate(across(domain:species, ~ na_if(. , ""))) %>%
  distinct() 

# setup new taxonomies
for(n in 1:nrow(new_taxonomies)){
  if (new_taxonomies[[n, "group_level"]] == "s"){
    new_taxonomies$species[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_genus[n])){
      new_taxonomies$t_species[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_species[n] <- as.character(trailing_zeros(new_taxonomies$t_genus[n]) + new_taxonomies[[n, "group_num"]])
    }
  } else if (new_taxonomies[[n, "group_level"]] == "g"){
    new_taxonomies$genus[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_family[n])){
      new_taxonomies$t_genus[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_genus[n] <- as.character(trailing_zeros(new_taxonomies$t_family[n]) + new_taxonomies[[n, "group_num"]])
    }
    for (taxa in c("species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "f"){
    new_taxonomies$family[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_order[n])){
      new_taxonomies$t_family[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_family[n] <- as.character(trailing_zeros(new_taxonomies$t_order[n]) + new_taxonomies[[n, "group_num"]])
    }
    for (taxa in c("genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "o"){
    new_taxonomies$order[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_class[n])){
      new_taxonomies$t_order[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_order[n] <- as.character(trailing_zeros(new_taxonomies$t_class[n]) + new_taxonomies[[n, "group_num"]])
    }
    for (taxa in c("family","genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "c"){
    new_taxonomies$class[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_phylum[n])){
      new_taxonomies$t_class[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_class[n] <- as.character(trailing_zeros(new_taxonomies$t_phylum[n]) + new_taxonomies[[n, "group_num"]])
    }
    for (taxa in c("order", "family","genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "p"){
    new_taxonomies$phylum[n] <- paste(new_taxonomies[[n, "group_name"]], "spp", new_taxonomies[[n, "group"]], sep = "_")
    if(is.na(new_taxonomies$t_domain[n])){
      new_taxonomies$t_phylum[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_phylum[n] <- as.character(trailing_zeros(new_taxonomies$t_domain[n]) + new_taxonomies[[n, "group_num"]])
    }
    for (taxa in c("class", "order", "family","genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } 
}

new_taxonomies <- new_taxonomies %>%
  select(c(new_taxid, domain, phylum, class, order, family, genus, species,
           t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species)) %>%
  rename(tax_id=new_taxid)

write_tsv(new_taxonomies, paste0(path_to_data_dir, "grouped_taxonomies.tsv"))
