#!/usr/bin/env Rscript

library(readr)
library(tidyverse)
options("scipen"=999)

#script to take in the minimap2 file and identify and label reads where the species can't be differentiated
args <- commandArgs(trailingOnly = TRUE)
minimap_stats <- read_csv(args[1])

#read in taxonomy info since samtools can't handle full headers
taxonomy <- read.csv("taxonomy.tsv", sep = "\t") %>%
  distinct() %>%
  column_to_rownames(var = "tax_id")

minimap_stats <- minimap_stats %>%
  mutate(QSPECIES = taxonomy[as.character(QTAXID), "t_species"],
         RSPECIES = taxonomy[as.character(RTAXID), "t_species"]) %>%
  filter(QSPECIES != RSPECIES)


taxid_pairs <- distinct(select(minimap_stats, c(QSPECIES, RSPECIES)))
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
read_list <- minimap_stats[,c("QNAME", "QTAXID", "QSPECIES")] %>% 
  rename(seq_id=QNAME, taxid=QTAXID, t_species=QSPECIES) %>%
  rbind(rename(minimap_stats[,c("RNAME", "RTAXID", "RSPECIES")], c(seq_id=RNAME, taxid=RTAXID, t_species=RSPECIES))) %>%
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

taxonomy <- taxonomy %>% rownames_to_column(var = "taxid")

read_list_with_tax <- read_list %>%
  mutate(taxid = as.character(taxid)) %>%
  left_join(taxonomy)


#make group names based on taxonomy
rename_groups <- read_list_with_tax %>%
  group_by(group) %>%
  mutate(across(where(is.character), ~ na_if(. , ""))) %>%
  mutate(n_uniq_domain=n_distinct(domain),
         n_uniq_phylum=n_distinct(phylum),
         n_uniq_class=n_distinct(class),
         n_uniq_order=n_distinct(order),
         n_uniq_family=n_distinct(family),
         n_uniq_genus=n_distinct(genus),
         n_uniq_species=n_distinct(species)) %>%
  mutate(group_taxa = case_when(n_uniq_genus < 4 ~ genus,
                                n_uniq_family < 4 ~ family,
                                n_uniq_order < 4 ~ order,
                                n_uniq_class < 4 ~ class,
                                n_uniq_phylum < 4 ~ phylum),
         group_level = case_when(n_uniq_genus == 1 ~ "s",
                                 n_uniq_family == 1 ~ "g",
                                 n_uniq_order == 1 ~ "f",
                                 n_uniq_class == 1 ~ "o",
                                 n_uniq_phylum == 1 ~ "c",
                                 .default = "p")) %>%
  select(c(group, group_taxa, group_level)) %>%
  distinct() %>%
  group_by(group, group_level) %>%
  summarize(group_name = paste(group_taxa, collapse = "-")) %>%
  mutate(info_group_name = paste0(group,"_",group_level,"_", group_name))
  

read_list_with_new_group_name <- read_list_with_tax %>%
  left_join(rename_groups)



# make one file for changing the headers of the relevant reads
# another file with the modified taxonomy, whatever that ends up being
# third is a file with all of the groups and who is in which (needs to be human and machine readable, so prob a csv)


#tsv for changing taxids
old2newIds <- read_list_with_new_group_name %>%
  mutate(new_taxid = paste0(taxid, "m", str_extract(group, "[0-9]*?$")),
         new_id = seq_id) %>%
  separate_wider_delim(new_id, ":", names = c("old_tax", "new_db", "new_seq")) %>%
  unite("new_id", c(new_taxid, new_db, new_seq), sep = ":" ) %>%
  select(c(seq_id, new_id))

write_tsv(old2newIds, "grouped-seq-headers.tsv", col_name = FALSE, escape = "none")


#### new taxonomies ####
#taxonomies to append to file
new_taxonomies <- read_list_with_new_group_name %>%
  mutate(new_taxid = paste0(taxid, "m", str_extract(group, "[0-9]*?$"))) %>%
  select(c(taxid, domain, phylum, class, order, family, genus, species,
           t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species,
           new_taxid, group, group_level, group_name, info_group_name)) %>%
  mutate(across(domain:species, ~ na_if(. , ""))) %>%
  distinct() 

# setup new taxonomies
for(n in 1:nrow(new_taxonomies)){
  if (new_taxonomies[[n, "group_level"]] == "s"){
    new_taxonomies$species[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_genus[n])){
      new_taxonomies$t_species[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_species[n] <- paste0(new_taxonomies$t_genus[n],"m",str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
    }
  } else if (new_taxonomies[[n, "group_level"]] == "g"){
    new_taxonomies$genus[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_family[n])){
      new_taxonomies$t_genus[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_genus[n] <- paste0(new_taxonomies$t_family[n],"m", str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
    }
    for (taxa in c("species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "f"){
    new_taxonomies$family[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_order[n])){
      new_taxonomies$t_family[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_family[n] <- paste0(new_taxonomies$t_order[n],"m", str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
    }
    for (taxa in c("genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "o"){
    new_taxonomies$order[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_class[n])){
      new_taxonomies$t_order[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_order[n] <- paste0(new_taxonomies$t_class[n],"m", str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
    }
    for (taxa in c("family","genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "c"){
    new_taxonomies$class[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_phylum[n])){
      new_taxonomies$t_class[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_class[n] <- paste0(new_taxonomies$t_phylum[n],"m", str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
    }
    for (taxa in c("order", "family","genus", "species")) {
      new_taxonomies[[n, taxa]] <- NA
      new_taxonomies[[n, paste0("t_",taxa)]] <- NA
    }
  } else if (new_taxonomies[[n, "group_level"]] == "p"){
    new_taxonomies$phylum[n] <- paste(new_taxonomies[[n, "group"]], new_taxonomies[[n, "group_name"]], "sp", sep = "_")
    if(is.na(new_taxonomies$t_domain[n])){
      new_taxonomies$t_phylum[n] <- new_taxonomies$new_taxid[n]
    } else {
      new_taxonomies$t_phylum[n] <- paste0(new_taxonomies$t_domain[n],"m", str_extract(new_taxonomies[[n, "group"]], "[0-9]*?$"))
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

write_tsv(new_taxonomies, "grouped_taxonomies.tsv")



#### human readable table of who is in which groups ####

#seq_id to group
seq_id_group <- read_list_with_new_group_name %>%
  select(c(seq_id, species, group, group_name))

write_csv(seq_id_group, "seq_id_to_group_list.csv")

#groups with species
group_to_species <- read_list_with_new_group_name %>%
  select(c(group, group_name, species)) %>%
  distinct() 

write_csv(group_to_species, "group_to_species_list.csv")
