#!/usr/bin/env Rscript

library(tidyverse)
library(readr)
library(gsubfn)


raw_emu_header <- read.pattern("name_list.txt", pattern = "^(\\S+) +(.*)$", quote = "")
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", , na.strings=c("","NA"))

emu_headers <- raw_emu_header %>%
  separate_wider_delim(V1, ":", names = c("tax_id", "db", "emu_num"))
colnames(emu_headers)[which(names(emu_headers) == "V2")] <- "header"


cat_taxon <- taxonomy %>%
  mutate(
    domain = if_else(is.na(domain), NA, paste0("d:", domain, "_", t_domain)),
    phylum = if_else(is.na(phylum), NA, paste0("p:", phylum, "_", t_phylum)),
    class = if_else(is.na(class), NA, paste0("c:", class, "_", t_class)),
    order = if_else(is.na(order), NA, paste0("o:", order, "_", t_order)),
    family = if_else(is.na(family), NA, paste0("f:", family, "_", t_family)),
    genus = if_else(is.na(genus), NA, paste0("g:", genus, "_", t_genus)),
    species = if_else(is.na(species), NA, paste0("s:", species, "_", t_species))
    ) %>%
  unite("lineage", domain:species, sep = ",", na.rm = TRUE) %>%
  mutate(lineage = str_replace_all(lineage, " ", "_")) %>%
  distinct() %>%
  column_to_rownames(var = "tax_id")


fasta_lineage <- emu_headers %>%
  mutate(lineage = cat_taxon[as.character(tax_id), "lineage"])


old2new_header <- fasta_lineage %>%
  mutate(emu_id = paste(tax_id, db, emu_num, sep = ":"),
         sintax_id = paste(tax_id, db, emu_num, sep = "_")) %>%
  unite("sintax_header", c(sintax_id, lineage), sep = ";tax=") %>%
  mutate(sintax_header = paste0(sintax_header, ";")) %>%
  select(c(emu_id, sintax_header))

write_tsv(old2new_header, "emu2sintax-header.tsv", col_name = FALSE, escape = "none")
