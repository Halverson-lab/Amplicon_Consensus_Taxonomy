# Recreating the Emu database
The steps were briefly addressed on [github](https://github.com/treangenlab/emu/issues/15).

    1. Download NCBI taxonomy for the names.dmp file and the nodes.dmp file.
    2. Download rrn fasta and NCBI 16S fasta. 
    3. Create a seq2taxid map file. 
    4. Run the emu build-database function 

### Conda env

```
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - taxonkit
  - csvtk
  - seqkit
  - minimap2
  - samtools
  - entrez-direct
```

## Download all the files needed
Use NCBI [FTP](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) to download the latest taxdump file and save it somewhere you can find it, I saved it directly into my taxonomy environment. 

Download the latest [rrnDB](https://rrndb.umms.med.umich.edu/downloads/).

Downloading the RefSeq database is more of a pain because you have to do it through BLAST.

```bash
#download database
 update_blastdb.pl --decompress 16S_ribosomal_RNA

#retrieve fasta
blastdbcmd -entry all -db <path_to_db>/16S_ribosomal_RNA -out ncbi_16S.fasta
```

## Create the seq2tax file

NCBI RefSeq and rrnDB use different style headers, so I handled them separately before merging.

### rrnDB

```bash

# For later application we need fasta headers starting with useful sequence IDs
# get fasta headers
seqkit seq -n rrnDB.fasta > rrnDB_header_list.txt

# get the RefSeq numbers from each header and generate new headers
cat rrnDB_header_list.txt \
    | csvtk add-header -t -n fasta_header \
    | csvtk mutate -t -f 1 -p "(?:[^\|]*\|){2}(.*?)\|" -n refSeq \
    | csvtk mutate2 -t -n newHeader -e '$refSeq + " " + $fasta_header'> rrnDB_header_refseq_list.txt


# File of rrnDB header to new rrnDB header with RefSeq ID as the sequence ID
cat rrnDB_header_refseq_list.txt | csvtk cut -t -f 1,3 | csvtk del-header -t > old2newheader.tsv
# replace the headers
seqkit replace -p '(.+)' -r '{kv}' -k old2newheader.tsv rrnDB.fasta > rrnDB_renamed.fasta


# Use refSeq numbers to get taxid
cat rrnDB_header_refseq_list.txt | csvtk cut -t -f 2 | csvtk del-header -t | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > rrnDB_seq2taxid.txt

# Filter out any of the fasta seq that were removed
csvtk cut -t -f 1 rrnDB_seq2taxid.txt > seq2keep.txt
seqkit grep -f seq2keep.txt rrnDB_renamed.fasta -o rrnDB_{date}.fasta
```

`rrnDB_{date}.fasta` and `rrnDB_seq2taxid.txt` are the two files you need from this. 


### NCBI 16S

```bash
#get fasta headers
seqkit seq -n ncbi_16S_{date}.fasta > ncbi_16S_header.txt

#get the RefSeq numbers from each header
cat ncbi_16S_header.txt | csvtk add-header -t -n fasta_header | csvtk mutate -t -f 1 -p "([^\s]+)" -n refSeq > ncbi_16S_header_refseq.txt

#Use refSeq numbers to get taxid
cat ncbi_16S_header_refseq.txt | csvtk cut -t -f 2 | csvtk del-header -t | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > refseq2taxid.txt

csvtk add-header -t -n refSeq,taxid refseq2taxid.txt > refseq2taxid2.txt

# Any seqs where the RefSeq record was suppressed will be removed from the list (this shouldn't change anything as long as your fasta file is current)
csvtk join -t -f refSeq ncbi_16S_header_refseq.txt refseq2taxid2.txt | csvtk cut -t -f 1,3 |csvtk del-header -t > ncbi_seq2taxid.txt
```
`ncbi_16S_{date}.fasta` and `ncbi_seq2taxid.txt` are the next two files.


Combine the corresponding files for each database to get the final seq2taxid and fasta files.
```bash
cat rrnDB_seq2taxid.txt ncbi_seq2taxid.txt > seq2taxid.txt

cat rrnDB_{date}.fasta ncbi_16S_{date}.fasta > 16S_seq.fasta
```
## Create the taxonomy file

You can technically build the database from here with the NCBI taxonomy files, but I prefer to build my own taxonomy table so its formatted how I want. This will produce a taxonomy that goes 

    - kingdom
    - phylum
    - class
    - order
    - family
    - genus
    - species
    - kingdom_taxid
    - phylum_taxid
    - class_taxid
    - order_taxid
    - family_taxid
    - genus_taxid
    - species_taxid

Having the taxids helps a lot later.

```bash
# Tell taxonkit where to find taxdump
export TAXONKIT_DB=/work/LAS/larryh-lab/Microbiome_AP/taxonomy-env/taxdump

# Build taxonomy.tsv
csvtk cut -t -f 2 seq2taxid.txt \
    | taxonkit lineage \
    | taxonkit reformat -t -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
    | csvtk cut -t -f -2 \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species,t_kingdom,t_phylum,t_class,t_order,t_family,t_genus,t_species -o taxonomy.tsv
```

## Build the database
Finally you can build the Emu database
```bash
emu build-database <db_name> --sequences 16S_seq.fasta --seq2tax seq2taxid.txt --taxonomy-list taxonomy.tsv
```


## Create the sintax database
Convert the previous EMU database to Sintax format 

Retrieve all the sequence headers from the EMU database fasta file.
```bash
seqkit seq species_taxid.fasta -n > name_list.txt
```

Take the list of sequence headers and the taxonomy list and use them to generate new headers in sintax format. You can probably do this in python instead, but I'm most familiar with R.
```R
library(tidyverse)
library(readr)
library(gsubfn)

setwd("//las-dfs-01.las.iastate.edu/lss/research/larryh-lab/Ashley-P/Microbiome_sequencing/Database_20241101")

raw_emu_header <- read.pattern("name_list.txt", pattern = "^(\\S+) +(.*)$", quote = "")

emu_headers <- raw_emu_header %>%
  separate_wider_delim(V1, ":", names = c("tax_id", "db", "emu_num"))
colnames(emu_headers)[which(names(emu_headers) == "V2")] <- "header"

taxonomy <- read.csv("emuDB_20241101/taxonomy.tsv", sep = "\t")
colnames(taxonomy)[which(names(taxonomy) == "kingdom")] <- "domain"
colnames(taxonomy)[which(names(taxonomy) == "t_kingdom")] <- "t_domain"


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
  unite("lineage", domain:species, sep = ",") %>%
  mutate(lineage = str_replace_all(lineage, " ", "_")) %>%
  distinct() %>%
  column_to_rownames(var = "tax_id")


fasta_lineage <- emu_headers %>%
  mutate(lineage = cat_taxon[as.character(tax_id), "lineage"])


old2new_header <- fasta_lineage %>%
  mutate(emu_id = paste(tax_id, db, emu_num, sep = ":"),
         sintax_id = paste(tax_id, db, emu_num, sep = "_")) %>%
  unite("emu_header", c(emu_id, header), sep = " ") %>%
  unite("sintax_header", c(sintax_id, lineage), sep = ";tax=") %>%
  select(-c(tax_id, db, emu_num))

write_tsv(old2new_header, "emu2sintax-header.tsv", col_name = FALSE, escape = "none")
```

Use `emu2sintax-header.tsv` as a list of key-value pairs to rename the EMU fasta file with the new headers and the output is the sintax database.

```bash
seqkit replace -p "(.+)" -r '{kv}' -k emu2sintax-header.tsv emu_db.fasta > sintax_db.fasta
```

You don't have to do anything else for the sintax database from here, its just needs the fasta file with properly formatted headers.



