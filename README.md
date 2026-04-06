# The Amplicon Consensus Taxonomy (ACT) pipeline & similarity-aware reference database (ACT-DB)

## DESCRIPTION
The **ACT pipeline** integrates outputs from three tools - [Emu](https://github.com/treangenlab/emu), [Sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html), and [LACA](https://github.com/yanhui09/laca) - to generate consensus taxonomic assignments for ONT long-read amplicons. Specifically, ACT cross-validates taxonomic calls from Emu and Sintax to increase classification confidence for reads from known species while reducing overclassification and exclusion of reads from rare and unknown taxa, and LACA preserves fine-level detail for unknown species by clustering sequences into OTUs.  


**ACT-DB** is a sequence-similarity-aware reference database of bacterial and archaeal 16S gene sequences from the [NCBI 16S rRNA RefSeq database](https://www.ncbi.nlm.nih.gov/refseq/) and the prokaryotic [Ribosomal RNA Operon Copy Number Database](https://rrndb.umms.med.umich.edu/) (rrnDB). Unlike traditional databases, ACT-DB clusters 16S sequences sharing >99.5% identity into named multi-taxa groups that are assigned in place of individual species. This approach reduces misclassification when 16S sequences cannot reliably distinguish closely related taxa, while preserving cases where other 16S copies provide species-level resolution. 

Although originally created for profiling full-length 16S amplicons, ACT will work for any size 16S amplicon if the pipeline parameters are adjusted. However, profiling using partial amplicons should be conducted using a region-specific version of ACT-DB, as partial and full-length amplicons differ in information content and similarity profiles. ACT can also be used for amplicon profiling using other marker regions (e.g., fungal ITS) if an appropriately formatted database is supplied. Instructions for updating the default ACT-DB, modifying the database, or building a custom database (e.g., partial 16S amplicons, fungal ITS) are provided [here](#updating-and-modifying-act-db-and-building-custom-ambiguity-aware-databases).

**WARNING:** To analyze amplicons from multiple marker regions (e.g., 16S and ITS), the pipeline needs to be run independently for each region using the appropriate database.

## ACT SETUP
* ACT should be run on an HPC system for optimal performance, but it can also be executed locally if sufficient resources are available.
* **REQUIREMENT:** ACT execution scripts require that a conda-compatible package manager (conda, mamba, or micromamba) is installed and accessible in the user’s HPC environment (e.g., available on $PATH).
```bash
# verify accessibility by running a version check in your shell
conda --version
mamba --version
micromamba --version
```

### 1. Clone the ACT GitHub repository
Clone the entire repository into the directory where you want to perform the analysis (e.g., /work/user/). All ACT steps will be executed within this folder.

```bash
git clone https://github.com/Halverson-lab/Amplicon_Consensus_Taxonomy.git
cd ./Amplicon_Consensus_Taxonomy
```

### 2. Modify the necessary files
The `Amplicon_Consensus_Taxonomy/` directory contains several items requiring user modification:

**All users** need to modify the configuration file:
1. The `config.txt` file contains all the settings required to build the conda environments and run ACT.
It comes pre-populated with example values that you must replace with your own information.
The file must be named config.txt, as all scripts reference this exact filename for configuration.

```bash
# edit the config file with preferred text editor
vim config.txt
```

**Users with raw, multiplexed reads** need to address the following:

2. The `0_raw_reads` directory is where ACT will look for raw reads. 
* Move or copy your raw read files (.fastq.gz format) into this directory. 
* Raw read files should be prefixed with a numerical library ID starting with 1_ (e.g., 1_sample.fastq). The library ID is used to differentiate read files originating from different sequencing runs (i.e., independently prepared ONT libraries). All read files originating from the same sequencing run should have the same library ID.

3. The `barcodes.fasta` file contains a list of barcode sequences with corresponding barcode IDs to be used for read demultiplexing (e.g., sample-specific barcodes added during PCR). ACT uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing) to perform demultiplexing and barcode trimming. The barcodes.fasta file lists the barcodes as pairs of linked adapters, as detailed in the [cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters). The `--revcomp` flag is used by default to check for reads in both orientations. This setting is also compatible with single barcodes.

ACT uses only the barcodes specified by the user in the config.txt file. This allows a research group to maintain a single comprehensive barcode file that can be used in full or as a subset for specific runs. Alternatively, the file may contain only the barcodes being used. For ACT to run properly, barcode IDs must be in the format [zero or more letters][numbers] with the numbers in three-digit format if using less than 100 barcodes (e.g., 001, 002, …) and four-digit format if using 100 or more barcodes (e.g, 0100, 0101, ...), such that there is always a leading 0. For example, `0101`, `BR010`, and `experimentAbarcode003` are all valid IDs but `120`, `BR100`, and `barcode0133experiment2` are not.

```
>barcode001
ATCCGGTCGGAGA...TCTCCGACCGGAT
```
#### _Additional notes on barcodes:_
* If multiple libraries using different barcodes were sequenced on a single flow cell, the `barcodes.fasta` file should include all barcodes from those libraries to account for potential cross-barcode contamination during demultiplexing. For example, if  library 1 contains barcodes 001–036 and library 2 contains barcodes 037–072, and both are sequenced together on a multiplexed flow cell, demultiplexing should use all 72 barcodes. Any reads that demultiplex into barcodes that were not part of the original sample pool (e.g., barcode 015 in library 2) should be discarded.
* ACT was designed for matched dual-barcoded reads (identical barcodes on both ends) and retains any read containing at least one of those barcodes. Combinatorial dual barcodes [different barcodes used on each end] are also supported if the barcode.fasta and config.txt files are adjusted accordingly. To require that the barcode is detected on both ends rather than just one, set REQ_BOTH=TRUE in the config file. If single barcodes are used, the `barcodes.fasta` file should be modified to list individual sequences rather than pairs. 
* For other demultiplexing control options, the `barcodes.fasta` file can be further modified following [cutadapt's adapter search parameters](https://cutadapt.readthedocs.io/en/stable/guide.html#search-parameters).


**Users with reads that are already demultiplexed and trimmed of adapters/barcodes** need to address the following:

2. All reads should be placed in a single folder and the path to that folder should be added into the `config.txt` file under `DEMUX_OUT=`. Reads may be provided in fastq or fastq.gz format and must be named using the ‘libraryID_barcodeID’ format, with library IDs starting at 1 and barcode IDs following the convention [zero or more letters][numbers with a leading zero] (see earlier description).

### 3. Build conda environments
Execution of the `environment_setup.sh` script builds all of the necessary conda environments for running ACT. Most of these environments will be saved in the `envs` directory. You will need to input "Y" several times to confirm the installation. Users may receive a pip warning regarding the `LACA` installation, but this warning can be disregarded. 

In needed, the script can be rerun but will not reconstruct the environments unless executed with the flag `-r | --reconstruct`. The reconstruct flag will delete the existing ACT environments and remake them from scratch. You may be prompted to confirm the deletion or overwriting of environments.

**NOTE:** Activate conda/mamba/micromamba before executing the following commands (replacing 'conda' as needed)

```bash
# Go to work directory
cd $WORK_DIR/
chmod +x ./environment_setup.sh
./environment_setup.sh

conda activate ACT-env

#to reconstruct the environments
./environment_setup.sh -r
```

## PREPARATION OF SEQUENCING READS
### 4. Inspecting reads for optimizing filtering settings.
Filtering parameters should be tailored to your read set. Begin by inspecting your reads with your preferred tool. An optional script for running [NanoPlot](https://github.com/wdecoster/NanoPlot) is included, as that is our recommended tool.

Determine quality score and length cutoffs based on your read distribution and expected amplicon size. During initial pre-processing, loose quality and length cutoffs should retain most reads, and the maximum length should be at least twice the expected amplicon length to recover any reads that may have been concatenated during basecalling. Final filtering should use stricter parameters that align closely with expected amplicon length.

If your library includes amplicons of multiple lengths, such as when pooling 16S amplicons generated with different universal primer sets, you can specify unique filtering parameters for each set in the `config.txt` file. 

```bash
# run the helper script to generate and submit the slurm batch script
cd $WORK_DIR
nanoplot_helper.sh

# after Nanoplot finishes, manually review the read length and quality distributions and use this information to populate the QC and Demux sections of the config.txt file.
```

### 5. Read demultiplexing and filtering
The `QC_demultiplex.sh` script can perform the following three tasks:
1. Pre-processing using loose filtering parameters
2. Demultiplexing with user-defined barcodes
3. Strict filtering for quality and length

Users can choose to perform tasks individually or as a group. Ensure the `QC and Demux Parameters` section of the `config.txt` file is fully populated before submitting the `QC_demultiplex.sh` script. 

**Note regarding demultiplexing parameters:**  
The `BARCODE_OVERLAP` parameter specifies how much of the barcode must be present for the [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) tool to consider it valid. This value should be equal to or slightly less than the full barcode length. As a general rule: for barcodes longer than 10 bp, set BARCODE_OVERLAP to barcode length minus 1 (e.g., BARCODE_OVERLAP=12 for 13-bp barcode); for barcodes of 10 bp or less, use the full barcode length.

**Users with raw, multiplexed reads** should run all steps using the `--all` flag.
**Users with reads that are already demultiplexed and trimmed of adapters/barcodes** should only perform strict filtering using the `--qc2` flag.
```bash
# run the QC_demultiplex.sh script to generate and submit the associated slurm batch scripts.
# The slurm scripts in this portion of ACT use job dependencies so each script must run and complete before the next one starts.
# Options:
# -h, --help      Display this help message
# -1, --qc1       Run QC and filtering pre-demultiplexing
# -d, --demux     Run demultiplexing
# -2, --qc2       Run QC and filtering post-demultiplexing
# -a, --all       Run demultiplexing with pre- and post- filtering

# to run all of the steps
QC_demultiplex.sh --all

# for users with demultiplexed and trimmed reads
QC_demultiplex.sh --qc2
```

We recommend inspecting your reads again after this step to check for potential issues, such as excessive read loss due to filtering parameters or an unexpected number of reads per barcode. One option for this inspection is to run NanoPlot again.

## CLUSTERING AND TAXONOMIC ASSIGNMENT
Within ACT, Emu and Sintax assign taxonomy to each read, while LACA clusters reads into OTUs that are subsequently classified by Sintax. **These steps can run simultaneously for efficiency.**

### 6. Clustering 
 
 Reads are clustered into OTUs with LACA. This step can be time-consuming, depending on the size of the read set. Before proceeding, verify that the path to LACA in the `config.txt` file is correct and matches the location where the `environment_setup.sh` script installed LACA. 
 
 *Note: Some modifications are made to LACA to enable compatibility with ACT. These changes are stored in the 'laca_changes' directory and are automatically applied by the 'laca_setup.sh' script.*

#### 6a. LACA Setup
 LACA requires all reads to be organized within a specific directory structure and using a specific file naming convention. The `laca_setup.sh` script will generate a slurm script to rename the files and organize them into the required structure. If the `-r` flag is used, the script will automatically submit this job; if not, you must manually submit the slurm script before running LACA clustering. 
 
 The setup script also generates the LACA `config.yaml` file in the `5_laca` directory. This file must be manually reviewed to ensure accuracy. Specifically, you will need to verify/edit the following:
 1. Update the primers listed in the quality control section (near line 42) and in the phylogenetic tree section (near line 178). 
 2. Check the Medaka model (near line 144) and confirm it matches the basecalling model used for your reads. A typical Medaka model name looks like: `r1041_e82_400bps_sup_v4.1.0`.
 
```bash
# Use the laca_setup.sh script to up the files and folders for laca and generate the config file
# Omit the '-r' flag to regenerate the config file and slurm scripts without automatically submitting the jobs

laca_setup.sh -r

# ensure the accuracy of the laca config.yaml file (primers, Medaka model) using your preferred text editor

vim 5_laca/config.yaml
```

#### 6b. Running LACA
Once the `config.yaml` file is verified to be accurate, submit the `5_laca_run.sh` script to run LACA.
```
cd slurm_scripts
sbatch 5_laca_run.sh

# Clustering will take a while, so you can run Emu and Sintax on your samples while it is still in progress and then run Sintax on the OTUs once clustering is complete.
```

### 7. Taxonomic assignment

The taxonomic assignments are done with Sintax and Emu. This step can be started while LACA is still running by using flags to submit only specific scripts (see details below). 

ACT uses the Emu- and Sintax-formatted databases specified in the `config.txt` file. By default, ACT will use the provided sequence-similarity-aware ACT-DB 16S database. The FASTA files for ACT-DB are gzipped due to GitHub file size limits but will be automatically unzipped by the `taxonomy_assignment.sh` script. 

If using an updated, modified, or custom database (see [here](#updating-and-modifying-act-db-and-building-custom-ambiguity-aware-databases) for details), ensure it follows the same formatting as the provided files and that the FASTA files are unzipped, since Emu and Sintax cannot process gzipped FASTA files.


**WARNING:** The `-o` and `-a` flags can only be run once LACA clustering is complete, as these flags perform Sintax-based taxonomic assignment of the OTUs. If you want to start taxonomic assignments while LACA is still running, use the example command below. 

```bash
# Use the flags to choose to run the slurm scripts or not
#Options:
# -h, --help      Display this help message
# -a, --all       Run all taxonomic assignments (can only be performed after LACA is finished)
# -e, --emu       Run EMU on samples
# -s, --sintax    Run sintax on samples
# -o, --otu       Run sintax on the LACA OTUs (can only be performed after LACA is finished)

# If running while LACA clustering is ongoing
taxonomy_assignment.sh -e -s

# Then run taxonomic assignment of OTUs once clustering is complete
taxonomy_assignment.sh -o

# If running all taxonomic assignment after LACA clustering is complete
taxonomy_assignment.sh -a
```


### 8. Generate the consensus OTU tables

Once clustering and taxonomic assignment are complete, the final consensus OTU table is generated using the `generate_consensus_table.R` script. For convenience, a shell script `generate_consensus.sh` is provided to simplify this process. If the pipeline has been run up to this point, all required files should already be in the proper locations and the script will handle everything automatically. Samples IDs in the output will use the library ID combined with the numerical portion of the barcode ID.

**Note regarding 'NA' threshold for generating consensus taxonomic assignments:**  
A critical parameter of the ACT decision tree is the Sintax 'NA' threshold. If the Sintax confidence falls below the NA threshold, taxonomy for that rank is withheld to avoid overclassification. The NA threshold is defined by the user in the 'Generate Consensus Parameters' section of the `config.txt` file (`NA_THRESHOLD`), with accepted values between 0 and 1. A higher value requires greater Sintax confidence for classification. If NA_THRESHOLD is not specified, ACT will use a default value of 0.3.

The `generate_consensus.sh` script can be run in an interactive session using the `-r` flag or submitted as a slurm script using the `-s` flag.

```bash
generate_consensus.sh
# Options:
# -h, --help      Display this help message
# -r, --run       Run function in current session
# -s, --slurm     Generate and submit slurm script
```

The **final output** is two abundance tables, named `abundance_table_with_OTU.tsv` and `abundance_table.tsv`, and one taxonomy file. Both abundance tables have sample IDs as the column names and TaxId as the row names, the difference lies in what is included in said TaxId. The abundance table with OTU has the OTU appended to the TaxId where possible (e.g. TaxId_286_OTU_4 would designate _Pseudomonas_ OTU 4) while the plain abundance table does not include OTU. The last file is an updated taxonomy table that matches the TaxId_OTU designations included in `abundance_table_with_OTU.tsv`. 


## TROUBLESHOOTING TIPS

Most ACT scripts will generate and automatically submit slurm batch scripts. All generated scripts are stored in the `slurm_scripts` directory. These slurm scripts are configured to notify users via the email specified in the `config.txt` file when the job finishes or fails. If a job array exits with a mixed signal, you can use the following command to identify which tasks in the array failed.
```bash
# replace $JOBID with the relevant slurm job id
# the first command prints the job and status, the second removes any marked as completed 
sacct -X -j $JOBID -o jobid%20,state%20 | grep -v COMPLETED
```

## UPDATING AND MODIFYING ACT-DB AND BUILDING CUSTOM AMBIGUITY-AWARE DATABASES
ACT-DB integrates 16S sequences from the [NCBI 16S rRNA RefSeq database]((https://www.ncbi.nlm.nih.gov/refseq/)) and the [Ribosomal RNA Database](https://rrndb.umms.med.umich.edu/). Users are encouraged to update ACT-DB, add sequences to ACT-DB, or construct a custom database to use with ACT to fit different experimental systems. When working with known organisms (e.g., laboratory strains), it is recommended that the corresponding sequences be manually added to the database files. Users can also provide custom databases, following the guidelines for creating an [Emu](https://github.com/treangenlab/emu) or [Sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html) database. The ACT pipeline will also automatically create a taxonomy file formatted for [microeco](https://chiliubio.github.io/microeco_tutorial/) as one of the final steps.

## Script for working with ACT-DB and building custom databases
For convenience, a `database_builder.sh` script is provided to make working with ACT-DB and building custom databases as easy as possible. The specific process for updating ACT-DB, adding sequences to ACT-DB, and building a custom database are outlined below. Multiple options can be combined in a single step, for example `database_builder.sh -s -a -p -g` would build a database from  a sintax formatted database, add the provided sequences, perform _in silico_ PCR, and finally group the results. Alternatively, these steps can be performed individually but should follow the order of 'build, add, trim, group', with grouping always performed last.
```bash
database_builder.sh

# Options:
# -h, --help      Display this help message
# -b, --build     Build new database
# -d, --default   Rebuild the default ACT database from latest NCBI 16S RefSeq and rrnDB, grouped by default
# -s, --sintax    Build new database from a sintax/usearch formatted database
# -a, --add       Add user provided sequences to database
# -n, --ncbi      Add sequences to database using list of NCBI accessions
# -p, --pcr       Use AmpliconHunter2 to perform in silico PCR on the database
# -g, --group     Group and rename highly similar sequences in the database
```
### IMPORTANT NOTE ON GROUPING (`-g`, `--group`)
ACT-DB is distinct from standard databases in that it accounts for sequence ambiguity by grouping highly similar 16S sequences into named multi-taxa groups that are assigned in place of individual taxa. Sequence similarity is estimated by global pairwise alignment of all database sequences using Minimap2.

Grouping behavior:
* If all sequences in a group belong to the same species, they are not renamed.
* If sequences in a group span multiple taxa, they are renamed with a new group identifier.
    * Example 1: A group containing two or more *Pseudomonas* species will be renamed at the species rank as *Pseudomonas* group 1.
    * Example 2: A group containing reads from multiple orders will be renamed at the order rank, e.g., Enterobacterales-Moraxellales group 2.

### Removal of problematic reads
Before forming groups, sequences that are identified as highly similar to multiple other species are removed, but only if they represent less than 5% of the sequences for that species. The most common example of this is when a partial 16S sequence matches many other species; if this partial sequence belongs to a species that has at least 20 other sequences in the database, then this partial sequence will be removed. If that species has fewer than 20 representative sequences, then the partial sequence would not be removed, to prevent losing information on intra-species diversity.

After the groups are identified, taxonomic outliers are removed if they make up ≤ 1% of a group at the phylum, class, or order level. For example, in a group of 100 sequences, if 99 belong to the order *Burkholderiales* and the remaining 1 belongs to *Bacillales* then that single *Bacillales* read would be removed as being likely erroneous. This is only applied at the phylum, class, and order level and reduces the number of groups that cross multiple phyla, classes, or orders.

The default similarity threshold used for grouping is 0.005, equating to 99.5% sequence identity. **Users can adjust the similarity threshold** to align with the error rate of their sequencing technology or specific read set by modifying the `SIM_THRESH` field in the `config.txt` file. 

To accompany the grouped database, the grouping function generates the following three lists for user reference: 
* sequences removed (sequences_to_be_removed.csv)
* sequences in each group (seq_id_to_group_list.csv)
* species in each group (group_to_species_list.csv)

## Rebuilding (i.e., updating) ACT-DB
We have provided some ready-to-use databases on OSF (insert doi to OSF once available), including two versions of ACT-DB built in Nov-2024 and Sep-2025. Highly similar sequences were grouped based on a similarity threshold of 99.5%, which represents the higher end of ONT sequencing accuracy achievable when the database was built (between Q20 (99%) and Q30 (99.9%)). 

ACT-DB can be rebuilt to reflect updated RefSeq and rrnDB sequences as follows:

### STEP 1: Modify the `config.txt` file
* Specify the URL for the desired rrnDB version in the `config.txt` file.
* *OPTIONAL:* change the sequence similarity threshold used for grouping using the `SIM_THRESH` parameter in the `config.txt` file (default=0.005, equating to 99.5% identity)
>**NOTE 1:** The rrnDB does not have a static link for downloading the latest version. To obtain the correct URL, visit the [rrnDB download page](https://rrndb.umms.med.umich.edu/downloads/) and copy the link to the latest `rrnDB-*_16S_rRNA.fasta.zip` file.  
**NOTE 2:** If you want to keep the original version of ACT-DB, make a renamed copy of the `/Amplicon_Consensus_Taxonomy/taxonomy_databases/` directory before proceeding further.

### STEP 2: Build the updated ACT-DB database
```bash
# The`-d` or `--default` flags will build the database then group it to output an updated ACT-DB 
database_builder.sh -d
```

## Modifying the existing database
Sequences can be added to an existing database, such as the provided ACT-DB, using either a list of NCBI accessions or user-provided files. If using a grouped database then sequences must be added prior to or simulataneous with grouping, this can be easily done by combining them into a single command e.g. `database_builder.sh -a -g`. The order of the flags after the command does not matter, as the addition of sequences will always be performed before the grouping step.
1. **Adding sequences using NCBI accessions:**
* Provide a text file with a list of RefSeq assembly IDs or nucleotide IDs (one per line) to `ACC_LIST` in the `config.txt` file 
* Run `database_builder.sh -n`   
*This will automatically pull the listed genomes, extract all annotated 16S rRNA sequences as well as the assigned taxonomy, and add everything to the existing database. Users will be notified of any accessions that can't be found.*
* Run `database_builder.sh -g`  
*This will re-group the database to include the newly added sequences.*

>**NOTE 1**: The `database_builder.sh -n` command uses NCBI's E-utilities, which can be finicky about the types of IDs it accepts. If an accession is not found, try using a different identifier. For example, the genome for *Acidovorax facilis* R-28 can be added with either the RefSeq chromosome ID	NZ_CP183985.1 or NCBI RefSeq assembly ID GCF_049913225.1, but not the GenBank assembly ID GCA_049913225.1.   
**NOTE 2:** The user does not control the taxonomy of assemblies added in this manner; they inherit the taxonomic assignment from NCBI; additionally, all identified 16S rRNA genes will be added. If you need more control over sequences and taxonomies, use the `--add` flag as described below.

2. **Modifying the database with user-provided files and the `--add` flag**  
* This option provides the most control over the sequences and taxonomies added, as the user provides all necessary data. However, preparing the data can be time-consuming, which is why we also provide the option of using accession IDs. 

* The `--add` option requires the following 2 files (or possibly 3 files**) to be listed in the `config.txt` file. 

| Config setting | File description |
| :------------- | :--------------- |
| ADD_USER_SEQ       | a fasta file of sequences to be added |
| ADD_USER_SEQ2TAX   | two-column tab-delimited file of sequence headers and TaxIds |
| ADD_USER_TAX**       | tab-delimited taxonomy file for each new TaxId |

>**Taxonomy information in the `ADD_USER_TAX` file is only required if adding TaxIds that do not exist in NCBI (e.g., newly sequenced species, synthetic sequences). This taxonomy file must contain the columns tax_id, domain, phylum, class, order, family, genus, species, t_domain, t_phylum, t_class, t_order, t_family, t_genus, and t_species. The columns starting with "t_" should contain the TaxId of the corresponding phylogeny, and the TaxIds in columns 'tax_id' and 't_species' should match. Ensure that new TaxIds do not already exist. Existing TaxIds can be found on the [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/).

Example for adding a new species (TaxId 8000001) and a synthetic sequence (TaxId 9000005); the new species is assigned existing TaxIds for every rank except species, and the synthetic sequence is given new TaxIds for all ranks except domain. 

| tax_id | domain | phylum | class | order | family | genus | species | t_domain | t_phylum | t_class | t_order | t_family | t_genus | t_species |
| ------ | ------ | ------ | ----- | ----- | ------ | ----- | ------- | -------- | -------- | ------- | ------- | -------- | ------- | --------- |
| 8000001 | Bacteria | Bacteroidota | Flavobacteriia | Flavobacteriales | Weeksellaceae | Chryseobacterium | Chryseobacterium S02 | 2 | 976 | 117743 | 200644 | 2762318 | 59732 | 8000001 |
| 9000005 | Bacteria | Synthetic | Synthetic | Synthetic | Synthetic | Ec5001 | Ec5001 | 2 | 9000000 | 9000001 | 9000002 | 9000003 | 9000004 | 9000005

* Once all files are ready and listed in the `config.txt` file, run the following commands:
```bash
database_builder.sh --add

database_builder.sh -g
```

## Building a custom database
Users can build a **completely custom database** using the `-b` or `--build` flags and providing two files in the `config.txt` file:

| Config setting | File description |
| :------------- | :--------------- |
| BUILD_USER_SEQ       | a fasta file of marker gene sequences with sequence IDs as headers |
| BUILD_USER_SEQ2TAX   | two-column, headerless, tab-delimited list of sequence IDs and their corresponding TaxIds |

The sequence IDs need to match in both files and the TaxIds should correspond to NCBI TaxIds. Sequences that cannot be assigned NCBI TaxIds should be added separately using the `-a | --add` flags after building the initial database. 

Users can also **build an ACT-compatible database from a [sintax/usearch formatted](https://www.drive5.com/usearch/manual/tax_annot.html) database** using the `-s` or `--sintax` flags. Specify the database in the `BUILD_USER_SEQ` field in the `config.txt` file. This command can be used to convert databases like [UNITE](https://unite.ut.ee/repository.php) into an ACT-compatible format.

## Extracting amplicons from a database
The `-p` or `--pcr` flags can be used to perform _in silico_ PCR on an existing database to extract the amplicons that would be produced by your primers. This can be used to make a database that focuses on specific regions and can increase the precision of taxonomic classification by removing noise from the database. It uses [AmpliconHunter2](https://github.com/rhowardstone/AmpliconHunter2) to extract all possible amplicons using user provided primers. Primers must be listed in the `config.txt` as `F_PRIMER` and `R_PRIMER` for forward and reverse and must be in 5' -> 3' orientation. Additional settings for AmpliconHunter2, such as min and max amplicon size and number of mismatches allowed, can be modified if desired but will otherwise run with defaults. If you are building a grouped database with extracted amplicons (database_builder.sh --group --pcr), the amplicons will be extracted first and then grouped by similarity.

| Config setting | Description | 
|:-------------- | :---------- | 
| AMPLICON_MIN | Minimum amplicon length (default: 50) |
| AMPLICON_MAX | Maximum amplicon length (default: 5000) |
| N_MISMATCH | Maximum mismatches (default: 3) |
| CLAMP | 3' clamp size (default: 2) |

## Other ACT-compatible databases
Pre-built ACT-compatible databases for Silva, GreenGenes2, GTDB can be found [here](add link) but are not recommended, for reasons described in the ACT [paper](add biorxiv link). None of the pre-built databases are grouped unless otherwise specified in their [documentation](link to the osf wiki page for the project). Scripts used to generate these databases are stored [here](add relative link to X folder). For fungal amplicons we have provided two versions of the [Eukaryome database](https://doi.org/10.1093/database/baae043) pre-built for ACT. One of these is amplicon extracted version (made using the `-p` default settings), as we found the original database to contain many large (20 kb+) mitochondrion sequences that would not be amplified by our primers (ITS1 TCCGTAGGTGAACCTGC and TW13 GGTCCGTGTTTCAAGACG).
