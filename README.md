# Amplicon Consensus Taxonomy (ACT) pipeline

IMPORTANT! This pipeline should be run on one gene at a time. Reads should be all 16S or all ITS (or whatever gene you're interested), not a mix of both. You have to run the pipeline independently for each set.

Make sure you have some version of conda installed. The pipeline can use conda, mamba, or micromamba.

All of these scripts are written to produce and run slurm batch scripts. All scripts will be generated in the `slurm_scripts` folder. Most of them will submit the scripts for you, unless otherwise specified.

## Files needed

First set up the environment by cloning the github repository into your working directory, i.e. `/work/user`.

```bash
git clone https://github.com/Halverson-lab/Amplicon_Consensus_Taxonomy.git
```

From here forward everything should be done in this folder. Move your raw read files into folder `0_raw_reads`. Raw reads files should be numbered at the beginning, starting with "1_". This number will be referred to as the library number.

## Config file

The `config.txt` file has all of the settings for running the pipeline. It comes pre-filled out with example info, you need to delete that and fill it out with your info. Once you're done save the file into the scripts folder.

```bash
# make a copy of the config file
cp example_config.txt config.txt

# edit the config file with preferred text editor
vim config.txt
```


## Set up the environment

Run the environment setup script once. This will build all of your environments. Most of these environments will be saved in the `envs` folder. You will need to input "Y" several times to confirm the installation. Users may receive a pip warning regarding the laca installation, but this warning can be disregarded. This script can be re-run if needed. 

```bash
# Go to work directory
cd $WORK_DIR/
chmod +x ./environment_setup.sh
./environment_setup.sh

conda activate ACT-env
```

## Database builder

User can use the provided database or construct their own. If using the provided database then skip to the next section. 

```bash
database_builder.sh

# Options:
#  -h, --help      Display this help message
#  -b, --build     Build new database
#  -d, --default   Build new database from latest NCBI 16S RefSeq and rrnDB
#  -a, --add       Add user provided sequences to database
#  -n, --ncbi      Add sequences to database using list of NCBI accesions
```

## Examining reads and setting quality control parameters.

The parameters you set for filtering should be based on what your read set looks like. Inspect your reads with your preferred tool. I've included a optional script for running NanoPlot, since that is my preferred tool. 

Inspect your reads and identify quality score and length cut offs based on your read distribution and expected amplicon size. The loose filtering parameters should cover the majority of reads and the max should be at least double your amplicon length, to enable the recovery of reads that may have been accidentally concatenated during basecalling. The stricter parameters are mostly based on your expected amplicon length.

If you have multiple primers in a library and are expecting different sequence lengths then you can identify the parameters for each set and specify them in the config file.

```bash
# run the script to generate the slurm batch script and submit it
cd $WORK_DIR
nanoplot_helper.sh

# after it finishes use it to fill in the QC and Demux portion of the config file, based on your read distributions
```

## Quality control and demultiplexing

Ensure all of the `QC and Demux Parameters` are filled out in the config file. The `BARCODE_OVERLAP` tells `cutadapt` how much of your barcode must be present to be considered valid and should be equal to your barcode length or slightly less. For barcodes greater than 10 bases long the overlap should be length minus 1 (13 bp barcode -> BARCODE_OVERLAP=12). For barcodes shoert than 10 bp it should be the length of the barcode.

### A note on barcodes

This was written to take in a fasta file of all paired numbered barcodes and then choose which you use. This allows the lab to make a single barcode file that anyone can use as a whole or as a subset. If you prefer, you can make a file that only contains the barcodes you use and number them from 1 to X. 

```
>barcode001
ATCCGGTCGGAGA...TCTCCGACCGGAT
```

If you have multiple libraries that have different barcodes but were run on the same flow cell then you should include all possible barcodes when demultiplexing due to the chance of barcode contamination. For example, I pooled barcodes 1-36 into one library and barcodes 37-72 into a second and sent them for sequencing then I would demultiplex with all 72, just to be safe. Any sequences that get binned into barcodes that shouldn't be in that sample should be tossed.


### Running the script

```bash
# run the script to generate the slurm batch scripts and submit it
# The slurm scripts in this job are run with dependencies so they have to run and finish in order
cd $WORK_DIR
QC_demultiplex.sh
# If you need to re-run this because it failed at some step, just delete the failed folders and leave the successful ones,
# it will only re-run the steps with empty folders
```

I recommend you inspect your reads after this step to check for potential issues, such as losing too many reads to filtering parameters or too many/too few reads in a barcode. One option for inspecting them is by running `NanoPlot` again.

## Clustering 
 
 Reads are clustered using [LACA](https://github.com/yanhui09/laca). This step can be time consuming, depending on how large your read set is. This is a good time to ensure the path to `LACA` in the config.txt file is correct and matches where the environment_setup script installed `LACA`. `LACA` requires all of the reads copied into a specific file structure. The `laca_setup.sh` script will generate a slurm script to move all files, and will submit the script if the `-r` flag is used. If the `-r` flag is not used then users will need to manually submit the slurm script for re-organizing the files before `LACA` can be run. Some modifications had to be made to `LACA` in order for it to run, these changes are in the `laca_changes` directory and will be automatically applied by the laca setup script. 
 
 The setup script will generate the `LACA` `config.yaml` file in the `5_laca` directory. This config file needs to be manually checked to ensure its correct. You will need to edit the primers in the quality control section, around line 42, and in the phylogenetic tree section, around line 178. Additionally, you need to check the medaka version model, around line 144, and ensure it matches the basecalling model of your sequences. Once you have checked it, you can manually submit the `laca_run.sh` slurm script. 

 If the clustering fails or times out you can resubmit the `laca_run.sh` slurm and it will resume where it left off.

```bash
# Set up laca (following their tutorial for the config file)
# use your preferred conda
conda activate laca

# Set up the files and folders for laca and generate the config file
# Use -r to run the slurm script, leave off the flag to just regenerate the config file and slurm scripts
cd $WORK_DIR
laca_setup.sh -r

# go in and check the laca config file, especially the medaka version and the primers
# the medaka version should be something like r1041_e82_400bps_sup_v4.1.0
vim $WORK_DIR/5_laca/config.yaml

# once you've checked your config file you can submit the slurm script
cd $WORK_DIR/slurm_scripts
sbatch laca_run.sh

conda deactivate
# clustering will take awhile, you can run EMU and Sintax on your samples while it runs and then run sintax on the OTUs when they're finished
```

## Taxonomic assignments

The taxonomic assignments are done with `sintax` and `EMU`. This step can be started while `LACA` is still running, using the flags to submit only the selected scripts. The databases for 16S genes is provided, but there is a guide on how to regenerate them included. If working with known organisms, such as inoculating with a lab strain, then I recommend manually adding the corresponding sequences to the files. You can also provide custom databases, following the guidelines for creating an [EMU](https://github.com/treangenlab/emu) and [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html) database. The taxonomy databases folder also contains a taxonomy file formatted for [microeco](https://chiliubio.github.io/microeco_tutorial/), to make it easier to analyze data later.

The fasta files for the taxonomy databases has to be gzipped due to github file size limits. These files will be unzipped by the `taxonomy_assignment.sh` script. If you are providing your own files, make sure they follow the same formatting as they provided files and that the fasta files are unzipped, as EMU and sintax can't use gzipped fasta files.


The `-o` and `-a` flags can only be run once the laca clustering is finished, as they perform the taxonomic assignment of the OTUs. If you want to start runing the taxonomic assignments while LACA is running you can run it using the example below. 

```bash
# Use the flags to choose to run the slurm scripts or not
#Options:
# -h, --help      Display this help message
# -a, --all       Run all taxonomic assignments (can only be performed after LACA is finished)
# -e, --emu       Run EMU on samples
# -s, --sintax    Run sintax on samples
# -o, --otu       Run sintax on the LACA OTUs


# If running before the OTUs are done clustering
cd $WORK_DIR
taxonomy_assignment.sh -e -s

# Then run the OTUs once they're finished
taxonomy_assignment.sh -o

# or if running everything at once
taxonomy_assignment.sh -a
```


# Generate the OTU table

## Set up the files you need 

The final OTU table is generated using `generate_consensus_table.R`, but to simplify things there is a shell command `generate_consensus.sh` that simplifies things. If the pipeline has been run to this point then all of the files should be in the proper location and the command will handle everything. You do need to prepare a csv that lists each sample ID (library-barcode 1-1) and whether or not it is gnotobiotic (T or F), in the format shown below. This file path should be saved in the config.txt file.

```
barcode,gnotobiotic
1-1,T
1-2,T
1-3,F
```

If you're using a SynCom then prepare another file that lists the species name with taxid for each member (the names must be in the same format as the Emu database).

```
species
Pseudomonas_putida_303
Agrobacterium_tumefaciens_358
```

## Running the script

Once you have the specified files you can run the `generate_consensus.sh` script. You can run it in an interactive session with the `-r` flag or you can generate and submit a slurm script with the `-s` flag.

```bash
generate_consensus.sh
# Options:
# -h, --help      Display this help message
# -r, --run       Run function in current session
# -s, --slurm     Generate and submit slurm script

```
