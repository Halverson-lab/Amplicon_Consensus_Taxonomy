# Amplicon Consensus Taxonomy (ACT) pipeline

IMPORTANT! This pipeline should be run on one gene at a time. Reads should be all 16S or all ITS (or whatever gene you're interested), not a mix of both. You have to run the pipeline independently for each set.

Make sure you have some version of conda installed. The pipeline can use conda, mamba, or micromamba.

All of these scripts are written to produce and run slurm batch scripts. All scripts will be generated in the `slurm_scripts` folder. Most of them will submit the scripts for you, unless otherwise specified.

## Files needed

First set up the environment by cloning the github repository into your working directory, i.e. `/work/larryh/user`.

```shell
git clone https://github.com/ashleyp1/Amplicon_Consensus_Taxonomy.git
```

From here forward everything should be done in this folder. Move your raw read files into folder `0_raw_reads`. Raw reads files should be numbered at the beginning, starting with "1_". This number will be referred to as the library number.

## Config file

The `config.txt` file has all of the settings for running the pipeline. It comes pre-filled out with example info, you need to delete that and fill it out with your info. Do not move or rename this file, simply open and edit it.


## Set up the environment

Run the environment setup script once. This will build all of your environments. Most of these environments will be saved in the `envs` folder.

```shell
# Go to work directory
cd $WORK_DIR
./environment_setup.sh
```

## Examining reads and setting quality control parameters.

The parameters you set for filtering should be based on what your read set looks like. Inspect your reads with your preferred tool. I've included a optional script for running NanoPlot, since that is my preferred tool. 

Inspect your reads and identify quality score and length cut offs based on your read distribution and expected amplicon size. The loose filtering parameters should cover the majority of reads and the max should be at least double your amplicon length, to enable the recovery of reads that may have been accidentally concatenated during basecalling. The stricter parameters are mostly based on your expected amplicon length.

If you have multiple primers in a library and are expecting different sequence lengths then you can identify the parameters for each set and specify them in the config file.

```shell
# run the script to generate the slurm batch script and submit it
./nanoplot_helper.sh

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

```shell
# run the script to generate the slurm batch scripts and submit it
# The slurm scripts in this job are run with dependencies so they have to run and finish in order 
./QC_demultiplex.sh
# If you need to re-run this because it failed at some step, just delete the failed folders and leave the successful ones,
# it will only re-run the steps with empty folders
```

I recommend you inspect your reads after this step to check for potential issues, such as losing too many reads to filtering parameters or too many/too few reads in a barcode.

## Clustering 
 
 Reads are clustered using [LACA](https://github.com/yanhui09/laca). This step can be time consuming, depending on how large your read set is. `LACA` requires all of the reads copied into a specific file structure. The `laca_setup.sh` script will generate a slurm script to move all files, and will submit the script if the `-r` flag is used. The setup script will also generate the `LACA` config file. This config file needs to be manually checked to ensure its correct. Once you have checked it, you can manually submit the `laca_run.sh` slurm script. 

 If the clustering fails or times out you can resubmit the `laca_run.sh` slurm and it will resume where it left off.

```shell
# Set up laca (following their tutorial for the config file)
# use your preferred conda
conda activate laca

# Set up the files and folders for laca and generate the config file
# Use -r to run the slurm script, leave off the flag to just regenerate the config file and slurm scripts
./laca_setup.sh -r

# go in and check the laca config file, especially the medaka version and the primers
vim $WORK_DIR/5_laca/config.yaml

# once you've checked your config file you can submit the slurm script
cd slurm_scripts
sbatch laca_run.sh

conda deactivate
# clustering will take awhile, you can run EMU and Sintax on your samples while it runs and then run sintax on the OTUs when they're finished
```

## Taxonomic assignments

The taxonomic assignments are done with `sintax` and `EMU`. This step can be started while `LACA` is still running, using the flags to submit only the selected scripts. The databases for 16S genes is provided, but there is a guide on how to regenerate them included. If working with known organisms, such as inoculating with a lab strain, then I recommend manually adding the corresponding sequences to the files. You can also provide custom databases, following the guidelines for creating an [EMU](https://github.com/treangenlab/emu) and [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html) database.

The `-o` and `-a` flags can only be run once the laca clustering is finished, as they perform the taxonomic assignment of the OTUs. If you want to start runing the taxonomic assignments while LACA is running you can run it using the example below. 

```shell
# Use the flags to choose to run the slurm scripts or not
# script usage: 
#   -e to run EMU
#   -s to run sintax on samples
#   -o to run sintax on OTUs
#   -a to run all of the above
#   -h to print options

# If running before the OTUs are done clustering
./taxonomy_assignment.sh -e -s

# Then run the OTUs once they're finished
./taxonomy_assignment.sh -o

# or if running everything at once
./taxonomy_assignment.sh -a
```


# Generate the OTU table

## Set up the files you need 

The final OTU table is generated using R. To make things easier I recommend saving all the necessary files in one folder, otherwise you can hardcode the locations in the R file.

`seqID_to_otu.tsv` and `sintax_OTUs.tsv` from `5_LACA`

`taxonomy.tsv` from the EMU database you used for the taxonomy assignments

Make a folder of all of the EMU read assignments (*_read-assignment-distributions.tsv) and a folder of all the sintax outputs.

Prepare a csv that lists each sample ID (library-barcode 1-HL01) and whether or not it is gnotobiotic (T or F), in the format shown below.

```
barcode,gnotobiotic
1-HL01,T
1-HL02,T
```

if you're using a SynCom then prepare another file that lists the species name with taxid for each member (the names must be in the same format as the Emu database).

```
species
Pseudomonas_putida_303
Agrobacterium_tumefaciens_358
```
