#!/bin/bash 
#### Slurm scripts for EMU and Sintax

emu_flag=false
otu_flag=false
sintax_flag=false

usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help      Display this help message"
 echo " -a, --all       Run all taxonomic assignments (can only be performed after LACA is finished)"
 echo " -e, --emu       Run EMU on samples"
 echo " -s, --sintax    Run sintax on samples"
 echo " -o, --otu       Run sintax on the LACA OTUs"
 echo " -t, --test      Generate the slurm scripts but do not submit them"
}

if [ $# -eq 0 ]; then
  echo "Error: At least one flag is required."
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a | --all) 
      emu_flag=true
      sintax_flag=true
      otu_flag=true
      echo "Running all taxonomic assignments" >&2
      ;;
    -e | --emu)
      emu_flag=true
      echo "Running EMU on samples" >&2
      ;;
    -s | --sintax) 
      sintax_flag=true
      echo "Running sintax on samples" >&2
      ;;
    -o | --otu)
      otu_flag=true
      echo "Running sintax on LACA OTUs" >&2
      ;;
    -t | --test)
      emu_flag=false
      sintax_flag=false
      otu_flag=false
      echo "Generating slurm scripts" >&2
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "$1 is not recognized" >&2
      usage
      exit 1
      ;;
  esac
  shift
done


#user defined variables
source config.txt

#throw error if any variables are missing from the config files
[[ -z "$EMAIL" ]] && { echo "EMAIL is empty" ; exit 1; }
[[ -z "$WORK_DIR" ]] && { echo "WORK_DIR is empty" ; exit 1; }
[[ -z "$LIBRARY" ]] && { echo "LIBRARY is empty" ; exit 1; }
[[ -z "$BARCODES" ]] && { echo "BARCODES is empty" ; exit 1; }
[[ -z "$BARCODE_SEQUENCE" ]] && { echo "BARCODE_SEQUENCE is empty" ; exit 1; }

[[ -z "$SINTAX_DB" ]] && { echo "SINTAX_DB is empty" ; exit 1; }
[[ -z "$SINTAX_THREADS" ]] && { echo "SINTAX_THREADS is empty" ; exit 1; }
[[ -z "$EMU_DB" ]] && { echo "EMU_DB is empty" ; exit 1; }
[[ -z "$EMU_THREADS" ]] && { echo "EMU_THREADS is empty" ; exit 1; }


if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    source activate $ENV_DIR/taxonomy-env
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $ENV_DIR/taxonomy-env
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_DIR/taxonomy-env
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

cd $WORK_DIR

#### make folders for sintax and emu
[[ -z "$POST_DEMUX_QC_OUT" ]] && { POST_DEMUX_QC_OUT="$WORK_DIR"/4_chopper_2 ; }

[[ -z "$LACA_OUT" ]] && { LACA_OUT="$WORK_DIR"/5_laca ; }
[[ -z "$EMU_OUT" ]] && { EMU_OUT="$WORK_DIR"/6_emu ; echo "emu output directory not provided, using default 6_emu" ; }
[[ -z "$SINTAX_OUT" ]] && { SINTAX_OUT="$WORK_DIR"/7_sintax ; echo "sintax output directory not provided, using default 7_sintax" ; }
[[ ! -e $EMU_OUT ]] && { mkdir $EMU_OUT ; }
[[ ! -e $EMU_OUT/read_assignments ]] && { mkdir $EMU_OUT/read_assignments ; }
[[ ! -e $EMU_OUT/minimap2_aln_stats ]] && { mkdir $EMU_OUT/minimap2_aln_stats ; }
[[ ! -e $SINTAX_OUT ]] && { mkdir $SINTAX_OUT ; }

#### if using the provided databases and they are still zipped then unzip them 
[[ -e $DATABASE_DIR/sintax_db.fasta.gz ]] && { cd $DATABASE_DIR ; gunzip *.gz ; }

cd $WORK_DIR/slurm_scripts
### Write the necessary slurm scripts
#if job time is empty then default is 3 hours per library libraries
[[ -z "$EMU_JOB_TIME" ]] && { EMU_JOB_TIME=$(($LIBRARY * 3)) ; }
ARRAY_SEQUENCE=$(IFS=,; echo "${BARCODE_SEQUENCE[*]}")

################################################ EMU ################################################
cat << EOF > 6_emu_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$EMU_JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$EMU_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "emu"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL
#SBATCH --array=$ARRAY_SEQUENCE


cd $EMU_OUT
ENV_DIR=$ENV_DIR
EMU_OUT=$EMU_OUT
READ_DIR=$POST_DEMUX_QC_OUT
EMU_DB=$EMU_DB
THREADS=$EMU_THREADS
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> 6_emu_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> 6_emu_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> 6_emu_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> 6_emu_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> 6_emu_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> 6_emu_slurm.sh
fi

cat << 'EOF' >> 6_emu_slurm.sh
for read in  $READ_DIR/*0"${SLURM_ARRAY_TASK_ID}".fast*; do # for each read in the read directory
    if [[ -e $read ]] ; then #if the read exists
        # strip the excess info from the sequence headers, by saving it as a temporary file 
        seqkit seq -i $read -o tmp_"$(basename "${read%%.*}")".fastq.gz 
        read2=tmp_"$(basename "${read%%.*}")".fastq.gz

        #run emu with the stripped file
        emu abundance --type lr:hq \
            --keep-counts \
            --keep-files \
            --keep-read-assignments \
            --output-unclassified \
            --threads ${THREADS} \
            $read2 \
            --db $EMU_DB \
            --output-dir $(basename "${read%%.*}") \
            --output-basename $(basename "${read%%.*}") 
                    
        # save minimap alignment stats to csv file for ACT
        minimap_to_csv.py "$(basename "${read%%.*}")"/*.sam "$(basename "${read%%.*}")"/"$(basename "${read%%.*}")"_aln_stats.csv

        # copy the important files to the relevant folders for later
        cp "$(basename "${read%%.*}")"/*_read-assignment-distributions.tsv $EMU_OUT/read_assignments
        cp "$(basename "${read%%.*}")"/*_aln_stats.csv $EMU_OUT/minimap2_aln_stats

        # remove the temporary read file
        rm $read2
    fi
done
EOF


################################################ Sintax samples ################################################
[[ -z "$SINTAX_JOB_TIME" ]] && { SINTAX_JOB_TIME=$(($LIBRARY * 3)) ; }
cat << EOF > 7_sintax_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$SINTAX_JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$SINTAX_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "sintax"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL
#SBATCH --array=$ARRAY_SEQUENCE


cd $SINTAX_OUT
ENV_DIR=$ENV_DIR
SINTAX_OUT=$SINTAX_OUT
READ_DIR=$POST_DEMUX_QC_OUT
SINTAX_DB=$SINTAX_DB

EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> 7_sintax_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> 7_sintax_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> 7_sintax_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> 7_sintax_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> 7_sintax_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> 7_sintax_slurm.sh
fi


cat << 'EOF' >> 7_sintax_slurm.sh
for read in  $READ_DIR/*0"${SLURM_ARRAY_TASK_ID}".fast*; do # for each read in the read directory
    #get number of reads
    READ_COUNT=$(seqkit stats $read -T | csvtk -t cut -f 4 | csvtk del-header)
    if [[ ! $READ_COUNT -eq 0 ]]; then # if there are >0 reads
        # strip the excess info from the sequence headers, by saving it as a temporary file 
        seqkit seq -i $read -o tmp_"$(basename "${read%%.*}")".fastq.gz
        read2=tmp_"$(basename "${read%%.*}")".fastq.gz

        # run sintax with the temp file
        vsearch --sintax \
            $read2 \
            --db $SINTAX_DB \
            --tabbedout $SINTAX_OUT/"$(basename "${read%%.*}")"_sintax.tsv \
            --sintax_cutoff 0.5 \
            --strand both \
            -notrunclabels \
            --sintax_random

        # remove the temporary read file
        rm $read2
    fi
done
EOF




################################################ Sintax OTUs ################################################

cat << EOF > 7_sintax_OTU_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$SINTAX_JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$SINTAX_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "sintax"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL


cd $LACA_OUT
ENV_DIR=$ENV_DIR
READ_DIR=$LACA_OUT
SINTAX_DB=$SINTAX_DB
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> 7_sintax_OTU_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> 7_sintax_OTU_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> 7_sintax_OTU_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> 7_sintax_OTU_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> 7_sintax_OTU_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> 7_sintax_OTU_slurm.sh
fi


cat << 'EOF' >> 7_sintax_OTU_slurm.sh
vsearch --sintax \
    $READ_DIR/rep_seqs.fasta \
    --db $SINTAX_DB \
    --tabbedout $READ_DIR/OTU_sintax.tsv \
    --sintax_cutoff 0.7 \
    --strand both \
    -notrunclabels \
    --sintax_random
EOF



################################################ submit slurm scripts ################################################

if [[ $emu_flag == "true" ]]; then
    sbatch 6_emu_slurm.sh
fi

if [[ $sintax_flag == "true" ]]; then
    sbatch 7_sintax_slurm.sh
fi

if [[ $otu_flag == "true" ]]; then
    if [[ -f $LACA_OUT/rep_seqs.fasta ]]; then
        sbatch 7_sintax_OTU_slurm.sh
    else
        echo "Cannot classify OTUs because LACA rep_seqs.fasta does not exist." >&2
        echo "Please make sure LACA has finished running before using -a or -o." >&2
        exit 1
    fi
fi

if [ $CONDA == "conda" ]; then
    conda deactivate
elif [ $CONDA == "mamba" ]; then
    mamba deactivate
elif [ $CONDA == "micromamba" ]; then
    micromamba deactivate
fi
