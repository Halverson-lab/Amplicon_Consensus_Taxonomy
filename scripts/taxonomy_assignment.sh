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
    -h | --help)
      usage
      exit 0
      ;;
    \?)
      echo "Invalid option" >&2
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
if [[ ! -e 6_emu ]]; then
    mkdir 6_emu
fi

if [[ ! -e 6_emu/read_assignments ]]; then
    mkdir 6_emu/read_assignments
fi

if [[ ! -e 7_sintax ]]; then
    mkdir 7_sintax
fi

#### if using the provided databases and they are still zipped then unzip them 
if [[ -e sintax_db.fasta.gz ]]; then
    gunzip taxonomy_databases/*.gz
fi



cd slurm_scripts
### Write the necessary slurm scripts

JOB_TIME=$(($LIBRARY * 3))
ARRAY_SEQUENCE=$(IFS=,; echo "${BARCODE_SEQUENCE[*]}")

################################################ EMU ################################################
cat << EOF > EMU_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$EMU_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "emu"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=$ARRAY_SEQUENCE


cd $WORK_DIR/6_emu
ENV_DIR=$ENV_DIR
WORK_DIR=$WORK_DIR/6_emu
READ_DIR=$WORK_DIR/4_NanoFilt_2
EMU_DB=$EMU_DB
THREADS=$EMU_THREADS
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> EMU_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> EMU_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> EMU_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> EMU_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> EMU_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> EMU_slurm.sh
fi

cat << 'EOF' >> EMU_slurm.sh
for read in  $READ_DIR/*0"${SLURM_ARRAY_TASK_ID}".fastq.gz; do
    if [[ -e $read ]] ; then
        emu abundance --type lr:hq \
            --keep-counts \
            --keep-read-assignments \
            --output-unclassified \
            --threads ${THREADS} \
            $read \
            --db $EMU_DB \
            --output-dir $(basename "$read" .fastq.gz) \
            --output-basename $(basename "$read" .fastq.gz) 
            
        cp "$(basename "$read" .fastq.gz)"/*_read-assignment-distributions.tsv $WORK_DIR/read_assignments/
    fi
done
EOF


################################################ Sintax samples ################################################
cat << EOF > Sintax_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$SINTAX_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "sintax"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=$ARRAY_SEQUENCE


cd $WORK_DIR/7_sintax
ENV_DIR=$ENV_DIR
WORK_DIR=$WORK_DIR/7_sintax
READ_DIR=$WORK_DIR/4_NanoFilt_2
SINTAX_DB=$SINTAX_DB
SINTAX_CUTOFF=$SINTAX_CUTOFF

EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> Sintax_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> Sintax_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> Sintax_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> Sintax_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> Sintax_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> Sintax_slurm.sh
fi


cat << 'EOF' >> Sintax_slurm.sh
for read in  $READ_DIR/*0"${SLURM_ARRAY_TASK_ID}".fastq.gz; do
    #get number of reads
    READ_COUNT=$(seqkit stats $read -T | csvtk -t cut -f 4 | csvtk del-header)
    if [[ ! $READ_COUNT -eq 0 ]]; then
        vsearch --sintax \
            $read \
            --db $SINTAX_DB \
            --tabbedout $WORK_DIR/"$(basename "$read" .fastq.gz)"_sintax.tsv \
            --sintax_cutoff 0.5 \
            --strand both \
            -notrunclabels \
            --sintax_random
    fi
done
EOF




################################################ Sintax OTUs ################################################

cat << EOF > Sintax_OTU_slurm.sh
#!/bin/bash 

#SBATCH --time=0-$JOB_TIME:00:00  # max job runtime
#SBATCH --cpus-per-task=$SINTAX_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "sintax"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


cd $WORK_DIR/5_laca
ENV_DIR=$ENV_DIR
READ_DIR=$WORK_DIR/5_laca
SINTAX_DB=$SINTAX_DB
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo "source activate $ENV_DIR/taxonomy-env" >> Sintax_OTU_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo "mamba activate $ENV_DIR/taxonomy-env" >> Sintax_OTU_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo "micromamba activate $ENV_DIR/taxonomy-env" >> Sintax_OTU_slurm.sh
fi


cat << 'EOF' >> Sintax_OTU_slurm.sh
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
    sbatch EMU_slurm.sh
fi

if [[ $sintax_flag == "true" ]]; then
    sbatch Sintax_slurm.sh
fi

if [[ $otu_flag == "true" ]]; then
    if [[ -f $WORK_DIR/5_laca/rep_seqs.fasta ]]; then
        sbatch Sintax_OTU_slurm.sh
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
