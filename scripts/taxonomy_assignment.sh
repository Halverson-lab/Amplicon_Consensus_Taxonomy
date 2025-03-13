#!/bin/bash
#### Slurm scripts for EMU and Sintax

all_flag=0
emu_flag=0
otu_flag=0
sintax_flag=0

while getopts 'aeosh' OPTION; do
  case "$OPTION" in
    a) all_flag=1
      ;;
    e) emu_flag=1
      ;;
    o) otu_flag=1
      ;;
    s) sintax_flag=1
      ;;
    h)
      printf "script usage: -e to run EMU
      -s to run sintax on samples
      -o to run sintax on OTUs
      -a to run all of the above
      " >&2
      exit 1
      ;;
  esac
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

if [[ ! -e 7_sintax ]]; then
    mkdir 7_sintax
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
ENV_DIR=$WORK_DIR/envs
WORK_DIR=$WORK_DIR/6_emu
READ_DIR=$WORK_DIR/4_NanoFilt_2
EMU_DB=$EMU_DB
THREADS=$EMU_THREADS
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> EMU_slurm.sh
    echo 'source activate $ENV_DIR/taxonomy-env' >> EMU_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> EMU_slurm.sh
    echo 'mamba activate $ENV_DIR/taxonomy-env' >> EMU_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> EMU_slurm.sh
    echo 'micromamba activate $ENV_DIR/taxonomy-env' >> EMU_slurm.sh
fi

cat << EOF >> EMU_slurm.sh
for n in {1..$LIBRARY};
EOF

cat << 'EOF' >> EMU_slurm.sh
do
    emu abundance --type lr:hq \
        --keep-counts \
        --keep-read-assignments \
        --output-unclassified \
        --threads ${THREADS} \
        $READ_DIR/"${n}"-"${SLURM_ARRAY_TASK_ID}".fastq.gz \
        --db $EMU_DB \
        --output-dir "${n}"-"${SLURM_ARRAY_TASK_ID}" \
        --output-basename "${n}"-"${SLURM_ARRAY_TASK_ID}" \
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
ENV_DIR=$WORK_DIR/envs
WORK_DIR=$WORK_DIR/7_sintax
READ_DIR=$WORK_DIR/4_NanoFilt_2
SINTAX_DB=$SINTAX_DB
SINTAX_CUTOFF=$SINTAX_CUTOFF

EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> Sintax_slurm.sh
    echo 'source activate $ENV_DIR/taxonomy-env' >> Sintax_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> Sintax_slurm.sh
    echo 'mamba activate $ENV_DIR/taxonomy-env' >> Sintax_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> Sintax_slurm.sh
    echo 'micromamba activate $ENV_DIR/taxonomy-env' >> Sintax_slurm.sh
fi

cat << EOF >> Sintax_slurm.sh
for n in {1..$LIBRARY};
EOF

cat << 'EOF' >> Sintax_slurm.sh
do
    vsearch --sintax \
        $READ_DIR/"${n}"-"${SLURM_ARRAY_TASK_ID}".fastq.gz \
        --db $SINTAX_DB \
        --tabbedout $WORK_DIR/"${n}"-"${SLURM_ARRAY_TASK_ID}"_sintax.tsv \
        --sintax_cutoff 0.5 \
        --strand both \
        -notrunclabels \
        --sintax_random
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


cd $WORK_DIR/8_OTU_sintax
ENV_DIR=$WORK_DIR/envs
READ_DIR=$WORK_DIR/5_laca
SINTAX_DB=$SINTAX_DB
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo 'source activate $ENV_DIR/taxonomy-env' >> Sintax_OTU_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo 'mamba activate $ENV_DIR/taxonomy-env' >> Sintax_OTU_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> Sintax_OTU_slurm.sh
    echo 'micromamba activate $ENV_DIR/taxonomy-env' >> Sintax_OTU_slurm.sh
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

if [[ $all_flag -eq 1 ]]; then
    emu_flag=1
    sintax_flag=1
    otu_flag=1
fi

if [[ $emu_flag -eq 1 ]]; then
    sbatch EMU_slurm.sh
fi

if [[ $sintax_flag -eq 1 ]]; then
    sbatch Sintax_slurm.sh
fi

if [[ $otu_flag -eq 1 ]]; then
    if [[ -f $WORK_DIR/5_laca/rep_seqs.fasta ]]; then
        sbatch Sintax_OTU_slurm.sh
    else
        echo "Cannot classify OTUs because LACA rep_seqs.fasta does not exist."
        echo "Please make sure LACA has finished running before using -a or -o."
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