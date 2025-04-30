#!/bin/bash
#### LACA setup

run_flag=0

while getopts 'r' OPTION; do
  case "$OPTION" in
    r) run_flag=1
      printf "slurm script submitted" >&2
      ;;
    ?)
      printf "slurm script not submitted" >&2
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

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"

elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"

else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

cd $WORK_DIR

################################################ Prep LACA read_dir ################################################
if [[ ! -e 5_laca ]]; then
    mkdir 5_laca
fi
cd 5_laca

if [[ ! -e demultiplexed_reads ]]; then
    mkdir demultiplexed_reads
fi


cd $WORK_DIR/slurm_scripts

cat << EOF > laca_file_setup.sh
#!/bin/bash

#SBATCH --time=1-0:00:00  # max job runtime
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "laca_file_setup"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


cd $WORK_DIR/5_laca/demultiplexed_reads

WORK_DIR=$WORK_DIR
LIBRARY=$LIBRARY
BARCODES=(${BARCODE_SEQUENCE[@]})

EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> laca_file_setup.sh
    echo 'source activate laca' >> laca_file_setup.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> laca_file_setup.sh
    echo 'mamba activate laca' >> laca_file_setup.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> laca_file_setup.sh
    echo 'micromamba activate laca' >> laca_file_setup.sh
fi

cat << 'EOF' >> laca_file_setup.sh

for i in $(seq 1 $LIBRARY);
do
    for b in "${BARCODES[@]}";
        do
            #go to folder
            cd $WORK_DIR/5_laca/demultiplexed_reads
            
            # New sample names
            SAMPLE=$((i * 100 + b))
            
            # new folder for each sample
            mkdir sample"${SAMPLE}"
            
            # copy the filtered reads into the new folder and rename them
            cp $WORK_DIR/4_NanoFilt_2/"$i"-"$b".fastq.gz sample"${SAMPLE}"/sample"${SAMPLE}".fastq.gz
            
            # unzip the reads
            cd sample"${SAMPLE}" && gunzip sample"${SAMPLE}".fastq.gz
            
        done
done

EOF

if [[ $run_flag -eq 1 ]]; then
    sbatch laca_file_setup.sh
fi


################################################ Prep laca config file ################################################

cd $WORK_DIR/5_laca/

# creates the config file for laca
laca init --dbdir $LACA_DIR \
    --workdir $WORK_DIR/5_LACA \
    --demuxdir $WORK_DIR/5_LACA/demultiplexed_reads \
    --jobs-max $MAX_JOBS \
    --no-primer-check \
    --ont


################################################ Laca slurm script ################################################

cd $WORK_DIR/slurm_scripts

cat << EOF > laca_run.sh
#!/bin/bash

#SBATCH --time=$LACA_JOB_TIME-0:00:00  # max job runtime
#SBATCH --cpus-per-task=$SLURM_MAX_CPUS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=$SLURM_MAX_MEM  # max memory
#SBATCH -J "laca_run"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


cd $WORK_DIR/5_laca/
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> laca_run.sh
    echo 'source activate laca' >> laca_run.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> laca_run.sh
    echo 'mamba activate laca' >> laca_run.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> laca_run.sh
    echo 'micromamba activate laca' >> laca_run.sh
fi

cat << EOF >> laca_run.sh


laca run all -j $SLURM_MAX_CPUS

EOF

################################################ Fix LACA files ################################################


# Replace the clust.smk and quant.smk files
cp $WORK_DIR/laca_changes/clust.smk $LACA_DIR/laca/workflow/rules/
cp $WORK_DIR/laca_changes/quant.smk $LACA_DIR/laca/workflow/rules/

# In the /laca/laca/workflow/envs directory, modify the yacrd.yaml file to specify minimap2=2.18
# In the /laca/laca/workflow/envs directory, modify the medaka.yaml file to change the order of the channels to make conda-forge first, and specify samtools=1.11
cp $WORK_DIR/laca_changes/medaka.yaml $LACA_DIR/laca/workflow/envs/
cp $WORK_DIR/laca_changes/yacrd.yaml $LACA_DIR/laca/workflow/envs/

