#!/bin/bash 
#### LACA setup

#user defined variables
source config.txt

echo "Setting up LACA and generating slurm scripts" >&2

run_flag=false


usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help      Display this help message"
 echo " -r, --run       Submit the generated slurm script"
}


while [[ $# -gt 0 ]]; do
  case "$1" in
    -r | --run) 
      run_flag=true
      echo "The slurm script will be submitted for you" >&2
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


#throw error if any variables are missing from the config files
[[ -z "$EMAIL" ]] && { echo "EMAIL is empty" ; exit 1; }
[[ -z "$WORK_DIR" ]] && { echo "WORK_DIR is empty" ; exit 1; }
[[ -z "$LIBRARY" ]] && { echo "LIBRARY is empty" ; exit 1; }
[[ -z "$BARCODES" ]] && { echo "BARCODES is empty" ; exit 1; }
[[ -z "$BARCODE_SEQUENCE" ]] && { echo "BARCODE_SEQUENCE is empty" ; exit 1; }

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    conda activate laca
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate laca
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate laca
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

cd $WORK_DIR


################################################ Prep LACA read_dir ################################################
if [[ -z "$LACA_OUT" ]]; then
    echo "LACA output directory not provided, using default"
    LACA_OUT="$WORK_DIR"/5_laca
fi

if [[ ! -e $LACA_OUT ]]; then
    mkdir $LACA_OUT
fi
cd $LACA_OUT

if [[ ! -e demultiplexed_reads ]]; then
    mkdir demultiplexed_reads
fi


cd $WORK_DIR/slurm_scripts

cat << EOF > 5_laca_file_setup.sh
#!/bin/bash 

#SBATCH --time=1-0:00:00  # max job runtime
#SBATCH --cpus-per-task=16  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "laca_file_setup"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL


cd $LACA_OUT/demultiplexed_reads

WORK_DIR=$WORK_DIR
LACA_OUT=$LACA_OUT
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

cat << 'EOF' >> 5_laca_file_setup.sh



for i in $(seq 1 $LIBRARY);
do
    for b in "${BARCODES[@]}";
        do
        READ=$WORK_DIR/4_chopper_2/"$i"_*0"$b".fastq.gz
        if [[ -e $READ ]] ; then
            #go to folder
            cd $LACA_OUT/demultiplexed_reads
            
            # New sample names
            SAMPLE=$((i * 100 + b))
            
            # new folder for each sample
            mkdir sample"${SAMPLE}"
            
            # copy the filtered reads into the new folder and rename them
            cp $WORK_DIR/4_chopper_2/"$i"_*0"$b".fastq.gz sample"${SAMPLE}"/sample"${SAMPLE}".fastq.gz
            
            # unzip the reads
            cd sample"${SAMPLE}" && gunzip sample"${SAMPLE}".fastq.gz
        fi
        done
done

EOF

if [[ $run_flag == "true" ]]; then
    sbatch 5_laca_file_setup.sh
fi


################################################ Prep laca config file ################################################

cd $LACA_OUT

# creates the config file for laca
laca init --dbdir $LACA_DIR \
    --workdir $LACA_OUT \
    --demuxdir $LACA_OUT/demultiplexed_reads \
    --jobs-max $MAX_JOBS \
    --no-primer-check \
    --ont


################################################ Laca slurm script ################################################

cd $WORK_DIR/slurm_scripts

cat << EOF > 5_laca_run.sh
#!/bin/bash 

#SBATCH --time=$LACA_JOB_TIME-0:00:00  # max job runtime
#SBATCH --cpus-per-task=$SLURM_MAX_CPUS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=$SLURM_MAX_MEM  # max memory
#SBATCH -J "laca_run"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL


cd $LACA_OUT
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

cat << EOF >> 5_laca_run.sh


laca run all -j $SLURM_MAX_CPUS

EOF
