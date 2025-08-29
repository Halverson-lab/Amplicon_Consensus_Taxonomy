#!/bin/bash 

## generate consensus abundance table

#user defined variables
source config.txt

#throw error if any variables are missing from the config files
[[ -z "$EMAIL" ]] && { echo "EMAIL is empty" ; exit 1; }
[[ -z "$WORK_DIR" ]] && { echo "WORK_DIR is empty" ; exit 1; }
[[ -z "$LIBRARY" ]] && { echo "LIBRARY is empty" ; exit 1; }
[[ -z "$BARCODES" ]] && { echo "BARCODES is empty" ; exit 1; }

if [[ -z "$NA_THRESHOLD" ]]; then
    echo "NA_THRESHOLD is empty, using default of 0.3" 
    NA_THRESHOLD=0.3
    sed -i 's/NA_THRESHOLD=/NA_THRESHOLD=0.3/g' config.txt
elif (( NA_THRESHOLD < 0 )); then
    echo "NA_THRESHOLD must be greater than 0" 
    exit 1
elif (( NA_THRESHOLD > 1 )); then 
    echo "NA_THRESHOLD must be less than 1" 
    exit 1
fi

if [[ -z "$SINTAX_THRESHOLD" ]]; then
    echo "SINTAX_THRESHOLD is empty, using default of 0.7" 
    SINTAX_THRESHOLD=0.7
    sed -i 's/SINTAX_THRESHOLD=/SINTAX_THRESHOLD=0.7/g' config.txt
elif (( SINTAX_THRESHOLD < 0 )); then
    echo "SINTAX_THRESHOLD must be greater than 0" 
    exit 1
elif (( SINTAX_THRESHOLD > 1 )); then 
    echo "SINTAX_THRESHOLD must be less than 1" 
    exit 1
elif (( NA_THRESHOLD > SINTAX_THRESHOLD )); then
    echo "NA_THRESHOLD cannot be greater than SINTAX_THRESHOLD"
    exit 1;
fi

[[ -z "$(ls -A $EMU_OUT/read_assignments/)" ]] && { echo "${EMU_OUT}/read_assignments/ is empty"; exit 1; }
[[ -z "$(ls -A $SINTAX_OUT)" ]] && { echo "${SINTAX_OUT} is empty"; exit 1; }
[[ ! -e $LACA_OUT/quant/seqID_to_otu.tsv ]] && { echo "${LACA_OUT}/quant/seqID_to_otu.tsv does not exist"; exit 1; }
[[ ! -e $LACA_OUT/sintax_OTUs.tsv ]] && { echo "${LACA_OUT}/OTU_sintax.tsv does not exist"; exit 1; }


run_flag=false
slurm_flag=false


usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help      Display this help message"
 echo " -r, --run       Run function in current session"
 echo " -s, --slurm     Generate and submit slurm script"
}

if [ $# -eq 0 ]; then
  echo "Error: At least one flag is required."
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r | --run) 
      run_flag=true
      echo "Run function in current session" >&2
      ;;
    -s | --slurm)
      slurm_flag=true
      echo "Generating and submitting slurm script" >&2
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



if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    source activate $WORK_DIR/envs/consensus-env
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $WORK_DIR/envs/consensus-env
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $WORK_DIR/envs/consensus-env
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

#generate the list of barcode names for R
seqkit seq -n $BARCODE_FILE > $WORK_DIR/barcode_names.txt

#run the R script
if [[ $run_flag == "true" ]]; then
    cd $WORK_DIR
    generate_consensus_table.R
fi


################################################ slurm script for generating consensus ################################################
### Write the necessary slurm scripts
if [[ $slurm_flag == "true" ]]; then
    cd $WORK_DIR/slurm_scripts
    cat << EOF > generate_consensus_slurm.sh
#!/bin/bash 

#SBATCH --time=0-4:00:00  # max job runtime
#SBATCH --cpus-per-task=$QC_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "generate_consensus_otu_table"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL

WORK_DIR=$WORK_DIR

cd $WORK_DIR/
source config.txt
EOF

    if [ $CONDA == "conda" ]; then
        echo 'eval "$(conda shell hook --shell bash)"' >> generate_consensus_slurm.sh
        echo "source activate $ENV_DIR/consensus-env" >> generate_consensus_slurm.sh
    elif [ $CONDA == "mamba" ]; then
        echo 'eval "$(mamba shell hook --shell bash)"' >> generate_consensus_slurm.sh
        echo "mamba activate $ENV_DIR/consensus-env" >> generate_consensus_slurm.sh
    elif [ $CONDA == "micromamba" ]; then
        echo 'eval "$(micromamba shell hook --shell bash)"' >> generate_consensus_slurm.sh
        echo "micromamba activate $ENV_DIR/consensus-env" >> generate_consensus_slurm.sh
    fi

    echo 'generate_consensus_table.R' >> generate_consensus_slurm.sh
    
    sbatch generate_consensus_slurm.sh

fi



