#!/bin/bash 
#### Nanoplot helper

#user defined variables
source config.txt

#throw error if var is empty
[[ -z "$EMAIL" ]] && { echo "EMAIL is empty" ; exit 1; }
[[ -z "$WORK_DIR" ]] && { echo "WORK_DIR is empty" ; exit 1; }
[[ -z "$LIBRARY" ]] && { echo "LIBRARY is empty" ; exit 1; }
[[ -z "$BARCODES" ]] && { echo "BARCODES is empty" ; exit 1; }

cd $WORK_DIR


if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    source activate $ENV_DIR/nanoplot-env
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $ENV_DIR/nanoplot-env
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_DIR/nanoplot-env
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

cd $WORK_DIR

if [[ ! -e 1_NanoPlot_raw_reads ]]; then
    mkdir 1_NanoPlot_raw_reads
fi

if [[ ! -e slurm_scripts ]]; then
    mkdir slurm_scripts
fi

cd slurm_scripts

cat << EOF > NanoPlot_slurm.sh
#!/bin/bash 
#SBATCH --time=0-2:00:00  # max job runtime
#SBATCH --cpus-per-task=8  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "NanoPlot"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-$LIBRARY


WORK_DIR=$WORK_DIR
cd $WORK_DIR

EOF


if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> NanoPlot_slurm.sh
    echo "source activate $ENV_DIR/nanoplot-env" >> NanoPlot_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> NanoPlot_slurm.sh
    echo "mamba activate $ENV_DIR/nanoplot-env" >> NanoPlot_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> NanoPlot_slurm.sh
    echo "micromamba activate $ENV_DIR/nanoplot-env" >> NanoPlot_slurm.sh
fi

cat << 'EOF' >> NanoPlot_slurm.sh
#Run nanoplot on each file in 0_raw_reads and output into 1_NanoPlot_raw_reads/
NanoPlot --fastq 0_raw_reads/"${SLURM_ARRAY_TASK_ID}"* -o 1_NanoPlot_raw_reads/"${SLURM_ARRAY_TASK_ID}"_Nanoplot_1 --N50 --huge --threads 8
EOF

sbatch NanoPlot_slurm.sh

if [ $CONDA == "conda" ]; then
    conda deactivate
elif [ $CONDA == "mamba" ]; then
    mamba deactivate
elif [ $CONDA == "micromamba" ]; then
    micromamba deactivate
fi
