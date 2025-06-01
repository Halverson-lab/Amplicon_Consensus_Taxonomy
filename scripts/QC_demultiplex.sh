#!/bin/bash
#### filter by length, demultiplex with cutadapt, filter again

#user defined variables
source config.txt

#throw error if any variables are missing from the config files
[[ -z "$EMAIL" ]] && { echo "EMAIL is empty" ; exit 1; }
[[ -z "$WORK_DIR" ]] && { echo "WORK_DIR is empty" ; exit 1; }
[[ -z "$LIBRARY" ]] && { echo "LIBRARY is empty" ; exit 1; }
[[ -z "$BARCODES" ]] && { echo "BARCODES is empty" ; exit 1; }
[[ -z "$BARCODE_SEQUENCE" ]] && { echo "BARCODE_SEQUENCE is empty" ; exit 1; }

[[ -z "$LOOSE_MIN_LENGTH" ]] && { echo "LOOSE_MIN_LENGTH is empty" ; exit 1; }
[[ -z "$LOOSE_MAX_LENGTH" ]] && { echo "LOOSE_MAX_LENGTH is empty" ; exit 1; }
[[ -z "$LOOSE_MIN_Q" ]] && { echo "LOOSE_MIN_Q is empty" ; exit 1; }


if [ $MULTI_LENGTH == "FALSE" ]; then
    [[ -z "$STRICT_MIN_LENGTH" ]] && { echo "STRICT_MIN_LENGTH is empty" ; exit 1; }
    [[ -z "$STRICT_MAX_LENGTH" ]] && { echo "STRICT_MAX_LENGTH is empty" ; exit 1; }
    [[ -z "$STRICT_MIN_Q" ]] && { echo "STRICT_MIN_Q is empty" ; exit 1; }
elif [ $MULTI_LENGTH == "TRUE" ]; then
    [[ -z "$NUM_OF_AMPLICON_SETS" ]] && { echo "NUM_OF_AMPLICON_SETS is empty" ; exit 1; }
    [[ -z "$LIB_SET_1" ]] && { echo "Define LIB_SET for 1 or more amplicon sets" ; exit 1; }
    [[ -z "$BARCODE_SET_1" ]] && { echo "Define BARCODE_SET for 1 or more amplicon sets" ; exit 1; }
    [[ -z "$MIN_1" ]] && { echo "Set MIN length for 1 or more amplicon sets" ; exit 1; }
    [[ -z "$MAX_1" ]] && { echo "Set MAX length for 1 or more amplicon sets" ; exit 1; }
    [[ -z "$MINQ_1" ]] && { echo "Set MINQ for 1 or more amplicon sets" ; exit 1; }
else
    echo "MULTI_LENGTH must be set as TRUE or FALSE"
    exit 1
fi


# Set up environments
cd $WORK_DIR/envs

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    source activate $WORK_DIR/envs/cutadapt-env
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $WORK_DIR/envs/cutadapt-env
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $WORK_DIR/envs/cutadapt-env
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

#Make the relevant directories if they don't already exist

cd $WORK_DIR

if [[ ! -e 2_NanoFilt_1 ]]; then
    mkdir 2_NanoFilt_1
fi
if [[ ! -e 3_Demultiplex ]]; then
    mkdir 3_Demultiplex
fi
if [[ ! -e 4_NanoFilt_2 ]]; then
    mkdir 4_NanoFilt_2
fi
if [[ ! -e slurm_scripts ]]; then
    mkdir slurm_scripts
fi


### Write the necessary slurm scripts
cd $WORK_DIR/slurm_scripts


################################################ First round of filtering, pre-demux ################################################
cat << EOF > NanoFilt_1_slurm.sh
#!/bin/bash

#SBATCH --time=0-4:00:00  # max job runtime
#SBATCH --cpus-per-task=$QC_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "2_NanoFilt_1"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=1-$LIBRARY

WORK_DIR=$WORK_DIR

cd $WORK_DIR/2_NanoFilt_1
EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> NanoFilt_1_slurm.sh
    echo 'source activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_1_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> NanoFilt_1_slurm.sh
    echo 'mamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_1_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> NanoFilt_1_slurm.sh
    echo 'micromamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_1_slurm.sh
fi

cat << EOF >> NanoFilt_1_slurm.sh

LOOSE_MIN_LENGTH=$LOOSE_MIN_LENGTH
LOOSE_MAX_LENGTH=$LOOSE_MAX_LENGTH
LOOSE_MIN_Q=$LOOSE_MIN_Q

READ_DIR=$WORK_DIR/0_raw_reads

EOF

cat << 'EOF' >> NanoFilt_1_slurm.sh


if gunzip --test $READ_DIR/"${SLURM_ARRAY_TASK_ID}"_* 2>/dev/null 1>/dev/null; then
  gunzip -c $READ_DIR/"${SLURM_ARRAY_TASK_ID}"_* \
    | NanoFilt --length $LOOSE_MIN_LENGTH --maxlength $LOOSE_MAX_LENGTH -q $LOOSE_MIN_Q \
    | gzip > "${SLURM_ARRAY_TASK_ID}"_filt.fastq.gz
else
  gzip $READ_DIR/"${SLURM_ARRAY_TASK_ID}"_*
  gunzip -c $READ_DIR/"${SLURM_ARRAY_TASK_ID}"_* \
    | NanoFilt --length $LOOSE_MIN_LENGTH --maxlength $LOOSE_MAX_LENGTH -q $LOOSE_MIN_Q \
    | gzip > "${SLURM_ARRAY_TASK_ID}"_filt.fastq.gz
fi

EOF



################################################ Demultiplex with cutadapt ################################################
cat << EOF > Demux_slurm.sh
#!/bin/bash

#SBATCH --time=0-${DEMUX_JOB_TIME}:00:00  # max job runtime
#SBATCH --cpus-per-task=$QC_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "3_Demux"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=1-$LIBRARY

cd $WORK_DIR/3_Demultiplex

WORK_DIR=$WORK_DIR
READ_DIR=$WORK_DIR/2_NanoFilt_1
BARCODE_FILE=$BARCODE_FILE
BARCODE_OVERLAP=$BARCODE_OVERLAP

EOF

if [ $CONDA == "conda" ]; then
    echo 'eval "$(conda shell hook --shell bash)"' >> Demux_slurm.sh
    echo 'source activate $WORK_DIR/envs/cutadapt-env' >> Demux_slurm.sh
elif [ $CONDA == "mamba" ]; then
    echo 'eval "$(mamba shell hook --shell bash)"' >> Demux_slurm.sh
    echo 'mamba activate $WORK_DIR/envs/cutadapt-env' >> Demux_slurm.sh
elif [ $CONDA == "micromamba" ]; then
    echo 'eval "$(micromamba shell hook --shell bash)"' >> Demux_slurm.sh
    echo 'micromamba activate $WORK_DIR/envs/cutadapt-env' >> Demux_slurm.sh
fi

cat << 'EOF' >> Demux_slurm.sh
cutadapt --revcomp --overlap $BARCODE_OVERLAP -j 16 -e $BARCODE_ERROR -a file:"${BARCODE_FILE}" -o "${SLURM_ARRAY_TASK_ID}"_demux-{name}.fastq.gz $READ_DIR/"${SLURM_ARRAY_TASK_ID}"_filt.fastq.gz
EOF


################################################ Second round of filtering, post-demux ################################################

JOB_TIME=$(($LIBRARY * 1))

# Write one slurm script with the strict variables if there is only one set of primer lengths
if [ $MULTI_LENGTH == "FALSE" ]; then
    MIN_LENGTH=$STRICT_MIN_LENGTH
    MAX_LENGTH=$STRICT_MAX_LENGTH
    MIN_Q=$STRICT_MIN_Q
    #Convert the barcode sequence into slurm format
    ARRAY_SEQUENCE=$(IFS=,; echo "${BARCODE_SEQUENCE[*]}")
    
    cat << EOF > NanoFilt_2_slurm.sh
#!/bin/bash

#SBATCH --time=0-${JOB_TIME}:00:00  # max job runtime
#SBATCH --cpus-per-task=$QC_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "NanoFilt_2"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=$ARRAY_SEQUENCE

WORK_DIR=$WORK_DIR
cd $WORK_DIR/4_NanoFilt_2
EOF

    if [ $CONDA == "conda" ]; then
        echo 'eval "$(conda shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
        echo 'source activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
    elif [ $CONDA == "mamba" ]; then
        echo 'eval "$(mamba shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
        echo 'mamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
    elif [ $CONDA == "micromamba" ]; then
        echo 'eval "$(micromamba shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
        echo 'micromamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
    fi

    cat << EOF >> NanoFilt_2_slurm.sh
READ_DIR=$WORK_DIR/3_Demultiplex

STRICT_MIN_LENGTH=$MIN_LENGTH
STRICT_MAX_LENGTH=$MAX_LENGTH
STRICT_MIN_Q=$MIN_Q

for n in {1..$LIBRARY};
EOF

    cat << 'EOF' >> NanoFilt_2_slurm.sh
do
    gunzip -c $READ_DIR/"$n"_*0"${SLURM_ARRAY_TASK_ID}".fastq.gz \
        | NanoFilt --length $STRICT_MIN_LENGTH --maxlength $STRICT_MAX_LENGTH -q $STRICT_MIN_Q \
        | gzip > "${n}"-"${SLURM_ARRAY_TASK_ID}".fastq.gz
done
EOF

# If there are multiple primer sets/lengths then pull those variables from the config file
elif [ $MULTI_LENGTH == "TRUE" ]; then
    #Loop through the amplicon sets 
    for m in $(seq 1 $NUM_OF_AMPLICON_SETS); 
    do
        #Get the variables for that set, uses indirect variable calling
        MIN_NAME="MIN_$m"
        MIN_LENGTH="${!MIN_NAME}"
        
        MAX_NAME="MAX_$m"
        MAX_LENGTH="${!MAX_NAME}"
        
        MINQ_NAME="MINQ_$m"
        MIN_Q="${!MINQ_NAME}"
        
        LIB_SET_NAME="LIB_SET_$m"
        LIB_SET=$LIB_SET_NAME[@]
        # Don't forget the parentheses when passing an array
        LIB_SEQUENCE=(${!LIB_SET})
        
        BARCODE_SET_NAME="BARCODE_SET_$m"
        tmp=$BARCODE_SET_NAME[@]
        BARCODE_SET_SEQUENCE=(${!tmp})
        # Convert the array to a comma delim list for the slurm script
        ARRAY_SEQUENCE=$(IFS=,; echo "${BARCODE_SET_SEQUENCE[*]}")
        
        #Write a slurm script for each amplicon set
        cat << EOF > NanoFilt_set"$m"_slurm.sh
#!/bin/bash

#SBATCH --time=0-${JOB_TIME}:00:00  # max job runtime
#SBATCH --cpus-per-task=$QC_THREADS  # number of processor cores
#SBATCH --nodes=1  # number of nodes
#SBATCH --mem=200G  # max memory
#SBATCH -J "NanoFilt_2"  # job name
#SBATCH --mail-user=$EMAIL  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=$ARRAY_SEQUENCE

WORK_DIR=$WORK_DIR
cd $WORK_DIR/4_NanoFilt_2
EOF

        if [ $CONDA == "conda" ]; then
            echo 'eval "$(conda shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
            echo 'source activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
        elif [ $CONDA == "mamba" ]; then
            echo 'eval "$(mamba shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
            echo 'mamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
        elif [ $CONDA == "micromamba" ]; then
            echo 'eval "$(micromamba shell hook --shell bash)"' >> NanoFilt_2_slurm.sh
            echo 'micromamba activate $WORK_DIR/envs/cutadapt-env' >> NanoFilt_2_slurm.sh
        fi

        cat << EOF >> NanoFilt_set"$m"_slurm.sh
READ_DIR=$WORK_DIR/3_Demultiplex

STRICT_MIN_LENGTH=$MIN_LENGTH
STRICT_MAX_LENGTH=$MAX_LENGTH
STRICT_MIN_Q=$MIN_Q

for n in ${LIB_SEQUENCE[@]};
EOF

        cat << 'EOF' >> NanoFilt_set"$m"_slurm.sh
do
    gunzip -c $READ_DIR/"$n"_*0"${SLURM_ARRAY_TASK_ID}".fastq.gz \
        | NanoFilt --length $STRICT_MIN_LENGTH --maxlength $STRICT_MAX_LENGTH -q $STRICT_MIN_Q \
        | gzip > "${n}"-"${SLURM_ARRAY_TASK_ID}".fastq.gz
done
EOF

    done
fi




################################################ Submit the slurm with each job dependent on the previous finishing ################################################

#Check if the corresponding folder is empty then submit the relevant slurm scripts to populate that file

if [ $MULTI_LENGTH == "TRUE" ]; then
    if [ -z "$( ls -A $WORK_DIR/2_NanoFilt_1 )" ]; then
        JOBID1=$(sbatch --parsable NanoFilt_1_slurm.sh)
        JOBID2=$(sbatch --parsable --dependency=afterok:$JOBID1 Demux_slurm.sh)
        for m in {1..$NUM_OF_AMPLICON_SETS};do sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_set"$m"_slurm.sh; done
    elif [ -z "$( ls -A $WORK_DIR/3_Demultiplex )" ]; then
        JOBID2=$(sbatch --parsable Demux_slurm.sh)
        sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_2_slurm.sh
        for m in {1..$NUM_OF_AMPLICON_SETS};do sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_set"$m"_slurm.sh; done
    elif [ -z "$( ls -A $WORK_DIR/4_NanoFilt_2 )" ]; then
        for m in {1..$NUM_OF_AMPLICON_SETS};do sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_set"$m"_slurm.sh; done
    else
        echo "Folders 2-4 already contain files"
    fi
else
    if [ -z "$( ls -A $WORK_DIR/2_NanoFilt_1 )" ]; then
        JOBID1=$(sbatch --parsable NanoFilt_1_slurm.sh)
        JOBID2=$(sbatch --parsable --dependency=afterok:$JOBID1 Demux_slurm.sh)
        sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_2_slurm.sh
    elif [ -z "$( ls -A $WORK_DIR/3_Demultiplex )" ]; then
        JOBID2=$(sbatch --parsable Demux_slurm.sh)
        sbatch --parsable --dependency=afterok:$JOBID2 NanoFilt_2_slurm.sh
    elif [ -z "$( ls -A $WORK_DIR/4_NanoFilt_2 )" ]; then
        sbatch --parsable NanoFilt_2_slurm.sh
    else
        echo "Folders 2-4 already contain files"
    fi
fi

if [ $CONDA == "conda" ]; then
    conda deactivate
elif [ $CONDA == "mamba" ]; then
    mamba deactivate
elif [ $CONDA == "micromamba" ]; then
    micromamba deactivate
fi
