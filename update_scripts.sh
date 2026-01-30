#!/bin/bash 
#### update scripts

#user defined variables
source config.txt


# Replace the clust.smk and quant.smk files
cp $WORK_DIR/laca_changes/clust.smk $LACA_DIR/laca/workflow/rules/
cp $WORK_DIR/laca_changes/quant.smk $LACA_DIR/laca/workflow/rules/

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    #put the ACT scripts in the path
    conda activate ACT-env
    cd $WORK_DIR/scripts
    chmod +x `ls`
    cp ./* $CONDA_PREFIX/bin/
    
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    #put the ACT scripts in the path
    mamba activate ACT-env
    cd $WORK_DIR/scripts
    chmod +x `ls`
    cp ./* $CONDA_PREFIX/bin/
    
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    #put the ACT scripts in the path
    micromamba activate ACT-env
    cd $WORK_DIR/scripts
    chmod +x `ls`
    cp ./* $CONDA_PREFIX/bin/
    
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi

cd $WORK_DIR/scripts
cp ./* $ENV_DIR/nanoplot-env/bin/
cp ./* $ENV_DIR/cutadapt-env/bin/
cp ./* $ENV_DIR/taxonomy-env/bin/
cp ./* $ENV_DIR/database-env/bin/
cp ./* $ENV_DIR/consensus-env/bin/
echo "Environments are ready"
