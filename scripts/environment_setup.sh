#!/bin/bash
#### Set up the environments

#user defined variables
source config.txt

# Clone LACA
cd $WORK_DIR
git clone https://github.com/yanhui09/laca.git


# Set up environments
cd $WORK_DIR/envs

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    conda env create --prefix ./nanoplot-env -f nanoplot.yaml
    conda env create --prefix ./cutadapt-env -f cutadapt.yaml
    conda env create --prefix ./taxonomy-env -f taxonomy.yaml
    cd $WORK_DIR/laca
    conda env create -n laca -f env.yaml 
    conda activate laca
    pip install --editable .
    conda deactivate
    echo "Environments ready"
    
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba env create --prefix ./nanoplot-env -f nanoplot.yaml
    mamba env create --prefix ./cutadapt-env -f cutadapt.yaml
    mamba env create --prefix ./taxonomy-env -f taxonomy.yaml
    cd $WORK_DIR/laca
    mamba env create -n laca -f env.yaml 
    mamba activate laca
    pip install --editable .
    mamba deactivate
    echo "Environments ready"
    
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba env create --prefix ./nanoplot-env -f nanoplot.yaml
    micromamba env create --prefix ./cutadapt-env -f cutadapt.yaml
    micromamba env create --prefix ./taxonomy-env -f taxonomy.yaml
    cd $WORK_DIR/laca
    micromamba env create -n laca -f env.yaml 
    micromamba activate laca
    pip install --editable .
    micromamba deactivate
    echo "Environments ready"
    
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi
