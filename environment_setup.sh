#!/bin/bash 
#### Set up the environments

#user defined variables
source config.txt

reconstruct_flag=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r | --reconstruct) 
      reconstruct_flag=true
      echo "Removing and reconstructing ACT environments" >&2
      ;;
    \?)
      echo "Invalid option" >&2
      exit 1
      ;;
  esac
  shift
done

# Clone LACA
cd $WORK_DIR
if [ $reconstruct_flag == "true" ]; then
    rm -r $ENV_DIR/nanoplot-env
    rm -r $ENV_DIR/cutadapt-env
    rm -r $ENV_DIR/taxonomy-env
    rm -r $ENV_DIR/database-env
    rm -r $ENV_DIR/consensus-env
    rm -r $LACA_DIR
fi


# Set up environments
cd $ENV_DIR

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    
    if [ $reconstruct_flag == "true" ]; then
        conda env remove -n ACT-env
        conda env remove -n laca
    fi
    
    [[ -z "$(conda list --name ACT-env)" ]] && { conda env create -n ACT-env -f ACT.yaml; }
    [[ ! -e $ENV_DIR/nanoplot-env ]] && { conda env create --prefix ./nanoplot-env -f nanoplot.yaml; }
    [[ ! -e $ENV_DIR/cutadapt-env ]] && { conda env create --prefix ./cutadapt-env -f cutadapt.yaml; }
    [[ ! -e $ENV_DIR/taxonomy-env ]] && { conda env create --prefix ./taxonomy-env -f taxonomy.yaml; }
    [[ ! -e $ENV_DIR/database-env ]] && { conda env create --prefix ./database-env -f database.yaml; }
    [[ ! -e $ENV_DIR/consensus-env ]] && { conda env create --prefix ./consensus-env -f consensus.yaml; }
    
    conda activate ./database-env
    git clone https://github.com/rhowardstone/AmpliconHunter2.git
    cd AmpliconHunter2; make;
    cp amplicon_hunter $ENV_DIR/database-env/bin/
    cd $ENV_DIR
    conda deactivate

    #prep laca env
    if [ -z "$( ls -A $LACA_DIR )" ]; then
        cd "$(dirname "$LACA_DIR")"
        git clone https://github.com/yanhui09/laca.git
        
        cd $LACA_DIR
        conda env create -n laca -f env.yaml 
        conda activate laca
        pip install --editable .
        conda deactivate

        # Replace the clust.smk and quant.smk files
        cp $WORK_DIR/laca_changes/clust.smk $LACA_DIR/laca/workflow/rules/
        cp $WORK_DIR/laca_changes/quant.smk $LACA_DIR/laca/workflow/rules/
        
        # In the /laca/laca/workflow/envs directory, modify the yacrd.yaml file to specify minimap2=2.18
        # In the /laca/laca/workflow/envs directory, modify the medaka.yaml file to change the order of the channels to make conda-forge first, and specify samtools=1.11
        cp $WORK_DIR/laca_changes/medaka.yaml $LACA_DIR/laca/workflow/envs/
        cp $WORK_DIR/laca_changes/yacrd.yaml $LACA_DIR/laca/workflow/envs/
                
    fi

    #put the ACT scripts in the path
    conda activate ACT-env
    cd $WORK_DIR/scripts
    chmod +x `ls`
    cp ./* $CONDA_PREFIX/bin/
    
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    
    if [ $reconstruct_flag == "true" ]; then
        mamba env remove -n ACT-env
        mamba env remove -n laca
    fi
    
    [[ -z "$(mamba list --name ACT-env)" ]] && { mamba env create -n ACT-env -f ACT.yaml; }
    [[ ! -e $ENV_DIR/nanoplot-env ]] && { mamba env create --prefix ./nanoplot-env -f nanoplot.yaml; }
    [[ ! -e $ENV_DIR/cutadapt-env ]] && { mamba env create --prefix ./cutadapt-env -f cutadapt.yaml; }
    [[ ! -e $ENV_DIR/taxonomy-env ]] && { mamba env create --prefix ./taxonomy-env -f taxonomy.yaml; }
    [[ ! -e $ENV_DIR/database-env ]] && { mamba env create --prefix ./database-env -f database.yaml; }
    [[ ! -e $ENV_DIR/consensus-env ]] && { mamba env create --prefix ./consensus-env -f consensus.yaml; }

    mamba activate ./database-env
    git clone https://github.com/rhowardstone/AmpliconHunter2.git
    cd AmpliconHunter2; make;
    cp amplicon_hunter $ENV_DIR/database-env/bin/
    cd $ENV_DIR
    mamba deactivate
    
    #prep laca env
    if [ -z "$( ls -A $LACA_DIR )" ]; then
        cd "$(dirname "$LACA_DIR")"
        git clone https://github.com/yanhui09/laca.git
        
        cd $LACA_DIR
        mamba env create -n laca -f env.yaml 
        mamba activate laca
        pip install --editable .
        mamba deactivate
        
        # Replace the clust.smk and quant.smk files
        cp $WORK_DIR/laca_changes/clust.smk $LACA_DIR/laca/workflow/rules/
        cp $WORK_DIR/laca_changes/quant.smk $LACA_DIR/laca/workflow/rules/
        
        # In the /laca/laca/workflow/envs directory, modify the yacrd.yaml file to specify minimap2=2.18
        # In the /laca/laca/workflow/envs directory, modify the medaka.yaml file to change the order of the channels to make conda-forge first, and specify samtools=1.11
        cp $WORK_DIR/laca_changes/medaka.yaml $LACA_DIR/laca/workflow/envs/
        cp $WORK_DIR/laca_changes/yacrd.yaml $LACA_DIR/laca/workflow/envs/
                
    fi
    
    #put the ACT scripts in the path
    mamba activate ACT-env
    cd $WORK_DIR/scripts
    chmod +x `ls`
    cp ./* $CONDA_PREFIX/bin/
    
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"

    if [ $reconstruct_flag == "true" ]; then
        micromamba env remove -n ACT-env
        micromamba env remove -n laca
    fi
    
    [[ -z "$(micromamba list --name ACT-env)" ]] && { micromamba env create -n ACT-env -f ACT.yaml; }
    [[ ! -e $ENV_DIR/nanoplot-env ]] && { micromamba env create --prefix ./nanoplot-env -f nanoplot.yaml; }
    [[ ! -e $ENV_DIR/cutadapt-env ]] && { micromamba env create --prefix ./cutadapt-env -f cutadapt.yaml; }
    [[ ! -e $ENV_DIR/taxonomy-env ]] && { micromamba env create --prefix ./taxonomy-env -f taxonomy.yaml; }
    [[ ! -e $ENV_DIR/database-env ]] && { micromamba env create --prefix ./database-env -f database.yaml; }
    [[ ! -e $ENV_DIR/consensus-env ]] && { micromamba env create --prefix ./consensus-env -f consensus.yaml; }

    micromamba activate ./database-env
    git clone https://github.com/rhowardstone/AmpliconHunter2.git
    cd AmpliconHunter2; make;
    cp amplicon_hunter $ENV_DIR/database-env/bin/
    cd $ENV_DIR
    micromamba deactivate
    
    #prep laca env
    if [ -z "$( ls -A $LACA_DIR )" ]; then
        cd "$(dirname "$LACA_DIR")"
        git clone https://github.com/yanhui09/laca.git
        
        cd $LACA_DIR
        micromamba env create -n laca -f env.yaml 
        micromamba activate laca
        pip install --editable .
        micromamba deactivate

        # Replace the clust.smk and quant.smk files
        cp $WORK_DIR/laca_changes/clust.smk $LACA_DIR/laca/workflow/rules/
        cp $WORK_DIR/laca_changes/quant.smk $LACA_DIR/laca/workflow/rules/
        
        # In the /laca/laca/workflow/envs directory, modify the yacrd.yaml file to specify minimap2=2.18
        # In the /laca/laca/workflow/envs directory, modify the medaka.yaml file to change the order of the channels to make conda-forge first, and specify samtools=1.11
        cp $WORK_DIR/laca_changes/medaka.yaml $LACA_DIR/laca/workflow/envs/
        cp $WORK_DIR/laca_changes/yacrd.yaml $LACA_DIR/laca/workflow/envs/
    fi
    
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
