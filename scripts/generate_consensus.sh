#!/bin/bash

## generate consensus abundance table

#user defined variables
source config.txt

[[ -z "$(ls -A $WORK_DIR/6_emu/read_assignments/)" ]] && { echo "${WORK_DIR}/6_emu/read_assignments/ is empty"; exit 1; }
[[ -z "$(ls -A $WORK_DIR/7_sintax/)" ]] && { echo "${WORK_DIR}/7_sintax/ is empty"; exit 1; }
[[ ! -e $WORK_DIR/5_laca/quant/seqID_to_otu.tsv ]] && { echo "5_laca/quant/seqID_to_otu.tsv does not exist"; exit 1; }
[[ ! -e $WORK_DIR/5_laca/quant/seqID_to_otu.tsv ]] && { echo "5_laca/quant/seqID_to_otu.tsv does not exist"; exit 1; }
[[ ! -e $WORK_DIR/5_laca/sintax_OTUs.tsv ]] && { echo "5_laca/sintax_OTUs.tsv does not exist"; exit 1; }


#run the R script