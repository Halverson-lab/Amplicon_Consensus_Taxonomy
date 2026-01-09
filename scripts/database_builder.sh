#!/bin/bash
#database builder and editor


add_flag=false
ncbi_flag=false
build_flag=false
default_flag=false
sintax_flag=false
group_flag=false

usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help      Display this help message"
 echo " -b, --build     Build new database"
 echo " -d, --default   Build new database from latest NCBI 16S RefSeq and rrnDB"
 echo " -g, --group     Group and rename highly similar sequences in the database"
 echo " -s, --sintax    Build new database from a sintax/UNITE formatted database"
 echo " -a, --add       Add user provided sequences to database"
 echo " -n, --ncbi      Add sequences to database using list of NCBI accessions"
}

if [ $# -eq 0 ]; then
  echo "Error: At least one flag is required."
  usage
  exit 1
fi


while [[ $# -gt 0 ]]; do
  case "$1" in
    -b | --build) 
      build_flag=true
      echo "Building database from provided files" >&2
      ;;
    -d | --default)
      build_flag=true
      default_flag=true
      group_flag=true
      echo "Building database from the latest NCBI 16S RefSeq and rrnDB files" >&2
      ;;
    -a | --add) 
      add_flag=true
      echo "Adding new sequences to database" >&2
      ;;
    -n | --ncbi)
      ncbi_flag=true
      echo "Retrieving sequences from NCBI and adding to database" >&2
      ;;
    -s | --sintax)
      build_flag=true
      sintax_flag=true
      echo "Building database using a sintax formatted database" >&2
      ;;
    -g | --group)
      group_flag=true
      echo "Database sequences will be grouped by similarity" >&2
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

#user defined variables
source config.txt

[[ -z "$DATABASE_DIR" ]] && { echo "DATABASE_DIR is required" ; exit 1; }

#check for necessary files if building from provided files
if [[ $build_flag == "true" ]]; then
    if [[ $default_flag == "true" ]]; then
        [[ -z $RRNDB  && -e $DATABASE_DIR/rrnDB.fasta ]] && { echo "please provide the URL to the latest rrnDB or provide is in fasta format in the database directory" ; exit 1; }
    elif [[ $sintax_flag == "true" ]]; then
        [[ -z "$BUILD_USER_SEQ" ]] && { echo "USER_SEQ is required for --sintax" ; exit 1; }
    else 
        [[ -z "$BUILD_USER_SEQ" ]] && { echo "USER_SEQ is required if building a custom database, run database_builder.sh -d to build a default database" ; exit 1; }
        [[ -z "$BUILD_USER_SEQ2TAX" ]] && { echo "USER_SEQ2TAX is required if building a custom database, run database_builder.sh -d to build a default database" ; exit 1; }
    fi
fi

[[ $add_flag == "true" ]] && [[ -z "$ADD_USER_SEQ" ]] && { echo "USER_SEQ is required for --add" ; exit 1; }
[[ $add_flag == "true" ]] && [[ -z "$ADD_USER_SEQ2TAX" ]] && { echo "USER_SEQ2TAX is required for --add" ; exit 1; }
[[ $ncbi_flag == "true" ]] && [[ -z "$ACC_LIST" ]] && { echo "ACC_LIST is required for --ncbi" ; exit 1; }
[[ $group_flag == "true" ]] && [[ -z "$SIM_THRESH" ]] && { echo "SIM_THRESH is required for --group" ; exit 1; }


# Set up environments
cd $WORK_DIR

if [ $CONDA == "conda" ]; then
    eval "$(conda shell hook --shell bash)"
    source activate $ENV_DIR/database-env
elif [ $CONDA == "mamba" ]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $ENV_DIR/database-env
elif [ $CONDA == "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_DIR/database-env
else
    echo "CONDA can be conda, mamba, or micromamba" 
    exit 1
fi


cd $DATABASE_DIR
rm database_builder.log
exec 1> database_builder.log

if [[ $default_flag  == "true" || $sintax_flag == "true" || $ncbi_flag == "true" || $add_flag == "true" ]]; then
    if [[ $default_flag  == "true" ]]; then
        #build with latest NCBI RefSeq and the rrnDB
    
        #download database
        [[ ! -e blast_16S_DB ]] && { mkdir blast_16S_DB ; } 
        cd blast_16S_DB 
        update_blastdb.pl --decompress 16S_ribosomal_RNA
        #retrieve fasta
        blastdbcmd -entry all -db 16S_ribosomal_RNA -out ../ncbi_16S.fasta
        cd $DATABASE_DIR
    
        #if the fasta file was not provided then download a new one with the provided URL
        if [[ ! -e $DATABASE_DIR/rrnDB.fasta ]]; then
            wget $RRNDB
            unzip rrnDB*.zip && rm -r rrnDB*.zip
            mv rrnDB*.fasta ./rrnDB.fasta 
        fi
    
        ## re-format rrnDB formatted files
        # For later application we need fasta headers starting with useful sequence IDs
        # get fasta headers
        seqkit seq -n rrnDB.fasta > rrnDB_header_list.txt
    
        # get the RefSeq numbers from each header and generate new headers
        cat rrnDB_header_list.txt \
            | csvtk add-header -t -n fasta_header \
            | csvtk mutate -t -f 1 -p "(?:[^\|]*\|){2}(.*?)\|" -n refSeq \
            | csvtk mutate2 -t -n newHeader -e '$refSeq + " " + $fasta_header'> rrnDB_header_refseq_list.txt 
    
        # File of rrnDB header to new rrnDB header with RefSeq ID as the sequence ID
        cat rrnDB_header_refseq_list.txt | csvtk cut -t -f 1,3 | csvtk del-header -t > old2newheader.tsv
        # replace the headers
        seqkit replace -p '(.+)' -r '{kv}' -k old2newheader.tsv rrnDB.fasta > rrnDB_renamed.fasta
    
        # Use refSeq numbers to get taxid
        cat rrnDB_header_refseq_list.txt | csvtk cut -t -f 2 | csvtk del-header -t | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > rrnDB_seq2taxid.txt
    
        # Filter out any of the fasta seq that were removed
        csvtk cut -t -f 1 rrnDB_seq2taxid.txt > seq2keep.txt
        seqkit grep -f seq2keep.txt rrnDB_renamed.fasta -o ACT_rrnDB.fasta
    
        #ACT_rrnDB.fasta and rrnDB_seq2taxid.txt are the only ones we need to keep
        rm rrnDB_header_list.txt rrnDB_header_refseq_list.txt old2newheader.tsv rrnDB_renamed.fasta seq2keep.txt
    
        ### re-format NCBI formatted files
        #get fasta headers
        seqkit seq -n ncbi_16S.fasta > ncbi_16S_header.txt
    
        #get the RefSeq numbers from each header
        cat ncbi_16S_header.txt | csvtk add-header -t -n fasta_header | csvtk mutate -t -f 1 -p "([^\s]+)" -n refSeq > ncbi_16S_header_refseq.txt
    
        #Use refSeq numbers to get taxid
        cat ncbi_16S_header_refseq.txt | csvtk cut -t -f 2 | csvtk del-header -t | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > refseq2taxid.txt
    
        csvtk add-header -t -n refSeq,taxid refseq2taxid.txt > refseq2taxid2.txt
    
        # Any seqs where the RefSeq record was suppressed will be removed from the list (this shouldn't change anything as long as your fasta file is current)
        csvtk join -t -f refSeq ncbi_16S_header_refseq.txt refseq2taxid2.txt | csvtk cut -t -f 2,3 |csvtk del-header -t > ncbi_seq2taxid.txt
    
        # ncbi_16S.fasta and ncbi_seq2taxid.txt are the ones to keep
        rm ncbi_16S_header.txt ncbi_16S_header_refseq.txt refseq2taxid.txt refseq2taxid2.txt
    
        # combine them into one 
        cat rrnDB_seq2taxid.txt ncbi_seq2taxid.txt > seq2taxid.txt
        cat ACT_rrnDB.fasta ncbi_16S.fasta > sequences.fasta
    
    
        # Build the taxonomy table
        # getting the NCBI taxonomy
        wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        mkdir taxdump
        tar -xzf taxdump.tar.gz -C taxdump && rm taxdump.tar.gz
        # Tell taxonkit where to find taxdump
        export TAXONKIT_DB=$DATABASE_DIR/taxdump
    
        # Build taxonomy.tsv
        csvtk cut -t -f 2 seq2taxid.txt \
        | taxonkit lineage \
        | taxonkit reformat2 -t -f "{domain|superkingdom|kingdom|acellular root}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}" \
        | csvtk cut -t -f -2 \
        | csvtk add-header -t -n taxid,domain,phylum,class,order,family,genus,species,t_domain,t_phylum,t_class,t_order,t_family,t_genus,t_species \
        | csvtk uniq -t -f taxid -o taxonomy.tsv
    fi
    
    if [[ $sintax_flag == "true" ]]; then
    
        seqkit seq -n $BUILD_USER_SEQ | csvtk add-header -t -n seqid -o read_taxon_list.tsv

        cat $BUILD_USER_SEQ  | seqkit replace -p ";.+" > sequences.fasta
        
        cat read_taxon_list.tsv \
            | csvtk sep -t -f 1 -s ";" -n seq_id,taxon,tmp \
            | csvtk cut -t -f 2,3 \
            | csvtk mutate -t -f taxon -n kingdom -p "k:(.*)" --na \
            | csvtk mutate -t -f taxon -n domain -p "d:(.*)" --na \
            | csvtk mutate -t -f taxon -n phylum -p "p:(.*)" --na \
            | csvtk mutate -t -f taxon -n class -p "c:(.*)" --na \
            | csvtk mutate -t -f taxon -n order -p "o:(.*)" --na \
            | csvtk mutate -t -f taxon -n family -p "f:(.*)" --na \
            | csvtk mutate -t -f taxon -n genus -p "g:(.*)" --na \
            | csvtk mutate -t -f taxon -n species -p "s:(.*)" --na \
            | csvtk replace -t -f kingdom -p ",.*" -r "" \
            | csvtk replace -t -f domain -p ",.*" -r "" \
            | csvtk replace -t -f phylum -p ",.*" -r "" \
            | csvtk replace -t -f class -p ",.*" -r "" \
            | csvtk replace -t -f order -p ",.*" -r "" \
            | csvtk replace -t -f family -p ",.*" -r "" \
            | csvtk replace -t -f genus -p ",.*" -r "" \
            | csvtk replace -t -f species -p ",.*" -r "" \
            | csvtk cut -t -f -2 > taxonomy_for_taxdump.tsv
            
        taxonkit create-taxdump -A 1 taxonomy_for_taxdump.tsv --out-dir taxdump
        
        rm read_taxon_list.tsv taxonomy_for_taxdump.tsv
        
        #tell taxonkit where to find the new taxdump files
        export TAXONKIT_DB=taxdump/
        
        #cp the taxid.map to make it easier to find
        cp taxdump/taxid.map ./seq2taxid.txt
        
        csvtk cut -t -f 2 seq2taxid.txt \
            | taxonkit lineage \
            | taxonkit reformat2 -t -f "{domain|superkingdom|kingdom|acellular root}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}" \
            | csvtk cut -t -f -2 \
            | csvtk add-header -t -n taxid,domain,phylum,class,order,family,genus,species,t_domain,t_phylum,t_class,t_order,t_family,t_genus,t_species \
            | csvtk uniq -t -f taxid -o taxonomy.tsv
    
    fi
    
    if [[ $ncbi_flag == "true" ]]; then
        #retieve based on list of accessions
        #/db_xref="RFAM:RF00177"
        ## takes in an accession list and spits out fasta files with the header locus_tag, record_name, and taxon:taxid  
        for ACC in $(cat $ACC_LIST); do 
            #download the genbank file
            esearch -db nuccore -query "${ACC}" | efetch -format gbwithparts > tmp.gbk
    
            if [ ! -s "tmp.gbk" ]; then 
                echo "Accession ${ACC} is invalid, please double check the input" >&2
            
            else
                #script to extract 16S rRNA sequences from the genbank file as fastas
                extract_16S_from_genbank_to_fasta.py tmp.gbk "${ACC%.*}".fasta
    
                #Append the fastas onto the fasta file
                cat "${ACC%.*}".fasta >> sequences.fasta
    
                #add seqid to taxid
                #get fasta headers
                seqkit seq -n "${ACC%.*}".fasta \
                | csvtk add-header -t -n fasta_header \
                | csvtk mutate -t -f 1 -p "taxon:(\d+)" -n taxid \
                | csvtk mutate -t -f 1 -p "([^\s]+)" -n ID \
                | csvtk cut -t -f 3,2 \
                | csvtk del-header -t > "${ACC%.*}"_seq2taxid.txt
                cat "${ACC%.*}"_seq2taxid.txt >> seq2taxid.txt
                rm "${ACC%.*}"*
            fi
    
        done
    
        #retrieve lineage info
        export TAXONKIT_DB=$DATABASE_DIR/taxdump
        csvtk cut -t -f 2 seq2taxid.txt \
            | taxonkit lineage \
            | taxonkit reformat2 -t -f "{domain|superkingdom|kingdom|acellular root}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}" \
            | csvtk cut -t -f -2 \
            | csvtk add-header -t -n taxid,domain,phylum,class,order,family,genus,species,t_domain,t_phylum,t_class,t_order,t_family,t_genus,t_species \
            | csvtk uniq -t -f taxid -o taxonomy.tsv
    fi
    
    if [[ $add_flag == "true" ]]; then
        # user provided
        cat $ADD_USER_SEQ >> sequences.fasta
        cat $ADD_USER_SEQ2TAX >> seq2taxid.txt
        if [[ -z "$ADD_USER_TAX" ]]; then
            cat $ADD_USER_TAX >> taxonomy.tsv
        fi 
    fi
    
    #remove special characters from the taxonomy
    cat taxonomy.tsv | tr -d ';#&%' | tr -d "\'" > safe_taxonomy.tsv && mv ./safe_taxonomy.tsv ./taxonomy.tsv
    
    # Build EMU formatted database
    emu build-database ACT_DB --sequences sequences.fasta --seq2tax seq2taxid.txt --taxonomy-list taxonomy.tsv
    mv ACT_DB/* ./ 
    rm -r ACT_DB
fi

# Directory should now contain species_taxid.fasta and taxonomy.tsv, whether just generated or user provided
[[ -z species_taxid.fasta ]] && { echo "species_taxid.fasta does not exist, please provide database files or run database_builder.sh -d to build a default database" >&2; exit 1; }
[[ -z taxonomy.tsv ]] && { echo "taxonomy.tsv does not exist, please provide database files or run database_builder.sh -d to build a default database" >&2; exit 1; }


### dereplicate the database
if [[ $group_flag == "true" ]]; then
    # abbreviate names because samtools has a maximum limit
    seqkit seq -i species_taxid.fasta > abbrev_species_taxid.fasta

    # minimap all against all and remove the abbreviated name file once done
    minimap2 -ax lr:hq -t $QC_THREADS -K 500M -N 50 abbrev_species_taxid.fasta abbrev_species_taxid.fasta -o database_alignment.sam
    rm abbrev_species_taxid.fasta

    #script to convert the sam file into a csv containing the relevant info
    minimap_to_csv.py database_alignment.sam database_alignment_stats.csv

    # filter out self to self and matches in the same taxid, then calculate the percent match
    # final output is a file of read pairs that failed the similarity threshold
    cat database_alignment_stats.csv \
        | csvtk filter2 -f '$QNAME != $RNAME' \
        | csvtk mutate --after QNAME -f QNAME -n QTAXID -p "^([0-9].*?):" \
        | csvtk mutate --after RNAME -f RNAME -n RTAXID -p "^([0-9].*?):" \
        | csvtk filter2 -f '$QTAXID != $RTAXID' \
        | csvtk mutate2 --after RLENGTH -n MINLENGTH -e '$QLENGTH < $RLENGTH ? $QLENGTH : $RLENGTH' \
        | csvtk mutate2 --after NM -n PERCENTMATCH -e '$NM / $MINLENGTH' -w 5 \
        | csvtk filter -f "PERCENTMATCH<$SIM_THRESH"  -o database_alignment_failed_reads.csv

    seqkit seq -n -i species_taxid.fasta > seq_id_list.txt
    # take in the config file and failed read alignments and output group notations, updated fasta headers, and updated taxonomy
    minimap_to_group.R $WORK_DIR/config.txt database_alignment_failed_reads.csv
    rm seq_id_list.txt

    # add the updated files to the existing db
    mv taxonomy.tsv ./ungrouped_taxonomy.tsv
    csvtk concat -t ungrouped_taxonomy.tsv grouped_taxonomies.tsv > taxonomy.tsv
    mv species_taxid.fasta ./ungrouped_species_taxid.fasta
    seqkit grep -v -f sequences_to_be_removed.csv ungrouped_species_taxid.fasta -o trimmed_species_taxid.fasta
    seqkit replace -p '^(\S+)(.+?)$' -r '{kv}$2' -k grouped-seq-headers.tsv trimmed_species_taxid.fasta --keep-key > species_taxid.fasta
    rm trimmed_species_taxid.fasta
fi

# make the name list for converting from EMU to Sintax
seqkit seq species_taxid.fasta -n > name_list.txt

# run the r script for converting EMU headers to Sintax headers
emu_to_sintax_db_converter.R
seqkit replace -p "\s.+" species_taxid.fasta | seqkit replace -p "(.+)" -r '{kv}' -k emu2sintax-header.tsv > sintax_db.fasta
rm name_list.txt
