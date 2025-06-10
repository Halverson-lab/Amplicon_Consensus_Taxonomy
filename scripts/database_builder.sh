#!/bin/bash
#database builder and editor


add_flag=false
ncbi_flag=false
build_flag=false
default_flag=false

usage() {
 echo "Usage: $0 [OPTIONS]"
 echo "Options:"
 echo " -h, --help      Display this help message"
 echo " -b, --build     Build new database"
 echo " -d, --default   Build new database from latest NCBI 16S RefSeq and rrnDB"
 echo " -a, --add       Add user provided sequences to database"
 echo " -n, --ncbi      Add sequences to database using list of NCBI accesions"
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
      echo "Building database from the latest NCBI 16S RefSeq and rrnDB files" >&2
      ;;
    -a | --add) 
      add_flag=true
      echo "adding new sequences to database" >&2
      ;;
    -n | --ncbi)
      ncbi_flag=true
      echo "retrieving sequences from NCBI and adding to database" >&2
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
    if [[ $default_flag == "false" ]]; then
        [[ -e $DATABASE_DIR/sequences.fasta ]] && { echo "sequences.fasta is not in the database directory, please provide database files or run database_builder.sh -bd to build a default database" ; exit 1; }
        [[ -e $DATABASE_DIR/seq2taxid.txt ]] && { echo "seq2taxid.txt  is not in the database directory, please provide database files or run database_builder.sh -bd to build a default database" ; exit 1; }
    fi

    if [[ $default_flag == "true" ]]; then
        [[ -z $RRNDB  && -e $DATABASE_DIR/rrnDB.fasta ]] && { echo "please provide the URL to the latest rrnDB or provide is in fasta format in the database directory" ; exit 1; }
    fi
    
fi

if [[ $build_flag == "true" && $default_flag == "false" ]]; then
    [[ -z $DATABASE_DIR/sequences.fasta ]] && { echo "sequences.fasta is not in the database directory, please provide database files or run database_builder.sh -bd to build a default database" ; exit 1; }
    [[ -z $DATABASE_DIR/seq2taxid.txt ]] && { echo "seq2taxid.txt  is not in the database directory, please provide database files or run database_builder.sh -bd to build a default database" ; exit 1; }
fi

#if adding to existing then need these params
if [[ $add_flag == "true" ]]; then
    [[ -z "$USER_SEQ" ]] && { echo "USER_SEQ is required for --add" ; exit 1; }
    [[ -z "$USER_SEQ2TAX" ]] && { echo "USER_SEQ2TAX is required for --add" ; exit 1; }
fi

if [[ $ncbi_flag == "true" ]]; then
    [[ -z "$ACC_LIST" ]] && { echo "ACC_LIST is required for --ncbi" ; exit 1; }
fi



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

if [[ $build_flag == "true" ]]; then
    #build with latest NCBI RefSeq and the rrnDB
    if [[ $default_flag == "true" ]]; then
        #download database
        mkdir blast_16S_DB && cd blast_16S_DB 
        update_blastdb.pl --decompress 16S_ribosomal_RNA
        #retrieve fasta
        blastdbcmd -entry all -db 16S_ribosomal_RNA -out ../ncbi_16S.fasta
        cd ../

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
    fi

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
    | taxonkit reformat -t -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
    | csvtk cut -t -f -2 \
    | csvtk add-header -t -n taxid,domain,phylum,class,order,family,genus,species,t_domain,t_phylum,t_class,t_order,t_family,t_genus,t_species -o taxonomy.tsv

    #remove special characters from the taxonomy
    cat taxonomy.tsv | tr -d ';#&%' | tr -d "\'" > safe_taxonomy.tsv && mv ./safe_taxonomy.tsv ./taxonomy.tsv

    # Build EMU formatted database
    emu build-database ACT_DB --sequences sequences.fasta --seq2tax seq2taxid.txt --taxonomy-list taxonomy.tsv
    mv ACT_DB/* ./ 
    rm -r ACT_DB

    # Directory should now contain species_taxid.fasta and taxonomy.tsv, whether just generated or user provided
    [[ -z species_taxid.fasta ]] && { echo "species_taxid.fasta does not exist, please provide database files or run database_builder.sh -bd to build a default database" >&2; exit 1; }
    [[ -z taxonomy.tsv ]] && { echo "taxonomy.tsv does not exist, please provide database files or run database_builder.sh -bd to build a default database" >&2; exit 1; }

    # make the name list for converting from EMU to Sintax
    seqkit seq species_taxid.fasta -n > name_list.txt
    
    # run the r script for converting EMU headers to Sintax headers
    emu_to_sintax_db_converter.R
    seqkit replace -p "\s.+" species_taxid.fasta | seqkit replace -p "(.+)" -r '{kv}' -k emu2sintax-header.tsv > sintax_db.fasta

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
        | taxonkit reformat -t -f "{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
        | csvtk cut -t -f -2 \
        | csvtk add-header -t -n taxid,domain,phylum,class,order,family,genus,species,t_domain,t_phylum,t_class,t_order,t_family,t_genus,t_species -o taxonomy.tsv

    #remove special characters from the taxonomy
    cat taxonomy.tsv | tr -d ';#&%' | tr -d "\'" > safe_taxonomy.tsv && mv ./safe_taxonomy.tsv ./taxonomy.tsv

    #rebuild the emu formatted database with the added sequences
    emu build-database ACT_DB --sequences sequences.fasta --seq2tax seq2taxid.txt --taxonomy-list taxonomy.tsv
    mv ACT_DB/* ./
    rm -r ACT_DB

    # make the name list for converting from EMU to Sintax
    seqkit seq species_taxid.fasta -n > name_list.txt

    # run the r script for converting EMU headers to Sintax headers
    emu_to_sintax_db_converter.R
    seqkit replace -p "\s.+" species_taxid.fasta | seqkit replace -p "(.+)" -r '{kv}' -k emu2sintax-header.tsv > sintax_db.fasta
fi

if [[ $add_flag == "true" ]]; then

    # user provided
    cat $USER_SEQ >> sequences.fasta
    cat $USER_SEQ2TAX >> seq2taxid.txt
    if [[ -z "$USER_TAX" ]]; then
        cat $USER_TAX >> taxonomy.tsv
        #remove special characters from the taxonomy
        cat taxonomy.tsv | tr -d ';#&%' | tr -d "\'" > safe_taxonomy.tsv && mv ./safe_taxonomy.tsv ./taxonomy.tsv

    fi 

    #rebuild the emu formatted database with the added sequences
    emu build-database ACT_DB --sequences sequences.fasta --seq2tax seq2taxid.txt --taxonomy-list taxonomy.tsv
    mv ACT_DB/* ./ 
    rm -r ACT_DB

    # make the name list for converting from EMU to Sintax
    seqkit seq species_taxid.fasta -n > name_list.txt

    # run the r script for converting EMU headers to Sintax headers
    emu_to_sintax_db_converter.R
    seqkit replace -p "\s.+" species_taxid.fasta | seqkit replace -p "(.+)" -r '{kv}' -k emu2sintax-header.tsv > sintax_db.fasta
fi

