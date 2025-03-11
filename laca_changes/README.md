# Changes to LACA snakemake files

Any edits I've made to `laca` for running my pipeline.

## clust.smk changes

File location     laca/laca/workflow/rules/clust.smk

During the `meshclust` step, `meshclust` will attempt to estimate a clustering threshold. This threshold has very strict limits and cannot be greater than 0.99 or `meshclust` throws an error that stops the pipeline. I added some lines to the `clust.smk` files that check if `meshclust` failed due to the estimated threshold being greater than 0.99 and then re-runs it with the threshold set to the max of 0.99. The changes are all in the first if statement of the shell block.

I posted this issue and the changes on github so it may be updated at some point.

### Original meshclust rule

```python 
rule meshclust:
    input: "clust/meshclust/{barcode}_{c1}_{c2}.fasta"
    output: "clust/meshclust/{barcode}_{c1}_{c2}.tsv"
    params:
        prefix = "clust/meshclust/{barcode}_{c1}_{c2}",
        t = "-t " + str(config["meshclust"]["t"]) if config["meshclust"]["t"] else "",
        max_batch_size = -1 if int(config["meshclust"]["max_batch_size"]) <= 0 else int(config["meshclust"]["max_batch_size"]) * 2,
    singularity: "docker://yanhui09/identity:latest"
    log: "logs/clust/meshclust/{barcode}_{c1}_{c2}.log"
    benchmark: "benchmarks/mechclust/{barcode}_{c1}_{c2}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell: 
        """
        if [ {params.max_batch_size} -eq -1 ]; then
          meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1
        else
          # number of reads in fastqs
          nlines=$(cat {input} | wc -l)
          if [ $nlines -le {params.max_batch_size} ]; then
            meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1
          else
            # determine the minimum partion size (number of batches), ceiling division
            min_part_size=$(((nlines + {params.max_batch_size} - 1) / {params.max_batch_size}))
            # determine the number of lines per batch, ceiling division, the nearest multiples of 2
            nlines_per_batch=$(((nlines + min_part_size * 2 - 1) / (min_part_size * 2) * 2))
            # if batches folder not exist, mkdir & split
            if [ ! -d {params.prefix} ]; then
              mkdir -p {params.prefix}
              split -l $nlines_per_batch -a3 -d --additional-suffix='.fasta' {input} {params.prefix}/b >{log} 2>&1
            fi
            for fa in {params.prefix}/b*.fasta; do
              batch_id=$(basename $fa | cut -d'.' -f1)
              if [ -f {params.prefix}/$batch_id.tsv ]; then
                continue
              fi
              meshclust -d $fa -o {params.prefix}/$batch_id.tsv -c {threads} {params.t} >> {log} 2>&1
              # add batchid after cluster (the first column)
              sed -i "s/\t/$batch_id\t/" {params.prefix}/$batch_id.tsv
              rm -f $fa
            done
            cat {params.prefix}/b*.tsv > {output}
            rm -rf {params.prefix} 
          fi 
        fi
        """
```

### Modified meschlust rule

```python
rule meshclust:
    input: "clust/meshclust/{barcode}_{c1}_{c2}.fasta"
    output: "clust/meshclust/{barcode}_{c1}_{c2}.tsv"
    params:
        prefix = "clust/meshclust/{barcode}_{c1}_{c2}",
        t = "-t " + str(config["meshclust"]["t"]) if config["meshclust"]["t"] else "",
        max_batch_size = -1 if int(config["meshclust"]["max_batch_size"]) <= 0 else int(config["meshclust"]["max_batch_size"]) * 2,
    singularity: "docker://yanhui09/identity:latest"
    log: "logs/clust/meshclust/{barcode}_{c1}_{c2}.log"
    benchmark: "benchmarks/mechclust/{barcode}_{c1}_{c2}.txt"
    threads: config["threads"]["large"]
    resources:
        mem = config["mem"]["large"],
        time = config["runtime"]["long"],
    shell: 
        """
        if [ {params.max_batch_size} -eq -1 ]; then
          meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1

          # Modification by AP
          # if meshclust fails due to threshold > 0.99 then re-run with threshold at 0.99
          if [! -f {output} ]; then
            a=$(awk ' /^Final threshold:/ {if ( $NF+0 > 0.99 ) print $NF } ' {log})
            if [ ! -z "$a" ]; then
              meshclust -d {input} -o {output} -c {threads} -t 0.99 > {log} 2>&1
            fi
          fi
        else
          # number of reads in fastqs
          nlines=$(cat {input} | wc -l)
          if [ $nlines -le {params.max_batch_size} ]; then
            meshclust -d {input} -o {output} -c {threads} {params.t} > {log} 2>&1
          else
            # determine the minimum partion size (number of batches), ceiling division
            min_part_size=$(((nlines + {params.max_batch_size} - 1) / {params.max_batch_size}))
            # determine the number of lines per batch, ceiling division, the nearest multiples of 2
            nlines_per_batch=$(((nlines + min_part_size * 2 - 1) / (min_part_size * 2) * 2))
            # if batches folder not exist, mkdir & split
            if [ ! -d {params.prefix} ]; then
              mkdir -p {params.prefix}
              split -l $nlines_per_batch -a3 -d --additional-suffix='.fasta' {input} {params.prefix}/b >{log} 2>&1
            fi
            for fa in {params.prefix}/b*.fasta; do
              batch_id=$(basename $fa | cut -d'.' -f1)
              if [ -f {params.prefix}/$batch_id.tsv ]; then
                continue
              fi
              meshclust -d $fa -o {params.prefix}/$batch_id.tsv -c {threads} {params.t} >> {log} 2>&1
              # add batchid after cluster (the first column)
              sed -i "s/\t/$batch_id\t/" {params.prefix}/$batch_id.tsv
              rm -f $fa
            done
            cat {params.prefix}/b*.tsv > {output}
            rm -rf {params.prefix} 
          fi 
        fi
        """
```


## quant.smk changes

File location     laca/laca/workflow/rules/quant.smk

My pipeline requires a file that lists every sequence ID and its OTU assignment, which `laca` doesn't generate by default. I have added two line to the `quant.smk` file that saves the necessaty table, since `laca` already generates it, one line to add the new output and one line to save a file in that output.

### Original rule

```python
# abudance matrix by read id
rule matrix_seqid: 
    input:
        rules.drep_consensus.output.tsv,
        rules.combine_cls.output,
        rules.drep_consensus.output.rep,
        fqs = lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_qced_barcodes(wc)),
    output: "quant/matrix_seqid.tsv"
    run:
        import pandas as pd
        # record the representative order from fasta file
        rep_order = []
        with open(input[2], "r") as rep:
            for line in rep:
                if line.startswith(">"):
                    rep_order.append(line.strip().split(">")[1])
        
        # OTU <- derep_cls -> cand_cls
        derep_cls = pd.read_csv(input[0], sep="\t", header=None)
        derep_cls.columns = ["rep_cls","cls"]
        # map to the 'rep_order' list as OTU, starting from 1
        derep_cls["OTU"] = derep_cls["rep_cls"].map(dict(zip(rep_order, range(1, len(rep_order)+1))))

        cand_cls = pd.read_csv(input[1], sep="\t", header=None)
        cand_cls.columns = ["seqid", "cls"]
        # merge on cls, left join, dummy files not included in {consensus}/{consensus}.fna 
        cls = derep_cls.merge(cand_cls, on="cls", how="left")
        # take header from fastq files, ^@, as a list
        for fq in input.fqs:
            barcode = os.path.basename(fq).split(".")[0]
            seqid = []
            with open(fq, "r") as f:
                for line in f:
                    if line.startswith("@"):
                        seqid.append(line.rstrip().split(" ")[0][1:])
            # if cls["seqid"] is in seqid, then 1, else 0
            cls[barcode] = cls["seqid"].isin(seqid).astype(int)
        
        # exclude rep_cls, cls, seqid
        # group by OTU, sum by barcode
        cls = cls.drop(["rep_cls", "cls", "seqid"], axis=1).groupby("OTU").sum()
        # append 'OTU_' after sorting
        cls.index = ["OTU_" + str(i) for i in cls.index.sort_values()]
        # rename OTU as "# OTU ID" for qiime2 import
        cls.index.name = "#OTU ID"
        # write to output
        cls.to_csv(output[0], sep="\t", header=True, index=True)
```

### Modified rule

```python
rule matrix_seqid: 
    input:
        rules.drep_consensus.output.tsv,
        rules.combine_cls.output,
        rules.drep_consensus.output.rep,
        fqs = lambda wc: expand("qc/qfilt/{barcode}.fastq", barcode=get_qced_barcodes(wc)),
    output: 
        "quant/matrix_seqid.tsv",
        "quant/seqID_to_otu.tsv"
    run:
        import pandas as pd
        # record the representative order from fasta file
        rep_order = []
        with open(input[2], "r") as rep:
            for line in rep:
                if line.startswith(">"):
                    rep_order.append(line.strip().split(">")[1])
        
        # OTU <- derep_cls -> cand_cls
        derep_cls = pd.read_csv(input[0], sep="\t", header=None)
        derep_cls.columns = ["rep_cls","cls"]
        # map to the 'rep_order' list as OTU, starting from 1
        derep_cls["OTU"] = derep_cls["rep_cls"].map(dict(zip(rep_order, range(1, len(rep_order)+1))))

        cand_cls = pd.read_csv(input[1], sep="\t", header=None)
        cand_cls.columns = ["seqid", "cls"]
        # merge on cls, left join, dummy files not included in {consensus}/{consensus}.fna 
        cls = derep_cls.merge(cand_cls, on="cls", how="left")

        # Modification by AP
        # save the seqID to OTU table to output
        cls.to_csv(output[1], sep="\t", header=True, index=True)

        # take header from fastq files, ^@, as a list
        for fq in input.fqs:
            barcode = os.path.basename(fq).split(".")[0]
            seqid = []
            with open(fq, "r") as f:
                for line in f:
                    if line.startswith("@"):
                        seqid.append(line.rstrip().split(" ")[0][1:])
            # if cls["seqid"] is in seqid, then 1, else 0
            cls[barcode] = cls["seqid"].isin(seqid).astype(int)
        
        # exclude rep_cls, cls, seqid
        # group by OTU, sum by barcode
        cls = cls.drop(["rep_cls", "cls", "seqid"], axis=1).groupby("OTU").sum()
        # append 'OTU_' after sorting
        cls.index = ["OTU_" + str(i) for i in cls.index.sort_values()]
        # rename OTU as "# OTU ID" for qiime2 import
        cls.index.name = "#OTU ID"
        # write to output
        cls.to_csv(output[0], sep="\t", header=True, index=True)
```