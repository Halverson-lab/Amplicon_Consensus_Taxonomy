#!/usr/bin/env python3

import sys
import pysam
import pandas as pd


#Take in two arguments, the read file and write file
#read sam file and make empty list
samfile = pysam.AlignmentFile(sys.argv[1])
d = []

#for each alignment, get the summary stats and save them to the list
for alignment in samfile.fetch():
    if alignment.reference_name:
        cigar_stats = alignment.get_cigar_stats()[0]
        n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    
        #query length
        q_length = int(alignment.query_length)
        if q_length==0:
            q_length = alignment.infer_query_length()

        #check that tag exists and save if it does
        if alignment.has_tag("AS"):
            AS_tag = alignment.get_tag("AS")
        else: 
            AS_tag = "NA"
        
        if alignment.has_tag("NM"):
            NM_tag = alignment.get_tag("NM")
        else: 
            NM_tag = "NA"

        d.append(
            {
                'QNAME': alignment.query_name,
                'QLENGTH': q_length,
                'RNAME': alignment.reference_name,
                'RLENGTH': alignment.reference_length,
                'AS': AS_tag,
                'NM': NM_tag,
                'CMATCH': cigar_stats[0],
                'CINS': cigar_stats[1],
                'CDEL': cigar_stats[2],
                'SOFT_CLIP': cigar_stats[4],
                'N_MISMATCH': n_mismatch
            }
        )

#Make the list into a pandas dataframe and then save as a csv
df = pd.DataFrame(d)
df.to_csv(sys.argv[2], index = False)
