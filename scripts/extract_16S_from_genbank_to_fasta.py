#!/usr/bin/env python

from Bio import SeqIO
import sys

#qualifier for 16S rRNA
qualifier = "RFAM:RF00177"

with open(sys.argv[2], 'w') as ofh:
    for seq_record in SeqIO.parse(sys.argv[1], 'genbank'):
            for seq_feature in seq_record.features:
                    if seq_feature.type=="source":
                        organism = str(seq_feature.qualifiers['organism'][0]).rsplit(" ", 1)[0]
                        taxon = str(seq_feature.qualifiers['db_xref'][0])

                    if seq_feature.type=="rRNA":
                        if seq_feature.qualifiers['db_xref'][0] == qualifier:
                            ofh.write(">%s %s %s\n%s\n" % (
                                seq_feature.qualifiers['locus_tag'][0],
                                seq_record.name,
                                taxon,
                                seq_feature.location.extract(seq_record).seq))
