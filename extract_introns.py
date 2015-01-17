#! /usr/bin/env python3

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# script to split intron/exons from Sanger sequences
# specify working directory and path to csv file
# modify csv file
# call script: python3 extract_introns.py

# script needs genes_introns.csv file in tab format names in csv file have to be identical with .nex file names
# format of extract_introns.csv: ortholog,start_intron[inclusive],end_intron[exclusive], only exon: ortholog
# all nexus files have to be in the same folder and as simplified nexus format (e.g. NOT formated with mesquite)

# define working directory
align_dir = "/Users/exon_intron_split"

# import csv file with information about intron positions
intron_file = open("/Users/exon_intron_split/extract_introns.csv")
for line in intron_file:
    parts = line.strip().split(',')
    gene = parts[0]
    bounds = ['1'] + parts[1:]
    
    # check if input file exists
    align_fn = "%s/%s.nex" % (align_dir, gene)
    if not os.path.exists(align_fn):
        continue
    
    # open input file
    records = list(SeqIO.parse(open(align_fn),"nexus"))
    
    # store exons and introns in lists
    features ={"whole": [], "coding": [], "noncoding": [[] for x in range(len(records))]}
    exon_count = 0
    intron_count = 0
    exon = True
    
    # split sequence into exons and introns
    for i in range(len(bounds)):
        first = (i == 0)
        last = (i == len(bounds)-1)
        start = int(bounds[i])
        sub_seqs = []
        # cut each sequence in alignment separately
        for j in range(len(records)):
            rec = records[j]
            end = int(bounds[i+1]) if not last else len(rec)
            feat_seq = rec.seq[start-1:end-1]
            id,name,desc = rec.id, rec.name, rec.description
            # introns get set in lower case (for visualization - but not working)
#             if not exon:
#                feat_seq = feat_seq.lower()
            feat_rec = SeqRecord(feat_seq, id=id, name=name, description=desc)
            sub_seqs.append(feat_rec)
            # create coding and edited sequence
            if first:
                features["whole"].append(feat_rec)
                if exon:
                    # add feature to coding sequences
                    features["coding"].append(feat_rec)
            else:
                features["whole"][j] += feat_rec
                if exon:
                    # add feature to coding sequences
                    features["coding"][j] += feat_rec
            if not exon:
                # add feature to noncoding sequences
                features["noncoding"][j].append(feat_rec)
        
        # add alignment block to features
        if exon:
            exon_count += 1
            features["exon%i" % exon_count] = sub_seqs
        else:
            intron_count += 1
            features["intron%i" % intron_count] = sub_seqs
        exon = not exon
    
    # write output to files
    for feat in features:
        if feat != "noncoding":
            outfile = open("%s/%s_%s.fasta" % (align_dir, gene, feat), "w")
            SeqIO.write(features[feat], outfile, "fasta")
            outfile.close()
        else:
            if sum([1 if len(l)>0 else 0 for l in features[feat]]) == 0:
                continue
            outfile = open("%s/%s_%s.fasta" % (align_dir, gene, feat), "w")
            for reclist in features[feat]:
                if len(reclist) > 0:
                    rec = reclist[0]
                    for i in range(1,len(reclist)):
                        rec += reclist[i]
                    SeqIO.write(rec, outfile, "fasta")
            outfile.close()

print("\nHave a nice day.\n")
