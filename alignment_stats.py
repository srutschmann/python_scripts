#!/usr/bin/env python

'''
script to calculate nucleotide diversity and segregating sites of sequence alignments

input: folder with FASTA files
note for larger data sets: tajimas'd and wattersons theta can be excluded

call script: python alignment_stats.py <folder> > alignment_stats.csv
'''

from __future__ import print_function
import dendropy
from dendropy.calculate import popgenstat
from glob import glob
import os
import sys

# define folder with input files
work_dir = sys.argv[1]
out_sep = ','

# print parameters
print(out_sep.join(['locus','seg_sites', 'nuc_div', 'taj_d', 'wat_theta']))

# calculate segregating_sites, nucleotide diversity, tajimas_d, watterson_theta
for filename in glob(os.path.join(work_dir, '*.fasta')):
    locus = os.path.basename(filename).split('.') [0]
    char_mat = dendropy.DnaCharacterMatrix.get_from_path(filename, 'fasta')

    nuc_div   = 0.0
    taj_d     = 0.0
    wat_theta = 0.0
    seg_sites = popgenstat.num_segregating_sites(char_mat)
    # calculate nucleotide diversity, tajimas'd, wattersons theta
    # (only makes sense if there are segregating sites)
    if seg_sites > 0:
        nuc_div   = popgenstat.nucleotide_diversity(char_mat)
        taj_d     = popgenstat.tajimas_d(char_mat)
        wat_theta = popgenstat.wattersons_theta(char_mat)

    print(out_sep.join([str(locus), str(seg_sites), str(nuc_div), str(taj_d), str(wat_theta)]))
