#! /usr/bin/env python
from __future__ import division, print_function
import sys
from Bio import AlignIO
from Bio.Alphabet.IUPAC import ambiguous_dna

# translate between degenerate DNA and Structure codes
base2geno = {
    'A': (0, 0),
    'C': (1, 1),
    'G': (2, 2),
    'T': (3, 3),
    'R': (0, 2),
    'Y': (1, 3),
    'S': (1, 2),
    'W': (0, 3),
    'K': (2, 3),
    'M': (0, 1),
    'N': (-9, -9),
    '?': (-9, -9),
    '-': (-9, -9)
}

if len(sys.argv) < 2:
  print("usage: %s alignment.nex" % sys.argv[0], file=sys.stderr)
  sys.exit(1)
alignment_fn = sys.argv[1]

# read alignment
# input can be stdin or command line param
align = AlignIO.read(alignment_fn, 'nexus')
align_len = align.get_alignment_length()

# get samples and their population from alignment
pop_cnt = 0
populations = {}
samples = {}
for rec in align:
  sample_id, pop = rec.id.rsplit('_',1)
  if pop not in populations:
    pop_cnt += 1
    populations[pop] = pop_cnt
  samples[rec.id] = populations[pop]
genotypes = {s: [] for s in samples}

# scan alignment for polymorphic loci and translate bases to numbers
for pos in range(align_len):
  polymorphic = False
  bases = set()
  loc = {rec.id:'-' for rec in align}
  
  for rec in align:
    base = rec.seq[pos].upper()
    bases.update(base)
    loc[rec.id] = base2geno[base]

  polymorphic = (len(bases - set('-N?')) > 1)
  if polymorphic:
    for s in loc:
      genotypes[s].append(loc[s])

# output STRUCTURE/STRUCTURAMA format
for s in sorted(genotypes.keys()):
  sys.stdout.write("%s\t%s" % (s, samples[s]))
  for a,b in genotypes[s]:
    # for STRUCTURAMA (diploid)
    # sys.stdout.write("\t(%i,%i)" % (a,b))
    # for STRUCTURE
    sys.stdout.write("\t%i\t%i" % (a,b))
  sys.stdout.write('\n')

print("Have a nice day.\n", file=sys.stderr)
