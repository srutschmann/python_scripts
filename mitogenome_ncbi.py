#!/usr/bin/python

# script to fetch annotated sequences within mitochondrial genomes
# !!!always keep it together with corresponding mitogenome_translation.csv file!!!

import sys
from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# get gene name translations
aliases = {}
f = open('translations.csv', 'rU')
for line in f:
    splitarr = line.strip().split(',')
    aliases[splitarr[0]] = splitarr[1]
f.close()

# ask user for specific options
organism = raw_input("For which taxon would you like to query? [Ephemeroptera]\n")
if len(organism)==0:
    organism = 'Ephemeroptera'
export_CDS = True if raw_input("CDS export? [y/n] ") == 'y' else False
if export_CDS:
    export_CDS_prot = True if raw_input("  export CDS as protein sequence? [y/n] ") == 'y' else False
export_tRNA = True if raw_input("tRNA export? [y/n] ") == 'y' else False
export_rRNA = True if raw_input("rRNA export? [y/n] ") == 'y' else False
export = {'CDS': export_CDS, 'tRNA': export_tRNA, 'rRNA': export_rRNA}
species = {}
records = []
genes = {}

# NCBI Entrez guidelines require that you provide an email address with each query
Entrez.email = "your.name@your-email.com"
if not Entrez.email:
    print "you must add your email address"
    sys.exit(2)

# search for mitochondrial sequences
search_term = "%s[Orgn] AND (mitochondrion[Titl] OR mitochondrial[Titl]) AND genome[Titl]" % organism
print "\nQuerying NCBI nucleotide database for\n\t'%s'" % search_term
handle = Entrez.esearch(db='nucleotide', term=search_term, usehistory='y')  #retmax=200)
search_results = Entrez.read(handle, validate=False)
handle.close()
count = int(search_results['Count'])
webenv = search_results['WebEnv']
query_key = search_results['QueryKey']
print "%i entries found." % count

# fetch actual sequences from GenBank
print "\nFetching %i records from GenBank..." % count
handle = Entrez.efetch(db="nucleotide", rettype="gb", webenv=webenv, query_key=query_key)
for rec in SeqIO.parse(handle, 'gb'):
    org_name = rec.annotations['organism']
    if org_name in species:
        species[org_name].append(rec)
    else:
        species[org_name] = [rec]
handle.close()
# select unique record for each species ("NC_*" records prefered)
for org in sorted(species, key=str.lower):
    if any(r.id.startswith('NC_') for r in species[org]):
        records += [r for r in species[org] if r.id.startswith('NC_')]
    else:
        records.append(species[org][0])
# collect genes from species records
for rec in records:
    accession = rec.annotations['accessions'][0]
    species = rec.annotations['organism']
    for feat in rec.features:
        if feat.type in export and export[feat.type]:
            gene_name = feat.qualifiers['product'][0] if 'product' in feat.qualifiers else feat.qualifiers['note'][0] #feat.qualifiers['gene'][0]
            gene_key = aliases[gene_name.lower()] if gene_name.lower() in aliases else gene_name
            # for CDS, ckeck wether protein or nucleotide sequence should be output
            if feat.type == 'CDS' and export_CDS_prot:
                # fasta_id = ("%s %s" % (species, organism)).replace(' ', '_')
                # print gene name in fasta header
                fasta_id = ("%s %s %s" % (species, organism, feat.qualifiers['protein_id'][0])).replace(' ', '_')
                seq = Seq(feat.qualifiers['translation'][0], alphabet=generic_protein)
            else:
                fasta_id = ("%s %s %s" % (species, organism, accession)).replace(' ', '_')
                seq = feat.extract(rec.seq)
            outseq = SeqRecord(seq, id=fasta_id, description=gene_key)
            if gene_key not in genes:
                genes[gene_key] = []                
            genes[gene_key].append(outseq)

# split duplicated genes in two instances
for g in genes.keys():
    if len(genes[g]) > len(set([s.id for s in genes[g]])):
        one = []
        two = []
        for seq in genes[g]:
            if any([s.id==seq.id for s in one]):
                two.append(seq)
            else:
                one.append(seq)
        genes.pop(g)
        genes['%s1' % g] = list(one)
        genes['%s2' % g] = list(two)

# write FASTA file for each gene
print "\nWriting FASTA output..."
for gene in genes:
    filename = "%s_%s.fasta" % (organism, gene.replace(' ', '-'))
    print "\t%s (%i sequences)" % (filename, len(genes[gene]))
    handle = open(filename, 'w')
    sorted_seqs = sorted(genes[gene], key=lambda s: s.id) # sort sequences by id
    SeqIO.write(sorted_seqs, handle, 'fasta')
    handle.close()

print '''
Have a nice day.
'''

