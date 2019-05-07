#!/usr/bin/python

'''
Query sequences from NCBI

Call script: ./ncbi_organismgene_query.py --organism[NAME] --gene[NAME] > outfile.fasta

Last modified: 2019/01
'''

# import libraries
from __future__ import division, print_function
import sys, os, re
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

batch_size = 5000 # record download limit
out_format = 'fasta'

# define and parse command line parameters
def parse_args():
  parser = argparse.ArgumentParser(description='Download sequences from NCBI nucleotide DB.')
  parser.add_argument('--organism', help='Organism name or taxidNNN')
  parser.add_argument('--gene', nargs='+', help='Gene name(s)')
  
  args = parser.parse_args()
  return args

def construct_query(args):
  queries = []
  if args.organism:
    queries.append(args.organism + '[Organism]')
  if args.gene:
    gene_query = '(' + ' OR '.join([gn+'[Gene Name]' for gn in args.gene]) + ')'
    queries.append(gene_query)
    
  query = ' AND '.join(queries)
  return query

def main(query, args):
  print('Using query: ' + query, file = sys.stderr)
  
  # specify user and query search
  Entrez.email = 'sereina.rutschmann@gmail.com'
  search_handle = Entrez.esearch(db='nucleotide', term=query, retmax=batch_size, usehistory='y')
  search_result = Entrez.read(search_handle)
  search_handle.close()

  # remember result details
  record_count = int(search_result['Count'])
  web_env = search_result['WebEnv']
  query_key = search_result['QueryKey']

  # download records
  print('Found %d records.' % record_count, file=sys.stderr)
  for start_index in range(0, record_count, batch_size):
    end_index = min(record_count, start_index+batch_size)
    print('Going to download record %i to %i' % (start_index+1, end_index), file=sys.stderr)
    fetch_handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text',
                                 retstart=start_index, retmax=batch_size,
                                 webenv=web_env, query_key=query_key)
    records = SeqIO.parse(fetch_handle, 'gb')
    for rec in records:
      # extract elements from GenBank record
      org_name = rec.annotations['organism']
      accession = rec.annotations['accessions'][0]
      seq_id = '_'.join([org_name, args.organism, accession])
      seq_rec = SeqRecord(rec.seq, id = seq_id, description='')
      
      # write output files
      sys.stdout.write(seq_rec.format(out_format))

if __name__ == '__main__':
  args = parse_args()
  query = construct_query(args)
  main(query, args)
