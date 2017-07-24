#Analysis from blasted sequences
#2017/7/19
#Clara Liu

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from test import *
import os
from shutil import copyfile


INCOMPLETE_CUTOFF = 0

def incomplete(seq):
  l = 0
  r = len(seq)-1
  while (seq[l]=='-' or seq[r]=='-'):
    if seq[l] == '-':
      l += 1
    else:
      r -= 1

  return l+len(seq)-1-r


class BlastedSeq:
  def __init__(self, seq, id, genome_name):
    self.seq = seq.upper()
    self.id = id
    self.gn = genome_name
    self.complete = incomplete(seq)<=INCOMPLETE_CUTOFF


def txt_to_seq(txt):
  file = open(txt)
  str_file = file.read()
  newline = str_file.find('\r\n')
  id = str_file[str_file.find('>')+1:newline]
  str_seq = ''
  newline += 2
  while str_file.find('\r\n',newline)>=0:
    r = str_file.find('\r\n',newline)
    str_seq += str_file[newline:r]
    newline = r+2

  #seq = SeqIO.read(str_seq,"fasta",IUPAC.unambiguous_dna)
  seq = Bio.Seq.Seq(str_seq,IUPAC.unambiguous_dna)

  genome_name = txt[txt.rfind('/')+1:txt.rfind('.fasta')]

  bSeq = BlastedSeq(seq,id,genome_name)

  return bSeq


def consensus_seq(bseq_names,gene,ref_len):
  cwd = os.getcwd()
  os.chdir('Mycobacterium_tuberculosis/output_seq/%s/' %gene)
  seqs = []
  for txt in bseq_names:
    bSeq = txt_to_seq(txt)
    if bSeq.complete and bSeq.id!='No alignment':
      seqs += [bSeq]
  os.chdir(cwd)

  print len(seqs)

  #print len(seqs)
  #print '' %counter

  counter = 0
  for bSeq in seqs:
    if len(bSeq.seq) == ref_len:
      counter += 1
      #print '%s has a different size of gene' %bSeq.gn

  print counter
  '''
  result = []
  for i in range(ref_len):
    atgc = [0,0,0,0]
    #print '%d/%d' %(i,length)
    for j in range(len(seqs)):
      if (seqs[j].seq)[i] == 'A':
        atgc[0] += 1
      elif (seqs[j].seq)[i] == 'T':
        atgc[1] += 1
      elif (seqs[j].seq)[i] == 'G':
        atgc[2] += 1
      elif (seqs[j].seq)[i] == 'C':
        atgc[3] += 1
    result += [atgc]

  consensus = ''
  for i in range(ref_len):
    if max(result[i]) == result[i][0]:
      consensus += 'A'
    elif max(result[i]) == result[i][1]:
      consensus += 'T'
    elif max(result[i]) == result[i][2]:
      consensus += 'G'
    elif max(result[i]) == result[i][3]:
      consensus += 'C'
  '''

  consensus=''
  return consensus

#two seq has to have the same length, which is divisble by 3
def count_diff(seq,ref_seq,txt):
  global identical
  seq = ''.join(str(seq).split('-'))
  prot = Bio.Seq.Seq(seq,IUPAC.unambiguous_dna).translate(to_stop=True)
  ref_prot = Bio.Seq.Seq(str(ref_seq),IUPAC.unambiguous_dna).translate(to_stop=True)
  pointer = 0
  silent = 0
  missense = 0
  aa = 0
  for pointer in range(len(seq)/3-1):
    nuc0 = (seq[pointer*3] == ref_seq[pointer*3])
    nuc1 = (seq[pointer*3+1] == ref_seq[pointer*3+1])
    nuc2 = (seq[pointer*3+2] == ref_seq[pointer*3+2])
    if pointer >= len(prot):
      continue
    if prot[pointer] == ref_prot[pointer]:
      silent += [nuc0,nuc1,nuc2].count(False)
    else:
      missense += [nuc0,nuc1,nuc2].count(False)
      aa += 1

  if seq == ref_seq:
    identical = 1
  elif len(prot)>len(ref_prot):
    identical = -1
  else:
    identical = 0

  #if len(prot)<len(ref_prot):
    #print 'Early stop codon: %s' %txt[:-6]

  return [silent,missense,aa,identical]

def count_diff2(seq,ref_seq):
  diff = 0
  for k in range(len(ref_seq)):
    if seq[k] != ref_seq[k]:
      diff += 1
  return diff

#os.chdir('Mycobacterium_tuberculosis/output_seq/Rv0001/')
#txt1 = '1423443.3.fasta'
#print txt_to_seq(txt1).seq

def run_consensus(gene_names):
  for gene in gene_names:
    f = file_names('Mycobacterium_tuberculosis/output_seq/%s/' %gene)
    print ('-- Gene: %s --' %gene)

    ref_gene = SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %gene),"fasta",IUPAC.unambiguous_dna).seq.upper()
    cons = consensus_seq(f,gene,len(ref_gene))
    for i in range(len(cons)):
      if cons[i]!=ref_gene[i]:
        print '%d: %s -> %s' %(i,cons[i],ref_gene[i])
    print '\n'

"""
Pick out the blasted sequences that meet the following critiria:
1. Consist of only 'ATGC' or '-' (no 'N' allowed)
2. Fully complete (allow deletion in the middle)
3. Same length as the reference gene (allow only deletion)
   Sequences with insertion will be picked out to folder - Insertion
"""

def only_atgc(str):
  a = str.count('A')+str.count('a')
  t = str.count('T')+str.count('t')
  g = str.count('G')+str.count('g')
  c = str.count('C')+str.count('c')
  return len(str)==a+t+g+c

def qc(bseq_names,ref_gene):
  #cwd = os.getcwd()
  output_seq_dir = 'Mycobacterium_tuberculosis/output_seq/%s/' %ref_gene
  ref_gene_seq = str(SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %ref_gene),"fasta",IUPAC.unambiguous_dna).seq.upper())
  #makemydir('Insertion/')
  #os.chdir(output_seq_dir)
  qc_dir = 'Mycobacterium_tuberculosis/QC_seq/%s/' %ref_gene
  makemydir(qc_dir)
  for txt in bseq_names:
    bSeq = txt_to_seq(output_seq_dir+txt)
    if only_atgc(bSeq.seq) and bSeq.id!='No alignment':
      if bSeq.complete:
        if len(bSeq.seq)==len(ref_gene_seq):
          copyfile('%s%s.fasta' %(output_seq_dir,bSeq.gn), '%s%s.fasta' %(qc_dir,bSeq.gn))
        elif len(bSeq.seq)>len(ref_gene_seq):
          makemydir('Mycobacterium_tuberculosis/Insertion/%s/' %ref_gene)
          copyfile('%s%s.fasta' %(output_seq_dir,bSeq.gn), 'Mycobacterium_tuberculosis/Insertion/%s/%s.fasta' %(ref_gene,bSeq.gn))

  #os.chdir(cwd)

def run_count_diff(bseq_names,ref_gene):
  ref_gene_seq = str(SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %ref_gene),"fasta",IUPAC.unambiguous_dna).seq.upper())
  seq_dir = 'Mycobacterium_tuberculosis/QC_seq/%s/' %ref_gene
  silent = 0
  missense = 0
  aa = 0
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    [s,m,a,i] = count_diff(bseq.seq,ref_gene_seq,txt)
    silent += s
    missense += m
    aa += a
  return [len(bseq_names),len(ref_gene_seq),silent,missense,aa]

def run_count_diff2(bseq_names,ref_gene):
  ref_gene_seq = str(SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %ref_gene),"fasta",IUPAC.unambiguous_dna).seq.upper())
  seq_dir = 'Mycobacterium_tuberculosis/QC_seq/%s/' %ref_gene

  dna = 0
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    d = count_diff2(bseq.seq,ref_gene_seq)
    dna += d
  return dna


g = file_names('Mycobacterium_tuberculosis/ref_genes')
genes = []
for gene in g:
  genes += [gene[:-6]]
with open('Mycobacterium_tuberculosis/count.txt','w') as txt:
  txt.write('Counts\nGene Name | Gene Length | Silent | Missense | AA\n')
for gene in genes:
  bseq_names = file_names('Mycobacterium_tuberculosis/QC_seq/%s/' %gene)
  [n,l,s,m,a] = run_count_diff(bseq_names,gene)
  #dna = run_count_diff2(bseq_names,gene)
  with open('Mycobacterium_tuberculosis/count.txt','a') as txt:
    txt.write('%s %d %d %d %d %d\n' %(gene,n,l,s,m,a))
  print '%s: n= %d len= %d silent= %d missense= %d aa= %d' %(gene,n,l,s,m,a)
 # print '%s %d' %(gene,dna)