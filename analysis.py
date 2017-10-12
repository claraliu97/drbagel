#Analysis from blasted sequences
#2017/7/19
#Clara Liu

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from drbagel import *
import os
from shutil import copyfile
import parameters
from run import run
import xlsxwriter


INCOMPLETE_CUTOFF = 0

def incomplete(seq):
  if len(seq)==0:
    return 0
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
  sep = ''
  if txt[-4:] == '.txt':
    sep = '\r\n'
  elif txt[-6:] == '.fasta':
    sep = '\n'
  else:
    raise ValueError
  newline = str_file.find(sep)
  id = str_file[str_file.find('>')+1:newline]
  str_seq = ''
  newline += len(sep)
  while str_file.find(sep,newline)>=0:
    r = str_file.find(sep,newline)
    str_seq += str_file[newline:r]
    newline = r+len(sep)

  #seq = SeqIO.read(str_seq,"fasta",IUPAC.unambiguous_dna)
  seq = Bio.Seq.Seq(str_seq,IUPAC.unambiguous_dna)

  genome_name = txt[txt.rfind('/')+1:txt.rfind('.fasta')]


  bSeq = BlastedSeq(str_seq,id,genome_name)

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
def count_diff(seq_str,ref_seq,txt):
  global identical
  seq = ''.join(str(seq_str).split('-'))
  prot = Bio.Seq.Seq(seq,IUPAC.unambiguous_dna).translate(to_stop=True)
  ref_prot = Bio.Seq.Seq(str(ref_seq),IUPAC.unambiguous_dna).translate(to_stop=True)
  pointer = 0
  silent = 0
  missense = 0
  aa = 0
  log = {'ref':ref_seq,'match':'','query':seq_str}
  for pointer in range(len(seq_str)/3-1):
    nuc0 = (seq_str[pointer*3] == ref_seq[pointer*3])
    nuc1 = (seq_str[pointer*3+1] == ref_seq[pointer*3+1])
    nuc2 = (seq_str[pointer*3+2] == ref_seq[pointer*3+2])
    if pointer >= len(prot):
      continue
    if prot[pointer] == ref_prot[pointer]:
      silent += [nuc0,nuc1,nuc2].count(False)
      if nuc0:
        log['match'] += ' '
      else:
        log['match'] += 'x'
      if nuc1:
        log['match'] += ' '
      else:
        log['match'] += 'x'
      if nuc2:
        log['match'] += ' '
      else:
        log['match'] += 'x'
    else:
      missense += [nuc0,nuc1,nuc2].count(False)
      aa += 1
      if nuc0:
        log['match'] += ' '
      else:
        log['match'] += '*'
      if nuc1:
        log['match'] += ' '
      else:
        log['match'] += '*'
      if nuc2:
        log['match'] += ' '
      else:
        log['match'] += '*'

  if seq == ref_seq:
    identical = 1
  elif len(prot)>len(ref_prot):
    identical = -1
  else:
    identical = 0

  #if len(prot)<len(ref_prot):
    #print 'Early stop codon: %s' %txt[:-6]

  return [silent,missense,aa,identical,log]

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
  dash = str.count('-')
  return len(str)==a+t+g+c+dash

def qc(bseq_names,ref_gene):
  d = os.getcwd()
  output_seq_dir = '%s/output_seq/%s/DNA/' %(parameters.species,ref_gene)
  ref_gene_seq = str(SeqIO.read(('%s/ref_genes/%s.fasta' %(parameters.species,ref_gene)),"fasta",IUPAC.unambiguous_dna).seq.upper())
  #makemydir('Insertion/')
  #os.chdir(output_seq_dir)
  qc_dir = '%s/QC_seq/%s/' %(parameters.species,ref_gene)
  makemydir(qc_dir)
  for txt in bseq_names:
    bSeq = txt_to_seq(output_seq_dir+txt)
    if only_atgc(bSeq.seq) and bSeq.id!='No alignment':
      if bSeq.complete:
        if len(bSeq.seq)==len(ref_gene_seq):
          copyfile(output_seq_dir+txt, qc_dir+txt)
          #copyfile('%s%s.fasta' %(output_seq_dir,bSeq.gn), '%s%s.fasta' %(qc_dir,bSeq.gn))
        elif len(bSeq.seq)>len(ref_gene_seq):
          with open('%s/RESULT/%s_log.txt' %(parameters.species,ref_gene),'a') as txt:
            txt.write('**** Insertion %s\n' %bSeq.gn)
          makemydir('%s/Insertion/%s/' %(parameters.species,ref_gene))
          copyfile('%s%s.fasta' %(output_seq_dir,bSeq.gn), '%s/Insertion/%s/%s.fasta' %(parameters.species,ref_gene,bSeq.gn))
        else:
          with open('%s/RESULT/%s_log.txt' %(parameters.species,ref_gene),'a') as txt:
            txt.write('**** Incomplete %s\n' %bSeq.gn)

  #os.chdir(cwd)

def run_count_diff(bseq_names,ref_gene):
  ref_gene_seq = str(SeqIO.read(('%s/ref_genes/%s.fasta' %(parameters.species,ref_gene)),"fasta",IUPAC.unambiguous_dna).seq.upper())
  seq_dir = '%s/QC_seq/%s/' %(parameters.species,ref_gene)
  silent = 0
  missense = 0
  aa = 0
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    [s,m,a,i,log] = count_diff(bseq.seq,ref_gene_seq,txt)
    silent += s
    missense += m
    aa += a
    if s>0 or m>0:
      with open('%s/RESULT/%s_log.txt' %(parameters.species,ref_gene),'a') as logtxt:
        logtxt.write('%s\n>%s\n\n' %('#'*75,txt))
        n = 0
        while n<len(log['match']):
          logtxt.write('%s\n' %log['ref'][n:n+75])
          logtxt.write('%s\n' %log['match'][n:n+75])
          logtxt.write('%s\n' %log['query'][n:n+75])
          n += 75
        logtxt.write('\n')

  return [len(bseq_names),len(ref_gene_seq),silent,missense,aa]

def run_count_seq(bseq_names,ref_gene,seqs):
  ref_gene_seq = str(SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %ref_gene),"fasta",IUPAC.unambiguous_dna).seq.upper())
  seq_dir = 'Mycobacterium_tuberculosis/Insertion/%s/' %ref_gene
  count = [0]*len(seqs)
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    for n in range(len(seqs)):
      if seqs[n] in bseq.seq:
        count[n] += 1
  for n in range(len(seqs)):
    print '%s: %d/%d' %(seqs[n],count[n],len(bseq_names))

def run_count_diff2(bseq_names,ref_gene):
  ref_gene_seq = str(SeqIO.read(('Mycobacterium_tuberculosis/ref_genes/%s.fasta' %ref_gene),"fasta",IUPAC.unambiguous_dna).seq.upper())
  seq_dir = 'Mycobacterium_tuberculosis/QC_seq/%s/' %ref_gene

  dna = 0
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    d = count_diff2(bseq.seq,ref_gene_seq)
    dna += d
  return dna


######### Display the mutation at each AA position #########

aa = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z','-']

def count_by_aa(seq_str,ref_seq,txt,aa_dict_lst):
  seq = ''.join(str(seq_str).split('-'))
  prot = Bio.Seq.Seq(seq,IUPAC.unambiguous_dna).translate(to_stop=True)
  ref_prot = Bio.Seq.Seq(str(ref_seq),IUPAC.unambiguous_dna).translate(to_stop=True)
  alignments = Bio.pairwise2.align.globalxs(ref_prot, prot,-10,-10)
  ref_aligned = alignments[0][0]
  prot_aligned = alignments[0][1]
  if ref_aligned!=ref_prot:
    if len(ref_aligned)>len(ref_prot):
      print "WRONG: %s" %txt
      return aa_dict_lst
    else:
      print "SUPER WRONG: %s" %txt
      return aa_dict_lst

  for pointer in range(len(ref_aligned)):
    if prot_aligned[pointer] != ref_prot[pointer]:
      aa_dict_lst[pointer][prot_aligned[pointer]] += 1

  return aa_dict_lst

def run_count_by_aa(bseq_names,ref_gene):
  ref_gene_seq = str(SeqIO.read(('%s/ref_genes/%s.fasta' %(parameters.species,ref_gene)),"fasta",IUPAC.unambiguous_dna).seq.upper())
  ref_prot = Bio.Seq.Seq(str(ref_gene_seq),IUPAC.unambiguous_dna).translate(to_stop=True)
  seq_dir = '%s/QC_seq/%s/' %(parameters.species,ref_gene)
  aa_dict_lst = [{} for i in ref_prot]
  for i in range(len(ref_prot)):
    for j in aa:
      aa_dict_lst[i][j] = 0
  for i in range(len(ref_prot)):
    aa_dict_lst[i][ref_prot[i]] = -1
  for txt in bseq_names:
    bseq = txt_to_seq(seq_dir+txt)
    aa_dict_lst = count_by_aa(bseq.seq,ref_gene_seq,txt,aa_dict_lst)
  return (ref_prot,aa_dict_lst)

def has_mut(aa_dict):
  mut = 0
  for k in aa_dict:
    if aa_dict[k] != -1:
      mut += aa_dict[k]
  return mut>0


curdir = os.getcwd()
run()
os.chdir(curdir)
g = file_names('%s/ref_genes/' %parameters.species)
genes = []
for gene in g:
  genes += [gene[:-6]]
#os.chdir(os.getcwd()[:os.getcwd().rfind('/')])
if True:
  makemydir('%s/RESULT/' %parameters.species)
  with open('%s/RESULT/count.txt' %parameters.species,'a') as txt:
    txt.write('Counts\nGene Name | # of genomes counted | Gene Length | Silent | Missense | AA\n')
  for gene in genes:
    #rename('Mycobacterium_tuberculosis/output_seq/%s/' %gene,'.txt','.fasta')
    txt_names = file_names('%s/output_seq/%s/DNA/' %(parameters.species,gene))
    with open('%s/RESULT/%s_log.txt' %(parameters.species,gene),'w') as logtxt:
      logtxt.write('Print Out: %s\n(Each entry: ref-query, x=silent, *=missense)\n\n' %gene)
    qc(txt_names,gene)
    bseq_names = file_names('%s/QC_seq/%s/' %(parameters.species,gene))
    [n,l,s,m,a] = run_count_diff(bseq_names,gene)
    #dna = run_count_diff2(bseq_names,gene)
    with open('%s/RESULT/count.txt' %parameters.species,'a') as txt:
      txt.write('%s %d %d %d %d %d\n' %(gene,n,l,s,m,a))
    print '%s: n= %d len= %d silent= %d missense= %d aa= %d' %(gene,n,l,s,m,a)
   # print '%s %d' %(gene,dna)

for gene in genes:
  #rename('Mycobacterium_tuberculosis/output_seq/%s/' %gene,'.txt','.fasta')
  #txt_names = file_names('%s/output_seq/%s/DNA/' %(parameters.species,gene))
  #qc(txt_names,gene)
  print gene
  bseq_names = file_names('%s/QC_seq/%s/' %(parameters.species,gene))
  (ref_prot,aa_dict_lst) = run_count_by_aa(bseq_names,gene)
  #dna = run_count_diff2(bseq_names,gene)
  workbook = xlsxwriter.Workbook(parameters.species+'/RESULT/'+gene+'.xlsx')
  ws0 = workbook.add_worksheet('by_aa')
  col = 1
  for j in range(len(aa)):
    ws0.write(j+2,0,aa[j])
  for i in range(len(aa_dict_lst)):
    if has_mut(aa_dict_lst[i]):
      ws0.write(0,col,i+1)
      ws0.write(1,col,ref_prot[i])
      for j in range(len(aa)):
        if aa_dict_lst[i][aa[j]]>0:
          ws0.write(j+2,col,aa_dict_lst[i][aa[j]])
      col += 1
  """
  for i in range(len(aa_dict_lst)):
    ws0.write(1,i+1,ref_prot[i])
    ws0.write(0,i+1,i+1)
  for j in range(len(aa)):
    ws0.write(j+2,0,aa[j])
  for i in range(len(aa_dict_lst)):
    for j in range(len(aa)):
      ws0.write(j+2,i+1,aa_dict_lst[i][aa[j]])
  """
  workbook.close()


"""
gene = 'Rv0050'
bseq_names = file_names('Mycobacterium_tuberculosis/Insertion/%s/' %gene)
seqs = ['GCCGCCTCCGCCGCCGCCGTCGG','CGCCGCCGCCGTCGGAGG','GCCGCCGCCGTCGGAGGT','CCGCCGCCGTCGGAGGTACCAC','GCCGCCGCCGTCGGAGGTAC']
run_count_seq(bseq_names,gene,seqs)
"""