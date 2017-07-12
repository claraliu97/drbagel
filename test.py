#Dr.Bagel Blast and Analysis Script
#Clara Liu
#2017/06/12

import Bio
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from parameters import *
import os

counter = 0
#return all the name of files in the given dir
def file_names(dir):
  f = []
  for (dirpath, dirnames, filenames) in os.walk(dir):
    f.extend(filenames)
    break
  return f

def folder_names(dir):
  f = []
  for (dirpath, dirnames, filenames) in os.walk(dir):
    f.extend(dirnames)
    break
  return f


def makemydir(dir):
  try:
    os.makedirs(dir)
  except OSError:
    pass

def create_db(input_dir,file_list):
  makemydir("db_%s/" %extension)
  #os.system("mkdir -p temp")
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    #make_db = "makeblastdb -in %s%s.%s -dbtype nucl -parse_seqids \

    make_db = "makeblastdb -in %s%s.%s -dbtype nucl \
    -out db_%s/%s/%s -title '%s'\n" %(input_dir,genome,extension,extension,genome,genome,genome)

    
    #with open("temp/test.sh", "w") as text_file:
      #windows
      #text_file.write("#!/bin/sh\n")
      #text_file.write(make_db)

    os.system(make_db)
    #os.system("sh temp/test.sh")
  #os.system("rm -r temp")


def create_blast(input_dir,file_list,ref_gene_name):
  global counter
  counter = 0
  print ref_gene_name

  outdir = "blast_%s/%s/" %(extension,ref_gene_name)
  #os.system("mkdir -p temp")
  makemydir(outdir)
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    if counter%73 ==0:
      print "blast: %d/100" %(counter/73)
    counter += 1

    make_blast = "blastn -query ref_genes/%s.fasta -task megablast -db db_%s/%s/%s \
    -outfmt 5 -out %s%s.xml -num_threads 4" %(ref_gene_name,extension,genome,genome,outdir,genome)
    
    #with open("temp/test.sh", "w") as text_file:
      #text_file.write(make_blast)

    #windows:
    os.system(make_blast)
    #macos:
    #os.system("sh temp/test.sh")
  #os.system("rm -r temp")

def dash_at_beginning(str):
  ls = 0
  le = 0
  rs = len(str)-1
  re = len(str)-1
  if (len(str)==0):
    return 0
  while (str[le]=='-' or str[rs]=='-'):
    if (str[le]=='-'):
      le += 1
    else:
      rs -= 1
  return le-ls+re-rs

def dash_at_beginning2(str):
  ls = 0
  le = 0
  rs = len(str)-1
  re = len(str)-1
  if (len(str)==0):
    return 0
  while (str[le]=='=' or str[rs]=='='):
    if (str[le]=='='):
      le += 1
    else:
      rs -= 1
  return le-ls+re-rs

def parse_single_blast(result_handle,ref_gene,genome_records,txt):

  global counter
  if counter%73 ==0:
    print "parse: %d/100" %(counter/73)
  counter += 1

  ref_gene = ref_gene.seq
  ref_protein = ref_gene.translate()

  blast_record = NCBIXML.read(result_handle)
  #find the best hit (filter the off-target hits)
  try:
    alignment = sorted(blast_record.alignments,cmp=lambda x,y: x.length-y.length)[-1]
    id_l = alignment.title.find(id_spliter)
    id_r = alignment.title.find(id_spliter,id_l+1)
    id = alignment.title[id_l+1:id_r]
    sbjct_gene = genome_records[id].seq
    #print id
    hsp = sorted(alignment.hsps,cmp=lambda x,y: x.identities-y.identities)[-1]

    #query = ref_gene.seq.upper()
    query = hsp.query.upper()
    len_add = hsp.query.count('-')+hsp.sbjct.count('-')

    len_upstream = hsp.query_start - 1
    len_downstream = len(ref_gene)-(hsp.query_end-hsp.query_start+1)-len_upstream
    #print (hsp.sbjct_start,hsp.sbjct_end)

    if hsp.sbjct_start<hsp.sbjct_end:
      sbjct_left_s = hsp.sbjct_start-len_upstream
      sbjct_left_e = hsp.sbjct_start
      sbjct_left = (sbjct_gene)[sbjct_left_s-1:sbjct_left_e-1].upper()

      sbjct_right_s = hsp.sbjct_end+1
      sbjct_right_e = hsp.sbjct_end+1+len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_s-1:sbjct_right_e-1].upper()

      #sbjct_from = hsp.sbjct_start-len_upstream
      #sbjct_to = sbjct_from+len(query)-1
      #sbjct = (sbjct_gene)[sbjct_from-1:sbjct_to].upper()
    else:
      sbjct_left_s = hsp.sbjct_end+1+len_upstream
      sbjct_left_e = hsp.sbjct_end+1
      sbjct_left = (sbjct_gene)[sbjct_left_e-1:sbjct_left_s-1].reverse_complement().upper()

      sbjct_right_s = hsp.sbjct_start
      sbjct_right_e = hsp.sbjct_start-len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_e-1:sbjct_right_s-1].reverse_complement().upper()

    if len(sbjct_left)<len_upstream:
      sbjct_left = '-'*(len_upstream-len(sbjct_left))+sbjct_left
    if len(sbjct_right)<len_downstream:
      sbjct_right = sbjct_right+'-'*(len_downstream-len(sbjct_right))

    sbjct = sbjct_left+hsp.sbjct.upper()+sbjct_right

      #sbjct_to = hsp.sbjct_start+len_upstream
      #sbjct_from = sbjct_to-len(query)+1
      #print "sbjct_to,sbjct_from: %d %d" %(sbjct_to,sbjct_from)
      #sbjct = (sbjct_gene)[sbjct_from-1:sbjct_to].reverse_complement().upper()


    query = (ref_gene[:hsp.query_start-1] + query + ref_gene[hsp.query_end:]).upper()

   # print (id,len(query),len(sbjct),counter)

    mismatch_g = []
    match_g = ''
    dash_count_g = dash_at_beginning(sbjct)
    for n in range(len(query)):
      if query[n]!=sbjct[n]:
        match_g += 'x'
        mismatch_g += [(n+1,query[n],sbjct[n])]
      else:
        match_g += ' '

    with open(txt, "a") as text_file:
      text_file.write('**** DNA Alignment ****\n')
      text_file.write('Sequence:\n%s\n' %alignment.title)
      text_file.write('Identities: %d/%d\n' %(match_g.count(' '),len(ref_gene)))
      #text_file.write('Positives: %d/%d\n' %(hsp.positives,len(hsp.query)))
      text_file.write('Mismatches:\n')
      print_parsed_blast([query,match_g,sbjct,mismatch_g],text_file)
      if dash_count_g>0.1*len(ref_gene):
        text_file.write('##### Invalid DNA > 10% #####\n')
        mismatch_g = []
        dash_count_g = 0

    sbjct_dna = ''.join(str(sbjct) .split('-'))

    protein = Seq(sbjct_dna, IUPAC.unambiguous_dna).translate()

    alignments = pairwise2.align.globalms(ref_protein, protein,2,0,-200,-100)
    align = pairwise2.format_alignment(*alignments[0])
    lst = align.split('\n')
    match_p = ''
    mismatch_p = []
    dash_count_p = dash_at_beginning(lst[2])
    for n in range(len(lst[0])):
      if lst[0][n] == lst[2][n]:
        match_p += ' '
      else:
        match_p += 'x'
        mismatch_p += [(n+1,lst[0][n],lst[2][n])]

    with open(txt, "a") as text_file:
      text_file.write('**** Protein Alignment ****\n')
      text_file.write('Identities: %d/%d\n' %(match_p.count(' '),len(ref_protein)))
      if dash_count_g>0.1*len(ref_gene):
        mismatch_p = []
        dash_count_p = 0
      else:
        print_parsed_blast([lst[0],match_p,lst[2],mismatch_p],text_file)

      text_file.write(' \n')

    return [len(mismatch_g),dash_count_g,len(mismatch_p),dash_count_p]

  except Bio.Data.CodonTable.TranslationError:
    open(txt,"a").write("Invalid DNA seq. Fail to translate\n \n")
    return [0,0,0,0]
  except IndexError:
    open(txt,"a").write("No alignment\n \n")
    return [0,0,0,0]


def count_single_invalid(result_handle,ref_gene,genome_records,txt):

  global counter
  if counter%73 ==0:
    print "parse: %d/100" %(counter/73)
  counter += 1

  ref_gene = ref_gene.seq
  ref_protein = ref_gene.translate()

  blast_record = NCBIXML.read(result_handle)
  #find the best hit (filter the off-target hits)
  try:
    alignment = sorted(blast_record.alignments,cmp=lambda x,y: x.length-y.length)[-1]
    id_l = alignment.title.find(id_spliter)
    id_r = alignment.title.find(id_spliter,id_l+1)
    id = alignment.title[id_l+1:id_r]
    sbjct_gene = genome_records[id].seq
    #print id
    hsp = sorted(alignment.hsps,cmp=lambda x,y: x.identities-y.identities)[-1]

    #query = ref_gene.seq.upper()
    query = hsp.query.upper()
    len_add = hsp.query.count('-')+hsp.sbjct.count('-')

    len_upstream = hsp.query_start - 1
    len_downstream = len(ref_gene)-(hsp.query_end-hsp.query_start+1)-len_upstream
    #print (hsp.sbjct_start,hsp.sbjct_end)

    if hsp.sbjct_start<hsp.sbjct_end:
      sbjct_left_s = hsp.sbjct_start-len_upstream
      sbjct_left_e = hsp.sbjct_start
      sbjct_left = (sbjct_gene)[sbjct_left_s-1:sbjct_left_e-1].upper()

      sbjct_right_s = hsp.sbjct_end+1
      sbjct_right_e = hsp.sbjct_end+1+len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_s-1:sbjct_right_e-1].upper()

      #sbjct_from = hsp.sbjct_start-len_upstream
      #sbjct_to = sbjct_from+len(query)-1
      #sbjct = (sbjct_gene)[sbjct_from-1:sbjct_to].upper()
    else:
      sbjct_left_s = hsp.sbjct_end+1+len_upstream
      sbjct_left_e = hsp.sbjct_end+1
      sbjct_left = (sbjct_gene)[sbjct_left_e-1:sbjct_left_s-1].reverse_complement().upper()

      sbjct_right_s = hsp.sbjct_start
      sbjct_right_e = hsp.sbjct_start-len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_e-1:sbjct_right_s-1].reverse_complement().upper()

    if len(sbjct_left)<len_upstream:
      sbjct_left = '-'*(len_upstream-len(sbjct_left))+sbjct_left
    if len(sbjct_right)<len_downstream:
      sbjct_right = sbjct_right+'-'*(len_downstream-len(sbjct_right))

    sbjct = sbjct_left+hsp.sbjct.upper()+sbjct_right

    query = (ref_gene[:hsp.query_start-1] + query + ref_gene[hsp.query_end:]).upper()


    dash_count_g = dash_at_beginning(sbjct)

    if dash_count_g>0.1*len(ref_gene):
      return 1
    else:
      return 0
  except Bio.Data.CodonTable.TranslationError:
    return 1
  except IndexError:
    return 1


def combine_parse_single_blast(result_handle,ref_gene,genome_records,txt):

  global counter
  if counter%73 ==0:
    print "parse: %d/100" %(counter/73)
  counter += 1

  ref_gene = ref_gene.seq
  ref_protein = ref_gene.translate()

  blast_record = NCBIXML.read(result_handle)
  #find the best hit (filter the off-target hits)
  try:
    alignment = sorted(blast_record.alignments,cmp=lambda x,y: x.length-y.length)[-1]
    id_l = alignment.title.find(id_spliter)
    id_r = alignment.title.find(id_spliter,id_l+1)
    id = alignment.title[id_l+1:id_r]
    sbjct_gene = genome_records[id].seq
    #print id
    hsp = sorted(alignment.hsps,cmp=lambda x,y: x.identities-y.identities)[-1]

    #query = ref_gene.seq.upper()
    query = hsp.query.upper()
    len_add = hsp.query.count('-')+hsp.sbjct.count('-')

    len_upstream = hsp.query_start - 1
    len_downstream = len(ref_gene)-(hsp.query_end-hsp.query_start+1)-len_upstream
    #print (hsp.sbjct_start,hsp.sbjct_end)

    if hsp.sbjct_start<hsp.sbjct_end:
      sbjct_left_s = hsp.sbjct_start-len_upstream
      sbjct_left_e = hsp.sbjct_start
      sbjct_left = (sbjct_gene)[sbjct_left_s-1:sbjct_left_e-1].upper()

      sbjct_right_s = hsp.sbjct_end+1
      sbjct_right_e = hsp.sbjct_end+1+len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_s-1:sbjct_right_e-1].upper()

      #sbjct_from = hsp.sbjct_start-len_upstream
      #sbjct_to = sbjct_from+len(query)-1
      #sbjct = (sbjct_gene)[sbjct_from-1:sbjct_to].upper()
    else:
      sbjct_left_s = hsp.sbjct_end+1+len_upstream
      sbjct_left_e = hsp.sbjct_end+1
      sbjct_left = (sbjct_gene)[sbjct_left_e-1:sbjct_left_s-1].reverse_complement().upper()

      sbjct_right_s = hsp.sbjct_start
      sbjct_right_e = hsp.sbjct_start-len_downstream
      sbjct_right = (sbjct_gene)[sbjct_right_e-1:sbjct_right_s-1].reverse_complement().upper()

    if len(sbjct_left)<len_upstream:
      sbjct_left = '='*(len_upstream-len(sbjct_left))+sbjct_left
    if len(sbjct_right)<len_downstream:
      sbjct_right = sbjct_right+'='*(len_downstream-len(sbjct_right))

    sbjct = sbjct_left+hsp.sbjct.upper()+sbjct_right

      #sbjct_to = hsp.sbjct_start+len_upstream
      #sbjct_from = sbjct_to-len(query)+1
      #print "sbjct_to,sbjct_from: %d %d" %(sbjct_to,sbjct_from)
      #sbjct = (sbjct_gene)[sbjct_from-1:sbjct_to].reverse_complement().upper()


    query = (ref_gene[:hsp.query_start-1] + query + ref_gene[hsp.query_end:]).upper()

   # print (id,len(query),len(sbjct),counter)

    dash_count_g = dash_at_beginning2(sbjct)
    g = ''
    for n in range(len(query)):
      if query[n]!=sbjct[n]:
        g += sbjct[n]
      else:
        g += ' '

    if dash_count_g>0.1*len(ref_gene):
      g = ''

    sbjct_dna = ''.join(str(sbjct) .split('='))
    sbjct_dna = ''.join(str(sbjct) .split('-'))

    protein = Seq(sbjct_dna, IUPAC.unambiguous_dna).translate()
    protein = protein[:protein.find('*')+1]

    alignments = pairwise2.align.globalms(ref_protein, protein,2,0,-200,-100)
    align = pairwise2.format_alignment(*alignments[0])
    lst = align.split('\n')
    p = ''
    dash_count_p = dash_at_beginning(lst[2])
    for n in range(len(lst[0])):
      if lst[0][n] == lst[2][n]:
        p += ' '
      else:
        p += lst[2][n]

    if dash_count_g>0.1*len(ref_gene):
      p = ''

    return [g,query,p,lst[0]]

  except Bio.Data.CodonTable.TranslationError:
    return ['']*4
  except IndexError:
    return ['']*4




def print_parsed_blast(lst,text_file):
  query = lst[0]
  match = lst[1]
  sbjct = lst[2]
  mismatch = lst[3]

  #for (n,q,s) in mismatch:
    #text_file.write('%d: %s->%s\n' %(n,q,s))

  if len(mismatch)>0:
    n = 0
    while (n<len(query)):
      text_file.write('%s\n' %query[n:n+75])
      text_file.write('%s\n' %match[n:n+75])
      text_file.write('%s\n' %sbjct[n:n+75])
      n += 75


def parse_blast(file_list,ref_gene_name):
  #E_VALUE_THRESH = 0.001
  global counter
  
  ref_gene = SeqIO.read('ref_genes/%s.fasta' %ref_gene_name,"fasta",IUPAC.unambiguous_dna)

  output_name = "blast_%s/%s.txt" %(extension,ref_gene_name)
  with open(output_name, "w") as f:
    f.write("%s %s\n" %(species, ref_gene_name))

  mut_sum = [0,0,0,0]
  #seq_cum = [[0]*len(ref_gene)]*5 #A,T,G,C

  counter = 0

  #try:
  file_name = "no genome input"
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]
    genome_records = SeqIO.index('%s/%s' %(extension,file_name),"fasta",IUPAC.unambiguous_dna)
    #print records["JLSA01000036"]
    result_handle = open("blast_%s/%s/%s.xml" %(extension,ref_gene_name,genome))

    [g,dg,p,dp] = parse_single_blast(result_handle,ref_gene,genome_records,output_name)
    mut_sum[0] += g
    mut_sum[1] += dg
    mut_sum[2] += p
    mut_sum[3] += dp

  return mut_sum

  #except:
    #print "fail at %s" %file_name

def count_invalid(file_list,ref_gene_name):
  #E_VALUE_THRESH = 0.001
  global counter
  
  ref_gene = SeqIO.read('ref_genes/%s.fasta' %ref_gene_name,"fasta",IUPAC.unambiguous_dna)

  output_name = "blast_%s/%s.txt" %(extension,ref_gene_name)
  with open(output_name, "w") as f:
    f.write("%s %s\n" %(species, ref_gene_name))

  invalid = 0
  #seq_cum = [[0]*len(ref_gene)]*5 #A,T,G,C

  counter = 0

  file_name = "no genome input"
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]
    genome_records = SeqIO.index('%s/%s' %(extension,file_name),"fasta",IUPAC.unambiguous_dna)
    result_handle = open("blast_%s/%s/%s.xml" %(extension,ref_gene_name,genome))

    invalid += count_single_invalid(result_handle,ref_gene,genome_records,output_name)

  return invalid

def combine_parse_blast(file_list,ref_gene_name):
  global counter
  
  ref_gene = SeqIO.read('ref_genes/%s.fasta' %ref_gene_name,"fasta",IUPAC.unambiguous_dna)
  ref_protein = ref_gene.seq.translate()

  len_g = len(ref_gene.seq)
  len_p = len_g/3

  output_name_g = "blast_%s/combine_g_%s.txt" %(extension,ref_gene_name)
  output_name_p = "blast_%s/combine_p_%s.txt" %(extension,ref_gene_name)
  with open(output_name_g, "w") as g:
    g.write("---- %s %s ----\n" %(species, ref_gene_name))
    g.write(str(ref_gene.seq)+'\n')
  with open(output_name_p, "w") as p:
    p.write("---- %s %s ----\n" %(species, ref_gene_name))
    p.write(str(ref_protein)+'\n')

  counter = 0
  to_handle = []

  #try:
  file_name = "no genome input"
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]
    genome_records = SeqIO.index('%s/%s' %(extension,file_name),"fasta",IUPAC.unambiguous_dna)
    #print records["JLSA01000036"]
    result_handle = open("blast_%s/%s/%s.xml" %(extension,ref_gene_name,genome))

    [a,b,c,d] = combine_parse_single_blast(result_handle,ref_gene,genome_records,output_name_g)
    if b == '' or d == '':
      continue
    elif b==len_g and d == len_p:
      with open(output_name_g,"a") as g:
        g.write(a+'\n')
      with open(output_name_p,"a") as p:
        p.write(c+'\n')
    else:
      to_handle += [[a,b,c,d]]

  for item in to_handle:
    with open(output_name_g,"a") as g:
      g.write(item[1]+'\n')
      n = 0
      while n<len(item[1]):
        if item[0][n] == itemp[1][n]:
          g.write(' ')
        else:
          g.write('x')
        n += 1
      g.write('\n')
      g.write(item[0]+'\n \n')

    with open(output_name_p,"a") as p:
      p.write(item[3]+'\n')
      n = 0
      while n<len(item[3]):
        if item[2][n] == itemp[3][n]:
          p.write(' ')
        else:
          p.write('x')
        n += 1
      p.write('\n')
      p.write(item[2]+'\n \n')