#drbagel.py
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
  for files in f:
    if 'DS' in files:
      f.remove(files)
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

def rename(dir,old_ext,new_ext):
  f = []
  for (dirpath, dirnames, filenames) in os.walk(dir):
    f.extend(filenames)
    break
  for file_name in f:
    if file_name[-len(old_ext):]==old_ext:
      file = dir+"/"+file_name
      base = os.path.splitext(file)[0]
      os.rename(file, base.lower() + new_ext)

def create_db_win(input_dir,file_list):
  makemydir("db_%s/" %extension)
  makemydir("temp/")

  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    make_db = "makeblastdb -in %s%s.%s -dbtype nucl \
    -out db_%s/%s/%s -title '%s'\n" %(input_dir,genome,extension,extension,genome,genome,genome)

    #with open("temp/test.sh", "w") as text_file:
      #text_file.write(make_db)

    os.system("sh %s" %make_db)
    #os.system("sh temp/test.sh")
  #os.system("rm -r temp/test.sh")

def create_db_mac(input_dir,file_list):
  makemydir("db_%s/" %extension)

  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    make_db = "makeblastdb -in %s%s.%s -dbtype nucl \
    -out db_%s/%s/%s -title '%s'\n" %(input_dir,genome,extension,extension,genome,genome,genome)

    os.system(make_db)


def create_blast_win(input_dir,file_list,ref_gene_name):
  counter = 0
  if len(file_list)/100>0:
    num_d_100 = len(file_list)/100
  else:
    num_d_100 = 1
  #print ref_gene_name

  outdir = "blast_%s/%s/" %(extension,ref_gene_name)
  makemydir(outdir)
  makemydir("temp/")
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    if counter%num_d_100 ==0:
      print "blast: %d/100" %(counter/num_d_100)
    counter += 1

    make_blast = "blastn -query ref_genes/%s.fasta -task megablast -db db_%s/%s/%s \
    -outfmt 5 -out %s%s.xml -num_threads 4" %(ref_gene_name,extension,genome,genome,outdir,genome)

    with open("temp/test.sh", "w") as text_file:
      text_file.write(make_blast)

    os.system("sh temp/test.sh")
  os.system("rm -r temp")

def create_blast_mac(input_dir,file_list,ref_gene_name):
  counter = 0
  if len(file_list)/100>0:
    num_d_100 = len(file_list)/100
  else:
    num_d_100 = 1
  print ref_gene_name

  outdir = "blast_%s/%s/" %(extension,ref_gene_name)
  makemydir(outdir)
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    if counter%num_d_100 ==0:
      print "blast: %d/100" %(counter/num_d_100)
    counter += 1

    make_blast = "blastn -query ref_genes/%s.fasta -task megablast -db db_%s/%s/%s \
    -outfmt 5 -out %s%s.xml -num_threads 4" %(ref_gene_name,extension,genome,genome,outdir,genome)

    os.system(make_blast)

def process_result_handle(genome,ref_gene_name):
  ref_gene = SeqIO.read('ref_genes/%s.fasta' %ref_gene_name,"fasta",IUPAC.unambiguous_dna)

  genome_records = SeqIO.index('%s/%s' %(extension,genome+'.'+extension),"fasta",IUPAC.unambiguous_dna)
  result_handle = open("blast_%s/%s/%s.xml" %(extension,ref_gene_name,genome))
  blast_record = NCBIXML.read(result_handle)
  #find the best hit (filter the off-target hits)
  try:
    alignment = sorted(blast_record.alignments,cmp=lambda x,y: x.length-y.length)[-1]
    id_l = alignment.title.find(id_spliter)
    id_r = alignment.title.find(id_spliter,id_l+1)
    id = alignment.title[id_l+1:id_r]
    sbjct_gene = genome_records[id].seq
    title = alignment.title
    hsp = sorted(alignment.hsps,cmp=lambda x,y: x.identities-y.identities)[-1]

    query = hsp.query.upper()

    len_upstream = hsp.query_start - 1
    len_downstream = len(ref_gene)-hsp.query_end
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

    sbjct = sbjct_left+hsp.sbjct.upper()+sbjct_right

    query = (ref_gene[:hsp.query_start-1] + query + ref_gene[hsp.query_end:]).upper()

    return [str(sbjct),str(query),len_upstream-len(sbjct_left),len_downstream-len(sbjct_right),title]

  except IndexError:
    return ['','',0,0,'']

class Blasted:

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



  def __init__(self,genome,ref_gene_name,option):
    self.ref_gene_name = ref_gene_name
    self.genome = genome
    self.left_miss = 0
    self.right_miss = 0
    self.sbjct = ''
    if option=='xml':
      [s,q,l,r,t] = process_result_handle(genome,ref_gene_name)
      #complete sbjct sequence (no dashes added even if incomplete)
      self.sbjct = s
      #complete query sequence with possible dash in the middle
      #self.query = q
      self.title = t
      self.left_miss = l
      self.right_miss = r
    elif option=='fasta':
      txt = 'output_seq/%s/%s' %(ref_gene_name,genome)
      file = open(txt)
      str_file = file.read()
      newline = str_file.find('\r\n')
      self.title = str_file[str_file.find('>')+1:newline]
      str_seq = ''
      newline += 2
      while str_file.find('\r\n',newline)>=0:
        r = str_file.find('\r\n',newline)
        str_seq += str_file[newline:r]
        newline = r+2
      self.sbjct = str_seq.upper()
      l = 0
      while(self.sbjct[l]=='-'):
        l += 1
      r = len(self.sbjct)-1
      while(self.sbjct[r]=='-'):
        r -= 1
      self.left_miss = l
      self.right_miss = len(self.sbjct)-1-r

  def validify(self):
    if self.sbjct == '':
      return False

  def get_ref_gene(self):
    return SeqIO.read('ref_genes/%s.fasta' %self.ref_gene_name,"fasta",IUPAC.unambiguous_dna)

  def get_ref_protein(self):
    return self.get_ref_gene().translate(to_stop=True)

  def get_title(self):
    return self.title

  def get_seq(self):
    return ''.join(self.sbjct.split('-'))

  def get_print_seq(self):
    return '-'*self.left_miss + self.sbjct + '-'*self.right_miss

  def get_protein(self):
    try:
      return str(Seq(self.get_seq(), IUPAC.unambiguous_dna).translate(to_stop=True))
    except Bio.Data.CodonTable.TranslationError:
      return ''
    except:
      return ''

  # make sure dir 'output_seq/AA' and 'output_seq/DNA' exist
  def write_single_seq(self):
    dir = 'output_seq/%s/' %(self.ref_gene_name)
    txt = '%s.fasta' %(self.genome)
    blast_error = []
    translate_error = []

    protein = self.get_protein()
    sbjct = self.get_print_seq()
    if sbjct!='' :
      with open(dir+"DNA/"+txt, "w") as text_file:
        text_file.write('>%s\n' %self.genome)
        n = 0
        while n<len(sbjct):
          text_file.write(sbjct[n:n+100]+'\n')
          n += 100
    else:
      blast_error += [self.genome]

    protein = self.get_protein()

    if protein!='' :
      with open(dir+"AA/"+txt, "w") as text_file:
        text_file.write('>%s\n' %self.genome)
        n = 0
        while n<len(protein):
          text_file.write(protein[n:n+100]+'\n')
          n += 100
    else:
      translate_error += [self.genome]

    return [blast_error,translate_error]

  def attach_exact_alignment(self,txt):
    pass

def write_seq(file_list,ref_gene_name):
  counter = 0
  if len(file_list)/100>0:
    num_d_100 = len(file_list)/100
  else:
    num_d_100 = 1

  makemydir('output_seq/%s/' %ref_gene_name)
  makemydir('output_seq/%s/DNA/' %ref_gene_name)
  makemydir('output_seq/%s/AA/' %ref_gene_name)

  #print '%s:' %ref_gene_name
  for file_name in file_list:
    if "DS" in file_name:
      continue

    if counter%num_d_100 == 0:
      print "write: %d/100" %(counter/num_d_100)
    counter += 1

    genome = file_name[:file_name.find('.'+extension)]
    blasted = Blasted(genome,ref_gene_name,'xml')
    [bE,tE] = blasted.write_single_seq()
    if len(bE)!=0 or len(tE)!=0:
      with open('output_seq/%s/error_seqs.txt' %ref_gene_name,'a') as txt:
        txt.write('Seqs not writen for %s: \n' %ref_gene_name)
        txt.write('BLAST error (no alignment):\n')
        for item in bE:
          txt.write(item+'\n')
        txt.write('Translation error: \n')
        for item in tE:
          txt.write(item+'\n')

def setup(species):
  makemydir(species)