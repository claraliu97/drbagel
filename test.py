#Dr.Bagel Blast and Analysis Script
#Clara Liu
#2017/06/12

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os

#return all the name of files in the given dir
def file_names(dir):
  f = []
  for (dirpath, dirnames, filenames) in os.walk(dir):
    f.extend(filenames)
    break
  return f

def makemydir(dir):
  try:
    os.makedirs(dir)
  except OSError:
    pass

def create_db(input_dir,file_list):
  makemydir("db_%s/" %extension)
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    #make_db = "makeblastdb -in %s%s.%s -dbtype nucl -parse_seqids \

    make_db = "makeblastdb -in %s%s.%s -dbtype nucl \
    -out db_%s/%s/%s -title '%s'\n" %(input_dir,genome,extension,extension,genome,genome,genome)

    os.system("mkdir -p temp")
    with open("temp/test.sh", "w") as text_file:
      text_file.write(make_db)

    os.system("sh temp/test.sh")
  os.system("rm -r temp")


def create_blast(input_dir,file_list):
  outdir = "blast_%s/%s/" %(extension,ref_gene_name)
  makemydir(outdir)
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]

    make_blast = "blastn -query ref_genes/%s.fasta -task megablast -db db_%s/%s/%s \
    -outfmt 5 -out %s%s.xml -num_threads 4" %(ref_gene_name,extension,genome,genome,outdir,genome)

    os.system("mkdir -p temp")
    with open("temp/test.sh", "w") as text_file:
      text_file.write(make_blast)

    os.system("sh temp/test.sh")
  os.system("rm -r temp")

def parse_single_blast2(result_handle,ref_gene,txt):
  pass

def parse_single_blast(result_handle,ref_gene,genome_records,txt):
  global counter
  counter += 1

  ref_gene = ref_gene.seq
  ref_protein = ref_gene.translate()

  blast_record = NCBIXML.read(result_handle)
  #find the best hit (filter the off-target hits)
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
  for n in range(len(query)):
    if query[n]!=sbjct[n]:
      match_g += ' '
      mismatch_g += [(n+1,query[n],sbjct[n])]
    else:
      match_g += '|'


  with open(txt, "a") as text_file:
    text_file.write('****DNA Alignment****\n')
    text_file.write('Sequence:\n%s\n' %alignment.title)
    text_file.write('Identities: %d/%d\n' %(match_g.count('|'),len(ref_gene)))
    #text_file.write('Positives: %d/%d\n' %(hsp.positives,len(hsp.query)))
    text_file.write('Mismatches:\n')

    print_parsed_blast([query,match_g,sbjct,mismatch_g],text_file)

  sbjct_dna = ''.join(str(sbjct) .split('-'))
  try:
    protein = Seq(sbjct_dna, IUPAC.unambiguous_dna).translate()

    alignments = pairwise2.align.globalms(ref_protein, protein,2,0,-100,-100)
    align = pairwise2.format_alignment(*alignments[0])
    lst = align.split('\n')
    match_p = ''
    mismatch_p = []
    for n in range(len(lst[0])):
      if lst[0][n] == lst[2][n]:
        match_p += '|'
      else:
        match_p += ' '
        mismatch_p += [(n+1,lst[0][n],lst[2][n])]

    with open(txt, "a") as text_file:
      text_file.write('****Protein Alignment****\n')
      text_file.write('Identities: %d/%d\n' %(match_p.count('|'),len(ref_protein)))

      print_parsed_blast([lst[0],match_p,lst[2],mismatch_p],text_file)

      text_file.write(' \n')

    return [len(mismatch_g),len(mismatch_p)]

  except:
    open(txt,"a").write("Invalid DNA seq. Fail to translate\n \n")
    return [0,0]


def print_parsed_blast(lst,text_file):
  query = lst[0]
  match = lst[1]
  sbjct = lst[2]
  mismatch = lst[3]

  for (n,q,s) in mismatch:
    text_file.write('%d: %s->%s\n' %(n,q,s))

  if len(mismatch)>0:
    n = 0
    while (n<len(query)):
      text_file.write('%s\n' %query[n:n+75])
      text_file.write('%s\n' %match[n:n+75])
      text_file.write('%s\n' %sbjct[n:n+75])
      n += 75


def parse_blast(file_list):
  #E_VALUE_THRESH = 0.001
  ref_gene = SeqIO.read('ref_genes/%s.fasta' %ref_gene_name,"fasta",IUPAC.unambiguous_dna)

  output_name = "blast_%s/%s.txt" %(extension,ref_gene_name)
  with open(output_name, "w") as f:
    f.write("%s %s\n" %(species, ref_gene_name))

  mut_sum = [0,0]

  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name[:file_name.find('.'+extension)]
    genome_records = SeqIO.index('%s/%s' %(extension,file_name),"fasta",IUPAC.unambiguous_dna)
    #print records["JLSA01000036"]
    result_handle = open("blast_%s/%s/%s.xml" %(extension,ref_gene_name,genome))

    [g,p] = parse_single_blast(result_handle,ref_gene,genome_records,output_name)
    mut_sum[0] += g
    mut_sum[1] += p

  return mut_sum

def run():
  os.chdir(species)
  input = "%s/" %(extension)
  #print input
  f = file_names(input)
  #create_db(input,f)
  #create_blast(input,f)
  [a,b] = parse_blast(f)
  print("---END---")
  #print "%s: #DNA muts = %d, #AA muts = %d" %(ref_gene_name,a,b)


counter = 0
species = "Mycobacterium_tuberculosis"
#"Mycobacterium_tuberculosis"
ref_gene_name = "Rv0001"
extension = "ffn"
id_spliter = ' ' #'|' if fna
run()