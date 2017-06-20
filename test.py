#Dr.Bagel Blast and Analysis Script
#Clara Liu
#2017/06/12

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os

ref_gene_name = "Rv0001"
out_folder = "output/"+ref_gene_name+"/"
output_name = out_folder+"Output.txt"

#return all the name of folders in the given dir
def file_dirs(dir):
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

def create_blast(file_list):
  makemydir("output/"+ref_gene_name)
  for file_name in file_list:
    if "DS" in file_name:
      continue
    genome = file_name+"/"+file_name
    genome_name = genome+".fna"
    out_name = out_folder+file_name+".xml"

    make_db = "makeblastdb -in "+ "input/"+genome_name +" -dbtype nucl -parse_seqids \
    -out "+out_folder+genome+" -title '"+genome+"'"
    make_blast = "blastn -query "+ ref_gene_name +".fasta -task megablast -db output/"+genome+" \
    -outfmt 5 -out "+ out_name +" -num_threads 4"

    os.system("mkdir -p temp")
    with open("temp/test.sh", "w") as text_file:
      text_file.write(make_blast)

    os.system("sh temp/test.sh")
  os.system("rm -r temp")

#E_VALUE_THRESH = 0.001
ref_gene = SeqIO.read(ref_gene_name+'.fasta',"fasta",IUPAC.unambiguous_dna)
ref_protein = ref_gene.seq.translate()

open(output_name, "w").write(ref_gene_name+"\n")


def parse_blast(file_list):
  mut_sum = [0,0]

  for file_name in file_list:
    if "DS" in file_name:
      continue
    result_handle = open(out_folder+file_name+".xml")
    blast_record = NCBIXML.read(result_handle)
    with open(output_name, "a") as text_file:
      alignment = blast_record.alignments[0]
      for hsp in alignment.hsps:
        #if hsp.expect < E_VALUE_THRESH:
        text_file.write('****DNA Alignment****\n')
        text_file.write('Sequence:\n%s\n' %alignment.title)
        text_file.write('Identities: %d/%d\n' %(hsp.identities,len(hsp.query)))
        text_file.write('Positives: %d/%d\n' %(hsp.positives,len(hsp.query)))
        text_file.write('Mismatches:\n')
        mismatch = []
        for n in range(len(hsp.query)):
          if hsp.match[n]!='|':
            mismatch += [(n+1,hsp.query[n],hsp.sbjct[n])]
        for (n,q,s) in mismatch:
          text_file.write('%d: %s->%s\n' %(n,q,s))

        if hsp.identities<len(ref_gene) or len(hsp.query)>len(ref_gene):
          mut_sum[0] += len(ref_gene)-hsp.identities
          n = 0
          while (n<len(hsp.query)):
            text_file.write(hsp.query[n:n+75]+'\n')
            text_file.write(hsp.match[n:n+75]+'\n')
            text_file.write(hsp.sbjct[n:n+75]+'\n')
            n += 75

        mut = ''.join(hsp.sbjct.split('-'))
        mut_dna = Seq(mut, IUPAC.unambiguous_dna)
        protein = mut_dna.translate()

        alignments = pairwise2.align.globalms(ref_protein, protein,2,0,-100,-100)
        align = pairwise2.format_alignment(*alignments[0])
        lst = align.split('\n')
        match = ''
        mismatch = []
        for n in range(len(lst[0])):
          if lst[0][n] == lst[2][n]:
            match += '|'
          else:
            match += ' '
            mismatch += [(n+1,lst[0][n],lst[2][n])]
        lst[1] = match
        text_file.write('****Protein Alignment****\n')
        text_file.write('Identities: %d/%d\n' %(lst[1].count('|'),len(lst[0])))
        #text_file.write('Positives: %d/%d\n' %(hsp.positives,len(hsp.query)))
        for (n,q,s) in mismatch:
            text_file.write('%d: %s->%s\n' %(n,q,s))
        if lst[1].count('|')<len(ref_protein) or len(lst[0])>len(ref_protein):
          mut_sum[1] += len(ref_protein)-lst[1].count('|')
          n = 0
          while (n<len(lst[0])):
            text_file.write(lst[0][n:n+75]+'\n')
            text_file.write(lst[1][n:n+75]+'\n')
            text_file.write(lst[2][n:n+75]+'\n')
            n += 75
        text_file.write(' \n')
  return mut_sum

def run():
  f = file_dirs("input")
  #create_blast(f)
  [a,b] = parse_blast(f)
  print("---END---")
  print "%s: #DNA muts = %d, #AA muts = %d" %(ref_gene_name,a,b)

run()