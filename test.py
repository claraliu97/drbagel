from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os

ref_gene_name = "test_gene.fasta"

#rename the .xls file into .xlsx file if .xls exist
f = []
for (dirpath, dirnames, filenames) in os.walk(os.getcwd()):
  f.extend(filenames)
  break
for file_name in f:
  if "DS" in file_name:
    continue
  genome = file_name+"/"+file_name
  genome_name = genome+".fna"
  out_name = "output/"+genome+".xml"

  test_sh = "makeblastdb -in "+ genome_name +" -dbtype nucl -parse_seqids \
  -out output/"+genome+" -title '"+genome+"' \n\
  blastn -query "+ ref_gene_name +" -task megablast -db output/"+genome+" \
  -outfmt 5 -out "+ out_name +" -num_threads 4"

  os.system("mkdir -p temp")
  with open("temp/test.sh", "w") as text_file:
    text_file.write(test_sh)

  os.system("sh temp/test.sh")

  #E_VALUE_THRESH = 0.001
  ref_gene = SeqIO.read(ref_gene_name,"fasta",IUPAC.unambiguous_dna)
  ref_protein = ref_gene.seq.translate()

  result_handle = open(out_name)
  blast_record = NCBIXML.read(result_handle)
  with open("output/Output.txt", "w") as text_file:
    for alignment in blast_record.alignments:
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

        n = 0
        while (n<len(hsp.query)):
          text_file.write(hsp.query[n:n+75]+'\n')
          text_file.write(hsp.match[n:n+75]+'\n')
          text_file.write(hsp.sbjct[n:n+75]+'\n')
          n += 75

        dna = Seq(hsp.sbjct, IUPAC.unambiguous_dna)
        protein = dna.translate()

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
      n = 0
      while (n<len(lst[0])):
        text_file.write(lst[0][n:n+75]+'\n')
        text_file.write(lst[1][n:n+75]+'\n')
        text_file.write(lst[2][n:n+75]+'\n')
        n += 75



  os.system("rm -r temp")
print("---END---")