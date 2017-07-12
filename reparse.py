#Run the blast script
#Clara Liu
#2017/06/30
from test import *
from parameters import *
import os


os.chdir(species)
input = "%s/" %(extension)
f = file_names(input)
"""
exist = folder_names("db_ffn/")
new_f = []
for filename in f:
	if not filename[:-4] in exist:
		new_f += [filename]
f = new_f
"""

#create_db(input,f)
ref_gene_files = ["Rv0001.fasta","Rv0016c.fasta"]

for ref_gene in ref_gene_files:
  ref_gene_name = ref_gene[:-6]
  #create_blast(input,f,ref_gene_name)
  [a,b,c,d] = parse_blast(f,ref_gene_name)
  summary = "%s: #DNA muts = %d(invalid = %d), #AA muts = %d(invalid = %d)" %(ref_gene_name,a-b,b,c-d,d)
  print summary
print("---END---")