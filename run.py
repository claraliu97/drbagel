#Run the blast script
#Clara Liu
#2017/06/30
from test import *
from parameters import *
import os


os.chdir(species)
input = "%s/" %(extension)
f = file_names(input)
#create_db(input,f)
ref_gene_files = file_names("ref_genes/")
print ref_gene_files
for ref_gene in ref_gene_files:
  ref_gene_name = ref_gene[:-6]
  create_blast(input,f,ref_gene_name)
  [a,b,c,d] = parse_blast(f,ref_gene_name)
  print "%s: #DNA muts = %d(invalid = %d), #AA muts = %d(invalid = %d)" %(ref_gene_name,a-b,b,c-d,d)
print("---END---")

