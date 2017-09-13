#Run the blast script
#Clara Liu
#2017/06/30
from drbagel import *
from parameters import *
import os


"""
exist = folder_names("db_ffn/")
new_f = []
for filename in f:
	if not filename[:-4] in exist:
		new_f += [filename]
f = new_f
"""

#create_db(input,f)


"""
new_ref = []
for ref in ref_gene_files:
	if ref not in ["Rv0001.fasta","Rv0016c.fasta"]:
		new_ref += [ref]
ref_gene_files = new_ref
"""
#f = ['1304279.3.ffn']
#ref_gene_files = ['Rv0001.fasta']

#create_blast(input,['1262525.3.ffn'],'Rv0001')

def blast_and_parse():
  #create_db(input,f)
  for ref_gene in ref_gene_files:
    ref_gene_name = ref_gene[:-6]
    create_blast(input,f,ref_gene_name)
    #[a,b,c,d] = parse_blast(f,ref_gene_name)
    #write_seq(f,ref_gene_name)
    #summary = "%s: #DNA muts = %d(invalid = %d), #AA muts = %d(invalid = %d)" %(ref_gene_name,a-b,b,c-d,d)
    #print summary
    #with open("result.txt", "a") as text_file:
        #text_file.write(summary+"\n")

  ("---END---")

def count_all_invalid():
  for ref_gene in ref_gene_files:
    ref_gene_name = ref_gene[:-6]
    inv = count_invalid(f,ref_gene_name)
    summary = "%s %d %d" %(ref_gene_name,inv,7311-inv)
    #print summary
    with open("result.txt", "a") as text_file:
        text_file.write(summary+"\n")
  print("---END---")

def combine_parse():
  for ref_gene in ref_gene_files:
    ref_gene_name = ref_gene[:-6]
    combine_parse_blast(f,ref_gene_name)
  print("---END---")


def rename_outputs():
  curdir = os.getcwd()
  ref_genes = file_names(curdir+'/ref_genes/')
  f = []
  for r in ref_genes:
    f += [r[:-6]]
  print f
  for ref_gene in f:
    rename(curdir+'/output_seq/%s/AA' %(ref_gene),'.txt','.fasta')
    rename(curdir+'/output_seq/%s/DNA' %(ref_gene),'.txt','.fasta')

#combine_parse()
#rename("fasta/",".txt",".fasta")
#rename("ref_genes/",".txt",".fasta")

def run():
  os.chdir(species)
  input = "%s/" %(extension)
  f = file_names(input)
  ref_gene_files = file_names("ref_genes/")
  with open("result.txt", "w") as text_file:
    text_file.write('****Result****\n')
    text_file.write('Reference gene name | Invalid | Valid')

  if db:
    if current_os == 'WIN':
      create_db_win(input,f)
    elif current_os == 'MAC':
      create_db_mac(input,f)
    else:
      print 'Error from Clara: Invalid OS'

  if blast:
    for ref_gene in ref_gene_files:
      ref_gene_name = ref_gene[:-6]
      if current_os == 'WIN':
        create_blast_win(input,f,ref_gene_name)
      elif current_os == 'MAC':
        create_blast_mac(input,f,ref_gene_name)
      else:
        print 'Error from Clara: Invalid OS'

  if write_individual_fasta:
    for ref_gene in ref_gene_files:
      ref_gene_name = ref_gene[:-6]
      write_seq(f,ref_gene_name)

#rename_outputs()