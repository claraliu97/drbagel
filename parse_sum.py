#parse to get seqs
#Clara Liu
#2017/07/19

def parse_seq_sum(txt, num):
  file = open(txt)
  counter = 0
  seqs = []
  l = 0
  r = 0
  while (counter < num):
    seq = ''
    line = file.readline()
    while (line!='\n'):
      seq += line
      line = file.readline()
    seqs += [seq]
    counter += 1
    print seq
  return seqs

def