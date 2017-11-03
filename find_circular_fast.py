#usage: python find_circular_fast.py contigs.fasta 150 2 contigs
## the last argument is either 'ids' or 'contigs'
import sys
from Bio import SeqIO

try:
  import edlib
except:
  print("ERROR: PLEASE RUN pip install edlib")
  sys.exit(1)

#read_length is 2nd argument; default 150 bp
try:
  read_length = int(sys.argv[2])
except:
  read_length = 150

# Max # of mismatches is 3rd argument; default 3
try:
  mismatch = int(sys.argv[3])
except:
  mismatch = 3

#FASTA file is 1st argument
for rec in SeqIO.parse(sys.argv[1], 'fasta'):

  start = str(rec.seq[:read_length])
  end = str(rec.seq[len(rec.seq)-read_length:])



  if edlib.align(start,end)['editDistance'] <= mismatch:
    if sys.argv[4] == 'contigs':
      print(">" + str(rec.id))
      print(str(rec.seq))
    else:
      print(rec.id)
