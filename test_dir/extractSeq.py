from Bio import SeqIO
import sys, re, os

fasta1=SeqIO.parse(sys.argv[1],"fasta")
fasta2=SeqIO.parse(sys.argv[2],"fasta")
for i in fasta1 :
    print(i.seq[19047])
    print(i.seq[32667])
    print(i.seq[125489])
for i in fasta2 :
    print(i.seq[19260])
    print(i.seq[37306])
    print(i.seq[78150])
