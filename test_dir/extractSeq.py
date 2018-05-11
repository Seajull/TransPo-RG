from Bio import SeqIO
import sys, re, os

fasta1=SeqIO.parse(sys.argv[1],"fasta")
fasta2=SeqIO.parse(sys.argv[2],"fasta")
for i in fasta1 :
    #print(i.seq[0:100])
    print(i.seq[42343:42444])
    print(i.seq[42392:42395])
   # print(i.seq[19211:19311])
for i in fasta2 :
    #print(i.seq[84:87])
    print(i.seq[42135:42236])
    print(i.seq[42184:42187])
