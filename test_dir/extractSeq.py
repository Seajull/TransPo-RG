from Bio import SeqIO
import sys, re, os

fasta1=SeqIO.parse(sys.argv[1],"fasta")
fasta2=SeqIO.parse(sys.argv[2],"fasta")
print("")
for i in fasta1 :
    #print(i.seq[78078:78079])
    print(i.seq[78028:78129])
    fa1=(i.seq[78028:78129])
for i in fasta2 :
    #print(i.seq[77870:77871])
    print(i.seq[77820:77921])
    fa2=(i.seq[77820:77921])
print("")
print(len(fa1))
print("")
if fa1 == fa2 : 
    print("similaire")

for i in range (0,101) :
    if fa1[i] != fa2[i] :
        print(str(i)+"\t"+fa1[i]+"\t"+fa2[i])
