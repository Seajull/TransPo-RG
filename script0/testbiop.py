from Bio import SeqIO
from subprocess import call
import sys, re, os, timeit

def parseFa() :
    for seqF in SeqIO.parse(sys.argv[1],"fasta") :
        print(seqF.id)
        print(len(seqF))
    return

#print(timeit.timeit(parseFa,number=1))
fileOut=(sys.argv[1]).split(".")
print(fileOut)
sys.argv[1]=fileOut[0]+"_out."+fileOut[1]
print(sys.argv[1])

#call(["python" ,"bioflank.py", "data/tab/v1.bed", "data/fasta/fastav1.fasta"])
