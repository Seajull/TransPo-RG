from pybedtools import BedTool
from Bio import SeqIO
import re, sys, os

dico={}
with open(sys.argv[1],"r") as wtf :
    for line in wtf :
        if line.split("\t")[0][0] != "#" :
            if int((line.split("\t")[0]).split("_")[-1]) not in dico.keys() :
                dico[int((line.split("\t")[0]).split("_")[-1])] = 1
            else :
                dico[int((line.split("\t")[0]).split("_")[-1])] += 1
print(dico)


#def parseFa() :
#    """
#        Later, in order to get the flank region and not
#        going over last position of the chromosome, we need
#        to index the length of them.
#        This function parse the inputted file fasta1 and
#        write in a tempfile all the chromosome ID and length.
#    """
#    with open("lenChr","w") as lenC :
#        for seqF in SeqIO.parse(sys.argv[1],"fasta") :
#            lenC.write(seqF.id+"\t"+str(len(seqF))+"\n")
#    return


#parseFa()
