#####################################################
## Script permettant de modifier un fichier tabulé ##
## pour récupérer x nucléotides de chaque coté     ##
## du variant avant de lancer getfasta             ##
#####################################################

from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile

x=100
tabIn = open(sys.argv[1],"r")
tab=(sys.argv[1]).split(".")
tabOut = (tab[0]+"_out."+tab[1])
ref1 = open(sys.argv[2],"r")

#fileTab = tabIn.readlines()
fasta = ref1.readlines()
fileTab= BedTool(sys.argv[1])

# Méthode permettant de récupérer l'ID des chromosomes et leurs tailles

def parseFa() :
    lenChr={}
    for seqF in SeqIO.parse(sys.argv[1],"fasta") :
        lenChr[seqF.id]=len(seqF) 
    return lenChr

lenChr=parseFa()

# Test si index des chromosomes existant sinon on le réalise

try :
    ind=open(sys.argv[2]+".fai","r")
except :
    call(["samtools","faidx",sys.argv[2]])
    ind=open(sys.argv[1]+".fai","r")

# On utilise la méthode slop de BedTools pour récupérer les régions flanquantes (+ et - x)

if (fileTab.file_type != "vcf"):
    fileTab.slop(b=x, g=sys.argv[2]+".fai" ,output=tabOut)
    
# Pour les fichiers VCF, on doit rajouter une colonne. On recrée donc un objet BedTools à partir
# du fichier VCF original.

else :
    res=""
    for feature in fileTab :
        if feature.stop+x-1 > lenChr[feature.chrom]:
            stop=lenChr[feature.chrom]
        else :
            stop=feature.stop+x-1
        if feature.start-x < 0 :
            start=0
        else :
            start=feature.start-x
        res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
    BedTool(res, from_string=True).saveas(tabOut)
