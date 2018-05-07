###############################################################
## Script permettant de récupérer les séquences flanquantes  ##
## à partir d'un fichier tabulé (VCF/BED/GFF3) d'une         ##
## référence fasta et de les aligner sur une autre           ##
## référence fasta                                           ##
###############################################################

from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile, argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f1", "--fasta1", dest="fasta1", default=None, help="Input de la référence fasta 1 (originelle)")
parser.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="Input de la référence fasta 2 (nouvelle)")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input d'un fichier tabulé fasta 1")
parser.add_argument("-o", "--output", dest="out", default=None, help="Fichier de sortie (.sam)")
parser.add_argument("-b", dest="flank", type=int, default=50, help="Taille de la fenêtre à sélectionner")
args = parser.parse_args()

if args.fasta1 == None:
    sys.exit("L'argument --fasta1 (-f1) est manquant")
if args.fasta2 == None:
    sys.exit("L'argument --fasta2 (-f2) est manquant")
if args.tabinput == None:
    sys.exit("L'argument --tabinput (-ti) est manquant")            

res1 = re.search("/?(\w+)\.",args.fasta1)
res2 = re.search("/?(\w+)\.",args.fasta2)
if args.out == None:
    args.out = "aln_"+res1.group(1)+"_"+res2.group(1)+".sam"


tabIn = open(args.tabinput,"r")
tab=(args.tabinput).split(".")
tabOut = (tab[0]+"_out."+tab[1])
fa1 = open(args.fasta1,"r")

fasta = fa1.readlines()
fileTab= BedTool(args.tabinput)

# Méthode permettant de récupérer l'ID des chromosomes et leurs tailles

def parseFa() :
    lenChr={}
    for seqF in SeqIO.parse(args.fasta1,"fasta") :
        lenChr[seqF.id]=len(seqF) 
    return lenChr

lenChr=parseFa()

# Test si index des chromosomes existant sinon on le réalise

def getFlank() :
    try :
        ind=open(args.fasta1+".fai","r")
        ind.close()
    except :
        call(["samtools","faidx",args.fasta1])

        # On utilise la méthode slop de BedTools pour récupérer les régions flanquantes
        # de taille passer en argument (50 par défaut)
        
    if (fileTab.file_type != "vcf"):
        fileTab.slop(b=args.flank, g=args.fasta1+".fai" ,output=tabOut)
    
        # Pour les fichiers VCF, on doit rajouter une colonne. On recrée donc un objet BedTools à partir
        # du fichier VCF original.

    else :
        res=""
        for feature in fileTab :
            if feature.stop+args.flank-1 > lenChr[feature.chrom]:
                stop=lenChr[feature.chrom]
            else :
                stop=feature.stop+args.flank-1
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank
                    res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
            BedTool(res, from_string=True).saveas(tabOut)

# Fichier [nom fichier tabulé]_out.[ext] (nom stocké dans tabOut) en sortie 


fasta1Out=(args.fasta1).split(".")
ref1=fasta1Out[0]+"_selected."+fasta1Out[1]

print("Récupération des séquences flanquantes dans le fichier "+args.fasta1+".\n")

#call(["bedtools","getfasta","-fi", args.fasta1, "-bed",tabOut, "-fo", ref1])
# ref1 (fasta) en sortie 


ind = input("Voulez-vous indexer "+res2.group(1)+" ? (oui/non)\n")
rep =[["OUI","yes","Oui","oui","o"],["Non","NON","non","n","no"]]
while ind not in rep[0] and ind not in rep[1]:
    ind = input("\nVeuillez rentrer oui ou non :\n")
if ind in rep[0]:
    print("\nIndexation du fichier "+args.fasta2+" par BWA.\n")
    #call(["bwa","index",args.fasta2])

print("Stockage de l'alignement dans le fichier "+args.out+".\n")

   # call(["bwa","mem", args.fasta2, ref1],stdout=bOut)

# Fichier .sam ici, next ?




