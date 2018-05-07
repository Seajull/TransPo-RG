###############################################################
## Script permettant de récupérer les séquences flanquantes  ##
## à partir d'un fichier tabulé (VCF/BED/GFF3) d'une         ##
## référence fasta et de les aligner sur une autre           ##
## référence fasta                                           ##
###############################################################

from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile, argparse, warnings, io, shutil

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,help="Affiche ces messages help")
parser.add_argument("-f1", "--fasta1", dest="fasta1", default=None, help="Input de la référence fasta 1 (originelle).")
parser.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="Input de la référence fasta 2 (nouvelle).")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input du fichier tabulé associé à fasta 1.")
parser.add_argument("-o", "--output", dest="out", default=None, help="Fichier de sortie (.sam).")
parser.add_argument("-b", dest="flank", type=int, default=50, help="Taille de la fenêtre à sélectionner (par défaut : 50)")
parser.add_argument("-t", "--type",dest="typeF", default=None, help="Sélectionne uniquement les types souhaiter dans un gff3 (exemple de format de l'option : 'mRNA,exon,CDS').")
parser.add_argument("-i", "--index",dest="index", action="store_true", help="Force l'indexation du fichier fasta 2.")
parser.add_argument("-v", "--verbose",dest="verbose", action="store_true", help="Active l'affichage")
parser.add_argument("-w", "--warning",dest="warn", action="store_true", help="Désactive l'affichage des warnings.")

args = parser.parse_args()

if args.warn :
    warnings.filterwarnings("ignore")
if args.fasta1 == None:
    sys.exit("ERREUR : L'option --fasta1 (-f1) est manquant")
if args.tabinput == None:
    sys.exit("ERREUR : L'option --tabinput (-ti) est manquant")            
if args.fasta2 == None:
    sys.exit("ERREUR : L'option --fasta2 (-f2) est manquant")
res1 = re.search("/?(\w+)\.",args.fasta1)
if args.fasta2 :
    res2 = re.search("/?(\w+)\.",args.fasta2)
    if args.out == None:
        args.out = "data/aln_"+res1.group(1)+"_"+res2.group(1)+".sam"

        
tabIn = open(args.tabinput,"r")
tab=(args.tabinput).split(".")
tabOut = (tab[0]+"_out."+tab[1])
fa1 = open(args.fasta1,"r")

fasta = fa1.readlines()
fileTab= BedTool(args.tabinput)

if args.typeF != None and fileTab.file_type != "gff" :
    print("")
    warnings.warn("L'option --type (-t) est ignorée car le fichier tabulé n'est pas au format GFF.",Warning)

# Méthode permettant de récupérer l'ID des chromosomes et leurs tailles

def parseFa() :
    with open("data/fasta/lenChr","w") as lenChr : # A modifier au plus vite avec tempfile ou stringIO but meh
        for seqF in SeqIO.parse(args.fasta1,"fasta") :
            lenChr.write(seqF.id+"\t"+str(len(seqF)))
    return 

parseFa()
#for i in lenChrtmp:
#    print(i)

# Test si index des chromosomes existant sinon on le réalise

def getFlank() :
    lenChr= open("data/fasta/lenChr","r")
    if args.verbose :
        print("\n ----- Création du fichier " +tabOut +". ----- \n")
        
    if (fileTab.file_type != "vcf"):
        fileTab.slop(b=args.flank, g="data/fasta/lenChr" ,output=tabOut)
        # Pour les fichiers VCF, on doit rajouter une colonne. On recrée donc un objet BedTools à partir
        # du fichier VCF original.
        
    else :
        res=""
        for feature in fileTab :
            for line in lenChr :
                lenC=re.search(feature.chrom+"\t(\d+)",line)
                if lenC :
                    break
            if feature.stop+args.flank-1 > int(lenC.group(1)):
                stop=int(lenC.group(1))
            else :
                stop=feature.stop+args.flank-1
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank
                    res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
            BedTool(res, from_string=True).saveas(tabOut)
    lenChr.close()
    return 
    # Fichier [nom fichier tabulé]_out.[ext] (nom stocké dans tabOut) en sortie
    
def cutGff() :
    splitType = args.typeF.split(",")
    typeFclean=[]
    for t in splitType :
        t=t.lower()
        if t not in typeFclean :
            typeFclean.append(t)
    if args.verbose :
        print(" ----- Récupération des features de type : \""+",".join(typeFclean)+"\" dans le GFF. ----- \n")
    typeFclean="|".join(typeFclean)
    with open(tabOut,"w") as out :
        call(["awk","tolower($3) ~ /"+typeFclean+"/",args.tabinput],stdout=out)
    return 
    

def getFasta() :
    fasta1Out=(args.fasta1).split(".")
    ref1=fasta1Out[0]+"_selected."+fasta1Out[1]
    if args.verbose :
        print(" ----- Récupération des séquences flanquantes dans le fichier "+args.fasta1+". ----- \n")
    call(["bedtools","getfasta","-fi", args.fasta1, "-bed",tabOut, "-fo", ref1])
    # ref1 (fasta) en sortieq
    return ref1

def index() :
    if args.index:
        if args.verbose :
            print(" ----- Indexation du fichier "+args.fasta2+" par BWA. ----- \n")
        call(["bwa","index",args.fasta2])
        print("")
    return

def align():
    ref1=getFasta()
    if args.verbose :
        print(" ----- Réalisation et stockage de l'alignement dans le fichier "+args.out+". ----- \n")
    with open(args.out,"w") as out:
        if args.verbose :
            call(["bwa","mem", args.fasta2, ref1],stdout=out)
        else :
            call(["bwa","mem","-v","0", args.fasta2, ref1],stdout=out)
    #Fichier .sam ici, next ?
    return

def getPosCds(tab,flank=0) :
    #if tab.file_type()
    warnings.resetwarnings()
    warnings.filterwarnings("error")
    dicoPos={}
    posGene=()
    with open(tab,"r") as out :
        for line in out :
            lineSplit=line.split("\t")
            typeA = lineSplit[2].lower()
            start =lineSplit[3]
            stop =lineSplit[4]
            if typeA == "gene" :
                posGene=(int(start),int(stop))
                if posGene not in dicoPos.keys():
                    dicoPos[posGene]=[]
            if typeA == "cds" :
                cdsStart=int(start)-int(posGene[0])+flank
                cdsStop=int(stop)-int(posGene[0])-flank
                if cdsStart > cdsStop :
                    warnings.warn("Start > stop",Warning)
                    break
                else :
                    dicoPos[posGene].append([cdsStart,cdsStop])
    warnings.resetwarnings()
    return dicoPos

#Problème pos gene et CDS -> indel ? 
#dicoPos1=getPosCds(args.tabinput)
#dicoPos2=getPosCds(tabout,args.flank)


getFlank()

if args.typeF != None and fileTab.file_type == "gff" :
    cutGff()
index()
align()


tabIn.close()
fa1.close()
