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
parser.add_argument("-te", "--tempfile",dest="tempf", action="store_true", help="Force la création du fichier fasta des séquences sélectionnées et du ficher tabulé intermédiaire.")

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
    else :
        args.out = "data/"+args.out


lenChr = tempfile.NamedTemporaryFile()

if not args.tempf : 
    tabOut=tempfile.NamedTemporaryFile()
    tabO=tabOut.name
    selectS = tempfile.NamedTemporaryFile()
    selectedSeq = selectS.name
else :
    fasta1Out=(args.fasta1).split(".")
    selectedSeq=fasta1Out[0]+"_selected."+fasta1Out[1]
    tab=(args.tabinput).split(".")
    tabO = (tab[0]+"_out."+tab[1])

fileTab=BedTool(args.tabinput)

if args.typeF != None and fileTab.file_type != "gff" :
    print("")
    warnings.warn("L'option --type (-t) est ignorée car le fichier tabulé n'est pas au format GFF.",Warning)

# Méthode permettant de récupérer l'ID des chromosomes et leurs tailles

def parseFa() :
    with open(lenChr.name,"w") as lenC : 
        for seqF in SeqIO.parse(args.fasta1,"fasta") :
            lenC.write(seqF.id+"\t"+str(len(seqF)))
    return 

parseFa()

# Test si index des chromosomes existant sinon on le réalise

def getFlank() :
    if args.verbose and args.tempf:
        print("\n ----- Création du fichier " +tabO +". ----- ")
        
    if (fileTab.file_type != "vcf"):
        fileTab.slop(b=args.flank, g=lenChr.name ,output=tabO)
        # Pour les fichiers VCF, on doit rajouter une colonne. On recrée donc un objet BedTools à partir
        # du fichier VCF original.
        
    else :
        with open(lenChr.name,"r") as lenC :
            res=""
            for feature in fileTab :
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                if feature.stop+args.flank-1 > int(lenghtC.group(1)):
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank-1
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank
                res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
        BedTool(res, from_string=True).saveas(tabO)
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
        print("\n ----- Récupération des features de type : \""+",".join(typeFclean)+"\" dans le GFF. -----")
    typeFclean="|".join(typeFclean)
    with open(tabO,"w") as out :
        call(["awk","tolower($3) ~ /"+typeFclean+"/",args.tabinput],stdout=out)
    return 
    

def getFasta() :
    if args.verbose :
        print("\n ----- Récupération des séquences flanquantes dans le fichier "+args.fasta1+". -----")
    BedTool(tabO).sequence(fi=args.fasta1).save_seqs(selectedSeq)
    # Sauvegarde dans le fichier désigné par la variable selectedSeq les séquences fasta correspondantes au informations contenue dans le fichier tabulé
    return 

def index() :
    if args.index:
        if args.verbose :
            print("\n ----- Indexation du fichier "+args.fasta2+" par BWA. ----- \n")
        call(["bwa","index",args.fasta2])
        print("")
    return

def align():
    getFasta()
    if args.verbose :
        print("\n ----- Réalisation et stockage de l'alignement dans le fichier "+args.out+". ----- \n")
    with open(args.out,"w") as out:
        if args.verbose :
            call(["bwa","mem", args.fasta2, selectedSeq],stdout=out)
        else :
            call(["bwa","mem","-v","0", args.fasta2, selectedSeq],stdout=out)
    #Fichier .sam ici, next ?
    return

def samToTab() :
    pos=0
    tabou=""
    samf=BedTool(args.out)
    with open(args.tabinput,"r") as tabi :
            for i in tabi :
                line=i.split("\t")
                #print(line)
                if line[0][0] == "#" :
                    tabou+= i 
                else :
                    for f in samf :
                        res=re.search(":(\d+)-\d+",f[0])
                        if res :
                            if int(res.group(1)) == int(line[1])-50 :
                                pos=int(f[3])+49
                                #print(pos)
                                break
                    if fileTab.file_type == "vcf" and f[5]=="100M":
                        tabou+=f[2]+" "+ str(pos) +" "+ " ".join(line[2:11])
                    elif fileTab.file_type == "bed" :
                        tabou+=f[2]+" "+ str(pos) +" "+ str(pos) +" "+ " ".join(line[2:3])
                    elif fileTab.file_type == "gff" :
                        pass            
            BedTool(tabou,from_string=True).saveas("out.vcf")
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
if fileTab.file_type == "vcf":
    samToTab()