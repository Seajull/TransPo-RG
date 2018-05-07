from pybedtools import BedTool
from subprocess import call
import sys, pprint, tempfile, warnings, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--type",dest="typeF", default=None, help="Sélectionne uniquement les types souhaiter dans un gff3 (exemple de format de l'option : 'mRNA,exon,CDS').")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input du fichier tabulé associé à fasta 1.")
parser.add_argument("-i", "--ignore",dest="ignore", action="store_true", help="Active l'affichage")
args = parser.parse_args()

obj=BedTool(args.tabinput)
ext=args.tabinput.split(".")
sortie="cdsPos."+ext[1]
#if (obj.file_type =="bed"):
warnings.filterwarnings("ignore")

if args.typeF != None and obj.file_type != "gff" :
        warnings.warn("L'option --type (-t) est ignorée car le fichier tabulé n'est pas au format GFF",Warning)

if obj.file_type == "gff" :
    pass
else :
    warnings.warn("notGFF",Warning)
    
#for a in obj :
    #print(a)
    #print(a.fields[2])
    #print(a.attrs.items())
   # print("\n")

def cutGff() :
    splitType = args.typeF.split(",")
    typeFclean=[]
    for t in splitType :
        t=t.lower()
        if t not in typeFclean :
            typeFclean.append(t)
    typeFclean="|".join(typeFclean)
    with open(sortie,"w") as out :
        call(["awk","tolower($3) ~ /"+typeFclean+"/",args.tabinput],stdout=out)
    return

if not args.ignore :
    cutGff()

def getPosCds(tab,flank=0) :
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
                geneStart=int(start)+flank
                geneStop=int(stop)-flank
                if geneStart > geneStop :
                    warnings.warn("Start > stop",Warning)                
                posGene=(geneStart,geneStop)
                if posGene not in dicoPos.keys():
                    dicoPos[posGene]=[]
            if typeA == "cds" :
                cdsStart=int(start)-int(posGene[0])+flank
                cdsStop=int(stop)-int(posGene[0])-flank
                if (cdsStart > cdsStop) or (cdsStart < posGene[0]) or (cdsStop > posGene[1]) :
                    warnings.warn("Start > stop",Warning)
                else :
                    dicoPos[posGene].append([cdsStart,cdsStop])
    warnings.resetwarnings()
    return dicoPos


dico=getPosCds(args.tabinput)

for i,j in dico.items() :
    print(i,j)


