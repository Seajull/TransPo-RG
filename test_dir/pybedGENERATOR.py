from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile, argparse, warnings,pybedtools

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", dest="out", default="aln_fastav1_fastav2.sam", help="Fichier de sortie (.sam).")
parser.add_argument("-e", "--verbose",dest="verbose", action="store_true", help="Active l'affichage")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input du fichier tabulé associé à fasta 1.")
args = parser.parse_args()

bed=BedTool(args.tabinput)
'''
vcf = BedTool("Chrx 130 . A T . . .", from_string=True)
print(vcf)

list_interv=[]
bed=BedTool(args.sam)
for a in bed :
    res=re.search("(\d+)-(\d+)",a.fields[0])
    sous=int(res.group(2))-int(res.group(1))
    list_interv.append([a[0],a[1],a[3]])
    #if sous != len(a) :
    #   print(a[0]+"\t"+a.fields[3]+"\t"+str(sous)+"\t" +str(len(a))+"\t"+a.fields[5]+"\t")
    

print(list_interv)
tabf=BedTool(list_interv)

if tabf == bed :
    print(tabf.file_type)

print(tabf)
#line1=tabf[0]
'''

def samToTab(dicoPos) :
    newTab=[]
    sam=BedTool(args.out)
    
    if bed.file_type=="gff" :
        for a in sam :
            if len(a) < 100 :
                continue
            res=re.search("(\d+)-(\d+)",a[0])
            if res :
                refPos=(int(res.group(1)),int(res.group(2)))           
                print(str(refPos)+"\t"+str(dicoPos.keys()))
                if refPos in dicoPos:
                    print("loooooooooooooooooool")
                
            newTab.append([a.chrom,a.start+50,a.stop-50,a[0]])
    else :
        for a in sam :
            newTab.append([a.chrom,a.start+50,a.stop-50,a[0]])       
    bedtoolOut=BedTool(newTab)
    return bedtoolOut 

#Créer un objet BedTool dans getPosCds
def getPosCds(tab,flank=0) :
    #warnings.resetwarnings()
    #warnings.filterwarnings("error")
    dicoPos={}
    posGene=()
    with open(tab,"r") as out :
        for line in out :
            lineSplit=line.split("\t")
            print(lineSplit)
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
    #warnings.resetwarnings()
    return dicoPos

dico=getPosCds(args.tabinput)
bed=samToTab(dico)
'''
for a in bed :
    print(a)
    #pass
'''
