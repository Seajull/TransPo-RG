#####################################################
## Script permettant de modifier un fichier tabulé ##
## pour récupérer x nucléotides de chaque coté     ##
## du variant avant de lancer getfasta             ##
#####################################################

import sys, re, os

tabIn = open(sys.argv[1],"r")
fileOut=(sys.argv[1]).split(".")
tabOut = open(fileOut[0]+"_out."+fileOut[1],"w")
ref1 = open(sys.argv[2],"r")

fileTab = tabIn.readlines()
fasta = ref1.readlines()

# Modifier x pour modifier la taille de la fênetres
# à récupérer (x avant le début du SNP/INDEL et x
# après la fin du SNP/INDEL)

x=50
lenChr={}

# Récupère les numéros des chromosomes du fichier tabulé

def getNumChr():
    numchromo=""
    for line in fileTab :
        numchr=re.search("^\w+(\d+)\t",line) ## utilisé if line[0:X] pour plus d'efficacité 
        if numchr:
            numchromo = numchr.group(1)
            if numchromo not in lenChr.keys():
                lenChr[numchromo]=0
    return

# Récupère la taille des séquences des chromosomes
# qui sont présent dans le fichier tabulé

def getLen() :
    count=False
    countNucl=0
    numchromo=""
    for line in fasta :
        if line[0] ==">" :
            if numchromo in lenChr.keys():
                lenChr[numchromo]=countNucl     
                countNucl=0
            numchr=re.search("chromosome_(\d+)",line)
            if numchr :
                numchromo = numchr.group(1)
                if numchromo in lenChr.keys():
                    count=True
                else :
                    count=False
        else :
            if count :
                for i in line :
                    if i != "\n":
                        countNucl+=1
    if numchromo in lenChr.keys():
        lenChr[numchromo]=countNucl      
    return

# Ecris dans un nouveau fichier tabulé les nouvelles position (+ et - x)

def newPos() :
    for line in fileTab :
        if line[0:11]=="chromosome_" :
            lineS=line.split("\t")
            if int(lineS[1]) < x :
                lineS[1] = 0
            else :
                lineS[1] = int(lineS[1])-x 
            if int(lineS[2])+x > lenChr[line[11]]:
                lineS[2] = lenChr[line[11]]
            else :
                lineS[2] = int(lineS[2])+x
            lineJ="\t".join(map(str,lineS))
            tabOut.write(lineJ)
        else :
            tabOut.write(line)
    return


getNumChr()
getLen()
newPos()
