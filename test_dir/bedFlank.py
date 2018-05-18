##################################################
# Script permettant de modifier le fichier BED   #
# pour récupérer x nucléotides de chaque coté    #
# du variant avant de lancer getfasta            #
##################################################

import sys, re, os

bedIn = open(sys.argv[1],"r")
fileOut=(sys.argv[1]).split(".")
bedOut = open(fileOut[0]+"_out."+fileOut[1],"w")
ref1 = open(sys.argv[2],"r")

bed = bedIn.readlines()
fasta = ref1.readlines()

# Modifier x pour modifier la tailler de la fênetres
# à récupérer (x avant le début du SNP/INDEL et x
# après la fin du SNP/INDEL)

x=50
lenChr={}

count=False
countNucl=0

# Récupère les numéros des chromosomes du fichier tabulé

for line in bed :
    numchr=re.search("^\w+(\d+)\t",line)
    if numchr:
        numchromo = numchr.group(1)
        if numchromo not in lenChr.keys():
            lenChr[numchromo]=0



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

for line in bed :
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
        bedOut.write(lineJ)
    else :
        bedOut.write(line)

