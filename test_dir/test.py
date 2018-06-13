#!/usr/bin/env/ python3

import sys, os, re

ref1 = open(sys.argv[1],"r") 
tab = open(sys.argv[2],"r")
output = open(sys.argv[3],"w")

fasta1 = ref1.read()
tab1 = tab.read()
chromo=""
pos={}
countvar={}

''' Nombre de nucléotide séléctionné avant et après le SNP '''
x=200

''' Récupération des positions des SNP '''
for line in tab1 :
    numchr = re.search("chr(\d+)\s",line)
    if numchr :
        numchromo = numchr.group(1)
        posvar = re.search("(\d+)\s",line)
        if posvar :
            if numchromo in pos.keys():
                pos[numchromo].append(posvar.group(1))
            else :
                pos[numchromo]=posvar.group(1)

''' Récupération des séquences autour des SNP '''
for line in fasta1 :
    chevron = re.search(">")
    if chevron :
        countnucl=0
        numchrfa = re.search("chromosome\s(\d+)",line)
        if numchrfa :
            numchromofa = numchrfa.group(1)
        
    else :
        if numchromofa in pos.keys() :
            for i in line :
                if countnucl in pos[numchromofa].values():
                    if numchromofa in countvar.keys() :
                        countvar[numchromofa]+=1
                    else :
                        countvar[numchromofa]=1
                    output.write(">chromosome "+str(numchromofa)+" variant " + str(countvar[numchromofa])+"\n")
                    if countnucl<x :
                        pos1=0
                    else :
                        pos1=countnucl-x
                    #if countnucl+oopsie
                    
                    for j in range(countnucl-x,countnucl+x+1):
                        output.write(linej)
                        
                countnucl+=1
            





            
