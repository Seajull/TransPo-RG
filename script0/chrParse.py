##########################################################
# Script permettant de modifier des entêtes fasta        #
# pour permettre d'utiliser getfasta de Bedtools         #
# (format d'entête : chrX avec X le numéro du chromosome #
##########################################################

import sys, os, re

fa = open(sys.argv[1],"r")
out=(sys.argv[1]).split(".")[0]+"_out.fasta"
output = open(out,"w")
fasta = fa.readlines()


for line in fasta :
    if line[0] == ">" :
        chromo = re.search("chromosome\s(\d+)",line)
        # Si l'on retrouve le numéro du chromosome dans l'entête, 
        # on modifie l'entête par ">chrX"
        if chromo :
            output.write(">Chr"+chromo.group(1)+"\n")
        # Sinon on garde l'entête intact
        # On peut décider d'enlever toute la séquence pour
        # un gain de mémoire / temps d'écriture
        else :
            output.write(line)
    else :
        output.write(line)
        
            
