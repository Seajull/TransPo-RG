from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile, argparse, warnings, io, shutil

tab=open(sys.argv[1],"r").readlines()
sam=open(sys.argv[2],"r").readlines()
tabl=BedTool(tab)
samB = BedTool(sam)
taout="out.vcf"
'''
for f in samB :
    #if f[5] == "100M" :
    #print(f[0:7]) de 0 a 6   
    tabou+=f.chrom+" "+str(f.start+51)+" "+str(f[2:11])+"\n"
    print(tabou)
''' 
def samToTab() :
    pos=0
    tabou=""
    for i in tab :
        line=i.split("\t")
        #print(line)
        if line[0][0] == "#" :
            tabou+= i 
        else :
            for f in samB :
                res=re.search(":(\d+)-\d+",f[0])
                if res :
                    if int(res.group(1)) == int(line[1])-50 :
                        pos=int(f[3])+49
                        print(pos)
                        break
            if tabl.file_type == "vcf" :
                tabou+=f[2]+" "+ str(pos) +" "+ " ".join(line[2:11])
            elif tabl.file_type == "bed" :
                tabou+=f[2]+" "+ str(pos) +" "+ str(pos) +" "+ " ".join(line[2:3])
            elif tabl.file_type == "gff" :
                pass
            
    BedTool(tabou,from_string=True).saveas(taout)
    return

samToTab()
 
