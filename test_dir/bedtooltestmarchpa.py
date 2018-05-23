from pybedtools import BedTool
import sys,re,os

vcf=sys.argv[1]

tabVcf=BedTool(vcf)

for feature in tabVcf :
    if feature.start+1 != feature.stop :
        print(str(feature.start)+"    "+str(feature.stop))

