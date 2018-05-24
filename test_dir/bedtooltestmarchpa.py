from pybedtools import BedTool
import sys,re,os

vcf=sys.argv[1]

tabVcf=BedTool(vcf)

for feature in tabVcf :
    if len(feature[4]) > 1:
        print(len(feature[3]))
