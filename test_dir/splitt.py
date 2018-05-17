import os, re, sys


strvcf="lol.vcf"
strbed="mdr.bed"
strrien="ptdr"

def spli(st) :
    return (st.split("."))

print(spli(strvcf)[0])
print(spli(strbed)[0])

print(spli(strrien)[0])


