import os, re, sys


strvcf="lol.vcf"
strbed="mdr.bed"
strrien="ptdr"
chaine="pp 1 p p 2p pp p3\tp4p5ppppp"
chaini="1 2 3\t4\t5\t 6"
def spli(st) :
    return (st.split("."))

#print(spli(strvcf)[0])
#print(spli(strbed)[0])
#
#print(spli(strrien)[0])
print(chaini.split("\t"))
