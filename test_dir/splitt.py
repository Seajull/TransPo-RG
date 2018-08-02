import os, re, sys


strvcf="lol.vcf"
strbed="mdr.bed"
strrien="ptdr"
chaine="pp 1 p p 2p pp p3\tp4p5ppppp"
chaini="1 2 3\t4\t5\t 6"
lastChar="1.2 .3  6 7 8 9"
lm="1.2.3.4.5.6"
#print(spli(strvcf)[0])
#print(spli(strbed)[0])
#
#print(spli(strrien)[0])
print(lm)
print(lm.split("."))
print(".".join(lm.split(".")))
