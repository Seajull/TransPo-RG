import sys, os, re

i=0
alphabet = []
for letter in range(97,123):
        alphabet.append(chr(letter))
with open(sys.argv[1],"w") as a:
    while i<=int(sys.argv[2]) :
        a.write(alphabet[i]+" "+str(i)+"\n")
        i+=1
