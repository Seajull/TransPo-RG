import sys, os, re

dico={}
dico[(1,546431,6546846)]="lol"
dico[(2,54,65)]="lol"

chouette=list(dico)


for i in chouette :
    if 2 in i :
        print(i)

print ([i for i in chouette if 2 in i])
