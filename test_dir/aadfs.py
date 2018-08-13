import sys,re,os

totlost=0
with open(sys.argv[1]) as l :
    for i in l :
        line = i.split("\t")
        if line[0]=="ID" :
            continue
        print(line[3])
        totlost+=int(line[3])
print(totlost)
