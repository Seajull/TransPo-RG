import sys,re,os

pos1=int(sys.argv[2])
pos2=pos1+1
with open(sys.argv[1]) as f :
    for line in f :
        lines=line.split("\t")
        if lines[0][0]=="#":
            continue
        if int(lines[pos1])>=int(lines[pos2]) :
            print(str(lines[pos1]) +"\t"+ str(lines[pos2]))
