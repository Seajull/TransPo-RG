import sys, re, os

with open(sys.argv[1],"r") as sysy :
    lul=[]
    for line in sysy :
        res=re.search(":(\d+)-(\d+)", line)
        if res :
            lul.append(int(res.group(1)))
if sorted(lul)==lul :
    print("sorted")
else :
    print("no")
