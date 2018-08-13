import sys,os,re
count=0
with open(sys.argv[1]) as m :
    for i in m :
        line=i.split("\t")
        if line[2] =="gene" :
            count+=1
print(count)
