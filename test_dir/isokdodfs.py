import sys,os,re

with open(sys.argv[1]) as j :
    for i in j :
        line = i.split("\t")
        if line[0][0] =="@" :
            continue
        res=re.search("Chr\d+",line[0])
        if res :
            print (line)

