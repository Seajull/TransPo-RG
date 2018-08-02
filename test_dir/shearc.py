import sys,os,re


with open(sys.argv[1]) as p :
    for l in p :
        line=l.split("\t")
        if line[0][0]=="#" :
            continue
        res=re.search("Sb01g000910.2",line[-1])
        res2=re.search("871551",line[3])
        if res :
            #print(res.group(0))
            #print(l)
            pass
        if res2 :
            print(l)
