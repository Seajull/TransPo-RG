import sys,re,os


with open(sys.argv[1]) as tabular :
    for i in tabular :
        if i[0]=="#" :
            continue
        if i[0]=="@" :
            continue
        line=i.split("\t")
        res=re.search("Sb01g000390",i)
        if res :
            print (i)

