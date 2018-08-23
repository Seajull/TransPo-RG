import sys,re,os


#with open(sys.argv[1])as aln :
#    for i in aln :
#        if i[0]=="@" :
#            continue
#        line=i.split("\t")
#        res=re.search("Z:((\d?)+\D+(\w?)+)",line[12])
#        if res :
#            print(res.group(1))


def mdParser():
    num=0
    with open(sys.argv[1])as aln :
        for i in aln :
            if i[0]=="@" :
                continue
            line=i.split("\t")
            print(line[12])
            if line[12][0:2]!="MD" :
                continue
            md=line[12].split(":")[-1]
            count=re.findall("\d+",md)
            tot=0
            for c in count :
                tot+=int(c)
            print(tot)
    return()

mdParser()
