import sys,os,re
prit=2

with open(sys.argv[1]) as p :
    for l in p :
        line=l.split("\t")
        if line[0][0]=="#" :
            continue
        res=re.search("Sb02g008411.2",line[-1])
        res2=re.search("871551",line[3])
        if res2 :
            #print(res.group(0))
            #print(l)
            pass
        if res :
           prit=0
        if prit==1 :
            break
        if prit == 0 :
            print(l)
            prit+=1
       
