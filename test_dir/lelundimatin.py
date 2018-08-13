import sys,re,os
numM=""
cl=0
lm=0
with open(sys.argv[1]) as l :
    for i in l :
        line =i.split("\t")
        if line[2]=="mRNA" :
            res=re.search("ID=Sb02g008411.2",line[-1])
            if res:
                cl=20
        if lm<cl :
            numM+=i
            lm+=1

print(numM)
