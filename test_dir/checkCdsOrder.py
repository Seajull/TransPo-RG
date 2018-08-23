import sys,re,os
sli={}
numGENE=0
jj=0
cpk=0
lol=""
with open(sys.argv[1]) as j :
    for i in j :
        line=i.split("\t")
        if line[0][0]=="#" :
            continue
        if line[2]=="gene" :
            numGENE+=1
            sli[numGENE]=[]
            if numGENE == 17346 :
                yes=True
                cpk=10
        if jj<cpk :
            lol+=i
            jj+=1
        if line[2]=="CDS" :
            res=re.search("Note=(\d+)",line[-1])
            if res :
                sli[numGENE].append(res.group(1))
def isSosted(li) :
    if sorted(li,key=int)==li :
        return True
    else :
        return False

for k,v in sli.items():
    if (isSosted(list(v)))==False :
        print(k,v)
 

