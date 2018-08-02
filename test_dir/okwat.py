import sys,re,os


with open(sys.argv[1]) as t :
    dic={}
    for p in t :
        line=p.split("\t")
        if line[0][0]=="#":
            continue
        if (line[0],line[3],line[4],line[2]) not in dic.keys() :
            dic[(line[0],line[3],line[4],line[2])]=1
            #dic[(line[0],line[3],line[4],line[2])]=1
        elif (line[0],line[3],line[4],line[2]) in dic.keys() :
            print(dic[(line[0],line[3],line[4],line[2])])
            dic[(line[0],line[3],line[4],line[2])]+=1
u=0
tot=0
for k,v in dic.items():
    if v>1 :
        u+=1
        tot+=v
        print(k,v)


print(u)
print(tot)
print(tot-u)
