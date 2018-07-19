import sys, os, re, csv

#lol =[(5,12,4),(2,6,45),(18,1,7),(1,30,3),(1,5,15),(1,8,20),(1,5,26)]
#print(sorted(lol, key= lambda ptdr: (str(ptdr[0]),ptdr[1],ptdr[2])))
fil=[]
with open(sys.argv[1], "r") as samsoul, open("samsoul2","w") as soulsam :
    for i in samsoul :
        if i[0]!="@" :
            fil.append(i.split("\t"))
        else :
            soulsam.write(i)
    #for i in sorted(samsoul) :
    #    soulsam.write(i)
    for i in sorted(fil, key=lambda ok: (ok[0][0],int(ok[0].split("_")[1].split(":")[0]),int(ok[0].split(":")[1].split("-")[0]),int(ok[0].split(":")[1].split("-")[1]))) :
        soulsam.write("\t".join(i))
