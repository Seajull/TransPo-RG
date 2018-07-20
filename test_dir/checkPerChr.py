import sys, os, re, collections


def check(fil) :
    lol={}
    with open(fil,"r") as ok :
        for i in ok :
            tab=i.split("\t")
            if tab[0][0] == "#" :
                continue
            if tab[0] not in lol.keys() :
                lol[tab[0]]=1
            else :
                lol[tab[0]]+=1
    return(lol)
n=1
#while n<(len(sys.argv)) :
#    print("\nFichier : " + sys.argv[n]+"\n")
#    mlp= collections.OrderedDict(sorted(check(sys.argv[n]).items(), key=lambda jm: int(jm[0].split("_")[-1])))
#    for i,v in mlp.items() :
#        print(i, v)
#    print("\n ----------------------------------------- ")
#    n+=1

dicOri=collections.OrderedDict(sorted(check(sys.argv[1]).items(), key=lambda jm: (jm[0].split("_")[0],int(jm[0].split("_")[-1]))))
dicNew=collections.OrderedDict(sorted(check(sys.argv[2]).items(), key=lambda jm: (jm[0].split("_")[0],int(jm[0].split("_")[-1]))))

print("ID\tLOSS\t%")
for key in dicOri.keys() :
    if key in dicNew.keys():
        loss=(int(dicOri[key])-int(dicNew[key]))
        if loss < 0 :
            sys.exit("ERROR : More line than original file")
        print(key+"\t"+str(loss)+"\t"+str("%.2f" % round(loss*100/dicOri[key],2))+"%")
