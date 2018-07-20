import os, sys, re

#lol=[0,1,2,3,4,5,6,7,8,9]
with open(sys.argv[1],"r")as lol :
    for i in lol :
        ipr+=1
        if int(i[0]) < 5 :
            print("continue \t" +str(i))
            continue
        print(i)
        break
    lol.seek(0)
    for i in lol :
        print("ça marche ? \t" + str(i))
