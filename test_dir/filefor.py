import os, sys, re
coiunt=0
with open(sys.argv[1],"r") as a, open(sys.argv[2],"r") as b :
    for j in b:
        for i in a :
            coiunt+=1
            if i[0]=="h" :
                if int(i[1])==8 :
                    break
            elif i[0]=="d" and i[1]==4 :
                break
            else :
                pass
            print(i)
print(coiunt)
