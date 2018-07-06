import os, sys, re
with open(sys.argv[1],"r") as fai :
    for i in fai :
        lul=i.split("\t")
        try :
            int(lul[0])
        except :
            print("cool")
            print(lul[0])
