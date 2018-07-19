import os, sys, re
j=0
with open(sys.argv[1],"r") as fai :
    while j < 30 :
        print(j)
        j+=1
        for i in fai :
            print(i)
