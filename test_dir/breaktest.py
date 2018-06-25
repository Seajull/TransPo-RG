import os, sys, re

lis=[0,1,2,3,4,5,6,7,8,9]
ais=[10,11,12,13,14,15,16,17,18,19]
with open(sys.argv[1],"r") as filo :
    for i in lis :
        print(i)
        for line in filo :
            if int(line[2]) == 2 :
                break
            print(line)
