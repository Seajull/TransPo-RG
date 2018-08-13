import sys,os,re

with open(sys.argv[1]) as un, open(sys.argv[2]) as deux :
    for i in un :
        if i not in deux :
            print(i)
