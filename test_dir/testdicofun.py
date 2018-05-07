import sys,re,os

chro="chromosome"
with open(sys.argv[1]) as tab :
    for i in tab :
        res = re.search("("+chro+")"+"(_1)",i)
        if res :
            print(res.group(1))
