from pybedtools import BedTool
import sys, re, os


saF = open(sys.argv[1],"r")

sam=BedTool(saF)


def parseCigar () :
    lengy=[]
    for i in sam :
        leng=0
        res =re.findall("\d+\w{1}",i[5])
        for i in res :
            if i[-1] in ["M","=","X","D","N"] :
                leng+=int(i[:-1])
        lengy.append(leng)
    return lengy





print(parseCigar())
