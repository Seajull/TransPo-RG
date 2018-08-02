import sys,re,os,argparse

parser=argparse.ArgumentParser()
parser.add_argument("-d",dest="direc",default="result",help="lul")


args=parser.parse_args()
di=args.direc
if di[-1]=="/" :
    di=di[:-1]

try :
    os.mkdir(di)
except :
    pass



print(di)
