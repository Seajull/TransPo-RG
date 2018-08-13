import sys,os,re
if sys.argv[1][-1]=="/" :
    dirr=sys.argv[1][:-1]
else :
    dirr=sys.argv[1]
try :
    os.mkdir(dirr)
except :
    pass
