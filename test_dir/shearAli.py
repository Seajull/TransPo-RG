import sys,os,re


with open(sys.argv[1]) as al :
    for i in al :
        res = re.search("Sb04g024740\.1",i)
        if res :
            print(i)
