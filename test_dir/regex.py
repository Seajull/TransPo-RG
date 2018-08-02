import re, sys, os


name="10"
res=re.search("[^\d]?(\d+)$", name)
if res :
    print(res.group(1))
