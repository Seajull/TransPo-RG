import re, sys, os



lul="ID=Sb01g0020453shfgd.dsfomifgdh.1;Parent=lololo"
lol="ID=Sb5612G1253.df1.5;."
res=re.search("ID=(.+);",lul)
ros=re.search("Parent=(\w+((\.?\w+)?)+)",lol)
ree=re.search("(ID|Parent)=(\w+)",lol)
if res :
    print(res.group(1))
if ros :
    print(ros.group(1))
if ree :
    print(ree.group(2))
