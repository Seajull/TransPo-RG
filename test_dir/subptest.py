import subprocess, time, sys, tempfile, shutil
from subprocess import call
##fast=open(sys.argv[1],"r")
#paoui=tempfile.NamedTemporaryFile()
#
#r="oskouuuuur "
#lol=0
#with open(paoui.name,"w") as osko :
#    while lol<10:
#        osko.write(r)
#        lol+=1
#
#z="viiiite "
#
#with open(paoui.name,"a") as oskor :
#    while lol<10:
#        oskor.write(z)
#        lol+=1
#    oskor.write(z)
#     
#
#        
#with open(paoui.name,"r") as oksoo :
#    for i in oksoo :
#        print(i)
#


call(["modulecmd python load bioinfo/bedtools/2.24.0"],shell=True)
