from subprocess import call, Popen,PIPE
import sys, os, re, datetime




with open("wat","w") as wat :
    p = Popen(["bwa", "mem", sys.argv[1], sys.argv[2]], stdout=wat, stderr=PIPE)

    err=p.communicate()
    print(err[-1].decode("utf-8"))

