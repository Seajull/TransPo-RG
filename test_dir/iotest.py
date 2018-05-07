import sys, io, os, re, tempfile

out = io.StringIO()

with open(sys.argv[1],"r") as filo :
    for line in filo :
        out.write(line)

out.seek(0)
for i in out :
    print (i)


out.close()       
