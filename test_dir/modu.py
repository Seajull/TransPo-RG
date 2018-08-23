from subprocess import call, Popen, PIPE, check_output
import os,sys,re, tempfile

p=Popen(["module list"],shell=True, stderr=PIPE)
ok,err=p.communicate()

res=re.findall("\w+/\w+/(?:\d+\.?)+",str(err))
print(res)
