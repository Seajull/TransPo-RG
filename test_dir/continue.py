from subprocess import call, Popen,PIPE
import sys, os, re, datetime


lp=[0,1,2,3,4,5,6,7,8,9]

for i in lp :
    if i%2 == 0 :
        continue
    print(i)
