import sys, re, os
koi
bed = open(sys.argv[1],"r")

tabFile = bed.readlines()
countLine=0
for line in tabFile :
    countLine+=1
    res=re.search("\w+\t(\d+)\t(\d+)\t",line)
    #print("Taille variant : "+str(int(res.group(2))-int(res.group(1))))
    if int(res.group(2))-int(res.group(1))==1 :
        print("SNP à la ligne "+str(countLine))
