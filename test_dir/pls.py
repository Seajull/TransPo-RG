import sys, re, os



with open(sys.argv[1],"r") as pls :
    for line in pls :
        linep=line.split("\t")
        if linep[0][0]=="@" :
            continue
        res=re.search("(\d+)M",linep[5])
        if res :
            if int(res.group(1))>101 :
                print(res.group(1))
