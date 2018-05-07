import timeit, sys, re

ref1 = open(sys.argv[1],"r")
fasta = ref1.readlines()







def res() :
    count=0
    for line in fasta :
        if re.search("chromosome",line):
            count+=1
    print(count)
    return

def compare() :
    count=0
    for line in fasta :
        if line[0:10]=="chromosome" :
            count+=1
   # print(count)
    return


def compare2() :
    count=0
    for line in fasta :
       if line.startswith("chromosome") :
            count+=1
    #print(count)
    return



def find() :
    count=0
    for line in fasta :
        
        if "chromosome" in line :
            count+=1
   # print(count)
    return

print("res ")
print((timeit.timeit(res, number=100))*1000)
print("compare ")
print((timeit.timeit(compare, number=100))*1000)
print("compare2 ")
print((timeit.timeit(compare2, number=100))*1000)
print("find ")
print((timeit.timeit(find, number=100))*1000)



