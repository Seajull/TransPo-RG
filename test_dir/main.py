from subprocess import call
import sys, re, time, argparse



parser = argparse.ArgumentParser()

parser.add_argument("-f1", "--fasta1", dest="fasta1", default=None, help="work in progress")
parser.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="work in progress")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="work in progress")
parser.add_argument("-o", "--output", dest="out", default=None, help="work in progress")
args = parser.parse_args()

print(args)

if args.fasta1 == None:
    sys.exit("L'argument --fasta1 (-f1) est manquant")
if args.fasta2 == None:
    sys.exit("L'argument --fasta2 (-f2) est manquant")
if args.tabinput == None:
    sys.exit("L'argument --tabinput (-ti) est manquant")            
if args.out == None:
    args.out=="aln.sam"

    
fasta1 = args.fasta1
fasta2 = args.fasta2
tabFile = args.tabinput
res1 = re.search("/?(\w+)\.",fasta1)
res2 = re.search("/?(\w+)\.",fasta2)

print("Mise à jour du fichier "+ args.tabinput +".\n")
call(["python" ,"bioFlank.py",args.tabFile,args.fasta1])

tabOut=tabFile.split(".")
tabFile=tabOut[0]+"_out."+tabOut[1]

fasta1Out=fasta1.split(".")
ref1=fasta1Out[0]+"_selected."+fasta1Out[1]

print("Récupération des séquences flanquantes dans le fichier "+sys.argv[1]+".\n")

call(["bedtools","getfasta","-fi", fasta1, "-bed",tabFile, "-fo", ref1])


ind = input("Voulez-vous indexer "+res1.group(1)+" ? (oui/non)\n")
rep =[["OUI","yes","Oui","oui","o"],["Non","NON","non","n","no"]]
while ind not in rep[0] and ind not in rep[1]:
    ind = input("\nVeuillez rentrer oui ou non :\n")
    
if ind in rep[0]:
    print("\nIndexation du fichier "+sys.argv[2]+" par BWA.\n")
    call(["bwa","index",fasta2])

fasta2Out=fasta2.split(".") 

bwaOut="aln_"+res1.group(1)+"_"+res2.group(1)+".sam"

print("Stockage de l'alignement dans le fichier "+bwaOut+".\n")

with open(bwaOut,"w") as bOut:
    call(["bwa","mem", fasta2, ref1],stdout=bOut)



