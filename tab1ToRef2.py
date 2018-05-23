###############################################################
## Script permettant de récupérer les séquences flanquantes  ##
## à partir d'un fichier tabulé (VCF/BED/GFF3) d'une         ##
## référence fasta et de les aligner sur une autre           ##
## référence fasta                                           ##
###############################################################

from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
import sys, re, os, timeit, tempfile, argparse, warnings, io, shutil

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,help="Affiche ces messages help")
parser.add_argument("-f1", "--fasta1", dest="fasta1", default=None, help="Input de la référence fasta 1 (originelle).")
parser.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="Input de la référence fasta 2 (nouvelle).")
parser.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input du fichier tabulé associé à fasta 1.")
parser.add_argument("-o", "--output", dest="out", default=None, help="Fichier de sortie (format tabulé).")
parser.add_argument("-b", dest="flank", type=int, default=50, help="Taille de la fenêtre à sélectionner (par défaut : 50)")
parser.add_argument("-t", "--type",dest="typeF", default=None, help="Sélectionne uniquement les types souhaiter dans un gff3 (exemple de format de l'option : 'mRNA,exon,CDS').")
parser.add_argument("-i", "--index",dest="index", action="store_true", help="Force l'indexation du fichier fasta 2.")
parser.add_argument("-v", "--verbose",dest="verbose", action="store_true", help="Active l'affichage")
parser.add_argument("-w", "--warning",dest="warn", action="store_true", help="Désactive l'affichage des warnings.")
parser.add_argument("-c", "--cds",dest="cds", action="store_true", help="Active la vérification des positions des CDS par rapport aux gènes ou au mRNA")
parser.add_argument("-p", "--prefix",dest="prefix", action="store_true", help="Permet de modifier le préfixe du fichier tabulé pour qu'il corresponde à celui du fichier fasta")
parser.add_argument("-te", "--tempfile",dest="tempf", action="store_true", help="Force la création du fichier fasta des séquences sélectionnées, du ficher tabulé intermédiaire et du fichier d'alignement .sam.")

args = parser.parse_args()

if args.warn :
    warnings.filterwarnings("ignore")
if args.fasta1 == None:
    sys.exit("ERREUR : L'option --fasta1 (-f1) est manquant")
if args.tabinput == None:
    sys.exit("ERREUR : L'option --tabinput (-ti) est manquant")
if args.fasta2 == None:
    sys.exit("ERREUR : L'option --fasta2 (-f2) est manquant")

res2 = re.search("/?(\w+)\.",args.fasta2)

fileTab=BedTool(args.tabinput)
ext=fileTab.file_type
if ext =="gff" :
    ext="gff3"

if args.out == None:
    args.out = "resultat/"+res2.group(1)+"_out."+ext
else :
    args.out = "resultat/"+args.out.split(".")[0]+"."+ext

lenChr = tempfile.NamedTemporaryFile()

try :
    os.mkdir("./resultat")
    if args.verbose :
        print("\n ----- Création du dossier resultat. -----")
except :
    pass

tabOut=tempfile.NamedTemporaryFile()
tabO=tabOut.name

tabPre=tempfile.NamedTemporaryFile()
prefixTab=tabPre.name

mod=tempfile.NamedTemporaryFile()
listMod=mod.name
with open(listMod,"r+") as out :
    call(["module list"],shell=True,stderr=out)
    out.seek(0)
    for line in out :
        res=re.findall("\)\s([^\s]+)\s+",line)
        if "listM" in locals() :
            for i in res :
                listM.append(i)
        else :
            listM=res
mandatoryMod=["bioinfo/bwa/0.7.15","system/python/3.4.3","bioinfo/bedtools/2.24.0"]
for i in mandatoryMod:
    if i not in listM :
        print(i)

if not args.tempf :
    aln=tempfile.NamedTemporaryFile()
    alnN=aln.name
    tabOPut=tempfile.NamedTemporaryFile()
    tabOP=tabOPut.name
    selectS = tempfile.NamedTemporaryFile()
    selectedSeq = selectS.name
else :
    alnN="resultat/aln_out_"+ext+".sam"
    fasta1Out=(args.fasta1).split(".")
    selectedSeq="resultat/"+fasta1Out[0].split("/")[-1]+"_selected_"+ext+"."+fasta1Out[1]
    tab=(args.tabinput).split(".")
    tabOP="resultat/"+tab[0].split("/")[-1]+"_out."+tab[1]

if args.typeF != None and ext != "gff3" :
    print("")
    warnings.warn("L'option --type (-t) est ignorée car le fichier tabulé n'est pas au format GFF.",Warning)

if args.typeF != None and ext == "gff3" :
    splitType = args.typeF.split(",")
    typeFclean=[]
    for t in splitType :
        t=t.lower()
        if t not in typeFclean :
            typeFclean.append(t)

# Méthode permettant de récupérer l'ID des chromosomes et leurs tailles à l'aide de biopython

def parseFa() :
    with open(lenChr.name,"w") as lenC :
        for seqF in SeqIO.parse(args.fasta1,"fasta") :
            lenC.write(seqF.id+"\t"+str(len(seqF)))
    return

parseFa()

def prefix() :
    nul=""
    fasta=SeqIO.parse(args.fasta1,"fasta")
    first_seq=next(fasta)
    if (fileTab[0][0][0:-1])!= (first_seq.id[0:-1]) :
        with open(prefixTab,"w") as pre :
            for feat in fileTab :
                nul+= first_seq.id[0:-1]+feat[0][-1]+" "+(" ".join(feat[1:])+"\n")
    BedTool(nul, from_string=True).saveas(tabO)
    return

def cutGff() :
    global typeFclean
    if args.verbose :
        print("\n ----- Récupération des features de type : \""+",".join(typeFclean)+"\" dans le GFF. -----")
    typeFclean="|".join(typeFclean)
    with open(tabO,"w") as out :
        call(["awk","tolower($3) ~ /"+typeFclean+"/",args.tabinput],stdout=out)
    return


def getFlank() :
    fileTab2=BedTool(args.tabinput)
    if (args.typeF != None and ext== "gff3") or args.prefix :
        fileTab2=BedTool(tabO)
    if args.verbose and args.tempf:
        print("\n ----- Création du fichier " +tabOP +". ----- ")

    if (fileTab2.file_type != "vcf"):
        fileTab2.slop(b=args.flank, g=lenChr.name ,output=tabOP)
        # Pour les fichiers VCF, on doit rajouter une colonne. On recrée donc un objet BedTools à partir
        # du fichier VCF original.
    else :
        with open(lenChr.name,"r") as lenC :
            res=""
            for feature in fileTab2 :
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                if feature.stop+args.flank-1 > int(lenghtC.group(1)):
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank-1
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank-1
                res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
        BedTool(res, from_string=True).saveas(tabOP)
    return
    # Fichier [nom fichier tabulé]_out.[ext] (nom stocké dans tabOut) en sortie  

def getFasta():
    if args.verbose :
        print("\n ----- Récupération des séquences flanquantes dans le fichier "+args.fasta1+". -----")
    BedTool(tabOP).sequence(fi=args.fasta1).save_seqs(selectedSeq)
    # Sauvegarde dans le fichier désigné par la variable selectedSeq les séquences fasta correspondantes au informations contenue dans le fichier tabulé
    return

def index() :
    if args.verbose :
        print("\n ----- Indexation du fichier "+args.fasta2+" par BWA. ----- \n")
        call(["bwa","index",args.fasta2])
        print("")
    return

def align():
    getFasta()
    if args.verbose :
        print("\n ----- Réalisation de l'alignement. ----- \n")
    with open(alnN,"w") as out:
        if args.verbose :
            call(["bwa","mem", args.fasta2, selectedSeq],stdout=out)
        else :
            call(["bwa","mem","-v","0", args.fasta2, selectedSeq],stdout=out)
    #Fichier .sam ici, next ?
    return



def parseCigar(sam) :
    lenCig=[]
    for i in sam :
        if int(i[1]) != 0: # On ignore les match complémentaire (flag 2048)
            continue
        leng=0
        res =re.findall("\d+\w",i[5])
        for i in res :
            if i[-1] in ["M","=","X","I","S"] :
                leng+=int(i[:-1])
        lenCig.append(leng)
    return lenCig

# C'est le bazar, à simplifier ou au moins à ordonner pour pouvoir facilement rajouter des filtres.
def samToTab() :
    start=0
    stop=0
    tabou=""
    samf=BedTool(alnN)
    lengh=parseCigar(samf)
    countLine=0
    if (ext=="gff3" and args.typeF != None) or args.prefix:
        args.tabinput=tabO
    with open(args.tabinput,"r") as tabi :
        for i in tabi : # i parcours tabinput (ou tabO après cutGff)
            line=i.split("\t")
            if line[0][0] == "#" :
                tabou+= i
            else :
                for f in samf : # f parcours le fichier .sam
                    if int(f[1])!=0 : # On ignore les match complémentaire (flag 2048)
                        continue
                    res=re.search(":(\d+)-(\d+)",f[0])
                    if res :
                        if ext == "gff3" :
                            if int(res.group(1)) == int(line[3])-args.flank -1 :
                                start=int(f[3])+args.flank
                                stop=start+lengh[countLine]-(args.flank*2)-1
                                #print(pos)
                                countLine+=1
                                break
                        elif ext == "bed" :
                            if int(res.group(1)) == int(line[1])-args.flank :
                                start=int(f[3])+args.flank
                                stop=start+lengh[countLine]-(args.flank*2)-1
                                #print(pos)
                                countLine+=1
                                break
                        elif ext == "vcf" :
                            if int(res.group(1)) == int(line[1])-args.flank-1:
                                start=int(f[3])+args.flank
                                break
                        else :
                            countLine+=1
                if ext == "vcf" and f[5]==str(args.flank*2+1)+"M": # Match parfait uniquement pour un SNP
                    tabou+=f[2]+" "+ str(start) +" "+ " ".join(line[2:11])
                elif ext == "bed" :
                    tabou+=f[2]+" "+ str(start) +" "+ str(stop) +" "+line[3].replace(" ","_")+"\n"
                elif ext == "gff3" :
                    tabou+=f[2]+" "+line[1]+" "+line[2]+" "+ str(start) +" "+ str(stop) +" "+" ".join(line[5:8])+" "+line[8].replace(" ","_")+"\n"
        if ext == "vcf" :
            BedTool(tabou,from_string=True).saveas(args.out)
        elif ext == "bed" :
            BedTool(tabou,from_string=True).saveas(args.out)
        elif ext == "gff3" :
            BedTool(tabou,from_string=True).saveas(args.out)
    return

def getPosCds(tab) :
    #warnings.resetwarnings()
    #warnings.filterwarnings("error")
    dicoPos={}
    posGene=()
    global typeFclean
    if args.typeF == None :
        typeD = True
    else :
        if ("gene" in typeFclean and "mrna" in typeFclean) :
            typeD = True
        else :
            typeD=False
    with open(tab,"r") as out :
        numGene=0
        for line in out :
            lineSplit=line.split("\t")
            typeA = lineSplit[2].lower()
            if typeD and typeA == "gene" :
                continue
            start = lineSplit[3]
            stop = lineSplit[4]
            if typeA == "gene" or typeA == "mrna" :
                numGene+=1
                posGene=(numGene,int(start),int(stop))
                if posGene not in dicoPos.keys():
                    dicoPos[posGene]=[]
            if typeA == "cds":
                cdsStart=int(start)-int(posGene[1])
                cdsStop=int(stop)-int(posGene[1])
                if cdsStart > cdsStop :
                    warnings.warn("Start > stop",Warning)
                else :
                    dicoPos[posGene].append([cdsStart,cdsStop])
    #warnings.resetwarnings()
    return dicoPos

def isComplete() :
    if ext == "gff3" :
        dicoPos1=getPosCds(args.tabinput)
        dicoPos2=getPosCds(args.out)
        geneInt=[]
        lastG=0
        for key in dicoPos1.keys() :
            for keys in dicoPos2.keys() :
                if keys[0]>lastG :
                    lastG=keys[0]
                if keys[0]==key[0]:
                    if dicoPos1[key]==dicoPos2[keys]:
                        geneInt.append(keys[0])
                   # else :
                   #     print(key[0])
                   #     print(dicoPos1[key])
                   #     print(dicoPos2[keys])
        for l in range(1,lastG+1) :
            if l in sorted(geneInt) :
                print("Gène "+str(l)+" intègre.")
            else :
                print("Gène "+str(l)+" non intègre.")
    return

if args.typeF != None and ext == "gff3" :
    cutGff()
prefix()
getFlank()
if args.index:
    index()
align()
samToTab()
if args.cds and ext=="gff3" and (args.typeF==None or (("gene" in typeFclean or "mrna" in typeFclean) and "cds" in typeFclean)) :
    isComplete()
