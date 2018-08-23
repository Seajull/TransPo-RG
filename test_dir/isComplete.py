from pybedtools import BedTool
from subprocess import call
import re, os, sys, pprint, tempfile, warnings, argparse, datetime

parser = argparse.ArgumentParser()
parser.add_argument("-t1", "--tabinput1", dest="tabinput1", default=None, help="Input du fichier tabulé associé à fasta 1.")
parser.add_argument("-t2", "--tabinput2", dest="tabinput2", default=None, help="Input du fichier tabulé associé à fasta 1.")
args = parser.parse_args()
typeY="cds,gene"
ext="gff3"
typ=(re.findall("[a-zA-Z0-9]+",typeY))
typeAclean=[l.lower() for l in typ]

ts=datetime.datetime.now()
def cutGff() :
    global typeAclean
    print("\n ----- Extracting GFF's annotation which match '"+",".join(typeAclean)+"'. -----")
    with open(args.tabinput1,"r") as inpTab, open("tabO","w") as out:
        for line in inpTab :
            fields = line.strip().split("\t")
            if fields[2].lower() in typeAclean :
                out.write(("\t".join(fields))+"\n")
    return

def getPosCds(tab) :
    typeA="gene"
    dicoPos={}
    posGene=()
    global typeAclean
    with open(tab,"r") as out :
        numGene=0
        for line in out :
            if line[0]=="#" :
                continue
            lineSplit=line.split("\t")
            typeA = lineSplit[2].lower()
            start = lineSplit[3]
            stop = lineSplit[4]
            if typeA == "gene" or typeA == "mrna" :
                numGene+=1
                getags=lineSplit[-1]
                posGene=(getags,numGene,int(start),int(stop))
                if posGene not in dicoPos.keys():
                    dicoPos[posGene]=[]
            if typeA == "cds":
                cdsStart=int(start)-int(posGene[1])
                cdsStop=int(stop)-int(posGene[1])
                if cdsStart > cdsStop :
                    warnings.resetwarnings()
                    warnings.filterwarnings("error")
                    warnings.warn("Start > stop",Warning)
                resTag=re.search("ID=(\w+)",getags)
                resTagCds=re.search("Parent=(\w+)",lineSplit[-1])
                if resTag and resTagCds:
                    if str(resTag.group(1))==str(resTagCds.group(1)) :
                        dicoPos[posGene].append([cdsStart,cdsStop])
    return dicoPos

def isComplete() :
    """
        This function call the function getPosCds(tab) with
        both inputted tabbed file and newly generated tabbed
        file and check if CDS in the newly generated file
        are in the same position within the mRNA (or gene).
        It take the name of the tabbed file newly generated
        in argument.
    """
    if ext == "gff3" :
        dicoPos1=getPosCds("tabO")
        dicoPos2=getPosCds(args.tabinput2)
        print("getposCDS done")
        outTab ="ptrd"
        geneInt=[]
        #lastG=0
        geneId=""
        cdsId=""
        countG=0
        selectable=False
        filtered = "# File generated the "+datetime.datetime.now().strftime("%d %b %Y") + " with following command line : \n"+"# "+" ".join(sys.argv)+"\n"
        for key1 in dicoPos1.keys() :
            for key2 in dicoPos2.keys() :
                if key2[0]==key1[0] :
                    if len(dicoPos1[key1]) == len(dicoPos2[key2]) :
                        geneInt.append(key2[1])
                       # for v in range (0,len(dicoPos1[key1])) :
                       #     print(dicoPos1[key1][v])
                       #     if dicoPos1[key1][v] == dicoPos2[key2][v] : # TODO : c'est de la merde.
                       #         geneOk+=1
                       # print(geneOk)
                       # if geneOk >= len(dicoPos1[key1]) : # here we can add/rm condition to accept or not the mRNA/gene
                       #    # add the mRNA/gene number to the list of "acceptable mRNA/gene to select"
                       #     geneOk = 0
                       # else :
                       #     geneOk = 0

        tu=datetime.datetime.now()
        print(tu-ts)
        print(geneInt)
        if "gene" in typeAclean :
            typeC="gene"
        elif "mrna" in typeAclean :
            typeC="mrna"
        with open(args.tabinput2,"r") as tabou :
            for line in tabou :
                if line[0]=="#":
                    continue
                lineS=line.strip().split("\t")
                if lineS[2].lower() == typeC : # TODO : unreadable
                    resTag=re.search("ID=(\w+)",lineS[-1])
                    if resTag :
                        geneId=resTag.group(1)
                if lineS[2] =="CDS" :
                    resTagCds=re.search("Parent=(\w+)",lineS[-1])
                    if resTagCds :
                        cdsId=resTagCds.group(1)
                if lineS[2].lower() == typeC : # TODO : unreadable
                    countG+=1
                    if countG in geneInt :
                        selectable=True
                    else :
                        selectable=False
                if lineS[2] =="CDS" and geneId!=cdsId :
                    selectable=False
                if selectable :
                    filtered+=("\s".join(lineS))+"\n"
            countG=0
        tt=datetime.datetime.now()
        print(tt-tu)
        print(tt-ts)

        print(" ----- Generating filtered GFF file '"+"filtered_"+outTab+"'. -----\n")
        BedTool(filtered, from_string=True, deli="\s").saveas("filtered_"+outTab)
    return



cutGff()
isComplete()

dicoPos2=getPosCds(args.tabinput2)
for k,v in dicoPos2.items() :
    res = re.search("Sb01g022180",str(k))
    if res :
        print(k,v)
        le=len(v)
dicoPos1=getPosCds("tabO")
for k,v in dicoPos1.items() :
    res = re.search("Sb01g022180",str(k))
    if res :
        print(k,v)
        if le == len(v) :
           print("ok")
