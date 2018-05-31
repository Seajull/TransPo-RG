from Bio import SeqIO
from subprocess import call
from pybedtools import BedTool
from version import getVersion
import sys, os, re, tempfile, argparse, warnings, pydoc

parser = argparse.ArgumentParser(add_help=False)

required = parser.add_argument_group("Required arguments")
optional = parser.add_argument_group("Optional arguments")
required.add_argument("-f1", "--fasta1", dest="fasta1",default=None, help="Input of reference fasta1 (old) <fasta>.")
required.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="Input of reference fasta2 (new) <fasta>.")
required.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input of tabbed file related to fasta1 <bed/gff/vcf>")
optional.add_argument("-b", dest="flank", type=int, default=50, help="Size of flank region to extract from each side of the annotation (default : 50).")
optional.add_argument("-c", "--cds",dest="cds", action="store_true", help="Enable control of postions of CDS inside mRNA (or gene).")
optional.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,help="Show this help message then exit.")
optional.add_argument("-i", "--index",dest="index", action="store_true", help="Create the index of fasta2.")
optional.add_argument("-n", "--notempfile",dest="notempf", action="store_true", help="Create all file instead of using temporary file.")
optional.add_argument("-o", "--output", dest="out", default=None, help="Output file (same format of tabbed file input).")
optional.add_argument("-t", "--type",dest="typeA", default=None, help="Only extract annotation of specified type (gff file only) (example : \"mRNA exon cds\" (case insensitive)).")
optional.add_argument("-v", "--verbose",dest="verbose", action="store_true", help="Enable message.")
optional.add_argument("-ver", "--version",dest="version", action="store_true", help="Show version and date of last update then exit.")
optional.add_argument("-w", "--warning",dest="warn", action="store_true", help="Disable warnings.")

if __name__ == '__main__':
    args = parser.parse_args()
    if args.version :
        getVersion()
        sys.exit()

    if args.warn :
        warnings.filterwarnings("ignore")

    if args.fasta1 == None:
        sys.exit("ERROR : Argument --fasta1 (-f1) is missing.")
    if args.fasta2 == None:
        sys.exit("ERROR : Argument --fasta2 (-f2) is missing.")
    if args.tabinput == None:
        sys.exit("ERROR : Argument --tabinput (-ti) is missing.")

    global ext
    ext=BedTool(args.tabinput).file_type
    if ext == "gff" :
        ext ="gff3"

    try:
        os.mkdir("result")
        if args.verbose :
            print("\n ----- Creating directory 'result/'. -----")
    except :
        pass

    # Creating tempfile (or file if args --notempfile is enable)
    tabOut=tempfile.NamedTemporaryFile()
    tabO=tabOut.name
    if not args.notempf :
        aln=tempfile.NamedTemporaryFile()
        alnN=aln.name
        selectS = tempfile.NamedTemporaryFile()
        selectedSeq = selectS.name
    else :
        alnN="result/aln_out_"+ext+".sam"
        fasta1Out=(args.fasta1).split(".")
        selectedSeq="result/"+fasta1Out[0].split("/")[-1]+"_selected_"+ext+"."+fasta1Out[1]

    if args.typeA != None and ext != "gff3" :
        print("")
        warnings.warn("Argument --typeA is ignored because --tabinput format isn't GFF.",Warning)

    if args.typeA != None and ext == "gff3" :
        typ=(re.findall("[a-zA-Z0-9]+",args.typeA))
        typeAclean=[l.lower() for l in typ]

def checkDependency() :
    """
        Check if all third-party program are installed by
        making a system call with "module list" then
        comparing with a required module list. If at least
        one thrird-party program isn't installed, it quit.
        Currently, this function check for environment
        modules but it need further test or a more generic
        test.
    """
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
    mandatoryMod=["bioinfo/bwa/0.7.15","bioinfo/bedtools/2.24.0"]
    goInstall=""
    for i in mandatoryMod:
        if i not in listM :
            goInstall += ("/".join(i.split("/")[1:]))+"  "
    if goInstall :
        sys.exit("ERROR : please, install following tools : " + goInstall)
    return

def parseFa() :
    """
        Later, in order to get the flank region and not
        going over last position of the chromosome, we need
        to index the length of them.
        This function parse the inputted file fasta1 and
        write in a tempfile all the chromosome ID and length.
    """
    lenChr = tempfile.NamedTemporaryFile()
    with open(lenChr.name,"w") as lenC :
        for seqF in SeqIO.parse(args.fasta1,"fasta") :
            lenC.write(seqF.id+"\t"+str(len(seqF)))
    return (lenChr)

def prefix() :
    """
        In order to link the tabbed file and the fasta
        file, we need to get the same chromosome ID. This
        function check for this and change ID prefix in the
        tabbed file to the same as the ID prefix in the
        fasta file.
    """
    tabPre=tempfile.NamedTemporaryFile()
    prefixTab=tabPre.name
    global change
    change=False
    fileTab=BedTool(args.tabinput)
    nul=""
    fasta=SeqIO.parse(args.fasta1,"fasta")
    first_seq=next(fasta)
    if (fileTab[0][0][0:-1])!= (first_seq.id[0:-1]) :
        with open(prefixTab,"w") as pre :
            for feat in fileTab :
                nul+= first_seq.id[0:-1]+feat[0][-1]+" "+(" ".join(feat[1:])+"\n")
        BedTool(nul, from_string=True).saveas(tabO)
        change=True
    return

def cutGff() :
    """
        If the argument --typeA is specified, we need to cut
        the gff in order to keep only specified type. We use
        lower() function on each type in order to be case
        insensitive.
    """
    global typeAclean
    if args.verbose :
        print("\n ----- Extracting GFF's annotation which match '"+",".join(typeAclean)+"'. -----")
    with open(args.tabinput,"r") as inpTab, open(tabO,"w") as out:
        for line in inpTab :
            fields = line.strip().split("\t")
            fields[2]=fields[2].lower()
            if fields[2] in typeAclean :
                out.write(("\t".join(fields))+"\n")
    return

# TODO : add ENG comment and remove FR comment 
# Méthode permettant de modifier le fichier tabulé pour ajouter des régions flanquantes 
def getFlank() :
    """
        Change the position in tabbed file in order to
        include flank region. For vcf file, we need to
        generate a new temporary tabbed file which contain
        the ID, start position and stop position. If the
        argument --notempfile is passed, it create the file
        [tabbed_file_name]_out.[ext] in the directory
        ./result.
    """
    if not args.notempf :
        tabOPut=tempfile.NamedTemporaryFile()
        tabOP=tabOPut.name
    else :
        tab=(args.tabinput).split(".")
        tabOP="result/"+tab[0].split("/")[-1]+"_out."+tab[1]
    lenChr=parseFa()
    fileTab2=BedTool(args.tabinput)
    if (args.typeA != None and ext== "gff3") or change :
        fileTab2=BedTool(tabO)
    if args.verbose and args.notempf:
        print("\n ----- Creating file '"+tabOP +"'. ----- ")
    if (fileTab2.file_type != "vcf"):
        fileTab2.slop(b=args.flank, g=lenChr.name ,output=tabOP)
    else :
        with open(lenChr.name,"r") as lenC :
            res=""
            for feature in fileTab2 :
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                if feature.stop+args.flank-1+(len(feature[3])-1) > int(lenghtC.group(1)):
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank-1+(len(feature[3])-1)
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank-1
                res += feature.chrom +" "+str(start)+" "+ str(stop)+"\n"
        BedTool(res, from_string=True).saveas(tabOP)
    try :
        return (tabOPut)
    except :
        return (tabOP)

def getFasta(tabOP):
    """
        This function just extract the sequence defined in
        the tabbed file create by getFlank(). The output is
        a fasta file. If --notempfile is enable, it create the
        file [fasta1_name]_selected_[ext].fasta in the
        directory ./result.
    """
    if args.verbose :
        print("\n ----- Extracting flanking sequences in '"+args.fasta1+"'. -----")
    if not args.notempf :
        BedTool(tabOP.name).sequence(fi=args.fasta1).save_seqs(selectedSeq)
    else :
        BedTool(tabOP).sequence(fi=args.fasta1).save_seqs(selectedSeq)
    return

def index() :
    """
        This function generate the index of inputted fasta2
        file with a system call to bwa index.
    """
    if args.verbose :
        print("\n ----- Generating index of '"+args.fasta2+"' using 'bwa index'. ----- \n")
    #TODO : remove system call and add a check if index is already generate
    call(["bwa","index",args.fasta2])
    print("")
    return

def align(tabOp):
    """
        This function make a system call to bwa mem to align
        fasta2 file with selected sequence of fasta1. The
        output file is in sam format.
        If argument --notempfile is enable it create the
        file "aln_out_[ext].sam".
    """
    getFasta(tabOp)
    if args.verbose :
        print("\n ----- Generating alignement file using 'bwa mem'. ----- \n")
    with open(alnN,"w") as out:
        if args.verbose :
            call(["bwa","mem", args.fasta2, selectedSeq],stdout=out)
            print("")
        else :
            call(["bwa","mem","-v","0", args.fasta2, selectedSeq],stdout=out)
            print("")
    return

def parseCigar(sam) :
    """
        This function parse the CIGAR code in alignment file
        in order to get the length of the match.
    """
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

def samToTab() :
    """
        In order to get the result of alignment usable for
        further  use, we create a new tabbed file in the
        same format as the inputted tabbed file. This file
        is the final output of this program.
        First of all, we extract the start and stop position
        and remove the flank region to get the initial
        annotation length. Then we generate the output file
        with these new position and keep all the information
        of the original tabbed file.
        The name of output file is by defaut
        "[fasta2_name]_out.[ext]" or the name that user
        specify with argument --output (but still with [ext]
        as extension.
    """
    res2 = re.search("/?(\w+)\.",args.fasta2)
    if args.out == None:
        args.out = "result/"+res2.group(1)+"_out."+ext
    else :
        args.out = "result/"+args.out.split(".")[0]+"."+ext
    start=0
    stop=0
    tabou=""
    samf=BedTool(alnN)
    lengh=parseCigar(samf)
    countLine=0
    if (ext=="gff3" and args.typeA != None) or change:
        args.tabinput=tabO
    with open(args.tabinput,"r") as tabi :
        for i in tabi : # tabi = tabbed file after all modification 
            line=i.split("\t")
            if line[0][0] == "#" :
                tabou+= i
            else :
                for f in samf : # samf = alignment file
                    if int(f[1])!=0 : # ignoring complementary match (flag 2048)
                        continue
                    res=re.search(":(\d+)-(\d+)",f[0])
                    if res :
                        if ext == "gff3" :
                            if int(res.group(1)) == int(line[3])-args.flank -1 :
                                start=int(f[3])+args.flank
                                stop=start+lengh[countLine]-(args.flank*2)-1
                                countLine+=1
                                break
                        elif ext == "bed" :
                            if int(res.group(1)) == int(line[1])-args.flank :
                                start=int(f[3])+args.flank
                                stop=start+lengh[countLine]-(args.flank*2)-1
                                countLine+=1
                                break
                        elif ext == "vcf" :
                            if int(res.group(1)) == int(line[1])-args.flank-1:
                                start=int(f[3])+args.flank
                                break
                        else :
                            countLine+=1
                #if f[11][-1]!="0" and f[5]=="101M":     # show ID of sequence which contain a missmatch
                #    print(f[0].split(":")[1]+"\t"+f[12])
                if ext == "vcf" and f[5]==str(args.flank*2+1)+"M": # perfect match only for snp
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
    return (args.out)

# TODO : add ENG comment and remove FR comment 
# Méthode permettant de récupérer les postions relatives des CDS dans les ARNm (par défaut) ou dans les gènes.
def getPosCds(tab) :
    """
        This function extract position of mRNA (by defaut,
        else it's position of gene) and the relative
        position of CDS within.
        It take a file in gff format and return a
        dictionary with tuple of mRNA (or gene) in keys and
        list of list of positon of CDS in value.
    """
    dicoPos={}
    posGene=()
    global typeAclean
    if args.typeA == None :
        typeD = True
    else :
        if ("gene" in typeAclean and "mrna" in typeAclean) :
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
                    warnings.resetwarnings()
                    warnings.filterwarnings("error")
                    warnings.warn("Start > stop",Warning)
                else :
                    dicoPos[posGene].append([cdsStart,cdsStop])
    return dicoPos

def isComplete(samtotabOut) :
    """
        This function call the function getPosCds(tab) with
        both inputted tabbed file and newly generated tabbed
        file and check if CDS in the newly generated file
        are in the same position within the mRNA (or gene).
    """
    if ext == "gff3" :
        dicoPos1=getPosCds(args.tabinput)
        dicoPos2=getPosCds(samtotabOut)
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
                print("Gene "+str(l)+" complete.")
            else :
                print("Gene "+str(l)+" incomplete.")
    return

if __name__ == "__main__":
    checkDependency()
    if args.typeA != None and ext == "gff3" :
        cutGff()
    prefix()
    tabOp = getFlank()
    if args.index:
        index()
    align(tabOp)
    samtotabOut=samToTab()
    if (args.cds and ext=="gff3" and (args.typeA==None or (("gene" in typeAclean or "mrna" in typeAclean) and "cds" in typeAclean))) :
        isComplete(samtotabOut)
