from pybedtools import BedTool
import sys, os, re, datetime

alnN=sys.argv[1]
ext="vcf"
aout="lul"
flank=50
tabinput=sys.argv[2]
def parseCigar(sam) :
    """
        This function parse the CIGAR code in alignment file
        in order to get the length of the match.
    """
    lenCig=[]
    for i in sam :
        if int(i[1]) != 0: # Ignoring complementary match (flag 2048)
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
        as extension). It return the name of the file created.
    """
    start=0
    stop=0
    tabou=""
    lengh=parseCigar(BedTool(alnN))
    countLine=0
    okw=False
    with open(tabinput,"r") as tabi, open(alnN,"r") as sam :
        for i in tabi : # tabi = tabbed file after all modification 
            line=i.split("\t")
            if line[0][0] == "#" :
                tabou += i
                continue
            tabou += "# File generated the "+datetime.datetime.now().strftime("%d %b %Y") + " with following command line : \n"+"# "+" ".join(sys.argv)+"\n"
            break
        tabi.seek(0)
        for f in sam : # samf = alignment file inside Bedtools object
            samf=f.split("\t")
            try :
                int(samf[1])
            except :
                continue
            if int(samf[1])!=0 : # ignoring complementary match (flag 2048)
                continue
            res=re.search(":(\d+)-(\d+)",samf[0])
            for i in tabi :
                line=i.split("\t")
                if line[0][0]=="#" :
                    continue
                if res :
                    if ext == "gff3" :
                        if int(res.group(1)) == int(line[3])-flank -1 :
                            start=int(samf[3])+flank
                            stop=start+lengh[countLine]-(flank*2)-1
                            countLine+=1
                            okw=True
                            break
                        else :
                            countLine+=1
                    elif ext == "bed" :
                        if int(res.group(1)) == int(line[1])-flank :
                            start=int(samf[3])+flank
                            stop=start+lengh[countLine]-(flank*2)-1
                            countLine+=1
                            okw=True
                            break
                        else :
                            countLine+=1
                    elif ext == "vcf" :
                        if int(res.group(1)) == int(line[1])-flank-1:
                            start=int(samf[3])+flank
                            okw=True
                            break
            #if f[11][-1]!="0" and f[5]=="101M":     # show ID of sequence which contain a missmatch
            #    print(f[0].split(":")[1]+"\t"+f[12]) 
            if ext == "vcf" and samf[5]==samf[12].split(":")[-1]+"M" and okw : # perfect match only for snp
                tabou+=line[0]+"\s"+ str(start) +"\s"+ "\s".join(line[2:11])
            elif ext == "bed" and okw:
                tabou+=line[0]+"\s"+ str(start) +"\s"+ str(stop) +"\s"+line[3]+"\n"
            elif ext == "gff3" and okw:
                tabou+=line[0]+"\s"+line[1]+"\s"+line[2]+"\s"+ str(start) +"\s"+ str(stop) +"\s"+"\s".join(line[5:8])+"\s"+line[8]+"\n"
            okw=False
        print(" ----- Creating file '"+aout+"'. ----- \n")
        BedTool(tabou,from_string=True, deli="\s").saveas(aout)
    return (aout)

lol=samToTab()
