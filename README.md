# TransPo-RG

----------
TransPo-RG, *Transfer of Position to Resequenced Genome* (or TPRG), is a python
script for transfering position of annotation from one genome assembly to a another
(e.g. newer version for a reference genome). It extract sequences with flanking
nucleotide (50bp in each side by default) and map them on the target genome.
Then we generate the output file in the same format of the entry file with only
changing positions. 

### Usage :

It's a command-line program. It needs three files in entry : two genome file in
FASTA format and one tabbed file containing position of annotation. The tabbed
file can be VCF, GFF3 or BED. 

Example of a simple command for a reference gemone fasta1, a target genome fasta2 and a
VCF file v1.vcf :

```
python transpo-rg.py -f1 fasta1.fasta -f2 fasta2.fasta -ti v1.vcf
```

The output file will only be a VCF file named v1\_out.vcf. It will go in directory
named result/ in your current directory. 
Both name are configurable with option so you can redirect output file in whatever
directory you want.
If you want to keep all the intermediate file like alignment file (in SAM), you can
use the option -n (--notempfile).


### Requires :  

* python (2.7+)
* bwa mem (0.7.15)
* bedtools (2.24.0)
* lib python :
    * biopython (1.65+) 
    * bedtools (2.24.0)
    * pybedtools (0.7.10)

#### Warning ! 

For pybedtools, replace the bedtools.py by the file found in tools/

----------
```
Version : 0.7.0
Last update : 23 Aug 2018

usage: transpo-rg.py [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-b FLANK] [-c]
                     [-d DIRECTORY] [-h] [-i] [-l] [-n] [-m MISMATCH] [-o OUT]
                     [-p] [-t TYPEA] [-u] [-v {0,1,2}] [-ver] [-w]

Required arguments :
  -f1 FASTA1, --fasta1 FASTA1
                        Input of genome fasta1 (reference) <fasta>.
  -f2 FASTA2, --fasta2 FASTA2
                        Input of genome fasta2 (target) <fasta>.
  -ti TABINPUT, --tabinput TABINPUT
                        Input of tabbed file related to fasta1 <bed/gff/vcf>

Optional arguments :
  -b FLANK, --flank FLANK
                        Size of flank region to extract from each side of the
                        annotation (default : 50).
  -c, --cds             Enable control of number and postions of CDS inside
                        mRNA (or gene). It generate a new tab file.
  -d DIRECTORY, --directory DIRECTORY
                        Name of the directory where files are generated
                        (default : result).
  -h, --help            Show this help message then exit.
  -i, --index           Create the index of fasta2 if it doesn't exist (-ii
                        for forcing).
  -l, --loss            Enable the creation of 'StatsLoss' file which contain
                        percentage of loss.
  -n, --notempfile      Create all file instead of using temporary file.
  -m MISMATCH, --mismatch MISMATCH
                        Specify the maximum of mismatch allowed (add option -p
                        if you want a percentage).
  -o OUT, --output OUT  Output file (same format of tabbed file input).
  -p, --percentage      Allows use of percentage for option --mismatch
  -t TYPEA, --type TYPEA
                        Only extract annotation of specified type (gff file
                        only) (example : "mRNA exon cds" (case insensitive)).
  -u, --update          Update the .version file by checking the git log
  -v {0,1,2}, --verbose {0,1,2}bbbbbbbbbbbbbbbbb
                        Change verbosity level.
  -ver, --version       Show version and date of last update then exit.
  -w, --warning         Disable warnings.
```

### Options :   
<br>
* #### --fasta1
    * It's the input of the genome which already have a tab file annotated.

* #### --fasta2 
    * It's the input of the genome for which we want to acquire the tab file. It will be
index by bwa.

* #### --tabintput
    * It's the tab file annotated related to genome fasta1. The format can be VCF, GFF or
BED.

* #### --flank
    * This option allow to control the lenght of flanking region that we select. It's 50
nucleotide by default.

* #### --cds 
    * Enable control of number and postions of CDS inside mRNA (or gene). It check is CDS
number 1 in his gene is still at number 1 after transfer.
It generate a new tab file.

* #### --directory
    * Allow to create a directory where output file will be redirect.

* #### --help
    * Show help message.

* #### --index
    * If just one -i is specify, the script will check if the index of fasta2 already exist,
else, it will be created. If -ii is specify, the index will be generate (even if it
already exist)

* #### --loss
    * Generate a statistic file which show difference of line between the first tab file
(tabinput) and the output tab file.

* #### --notempfile
    * Force the script to not use tempfile for secondary file like alignment file. The extra
files will be located in the output directory.

* #### --mismatch
    * Specify the number of mismatch allowed. If option --percentage is specify, the number
pass will be a percentage, else it's in number of nucleotide. Example : -m 10 will
accept 10 nucleotide max mismatch and -m 10 -p will accept max 10 % of the lenght of
the match of mismatch.

* #### --output 
    * Allow to specify a name for the output tab file.

* #### --percentage
    * Convert the option --mismatch to percentage. See above at --mismatch.

* #### --type
    * Specify the type to extract from a GFF (Example : cds, gene ,mrna). The other one are
ignored.

* #### --verbose
    * Change verbosity level. By default, verbose is set at 1 (recommended). 0 will show
only warnings and errors and 2 will show index and alignment log.

* #### --version
    * Show version and date of last update then exit.

* #### --warning
    * Disable the print of warnings message.

### Example of use :

If you have a reference genome named fasta1.fa with his 3 tab file v1.bed, v1.vcf and v1.gff
and you want to transfer those data to a other genome named fasta2.fa, the recommended
command line is :

* For the bed file :

        python transpo-rg.py -f1 fasta1.fa -f2 fasta2.fa -ti v1.bed -i -l
If you want to look to intermediate file, add option -n

* For the VCF file :

        python transpo-rg.py -f1 fasta1.fa -f2 fasta2.fa -ti v1.vcf -i -l

* For the GFF file :

    * Full file :

            python transpo-rg.py -f1 fasta1.fa -f2 fasta2.fa -ti v1.gff3 -i -l

    * Only mRNA and CDS :

            python transpo-rg.py -f1 fasta1.fa -f2 fasta2.fa -ti v1.gff3 -i -l -c -t "mrna,cds"

### Data for test :
<br>
Data and result are available in the cluster at the following path : /gs7k1/projects/GenomeHarvest/TransPo-RG/example\_data/


It's the first version of Sorghum Bicolor (Sbicolor\_79) and his tab file (VCF, BED and GFF3 found in example\_data/tab/) with the third version of Sorghum Bicolor (Sbicolor\_313).


All the file found in /example\_data/result/genecutGFF/ were generated by the following command line : 

        python transpo-rg.py -f1 example_data/genome/Sbicolor_79.assembly.fna -f2 example_data/genome/Sbicolor_313_v3.1.assembly.fna -ti example_data/tab/Sbicolor_79.gff3 -n -l -d genecutGFF -v 2 -c -t cds,gene
