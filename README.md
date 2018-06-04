### Requires :  

* bwa mem (0.7.15)
* bedtools (2.24.0)
* python 3 (3.4.3)
    * biopython (1.65+) 
    * bedtools (2.24.0)
    * pybedtools (0.7.10)


----------


``` python3
usage: tab1ToRef2.py [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-b FLANK] [-c]
                     [-h] [-i] [-n] [-o OUT] [-t TYPEA] [-v] [-ver] [-w]

Required arguments:
  -f1 FASTA1, --fasta1 FASTA1
                        Input of reference fasta1 (old) <fasta>.
  -f2 FASTA2, --fasta2 FASTA2
                        Input of reference fasta2 (new) <fasta>.
  -ti TABINPUT, --tabinput TABINPUT
                        Input of tabbed file related to fasta1 <bed/gff/vcf>

Optional arguments:
  -b FLANK              Size of flank region to extract from each side of the
                        annotation (default : 50).
  -c, --cds             Enable control of postions of CDS inside mRNA (or
                        gene).
  -h, --help            Show this help message then exit.
  -i, --index           Create the index of fasta2.
  -n, --notempfile      Create all file instead of using temporary file.
  -o OUT, --output OUT  Output file (same format of tabbed file input).
  -t TYPEA, --type TYPEA
                        Only extract annotation of specified type (gff file
                        only) (example : "mRNA exon cDs" (case insensitive)).
  -v, --verbose         Enable message.
  -ver, --version       Show version and date of last update then exit.
  -w, --warning         Disable warnings.
```
