
#rule scriptEntete:
#    input :
#        "{sample1}.fasta"
#    shell :
#        "python chrParse.py {input}"

rule flankingReg :
    input :
        "{tab}.bed",
        "{sample1}.fasta"
    output :
        "{tab}_{sample1}.bed",
    shell :
        "python bedFlank.py {input} {output.tab}"
    #script :
    #    "bedFlank.py {input} {output}"
    
rule getfasta:
    input :
        "{sample1}.fasta"
    output :
        "{sample1}_selected.fasta"
    shell :
        "bedtools getfasta -fi {input} -bed rules.flankingReg.output -fo {output}"

#rule bwa_index :
    
#rule bwa_mem :
 #   input :
  #      "{sample1}_selected.fasta",
   #     "{sample2}.fasta"
   # output :
   #     "{sample1}.bam"
   # shell :
   #     "bwa 
