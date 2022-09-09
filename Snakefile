# Snakemake workflow for GATK4 snps g.vcf
# Guide - ## guide https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

# Samples 
SAMPLES = ["test1", "test3"]

# final rule which will be the merged SNPs file (the complete cohort)
rule all:
    input:
        "/Users/tsewell/gatk_snake/merged/cohort.g.vcf.gz"

'''
rule bwa_map:
    input:
        reads = expand("path/to/reads/{sample}.fastq, sample=SAMPLES"),
        ref = "/Users/tsewell/gatk_snake/ref/ref.fa"
    output:
        bam = "/Users/tsewell/gatk_snake/bams/{sample}.bam"
    conda:
        "/Users/tsewell/bwa.yaml"
    shell: 
        "bwa mem {input.reads} {input.ref} | samtools view -Sb > {output.bam}"

# INDEX FILES NEEDED - EITHER WE NEED A RULE HERE OR WE CAN INCLUDE IN HAPLOTYPE RULE
rule index_bams: 
    input: 
        bam = "/Users/tsewell/gatk_snake/bams/{sample}"
    output:
        index = "/Users/tsewell/gatk/bams/{sample}.bai"
    shell:
        "samtools index {input.bam} -o ${output.index}
'''

# Haplotype caller to call snps with GVCF flag (-ERC) 
rule haplotype:
    input:
        bam = "/Users/tsewell/gatk_snake/bams/{sample}.bam",
        ref = "/Users/tsewell/gatk_snake/ref/ref.fa"
    output:
        vcf = "/Users/tsewell/gatk_snake/vcfs/{sample}.g.vcf.gz"
    conda: 
        "/Users/tsewell/gatk_snake/gatk4.yaml"
    shell: 
        "gatk --java-options '-Xmx32g' HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} --tmp-dir /Users/tsewell/gatk_snake/tmp -ERC GVCF"

# merge files into a single cohort VCF
rule merge:
    input:
        vcfs = expand("/Users/tsewell/gatk_snake/vcfs/{sample}.g.vcf.gz", sample=SAMPLES),
        ref = "/Users/tsewell/gatk_snake/ref/ref.fa"
    output:
        glist = "/Users/tsewell/gatk_snake/merged/list/gvcfs.list",
        cohort = "/Users/tsewell/gatk_snake/merged/cohort.g.vcf.gz"
    conda: 
        "/Users/tsewell/gatk_snake/gatk4.yaml"
    shell:
        "ls /Users/tsewell/gatk_snake/vcfs/*.vcf.gz > /Users/tsewell/gatk_snake/merged/list/gvcfs.list | " 
        "gatk CombineGVCFs -R {input.ref} --variant {output.glist} --tmp-dir /Users/tsewell/gatk_snake/tmp -O {output.cohort}"

#ls gvcfs/*.vcf >gvcfs.list!!!!!!
#gatk --java-options "-Xmx180G -XX:ParallelGCThreads=36" CombineGVCFs -R $ref --variant gvcfs.list --dbsnp $DBSNP -O #combined_gvcf.vcf 

###### PROBLEM!!! THERE NEEDS TO BE A --variant for each vcf input files ####