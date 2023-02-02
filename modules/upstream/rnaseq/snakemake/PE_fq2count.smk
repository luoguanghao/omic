# RNAseq fq2count PE\
# nohup snakemake -p --keep-going --snakefile ./PE_fq2count.smk   --cores 16 > 1.log &
# snakemake -np --snakefile ./PE_fq2count.smk
# snakemake --snakefile ./PE_fq2count.smk --unlock
# 
## QC 要删掉
## hisat2
## samtools
## featureCount


configfile: "PE_fq2count.yaml"

#def get_fq_file_1(wildcards):
#    return config["fq_file_1"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]



rule all:
    input:
        expand("result/{sample}/QC/{sample}_1_fastqc.html", sample=config["samples"]), ### <---- 删掉
        expand("result/{sample}/QC/{sample}_2_fastqc.html", sample=config["samples"]), ### <---- 删掉
        expand("result/{sample}/{sample}.sorted.bai", sample=config["samples"]),
        expand("result/{sample}/{sample}.sorted.bam", sample=config["samples"]),
        featurecount_res='result/fc.txt'


rule fastqc:  ### <---- 删掉
    input:
        fq_file1=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_1.fq.gz",sample=wildcards.sample),
        fq_file2=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_2.fq.gz",sample=wildcards.sample),
    output:
        qc="result/{sample}/QC/{sample}_1_fastqc.html"
    params:
        index=config["hisat2_index"], #<<<
    log:
        "logs/QC/{sample}.log"        
    threads: 8
    shell:
        """fastqc -t 8 \
        {output.qc} \
        {input.fq_file1} {input.fq_file2} 1> {log} 2>&1"""


rule hisat2:
    input:
        fq_file1=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_1.fq.gz",sample=wildcards.sample),
        fq_file2=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_2.fq.gz",sample=wildcards.sample),
    output:
        sam_file=temp("result/{sample}/{sample}.sam")
    params:
        index=config["hisat2_index"], #<<<
    log:
        "logs/hisat2/{sample}.log"        
    threads: 8
    shell:
        """hisat2 -p 8 -t \
        -x {params.index} \
        -1 {input.fq_file1} -2 {input.fq_file2} \
        -S {output.sam_file} 1> {log} 2>&1"""

rule samtools:
    input:
        sam_file="result/{sample}/{sample}.sam"
    output:
        bam_file=temp("result/{sample}/{sample}.bam"),
        sorted_bam_file=temp("result/{sample}/{sample}.sorted.bam"),
        bai_file=temp("result/{sample}/{sample}.sorted.bai"),
    log:
        "logs/samtools/{sample}.log"
    threads: 8
    shell:
        """samtools view  -S {input.sam_file} -b > {output.bam_file}
        samtools sort {output.bam_file} -o {output.sorted_bam_file} -@ 8 1>> {log} 2>&1
        samtools index {output.sorted_bam_file} {output.bai_file} -@ 8 1>> {log} 2>&1"""
        
        
rule featurecount:
    input: # 批量给，参考PoN的代码
        sorted_bam_files = expand("result/{sample}/{sample}.sorted.bam", sample=config["samples"]),
        gtf=config["gtf"],
    output:
        fc_file='result/fc.txt',
    params:    
        " ".join(expand("result/{sample}/{sample}.sorted.bam", sample=config["samples"]))
    log:
        "logs/featurecount.log"
    threads: 8
    shell:
        """featureCounts \
        -T 8 -p -O -g gene_name \
        -a {input.gtf} \
        -o {output.fc_file} {params} 1> {log} 2>&1"""

        
        
