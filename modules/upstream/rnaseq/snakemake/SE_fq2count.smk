
# RNAseq fq2count SE 待补充，等分析到再说

configfile: "PE_fq2count.yaml"

#def get_fq_file_1(wildcards):
#    return config["fq_file_1"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]




rule all:
    input:
        expand("result/star_2pass/{sample}/{sample}_Aligned.sortedByCoord.out.bai", sample=config["samples"]),
        expand("result/VariantFiltration/{sample}.filter.vcf.gz", sample=config["samples"])



rule hisat2:
    input:
        fq_file1=lambda wildcards: expand("/gatk/top/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_1.fq.gz",sample=wildcards.sample),
        fq_file2=lambda wildcards: expand("/gatk/top/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_2.fq.gz",sample=wildcards.sample),
        index=config["hisat2_index"], #<<<
    output:
        sam_file="result/{sample}/{sample}.sam"
    log:
        "logs/hisat2/{sample}.log"        
    threads: 8
    shell:
        """hisat2 -p 8 -t \
        -x {input.index} \
        -1 {input.fq_file1} -2 {input.fq_file2} \
        -S {output.bam_file} 1> {log} 2>&1"""



rule samtools:
    input:
        sam_file="result/{sample}/{sample}.sam"
    output:
        bam_file="result/{sample}/{sample}.bam"
        sorted_bam_file="result/{sample}/{sample}.sorted.bam"
        bai_file="result/{sample}/{sample}.bai"
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
        "logs/samtools/{sample}.log"
    threads: 8
    shell:
        """featureCounts \
        -T 8 -p -O -g gene_name \
        -a {input.gtf} \
        -o {output.fc_file} \
        {params}
        """




