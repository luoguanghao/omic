# Star

configfile: "star_pe.yaml"

def get_fq_file_1(wildcards):
    return config["samples"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]

# nohup snakemake -p --keep-going --snakefile ./star_pe.smk  --cores 8 > 1.log &
# snakemake -np --snakefile ./star_pe.smk
# snakemake --snakefile ./star_pe.smk --unlock

rule all:
    input:
        expand("result/STAR/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"])

#rule STAR_genomeGenerate:
#    input:
#        fa_ref=config["ref"], #<<<
#        gtf_ref=config["gtf"],
#    output:
#        "star_index/"
#    log:
#        "logs/index_1pass.log"
#    shell:
#        """STAR --runMode genomeGenerate \
#        --runThreadN 4 \
#        --genomeDir ./star_index/  \
#        --genomeFastaFiles {input.fa_ref}  \
#        --sjdbGTFfile {input.gtf_ref} > {log} 2>&1"""


rule STAR:
    input:
        fq_file1=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_1.fq.gz",sample=wildcards.sample),
        fq_file2=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_2.fq.gz",sample=wildcards.sample),
        index_1pass=config["star_index"],
        gtf=config["gtf"]
    output:
        star_sam="result/STAR/{sample}_Aligned.sortedByCoord.out.bam",
    
    params:
        outFileNamePrefix_1pass="result/STAR/{sample}_",
    
    log:
        "logs/STAR/STAR_{sample}.log"
    threads: 4
    shell:
        """STAR --genomeDir {input.index_1pass} \
            --readFilesIn {input.fq_file_1} {input.fq_file_2} \
            --runThreadN 4 \
            --outFilterMultimapScoreRange 1 \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 10 \
            --alignIntronMax 500000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 2 \
            --alignSJDBoverhangMin 1 \
            --genomeLoad NoSharedMemory \
            --limitBAMsortRAM 70000000000 \
            --readFilesCommand zcat \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterScoreMinOverLread 0.33 \
            --sjdbOverhang 100 \
            --outSAMstrandField intronMotif \
            --outSAMattributes NH HI NM MD AS XS \
            --sjdbGTFfile {input.gtf} \
            --limitSjdbInsertNsj 2000000 \
            --outSAMunmapped None \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMheaderHD @HD VN:1.4 \
            --twopassMode Basic \
            --outSAMmultNmax 1 \
            --outFileNamePrefix {params.outFileNamePrefix_1pass} > {log} 2>&1"""




