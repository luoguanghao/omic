# Star

configfile: "star_se.yaml"

def get_fq_file_1(wildcards):
    return config["fq_file_1"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]

# nohup snakemake -p --keep-going --snakefile ./star_se.smk --rerun-incomplete  --cores 8 > 1.log &
# snakemake -np --snakefile ./star_se.smk
# snakemake --snakefile ./star_se.smk --unlock

rule all:
    input:
        expand("result/STAR/{sample}_Aligned.sortedByCoord.out.bam", sample=config["fq_file_1"])


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
        fq_file_1=get_fq_file_1,
        #fq_file_2=get_fq_file_2,
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
            --readFilesIn {input.fq_file_1} \
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
            --readFilesCommand cat \
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
            --readFilesCommand zcat \            
            --outFileNamePrefix {params.outFileNamePrefix_1pass} > {log} 2>&1"""





