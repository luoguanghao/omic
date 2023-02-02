# Star Fusion


# configfile: "rna_snv_upstream_hcaller.yaml"

def get_fq_file_1(wildcards):
    return config["fq_file_1"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]



rule all:
    input:
        "result/STAR/{sample}_Aligned.out.sam"


rule STAR_genomeGenerate:
    input:
        fa_ref=config["ref"], #<<<
        gtf_ref=config["gtf"],
    output:
        "star_index/"
    log:
        "logs/index_1pass.log"
    shell:
        """STAR --runMode genomeGenerate \
        --runThreadN 2 \
        --genomeDir ./star_index/  \
        --genomeFastaFiles {input.fa_ref}  \
        --sjdbGTFfile {input.gtf_ref} > {log} 2>&1"""


rule STAR:
    input:
        fq_file_1=get_fq_file_1,
        #fq_file_2=get_fq_file_2,
        index_1pass="star_index/"
    output:
        star_sam="result/STAR/{sample}_Aligned.out.sam",
    params:
        outFileNamePrefix_1pass="result/STAR/{sample}_",
    log:
        "logs/STAR/STAR_{sample}.log"
    shell:
        """STAR --runThreadN 4 \
        --genomeDir {input.index_1pass} \
        --readFilesIn {input.fq_file_1} \
        --readFilesCommand "gunzip -c" \
        --twopassMode Basic \
        --outReadsUnmapped None \
        --outSAMunmapped Within \
        --outSAMstrandField intronMotif \
        --outFileNamePrefix {params.outFileNamePrefix_1pass} \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:GRPundef \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 > {log} 2>&1"""




rule do_samtools:
    input:
        star_sam="result/STAR/{sample}_Aligned.out.sam",
    output:
        bai="result/STAR/{sample}.bai",
    params:
        bam="result/STAR/{sample}.bam",
        sort_bam="result/STAR/{sample}.bam"
    log:
        "logs/STAR/samtools_{sample}.log"
    shell:
        """samtools view -S {input.star_sam} -b -o {params.bam} \
        samtools sort -@ 4 -o {params.sort_bam} {params.bam} > {log} 2>&1 \
        samtools index {params.sort_bam} {output.bai} >> {log} 2>&1"""

# =================





