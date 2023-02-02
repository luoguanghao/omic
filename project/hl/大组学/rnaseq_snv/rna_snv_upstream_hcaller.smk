# RNAseq Call SNV
# 
#
configfile: "rnaseq_snv.yaml"

def get_fq_file_1(wildcards):
    return config["fq_file_1"][wildcards.fq_1]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]



rule all:
    input:
        expand("VariantFiltration/{sample}_filtered.vcf.gz", sample=config["fq_file_1"]),


rule STAR_genomeGenerate:
    input:
        fa_ref=config["ref"], #<<<
        gtf_ref=config["gtf"],
    output:
        "{sample}/{sample}.star.bam"
    log:
        "logs/{sample}/SplitNCigarReads_{sample}.log"
    shell:
        "STAR --runMode genomeGenerate 
        --runThreadN 2
        --genomeDir ./star_index/
        --genomeFastaFiles {input.fa_ref}
        --sjdbGTFfile {input.gtf_ref}"


rule STAR:
    input:
        fq_file_1=get_fq_file_1,
        #fq_file_2=get_fq_file_2,
        fa_ref=config["ref"], #<<<
        gtf_ref=config["gtf"],
    output:
        star_1="start_1pass/{sample}_Aligned.out.sam"
        genomeGenerate_2="star_index_2pass/"
        star_2="start_2pass/{sample}_Aligned.out.sam"
    log:
        "logs/{sample}/STAR_{sample}.log"
    shell:
        "STAR --runThreadN 2 \
        --genomeDir ./star_index/ \
        --readFilesIn {fq_file_1} \
        --readFilesCommand zcat \
        --outFileNamePrefix ./star_1pass/{sample}_
        STAR --runMode genomeGenerate \
        --runThreadN 2 \
        --genomeDir ./star_index_2pass/ \
        --genomeFastaFiles {input.fa_ref} \
        --sjdbFileChrStartEnd ./star_1pass/{sample}_SJ.out.tab
        STAR --runThreadN 2 \
        --genomeDir ./star_index_2pass/ \
        --readFilesIn {fq_file_1} \
        --readFilesCommand zcat \
        --outFileNamePrefix ./star_2pass/{sample}_"



# =================

rule AddOrReplaceReadGroups:
    input:
        sam_file="start_2pass/{sample}_Aligned.out.sam"
    output:
        "{sample}/{sample}.rg_added_sorted.bam"
    log:
        "logs/{sample}/AddOrReplaceReadGroups_{sample}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' AddOrReplaceReadGroups "
        "I = {sam_file} "
        "O = {output} "
        "SO = coordinate "
        "RGID = {sample} "
        "RGLB = rna "
        "RGPL = illumina "
        "RGPU = novo "
        "RGSM = {sample} "



# =================



rule MarkDuplicates:
    input:
        "{sample}/{sample}.rg_added_sorted.bam"
    output:
        bam = "{sample}/{sample}.markdup.bam",
        txt = "{sample}/{sample}.markdup_metrics.txt"
    params:
        aso = r"queryname",
        so = r"coordinate"
    shell:
        "gatk MarkDuplicates \
        --INPUT {input} \
        --OUTPUT /dev/stdout \
        --METRICS_FILE {output.txt} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER '{params.aso}' \
        --CREATE_MD5_FILE false \
        --REMOVE_DUPLICATES false | \
        gatk SortSam \
        --INPUT /dev/stdin \
        --OUTPUT {output.bam} \
        --SORT_ORDER '{params.so}' \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false"

# =================


rule SplitNCigarReads:
    input:
        bam_file="{sample}/{sample}.markdup.bam",
        ref=config["ref"],
    output:
        "{sample}/{sample}.SplitNCigarReads.bam"
    log:
        "logs/{sample}/SplitNCigarReads_{sample}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' SplitNCigarReads "
        "-R {input.ref}"
        "-I {input.bam_file} "
        "-O {output}"


# =================


rule BaseRecalibrator:
    input:
        "{sample}/{sample}.SplitNCigarReads.bam"
    output:
        "{sample}/{sample}.bqsr_data.table"
    shell:
        "gatk BaseRecalibrator \
        -R database/hg19.fa \
        -I {input} \
        --use-original-qualities \
        -O {output} \
        --known-sites /gatk/my_dir/reference/gatk_bundle/hg38/hapmap_3.3.hg19.sites.vcf \ ##################################
        --known-sites /gatk/my_dir/reference/gatk_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \  
        --known-sites /gatk/my_dir/reference/gatk_bundle/hg38/dbsnp_146.hg38.vcf.gz \
        --known-sites /gatk/my_dir/reference/gatk_bundle/hg38/1000G_phase1.indels.hg19.sites.vcf \   ##################################
        --known-sites /gatk/my_dir/reference/gatk_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"  ##################################


# =================


rule ApplyBQSR:
    input:
        bam_file="{sample}/{sample}.SplitNCigarReads.bam",
        bqsr_tbl="{sample}/{sample}.bqsr_data.table",
        ref=config["ref"],
    output:
        "{sample}/{sample}.BQSR.vcf.gz"
    log:
        "logs/ApplyBQSR/ApplyBQSR_{sample}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' ApplyBQSR "
        "--add-output-sam-program-record "
        "-R {input.ref} "
        "-I {input.bam_file} "
        "--use-original-qualities "
        "--bqsr-recal-file ${input.bqsr_tbl}"
        "--CREATE_INDEX true" # <<<<<<<<<<< index!!!!!!!
        "-O {output} "

# =================

rule HaplotypeCaller:
    input:
        bam_file="{sample}/{sample}.BQSR.vcf.gz",
        ref=config["ref"],
        snp=config["snp"],
    output:
        "HaplotypeCaller/{sample}.vcf.gz"
    log:
        "logs/HaplotypeCaller/HaplotypeCaller_{sample}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  HaplotypeCaller "
        "-I {input.bam_file} " # v
        "-R {input.ref} " # v
        "--dbsnp {input.snp} "
        #"--max-mnp-distance 0 "
        "-O {output} > {log} 2>&1"

# =================


rule MergeVcfs:
    input:
        # multi sample
        vcf_files = expand('HaplotypeCaller/{sample}.vcf.gz', sample=config["samples"]),
    output:
        merge_vcf="MergeVcfs/{sample}_merged.vcf.gz"
    params:
        " -I ".join(expand("HaplotypeCaller/{sample}.vcf.gz", sample=config["samples"]))
    log:
        "logs/merge_vcf/merge_vcfs_{sample}_merged.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' MergeVcfs "
        "-I {params} "
        "-O {output.merge_vcf}"


# =================


rule VariantFiltration:
    input:
        vcf="MergeVcfs/{sample}_merged.vcf.gz""
    output:
        fil_vcf="VariantFiltration/{sample}_filtered.vcf.gz"
    log:
        "logs/VariantFiltration/VariantFiltration_{sample}_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.vcf} "
        "--window 35 "
        '--filter-name "FS" '
        '--filter "FS > 30.0" '
        '--filter-name "FS" '
        '--filter "QD < 2.0" '        
        "-O {output.fil_vcf} "






