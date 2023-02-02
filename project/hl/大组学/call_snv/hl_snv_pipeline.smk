configfile: "bam2vcf.yaml"

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_tumor(wildcards):
    return config["samples_tumor"][wildcards.sample_2]


rule all:
    input:
        #'pon/gatk4_mutect2_pon.vcf.gz',
        #expand("MuTect2_tumor/{sample_2}.vcf.gz", sample_2=config["samples_tumor"]),
        #expand("contamination_tumor/{sample_2}.contamination.table", sample_2=config["samples_tumor"])
        #expand("FilterMutectCalls_tumor/{sample_2}.filtered.vcf.gz", sample_2=config["samples_tumor"]),
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg38_multianno.vcf", sample_2=config["samples_tumor"]),
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg38_multianno.txt", sample_2=config["samples_tumor"])


#####################################
# UpStream
#####################################

rule bwa_map:
    input:
        ref=config["ref"],
        "{sample}/clean_fq/{sample}.L1-B1.R1_val_1.fq.gz",
        "{sample}/clean_fq/{sample}.L1-B1.R2_val_2.fq.gz"
    output:
        "{sample}/{sample}.mem.bam"
    params:
        rg=r"@RG\tID:{sample}.L1-B1\tSM:{sample}\tLB:L1-B1\tPL:ILLUMINA"
    shell:
        "bwa mem \
        -M -Y \
        -R '{params.rg}' \
        -t 16 \
        {input} | \
        samtools view -1 - > {output}"



rule MarkDuplicates:
    input:
        "{sample}/{sample}.mem.bam"
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


rule BaseRecalibrator:
    input:
        "{sample}/{sample}.markdup.bam"
    output:
        "{sample}/{sample}.bqsr_data.table"
    shell:
        "gatk BaseRecalibrator \
        -R database/hg19.fa \
        -I {input} \
        --use-original-qualities \
        -O {output} \
        --known-sites database/hapmap_3.3.hg19.sites.vcf \
        --known-sites database/1000G_omni2.5.hg19.sites.vcf \
        --known-sites database/dbsnp_138.hg19.sites.vcf \
        --known-sites database/1000G_phase1.indels.hg19.sites.vcf \
        --known-sites database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"


rule ApplyBQSR:
    input:
        bam = "{sample}/{sample}.markdup.bam",
        table = "{sample}/{sample}.bqsr_data.table"
    output:
        "{sample}/{sample}.bqsr.bam"
    shell:
        "gatk ApplyBQSR \
        -R database/hg19.fa \
        -I {input.bam} \
        -O {output} \
        -bqsr {input.table} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --create-output-bam-index true \
        --create-output-bam-md5 true \
        --use-original-qualities"




$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation  -I ${sample}_marked.bam  -O ${sample}_marked_fixed.bam  -SO coordinate  1>${sample}_log.fix 2>&1


samtools index ${sample}_marked_fixed.bam





#####################################
# Call
#####################################


rule HaplotypeCaller_tumor:
    input:
        bam_file=get_bam_file_tumor,
        ref=config["ref"],
        snp=config["snp"],
    output:
        "HaplotypeCaller_tumor/{sample_2}.vcf.gz"
    log:
        "logs/HaplotypeCaller_tumor/HaplotypeCaller_tumor_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  HaplotypeCaller "
        "-I {input.bam_file} " # v
        "-R {input.ref} " # v
        "--dbsnp {input.snp} "
        #"--max-mnp-distance 0 "
        "-O {output} > {log} 2>&1"



rule SelectVariants_tumor_snp:
    input:
        vcf="HaplotypeCaller_tumor/{sample_2}.vcf.gz"
    output:
        sel_vcf="SelectVariants_tumor/{sample_2}_snp.vcf.gz"
    log:
        "logs/SelectVariants_tumor/SelectVariants_tumor_{sample_2}_snp.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type SNP "
        "-V {input.vcf} "
        "-O {output.sel_vcf} "


rule SelectVariants_tumor_indel:
    input:
        vcf="HaplotypeCaller_tumor/{sample_2}.vcf.gz"
    output:
        sel_vcf="SelectVariants_tumor/{sample_2}_indel.vcf.gz"
    log:
        "logs/SelectVariants_tumor/SelectVariants_tumor_{sample_2}_indel.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type INDEL "
        "-V {input.vcf} "
        "-O {output.sel_vcf} "

# ===

rule VariantFiltration_tumor_snp:
    input:
        sel_vcf="SelectVariants_tumor/{sample_2}_snp.vcf.gz"
    output:
        fil_vcf="SelectVariants_tumor/{sample_2}_snp_filtered.vcf.gz"
    log:
        "logs/SelectVariants_tumor/SelectVariants_tumor_{sample_2}_snp_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} "


rule VariantFiltration_tumor_indel:
    input:
        sel_vcf="SelectVariants_tumor/{sample_2}_indel.vcf.gz"
    output:
        fil_vcf="SelectVariants_tumor/{sample_2}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants_tumor/SelectVariants_tumor_{sample_2}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} "


rule MergeVcfs_tumor:
    input:
        fil_vcf_indel="SelectVariants_tumor/{sample_2}_indel_filtered.vcf.gz",
        fil_vcf_snp="SelectVariants_tumor/{sample_2}_snp_filtered.vcf.gz"
    output:
        merge_vcf="SelectVariants_tumor/{sample_2}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants_tumor/SelectVariants_tumor_{sample_2}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' MergeVcfs "
        "-I {input.fil_vcf_indel} "
        "-I {input.fil_vcf_snp} "
        "-O {output.merge_vcf}"



##############################
# ANNOVAR
##############################

rule ANNOVAR_tumor:
    input:
        vcf="SelectVariants_tumor/{sample_2}_indel_filtered.vcf.gz",
        # contamin='contamination_tumor/{sample_2}.contamination.table',   
    output:
        multianno_vcf='ANNOVAR_tumor/{sample_2}.filtered.hg38_multianno.vcf',
        multianno_txt='ANNOVAR_tumor/{sample_2}.filtered.hg38_multianno.txt'
    params:
        "ANNOVAR_tumor/{sample_2}.filtered"
    log:
        "logs/ANNOVAR_tumor/ANNOVAR_tumor_{sample_2}.log"
    shell:
        "/biosoft/annovar/table_annovar.pl "
        "{input.vcf} "
        "/biosoft/annovar/humandb/ "
        "--buildver hg38 "
        "--protocol refGene,exac03,gnomad_genome " # refGene
        "--operation g,f,f -vcfinput "
        "--outfile {params} > {log} 2>&1"

# 有些参数和师兄的不一样，再看看
#




#time /media/whl/CRISPRAnalyzeR/biology/gatk-4.0.8.1/gatk SelectVariants -select-type SNP -V ${sample}_raw.vcf -O ${sample}_raw.snp.vcf
#time /media/whl/CRISPRAnalyzeR/biology/gatk-4.0.8.1/gatk SelectVariants -select-type INDEL -V ${sample}_raw.vcf -O ${sample}_raw.indel.vcf


time /media/whl/CRISPRAnalyzeR/biology/gatk-4.0.8.1/gatk VariantFiltration \
    -V ${sample}_raw.indel.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ${sample}_raw.indel.filter.vcf

time /media/whl/CRISPRAnalyzeR/biology/gatk-4.0.8.1/gatk VariantFiltration \
    -V ${sample}_raw.snp.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ${sample}_raw.snp.filter.vcf

time /media/whl/CRISPRAnalyzeR/biology/gatk-4.0.8.1/gatk MergeVcfs \
    -I ${sample}_raw.snp.filter.vcf \
    -I ${sample}_raw.indel.filter.vcf \
    -O ${sample}.filter.vcf



















