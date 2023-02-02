# https://www.jianshu.com/p/7a3de6b8e503 质控流程
# https://blog.csdn.net/weixin_39691233/article/details/111703829 整个流程

rule trim_galore:
    input:
        "{sample}/{sample}.L1-B1.R1.fastq.gz",
        "{sample}/{sample}.L1-B1.R2.fastq.gz"
    output:
        "{sample}/clean_fq/{sample}.L1-B1.R1_val_1.fq.gz",
        "{sample}/clean_fq/{sample}.L1-B1.R2_val_2.fq.gz"
    shell:
        "trim_galore --paired --retain_unpaired -q 20 --phred33 --length 30 --stringency 3 --gzip --cores 4 -o {wildcards.sample}/clean_fq {input}"



rule bwa_map:
    input:
        "database/hg19.fa",
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












