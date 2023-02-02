# MuTect2
# Create PoN
# MuTect2 no matched mode
# 
# 在GenomicsDBImport的-L参数中，需要指定interval，这里指定 -L 1 -L 2 -L 3 -L 4 ...
# 
# gatk snv流程的软件基本都是单线程的
# 

configfile: "bam2vcf.yaml"

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_tumor(wildcards):
    return config["samples_tumor"][wildcards.sample_2]

# ====

rule all:
    input:
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.vcf", sample_2=config["samples_tumor"]),
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.txt", sample_2=config["samples_tumor"])



##############################
# DownsampleSam & Pooling
##############################

rule DownsampleSam:
    input:
        bam_file=get_bam_file_tumor,
    output:
        bam_file='downsampling_bam/{sample_1}.downsampling.bam'
    log:
        "logs/downsampling_bam/downsampling_bam_{sample_1}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  DownsampleSam  "
        "-I {input.bam_file} "
        "-O {output.bam_file} "
        "-P 0.5 > {log} 2>&1"


# /gatk/gatk DownsampleSam \
# I=/gatk/top/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/ctrl/bam/HW3.final.bam \
# O=HW3.downsampling.bam \
# P=0.5



rule MergeSamFiles:
    input:
        bam_file=expand('downsampling_bam/{sample_1}.downsampling.bam', sample_1=config["samples_ctrl"]),
        # bam_file='downsampling_bam/{sample_1}.downsampling.bam'
    output:
        merge_bam_file='MergeSamFiles/merge_sam.downsampling.bam',
        merge_bai_file='MergeSamFiles/merge_sam.downsampling.bai'
    log:
        "logs/MergeSamFiles.log"
    params:
        " -I ".join(expand('downsampling_bam/{sample_1}.downsampling.bam', sample_1=config["samples_ctrl"]))    
    shell:
        "/gatk/gatk --java-options '-Xmx40G' DownsampleSam "
        "-I {params} "
        "-MSD true "
        "-CREATE_INDEX true "
        "-O {output.merge_bam_file} > {log} 2>&1"








############################################
############################################

##############################
# MuTect2 with pooling match
##############################
# 需要把GRCh37.preprocessed.interval_list挪到目录下

# tumor sample
rule MuTect2_tumor_pooling:
    input:
        ref=config["ref"],
        bam_file=get_bam_file_tumor,
        interval_list="GRCh37.preprocessed.interval_list",  # <<<<<<<
        pooling_bam_file='MergeSamFiles/merge_sam.downsampling.bam'
        gr="/gatk/my_dir/Download/gatk_bundle/gn/af-only-gnomad.raw.sites.b37.vcf.gz"
    output:
        vcf="MuTect2_tumor/{sample_2}.vcf.gz",
        f1r2="MuTect2_tumor/{sample_2}.f1r2.tar.gz"
    log:
        "logs/MuTect2_tumor/MuTect2_tumor_{sample_2}.log"    
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  Mutect2 "
        "-R {input.ref} "
        "-L {input.interval_list} "
        "-I {input.bam_file} "
        "-I {input.pooling_bam_file} "
        "-normal  merge_sam.downsampling "
        "--germline-resource {input.gr} "
        #"--max-mnp-distance 0 "
        "--f1r2-tar-gz {output.f1r2} "
        "-O {output.vcf} > {log} 2>&1"




##############################
# Pileup & Contamination
##############################

rule GetPileupSummaries:
    input:
        bam_file=get_bam_file_tumor,
        v_file="/gatk/my_dir/Download/gatk_bundle/gn/small_exac_common_3_b37.vcf.gz", ##
    output:
        'pileups_tumor/{sample_2}.pileups.table'
    log:
        "logs/GetPileupSummaries/GetPileupSummaries_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  GetPileupSummaries  "
        "-I {input.bam_file} "
        "-V {input.v_file} "
        "-L {input.v_file} "
        "-O {output} > {log} 2>&1"


rule CalculateContamination:
    input:
        'pileups_tumor/{sample_2}.pileups.table'
    output:
        ct='contamination_tumor/{sample_2}.contamination.table',
        st='contamination_tumor/{sample_2}.segments.table'
    log:
        "logs/CalculateContamination/CalculateContamination_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' CalculateContamination "
        "-I {input} "
        "--tumor-segmentation {output.st} "
        "-O {output.ct} > {log} 2>&1"



##############################
# LearnReadOrientationModel
##############################
# gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

rule LearnReadOrientationModel:
    input:
        "MuTect2_tumor/{sample_2}.f1r2.tar.gz"
    output:
        'LearnReadOrientationModel/{sample_2}.read-orientation-model.tar.gz'
    log:
        "logs/LearnReadOrientationModel/LearnReadOrientationModel_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  LearnReadOrientationModel "
        "-I {input} "
        "-O {output} > {log} 2>&1"



##############################
# Filter
##############################


rule FilterMutectCalls_tumor:
    input:
        vcf="MuTect2_tumor/{sample_2}.vcf.gz",
        ref=config["ref"],
        contamin='contamination_tumor/{sample_2}.contamination.table',
        segment='contamination_tumor/{sample_2}.segments.table',
        orientation='LearnReadOrientationModel/{sample_2}.read-orientation-model.tar.gz'
    output:
        'FilterMutectCalls_tumor/{sample_2}.filtered.vcf.gz'
    log:
        "logs/FilterMutectCalls_tumor/FilterMutectCalls_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  FilterMutectCalls  "
        "-R {input.ref} "
        "--tumor-segmentation {input.segment} "
        "--contamination-table {input.contamin} "
        "--ob-priors {input.orientation} "
        "-V {input.vcf} "
        "-O {output} > {log} 2>&1"


##############################
# add chr with use of bcftools
##############################

rule add_chr_tumor:
    input:
        chr_name_change="chr_name_change.txt",
        vcf='FilterMutectCalls_tumor/{sample_2}.filtered.vcf.gz',
        # contamin='contamination_tumor/{sample_2}.contamination.table',   
    output:
        'add_chr_tumor/{sample_2}.addchr.filtered.vcf.gz'
    log:
        "logs/add_chr_tumor/addchr_{sample_2}.log"
    shell:
        "bcftools annotate "
        "--rename-chrs {input.chr_name_change} "
        "{input.vcf} | "
        "bgzip -c > {output} 2> {log}"

# bcftools annotate --rename-chrs chr_name_change.txt ./FilterMutectCalls/CHX20170528.filtered.vcf.gz | bgzip -c > CHX20170528.addchr.filtered.vcf.gz


##############################
# ANNOVAR
##############################

rule ANNOVAR_tumor:
    input:
        vcf_addchr='add_chr_tumor/{sample_2}.addchr.filtered.vcf.gz',
        # contamin='contamination_tumor/{sample_2}.contamination.table',   
    output:
        hg19_multianno_vcf='ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.vcf',
        hg19_multianno_txt='ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.txt'
    params:
        "ANNOVAR_tumor/{sample_2}.filtered"
    log:
        "logs/ANNOVAR_tumor/ANNOVAR_tumor_{sample_2}.log"
    shell:
        "/biosoft/annovar/table_annovar.pl "
        "-vcfinput {input.vcf_addchr} "
        "/biosoft/annovar/humandb/ "
        "--buildver hg19 "
        "--protocol refGene,knownGene,clinvar_20210501 " # refGene
        "--operation g,g,f "
        "--outfile {params} > {log} 2>&1"



