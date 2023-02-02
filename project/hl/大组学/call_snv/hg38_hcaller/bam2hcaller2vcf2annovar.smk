# HaplotypeCaller
# dbsnp_138.hg19.vcf.gz should be dropped the 'chr'
# should be index
## snp file: dbsnp_138.hg19.drop_chr.vcf.gz,dbsnp_138.hg19.drop_chr.vcf.gz.tbi


#configfile: "/gatk/my_dir/my_upstream_project/hl_call_cnv/gatk_sCNV/try_again_4/bam2PoN2Normalized2Seg.yaml"
#
# nohup snakemake -p --snakefile ./bam2hcaller2vcf2annovar.smk --cores 12 > all.log &
#
configfile: "bam2hcaller2vcf2annovar.yaml"

def get_bam_file(wildcards):
    return config["samples"][wildcards.sample]




rule all:
    input:
        expand("ANNOVAR/{sample}.filtered.hg38_multianno.vcf", sample=config["samples"]),
        expand("ANNOVAR/{sample}.filtered.hg38_multianno.txt", sample=config["samples"])



#####################################
# Call
#####################################


rule HaplotypeCallerl:
    input:
        bam_file=get_bam_file,
        ref=config["ref"],
        snp=config["snp"],
    output:
        "HaplotypeCaller/{sample}.vcf.gz"
    log:
        "logs/HaplotypeCaller/HaplotypeCaller_{sample}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  HaplotypeCaller "
        #"-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 "
        #"-L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY "
        "-I {input.bam_file} " # v
        "-R {input.ref} " # v
        "--dbsnp {input.snp} "
        #"--max-mnp-distance 0 "
        "-O {output} > {log} 2>&1"



rule SelectVariants_snp:
    input:
        vcf="HaplotypeCaller/{sample}.vcf.gz"
    output:
        sel_vcf="SelectVariants/{sample}_snp.vcf.gz"
    log:
        "logs/SelectVariants/SelectVariants_{sample}_snp.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type SNP "
        "-V {input.vcf} "
        "-O {output.sel_vcf} > {log} 2>&1"


rule SelectVariants_indel:
    input:
        vcf="HaplotypeCaller/{sample}.vcf.gz"
    output:
        sel_vcf="SelectVariants/{sample}_indel.vcf.gz"
    log:
        "logs/SelectVariants/SelectVariants_{sample}_indel.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type INDEL "
        "-V {input.vcf} "
        "-O {output.sel_vcf} > {log} 2>&1"

# ===

rule VariantFiltration_snp:
    input:
        sel_vcf="SelectVariants/{sample}_snp.vcf.gz"
    output:
        fil_vcf="SelectVariants/{sample}_snp_filtered.vcf.gz"
    log:
        "logs/SelectVariants/SelectVariants_{sample}_snp_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.sel_vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} > {log} 2>&1"


rule VariantFiltration_indel:
    input:
        sel_vcf="SelectVariants/{sample}_indel.vcf.gz"
    output:
        fil_vcf="SelectVariants/{sample}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants/SelectVariants_{sample}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.sel_vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} > {log} 2>&1"


rule MergeVcfs:
    input:
        fil_vcf_indel="SelectVariants/{sample}_indel_filtered.vcf.gz",
        fil_vcf_snp="SelectVariants/{sample}_snp_filtered.vcf.gz"
    output:
        merge_vcf="SelectVariants/{sample}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_{sample}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' MergeVcfs "
        "-I {input.fil_vcf_indel} "
        "-I {input.fil_vcf_snp} "
        "-O {output.merge_vcf} > {log} 2>&1"



##############################
# ANNOVAR
##############################

rule ANNOVAR:
    input:
        vcf="SelectVariants/{sample}_indel_filtered.vcf.gz",
        # contamin='contamination_normal/{sample_1}.contamination.table',   
    output:
        multianno_vcf='ANNOVAR/{sample}.filtered.hg38_multianno.vcf',
        multianno_txt='ANNOVAR/{sample}.filtered.hg38_multianno.txt'
    params:
        "ANNOVAR/{sample}.filtered"
    log:
        "logs/ANNOVAR/ANNOVAR_{sample}.log"
    shell:
        "/gatk/biosoft/annovar/table_annovar.pl "
        "{input.vcf} "
        "/gatk/biosoft/annovar/humandb/ "
        "--buildver hg38 "
        "--protocol refGene,exac03,gnomad211_exome " # refGene
        "--operation g,f,f -vcfinput "
        "--outfile {params} > {log} 2>&1"

# 有些参数和师兄的不一样，再看看
#


