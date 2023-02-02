

# MuTect2
# Create PoN
#configfile: "/gatk/my_dir/my_upstream_project/hl_call_cnv/gatk_sCNV/try_again_4/bam2PoN2Normalized2Seg.yaml"
configfile: "bam2vcf.yaml"

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_normal(wildcards):
    return config["samples_normal"][wildcards.sample_1]



rule all:
    input:
        expand("ANNOVAR_normal/{sample_1}.filtered.hg19_multianno.vcf", sample_1=config["samples_normal"]),
        expand("ANNOVAR_normal/{sample_1}.filtered.hg19_multianno.txt", sample_1=config["samples_normal"])



#####################################
# Call
#####################################


rule HaplotypeCaller_normal:
    input:
        bam_file=get_bam_file_normal,
        ref=config["ref"],
        snp=config["snp"],
    output:
        "HaplotypeCaller_normal/{sample_1}.vcf.gz"
    log:
        "logs/HaplotypeCaller_normal/HaplotypeCaller_normal_{sample_1}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  HaplotypeCaller "
        "-I {input.bam_file} " # v
        "-R {input.ref} " # v
        "--dbsnp {input.snp} "
        #"--max-mnp-distance 0 "
        "-O {output} > {log} 2>&1"



rule SelectVariants_normal_snp:
    input:
        vcf="HaplotypeCaller_normal/{sample_1}.vcf.gz"
    output:
        sel_vcf="SelectVariants_normal/{sample_1}_snp.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_normal_{sample_1}_snp.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type SNP "
        "-V {input.vcf} "
        "-O {output.sel_vcf} "


rule SelectVariants_normal_indel:
    input:
        vcf="HaplotypeCaller_normal/{sample_1}.vcf.gz"
    output:
        sel_vcf="SelectVariants_normal/{sample_1}_indel.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_normal_{sample_1}_indel.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  SelectVariants "
        "-select-type INDEL "
        "-V {input.vcf} "
        "-O {output.sel_vcf} "

# ===

rule VariantFiltration_normal_snp:
    input:
        sel_vcf="SelectVariants_normal/{sample_1}_snp.vcf.gz"
    output:
        fil_vcf="SelectVariants_normal/{sample_1}_snp_filtered.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_normal_{sample_1}_snp_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.sel_vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} "


rule VariantFiltration_normal_indel:
    input:
        sel_vcf="SelectVariants_normal/{sample_1}_indel.vcf.gz"
    output:
        fil_vcf="SelectVariants_normal/{sample_1}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_normal_{sample_1}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' VariantFiltration "
        "-V {input.sel_vcf} "
        '--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filter-name "PASS" '
        "-O {output.fil_vcf} "


rule MergeVcfs_normal:
    input:
        fil_vcf_indel="SelectVariants_normal/{sample_1}_indel_filtered.vcf.gz",
        fil_vcf_snp="SelectVariants_normal/{sample_1}_snp_filtered.vcf.gz"
    output:
        merge_vcf="SelectVariants_normal/{sample_1}_indel_filtered.vcf.gz"
    log:
        "logs/SelectVariants_normal/SelectVariants_normal_{sample_1}_indel_filtered.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G' MergeVcfs "
        "-I {input.fil_vcf_indel} "
        "-I {input.fil_vcf_snp} "
        "-O {output.merge_vcf}"



##############################
# ANNOVAR
##############################

rule ANNOVAR_normal:
    input:
        vcf="SelectVariants_normal/{sample_1}_indel_filtered.vcf.gz",
        # contamin='contamination_normal/{sample_1}.contamination.table',   
    output:
        multianno_vcf='ANNOVAR_normal/{sample_1}.filtered.hg38_multianno.vcf',
        multianno_txt='ANNOVAR_normal/{sample_1}.filtered.hg38_multianno.txt'
    params:
        "ANNOVAR_normal/{sample_1}.filtered"
    log:
        "logs/ANNOVAR_normal/ANNOVAR_normal_{sample_1}.log"
    shell:
        "/biosoft/annovar/table_annovar.pl "
        "{input.vcf} "
        "/biosoft/annovar/humandb/ "
        "--buildver hg19 "
        "--protocol refGene,exac03,gnomad211_exome " # refGene
        "--operation g,f,f -vcfinput "
        "--outfile {params} > {log} 2>&1"

# 有些参数和师兄的不一样，再看看
#

























































































