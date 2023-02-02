# MuTect2
# Create PoN
#configfile: "/gatk/my_dir/my_upstream_project/hl_call_cnv/gatk_sCNV/try_again_4/bam2PoN2Normalized2Seg.yaml"
configfile: "bam2vcf.yaml"

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_tumor(wildcards):
    return config["samples_tumor"][wildcards.sample_2]

# ====

rule all:
    input:
        #'pon/gatk4_mutect2_pon.vcf.gz',
        #expand("MuTect2_tumor/{sample_2}.vcf.gz", sample_2=config["samples_tumor"]),
        #expand("contamination_tumor/{sample_2}.contamination.table", sample_2=config["samples_tumor"])
        #expand("FilterMutectCalls_tumor/{sample_2}.filtered.vcf.gz", sample_2=config["samples_tumor"]),
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.vcf", sample_2=config["samples_tumor"]),
        expand("ANNOVAR_tumor/{sample_2}.filtered.hg19_multianno.txt", sample_2=config["samples_tumor"])


####################
# Create PoN
####################

# ctrl sample
rule MuTect2_ctrl:
    input:
        bam_file=get_bam_file_ctrl,
        interval_list="GRCh37.preprocessed.interval_list",  # <<<<<<<
        ref=config["ref"]
    output:
        "MuTect2_ctrl/{sample_1}.vcf.gz"
    log:
        "logs/MuTect2_ctrl/MuTect2_ctrl_{sample_1}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  Mutect2 "
        "-I {input.bam_file} "
        "-L {input.interval_list} "
        "-R {input.ref} "
        "--max-mnp-distance 0 "
        "-O {output} > {log} 2>&1"

        

# PoN
rule GenomicsDBImport:
    input:
        vcf_files = expand('MuTect2_ctrl/{sample_1}.vcf.gz', sample_1=config["samples_ctrl"]),
        interval_list="GRCh37.preprocessed.interval_list"  # <<<<<<<
        # Optional inputs. Omit if unused.

    output:
        # Required output PON file. MUST be hdf5 format.
        "pon_db"
    log:
        "logs/GenomicsDBImport.log"
    params:
        " -V ".join(expand("MuTect2_ctrl/{sample_1}.vcf.gz", sample_1=config["samples_ctrl"]))
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' GenomicsDBImport "
        "-V {params} "
        "-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y "
        "--genomicsdb-workspace-path {output} > {log} 2>&1  "
        # --reader-threads
        



rule CreateSomaticPanelOfNormals:
    input:
        ref=config["ref"],
        pon_db="pon_db"
    output:
        "pon/gatk4_mutect2_pon.vcf.gz"
    log:
        "logs/CreateSomaticPanelOfNormals.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' CreateSomaticPanelOfNormals "
        "-R {input.ref} "
        # "--germline-resource database/af-only-gnomad.raw.sites.hg19.vcf.gz "
        "-V gendb://{input.pon_db} "
        "-O {output} > {log} 2>&1"




##############################
# MuTect2 no matched mode
##############################


# tumor sample
rule MuTect2_tumor:
    input:
        bam_file=get_bam_file_tumor,
        interval_list="GRCh37.preprocessed.interval_list",  # <<<<<<<
        ref=config["ref"],
        pon="pon/gatk4_mutect2_pon.vcf.gz",
        gr="/gatk/my_dir/Download/gatk_bundle/gn/af-only-gnomad.raw.sites.b37.vcf.gz"
    output:
        vcf="MuTect2_tumor/{sample_2}.vcf.gz",
        f1r2="MuTect2_tumor/{sample_2}.f1r2.tar.gz"
    log:
        "logs/MuTect2_tumor/MuTect2_tumor_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  Mutect2 "
        "-I {input.bam_file} "
        "-L {input.interval_list} "
        "-R {input.ref} "
        "--panel-of-normals {input.pon} "
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
        "{input.vcf_addchr} "
        "/biosoft/annovar/humandb/ "
        "--buildver hg19 "
        "--protocol refGene,exac03,gnomad_genome " # refGene
        "--operation g,f,f -vcfinput "
        "--outfile {params} > {log} 2>&1"

'''
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
        "{input.vcf_addchr} "
        "/biosoft/annovar/humandb/ "
        "--buildver hg19 "
        "--protocol ensGene,refGene,knownGene,clivar_20210501 " # refGene
        "--operation g,g,g,f -vcfinput "
        "--outfile {params} > {log} 2>&1"
'''


##############################
# Funcotator
##############################

rule FilterMutectCalls_tumor:
    input:
        vcf="MuTect2_tumor/{sample_2}.vcf.gz",
        ref=config["ref"],
        # contamin='contamination_tumor/{sample_2}.contamination.table',   
    output:
        'FilterMutectCalls_tumor/{sample_2}.filtered.vcf.gz'
    log:
        "logs/FilterMutectCalls_tumor/FilterMutectCalls_{sample_2}.log"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  Funcotator "
        "--variant variants.vcf "
        "--reference Homo_sapiens_assembly19.fasta "
        "--ref-version hg19  "
        "--data-sources-path funcotator_dataSources.v1.2.20180329 "
        "--output variants.funcotated.maf "
        "--output-file-format MAF > {log} 2>&1"

#./gatk Funcotator \
#     --variant variants.vcf \
#     --reference Homo_sapiens_assembly19.fasta \
#     --ref-version hg19 \
#     --data-sources-path funcotator_dataSources.v1.2.20180329 \
#     --output variants.funcotated.maf \
#     --output-file-format MAF


























