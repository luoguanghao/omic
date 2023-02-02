# MuTect2
# Create PoN
#
# 在GenomicsDBImport的-L参数中，需要指定interval，这里指定 -L 1 -L 2 -L 3 -L 4 ...
#

configfile: "bam2vcf.yaml"

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_tumor(wildcards):
    return config["samples_tumor"][wildcards.sample_2]

# ====

rule all:
    input:
        expand("MuTect2_ctrl/{sample_1}.vcf.gz", sample_1=config["samples_ctrl"]),
        #expand("CollectReadCounts_tumor/{sample_2}.hdf5", sample_2=config["samples_tumor"]),
        'pon/gatk4_mutect2_pon.vcf.gz',
        #expand("denoise/{sample_2}.called.seg", sample_2=config["samples_tumor"])



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
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  Mutect2 "
        "-I {input.bam_file} "
        "-L {input.interval_list} "
        "-R {input.ref} "
        "--max-mnp-distance 0 "
        "-O {output}"




# PoN
rule GenomicsDBImport:
    input:
        vcf_files = expand('MuTect2_ctrl/{sample_1}.vcf.gz', sample_1=config["samples_ctrl"]),
        # interval_list="GRCh37.preprocessed.interval_list"  # <<<<<<<
        # Optional inputs. Omit if unused.

    output:
        # Required output PON file. MUST be hdf5 format.
        "pon_db"
    params:
        " -V ".join(expand("MuTect2_ctrl/{sample_1}.vcf.gz", sample_1=config["samples_ctrl"]))
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' GenomicsDBImport "
        "-V {params} "
        "-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y "
        "--genomicsdb-workspace-path {output}"
        



rule CreateSomaticPanelOfNormals:
    input:
        ref=config["ref"]
        pon_db="pon_db"
    output:
        "pon/gatk4_mutect2_pon.vcf.gz"
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' CreateSomaticPanelOfNormals "
        "-R {input.ref} "
        # "--germline-resource database/af-only-gnomad.raw.sites.hg19.vcf.gz "
        "-V gendb://{input.pon_db} "
        "-O {output}"




















































