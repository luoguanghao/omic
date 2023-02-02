
configfile: "bam2PoN.yaml"

rule all:
    input:
        expand("CollectReadCounts_ctrl/{sample_1}.hdf5", sample_1=config["samples_ctrl"]),
        expand("CollectReadCounts_tumor/{sample_2}.hdf5", sample_2=config["samples_tumor"]),
        'pon/pon.hdf5'

# ====

def get_bam_file_ctrl(wildcards):
    return config["samples_ctrl"][wildcards.sample_1]

def get_bam_file_tumor(wildcards):
    return config["samples_tumor"][wildcards.sample_2]



# ctrl sample
rule CollectReadCounts_ctrl:
    input:
        bam_file=get_bam_file_ctrl,
        interval_list="GRCh37.preprocessed.interval_list",
        ref=config["ref"]
    output:
        "CollectReadCounts_ctrl/{sample_1}.hdf5"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  CollectReadCounts "
        "-I {input.bam_file} "
        "-L {input.interval_list} "
        "-R {input.ref} "
        "--format HDF5 "
        "--interval-merging-rule OVERLAPPING_ONLY "
        "--output {output} "

        
# tumor
rule CollectReadCounts_tumor:
    input:
        bam_file=get_bam_file_tumor,
        interval_list="GRCh37.preprocessed.interval_list",
        ref=config["ref"]
    output:
        "CollectReadCounts_tumor/{sample_2}.hdf5"
    shell:
        "/gatk/gatk --java-options '-Xmx40G'  CollectReadCounts "
        "-I {input.bam_file} "
        "-L {input.interval_list} "
        "-R {input.ref} "
        "--format HDF5 "
        "--interval-merging-rule OVERLAPPING_ONLY "
        "--output {output} "

# PoN
rule create_read_count_panel_of_normals:
    input:
        count_files = expand('CollectReadCounts_ctrl/{sample_1}.hdf5', sample_1=config["samples_ctrl"])
        # Optional inputs. Omit if unused.

    output:
        # Required output PON file. MUST be hdf5 format.
        'pon/pon.hdf5'
    params:
        " -I ".join(expand("CollectReadCounts/{sample_1}.hdf5", sample_1=config["samples_ctrl"]))
    shell:
        "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' CreateReadCountPanelOfNormals "
        "{params} "
        "--minimum-interval-median-percentile 5.0 "
        "-O {output}"






