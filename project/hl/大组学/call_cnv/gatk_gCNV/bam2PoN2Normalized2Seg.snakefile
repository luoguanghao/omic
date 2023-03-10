
configfile: "bam2PoN2Normalized2Seg.yaml"

rule all:
    input:
        expand("CollectReadCounts_ctrl/{sample_1}.hdf5", sample_1=config["samples_ctrl"]),
        expand("CollectReadCounts_tumor/{sample_2}.hdf5", sample_2=config["samples_tumor"]),
        'pon/pon.hdf5',
        expand("denoise/{sample_2}.called.seg", sample_2=config["samples_tumor"])

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
        " -I ".join(expand("CollectReadCounts_ctrl/{sample_1}.hdf5", sample_1=config["samples_ctrl"]))
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' CreateReadCountPanelOfNormals "
        "-I {params} "
        "--minimum-interval-median-percentile 5.0 "
        "-O {output}"

##############################
# Denoise read count by PoN
##############################

rule denoise_read_counts:
    input:
        # Required inputs in either hdf5 or tsv format.
        # Inputs should be generated by running CollectReadCounts.
        count_file = 'CollectReadCounts_tumor/{sample_2}.hdf5',
        pon = 'pon/pon.hdf5',
    output:
        # Required output PON file. MUST be hdf5 format.
        standardized_copy_ratios = 'denoise/{sample_2}.standardized_copy_ratios.tsv',
        denoised_copy_ratios = 'denoise/{sample_2}.denoised_copy_ratios.tsv',
    params:
        # Optional parameters. Omit if unused.
        extra = '',
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' DenoiseReadCounts "
        "-I {input.count_file} "
        "--count-panel-of-normals {input.pon} "
        "--standardized-copy-ratios {output.standardized_copy_ratios} "
        "--denoised-copy-ratios {output.denoised_copy_ratios}"





rule plot_denoised_copy_ratios:
    input:
        # Required inputs in either hdf5 or tsv format.
        # Inputs should be generated by running CollectReadCounts.
        standardized_copy_ratios = 'denoise/{sample_2}.standardized_copy_ratios.tsv',
        denoised_copy_ratios = 'denoise/{sample_2}.denoised_copy_ratios.tsv',
        ref_dict = '/gatk/my_dir/reference/genome/hg19/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.dict'
    output:
        extra = '',
    params:
        # Optional parameters. Omit if unused.
        '{sample_2}'
    shell:
        "/gatk/gatk PlotDenoisedCopyRatios "
        "--standardized-copy-ratios {input.standardized_copy_ratios} "
        "--denoised-copy-ratios {input.denoised_copy_ratios} "
        "--sequence-dictionary {input.ref_dict} "
        "--minimum-contig-length 46709983 "
        "--output sandbox/plot/ "
        "--output-prefix {params}"



## ModelSegments

rule model_segment:
    input:
        # Required inputs in either hdf5 or tsv format.
        # Inputs should be generated by running CollectReadCounts.
        denoised_copy_ratios = 'denoise/{sample_2}.denoised_copy_ratios.tsv'
    output:
        cr_seg_file = 'denoise/{sample_2}.cr.seg',
    params:
        # Optional parameters. Omit if unused.
        '{sample_2}'
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' ModelSegments "
        "--denoised-copy-ratios {input.denoised_copy_ratios} "
        "--output denoise/ "
        "--output-prefix {params}"


rule call_copy_ratio_segments:
    input:
        cr_seg_file = 'denoise/{sample_2}.cr.seg'
    output:
        called_seg_file = 'denoise/{sample_2}.called.seg',
    shell:
        "/gatk/gatk --java-options '-Xmx40G -Djava.io.tmpdir=./' CallCopyRatioSegments "
        "--input {input.cr_seg_file} "
        "--output {output.called_seg_file}"
















