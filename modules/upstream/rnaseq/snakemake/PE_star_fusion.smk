# RNAseq star-fusion PE

configfile: "PE_star_fusion.yaml"

def get_fq_file_1(wildcards):
    return config["fq_file_1"][wildcards.sample]

#def get_fq_file_2(wildcards):
#    return config["fq_file_2"][wildcards.fq_2]

# star-fusion.fusion_predictions.tsv

rule all:
    input:
        expand("result/{sample}/star-fusion.fusion_predictions.tsv", sample=config["samples"]),


rule star_fusion:
    input:
        fq_file1=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_1.fq.gz",sample=wildcards.sample),
        fq_file2=lambda wildcards: expand("/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/rnaseq/raw_data/pair_end/{sample}/{sample}_2.fq.gz",sample=wildcards.sample),
        star_fusion_lib=config["star_fusion_lib"], #<<<
    output:
        fusion_pred="result/{sample}/star-fusion.fusion_predictions.tsv"
    params:
        "result/{sample}/"
    log:
        "logs/star_fusion/{sample}.log"        
    threads: 8
    shell:
        """STAR-Fusion -p 8 -t \
        --genome_lib_dir {input.star_fusion_lib} \
        --left_fq {input.fq_file1} --right_fq {input.fq_file2} \
        --output_dir {params} 1> {log} 2>&1"""
