######################
# 改用新的参考基因组
# 染色体名不要chr
# 抽空构建snakemake流程
#


GATK=/gatk/gatk

#############################
# prepare for interval_list
#############################
#bed=
#ref=
#dict=
#work_dir=

# bed=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/novo_wes_add_chr.bed
bed=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/novo_wes.bed
ref=/gatk/my_dir/reference/genome/hg19/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
dict=/gatk/my_dir/reference/genome/hg19/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.dict
work_dir=.

## bed to intervals_list
# 输入bed文件
# 输出 interval_list文件
$GATK BedToIntervalList -I ${bed} -O ${work_dir}/GRCh37.exon.interval_list -SD ${dict}
# $GATK BedToIntervalList --SORT -I ${bed} -O ${work_dir}/hg19.exon.interval_list -SD ${dict}


## Preprocess Intervals
## 输入 interval_list，参考基因组dict，参考基因组fa，一些 参数
## 输出 targets.preprocessed.interval_list
$GATK  PreprocessIntervals \
-L ${work_dir}/GRCh37.exon.interval_list \
--sequence-dictionary ${dict} \
--reference ${ref}  \
--padding 250 \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--output ${work_dir}/GRCh37.preprocessed.interval_list


######################
##  get read count
######################
ref=/gatk/my_dir/reference/genome/hg19/hg19.fa
bam_file=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/fq_file/tumor/FJ20190320_FDHE21H002338_bqsr.bam
readcount_file=FJ20190320_FDHE21H002338.hdf5

$GATK  --java-options "-Xmx40G"  CollectReadCounts \
  -I ${bam_file} \
  -L /gatk/my_dir/my_upstream_project/hl_call_cnv/data/targets.preprocessed.interval_list \
  -R ${ref} \
  --format HDF5 \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output sandbox/${readcount_file}
# ===================
ref=/gatk/my_dir/reference/genome/hg19/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
bam_file=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/bam_file/ctrl/FJ20190320.final.bam
readcount_file=FJ20190320.hdf5

$GATK  --java-options "-Xmx40G"  CollectReadCounts \
  -I ${bam_file} \
  -L ${work_dir}/GRCh37.preprocessed.interval_list \
  -R ${ref} \
  --format HDF5 \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output sandbox/${readcount_file}


######################
##  PoN
######################
#normal_readcount_files=
#pon_file=
normal_readcount_files=
pon_file=cnvponC.pon.hdf5

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" CreateReadCountPanelOfNormals \
    -I sandbox/HW3.hdf5 \
    -I sandbox/HW4.hdf5 \
    -I sandbox/LJ.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    -O sandbox/${pon_file}

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" CreateReadCountPanelOfNormals \
    -I sandbox/HW3.hdf5 \
    -I sandbox/HW4.hdf5 \
    -I sandbox/LJ.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    -O sandbox/cnvponC.pon.hdf5



###############################
# denoise & segmentation
###############################
# 
work_dir=
pon_file=
tumor_readcount_file= # <<<<<<<<<<<<<<
scr_file=
dcr_file=
plot_dir=
output_prefix=
## denoise

gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I sandbox/${tumor_readcount_file} \
    --count-panel-of-normals sandbox/${pon_file} \
    --standardized-copy-ratios sandbox/${scr_file} \
    --denoised-copy-ratios sandbox/${dcr_file}

gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I sandbox/FJ20190320_FDHE21H002338.hdf5 \
    --count-panel-of-normals sandbox/cnvponC.pon.hdf5 \
    --standardized-copy-ratios sandbox/FJ20190320_FDHE21H002338.standardizedCR.tsv \
    --denoised-copy-ratios sandbox/FJ20190320_FDHE21H002338.denoisedCR.tsv



## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios


# test tutorial file in 4.2.0

$GATK PlotDenoisedCopyRatios \
    --standardized-copy-ratios sandbox/${scr_file} \
    --denoised-copy-ratios sandbox/${dcr_file} \
    --sequence-dictionary ${dict} \
    --minimum-contig-length 46709983 \
    --output sandbox/${plot_dir} \
    --output-prefix ${output_prefix}




## ModelSegments


# no allelic
$GATK --java-options "-Xmx4g" ModelSegments \
    --denoised-copy-ratios sandbox/${dcr_file} \
    --output sandbox \
    --output-prefix ${output_prefix}

$GATK CallCopyRatioSegments \
    --input sandbox/${output_prefix}.cr.seg \
    --output sandbox/${output_prefix}.called.seg









