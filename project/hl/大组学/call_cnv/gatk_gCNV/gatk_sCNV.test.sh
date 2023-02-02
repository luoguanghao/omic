# 记得完善log输出
#
#

GATK=/gatk/gatk
ref=/gatk/my_dir/reference/genome/hg19/hg19.fa
dict=/gatk/my_dir/reference/genome/hg19/hg19.dict
#ref=/gatk/my_dir/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
#dict=/gatk/my_dir/reference/gatk_bundle/Homo_sapiens_assembly38.dict

bed=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/novo_wes_add_chr.bed

work_dir=.

## bed to intervals_list
# 输入bed文件
# 输出 interval_list文件
# $GATK BedToIntervalList -I ${bed} -O ~/wes_cancer/data/hg38.exon.interval_list -SD ${dict}
$GATK BedToIntervalList --SORT -I ${bed} -O ${work_dir}/hg19.exon.interval_list -SD ${dict}


## Preprocess Intervals
## 输入 interval_list，参考基因组dict，参考基因组fa，一些 参数
## 输出 targets.preprocessed.interval_list
$GATK  PreprocessIntervals \
-L ${work_dir}/hg19.exon.interval_list \
--sequence-dictionary ${dict} \
--reference ${ref}  \
--padding 250 \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--output ${work_dir}/targets.preprocessed.interval_list

# 4.0.1.1
/gatk/gatk PreprocessIntervals \
-L ./hg19.exon.interval_list \
--sequence-dictionary /gatk/my_dir/reference/genome/hg19/hg19.dict \
--reference /gatk/my_dir/reference/genome/hg19/hg19.fa  \
--padding 250 \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--output ./targets.preprocessed.interval.list

######################
##  get read count
######################

bam_dir=

#bam_file=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/fq_file_old/HW3_bqsr.bam
bam_file=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/bam_file/ctrl/HW3.final.chr.bam
output_file=HW3_2.hdf5
# new CollectReadCounts...  this is ok ##########
$GATK  --java-options "-Xmx40G"  CollectReadCounts \
  -I ${bam_file} \
  -L /gatk/my_dir/my_upstream_project/hl_call_cnv/data/targets.preprocessed.interval_list \
  -R ${ref} \
  --format HDF5  \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output sandbox/${output_file}

# 4.0.1.1
$GATK CollectFragmentCounts \
    -I ${bam_file} \
    -L targets.preprocessed.interval.list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${output_file}

# 4.0.1.1 demo ...
gatk CollectFragmentCounts \
    -I tumor.bam \
    -L sandbox/targets_C.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O sandbox/tumor.counts.hdf5

# test...
$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  CollectReadCounts \
  -I /gatk/my_dir/my_upstream_project/hl_call_cnv/data/bam_file/ctrl/HW3.final.chr.bam \
  -L /gatk/my_dir/my_upstream_project/hl_call_cnv/data/targets_C.preprocessed.interval_list \
  -R ${ref} \
  --format HDF5  \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output sandbox/tumor.counts.hdf5


targets_C.preprocessed.interval_list



# test tutorial file in 4.2.0
$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  CollectReadCounts \
  -I /gatk/my_dir/my_upstream_project/hl_call_cnv/data/gatk_tutorial/tutorial_11682/tumor.bam \
  -L /gatk/my_dir/my_upstream_project/hl_call_cnv/data/targets_C.preprocessed.interval_list \
  -R ${ref} \
  --format HDF5  \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output sandbox/tumor.counts.hdf5





######################
##  PoN
######################

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" CreateReadCountPanelOfNormals \
    -I HG00133.alt_bwamem_GRCh38DH.20150826.GBR.exome.counts.hdf5 \
    -I HG00733.alt_bwamem_GRCh38DH.20150826.PUR.exome.counts.hdf5 \
    -I NA19654.alt_bwamem_GRCh38DH.20150826.MXL.exome.counts.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    -O sandbox/cnvponC.pon.hdf5

$(for i in {1..6} ;do echo "--input ./8.cnv/gatk/counts/case${i}_germline.clean_counts.hdf5" ;done)


# test tutorial file in 4.2.0
tutor1=/gatk/my_dir/my_upstream_project/hl_call_cnv/data/gatk_tutorial/tutorial_11682

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" CreateReadCountPanelOfNormals \
    -I ${tutor1}/HG00133.alt_bwamem_GRCh38DH.20150826.GBR.exome.counts.hdf5 \
    -I ${tutor1}/HG00733.alt_bwamem_GRCh38DH.20150826.PUR.exome.counts.hdf5 \
    -I ${tutor1}/NA19654.alt_bwamem_GRCh38DH.20150826.MXL.exome.counts.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    -O sandbox/cnvponC.pon.hdf5



#########################
## denoise
#########################

gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I hcc1143_T_clean.counts.hdf5 \
    --count-panel-of-normals cnvponC.pon.hdf5 \
    --standardized-copy-ratios sandbox/hcc1143_T_clean.standardizedCR.tsv \
    --denoised-copy-ratios sandbox/hcc1143_T_clean.denoisedCR.tsv

# test tutorial file in 4.2.0

gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I ${tutor1}/hcc1143_T_clean.counts.hdf5 \
    --count-panel-of-normals sandbox/cnvponC.pon.hdf5 \
    --standardized-copy-ratios sandbox/hcc1143_T_clean.standardizedCR.tsv \
    --denoised-copy-ratios sandbox/hcc1143_T_clean.denoisedCR.tsv


############################
## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios
############################




# test tutorial file in 4.2.0

$GATK PlotDenoisedCopyRatios \
    --standardized-copy-ratios sandbox/hcc1143_T_clean.standardizedCR.tsv \
    --denoised-copy-ratios sandbox/hcc1143_T_clean.denoisedCR.tsv \
    --sequence-dictionary ${dict} \
    --minimum-contig-length 46709983 \
    --output sandbox/plots \
    --output-prefix hcc1143_T_clean






###################################################
# Count ref and alt alleles at common germline variant sites using CollectAllelicCounts
###################################################


# demo...
gatk --java-options "-Xmx4g" ModelSegments \
    --denoised-copy-ratios hcc1143_T_clean.denoisedCR.tsv \
    --allelic-counts hcc1143_T_clean.allelicCounts.tsv \ # 可选
    --normal-allelic-counts hcc1143_N_clean.allelicCounts.tsv \ # 可选
    --output sandbox \
    --output-prefix hcc1143_T_clean


####################################################################
# Call copy-neutral, amplified and deleted segments with CallCopyRatioSegments
####################################################################
# demo...
gatk CallCopyRatioSegments \
    --input hcc1143_T_clean.cr.seg \
    --output sandbox/hcc1143_T_clean.called.seg

# my_data
$GATK CallCopyRatioSegments \
    --input sandbox/FJ20190320_FDHE21H002338_T_clean.cr.seg \
    --output sandbox/FJ20190320_FDHE21H002338_T_clean.called.seg



###########################################################################
# Plot modeled copy ratio and allelic fraction segments with PlotModeledSegments
###########################################################################

gatk PlotModeledSegments \
    --denoised-copy-ratios hcc1143_T_clean.denoisedCR.tsv \
    --allelic-counts hcc1143_T_clean.hets.tsv \
    --segments hcc1143_T_clean.modelFinal.seg \
    --sequence-dictionary Homo_sapiens_assembly38.dict \
    --minimum-contig-length 46709983 \
    --output sandbox/plots \
    --output-prefix hcc1143_T_clean






















# ===================================================

gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' PreprocessIntervals \
    -L ../tutorial_11682/targets_C.interval_list \
    -R /gatk/ref/Homo_sapiens_assembly38.fasta \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O mytest/sandbox/targets_C.preprocessed.interval_list




gatk CollectFragmentCounts \
    -I tumor.bam \
    -L targets_C.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O sandbox/tumor.counts.hdf5






















