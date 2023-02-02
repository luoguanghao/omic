############################
# 流程正在测试
# WES upstream pipeline
# from fastq to bam&bai
#
# func #
## bwa
## samtools view&sort
## MarkDuplicates
## FixMateInformation
## samtools index ***
## BaseRecalibrator
## ApplyBQSR
#
# output file #
## ${sample}_sorted_marked_fixed.bam.bai
## ${sample}_sorted_marked_fixed.bam
## ${sample}_bqsr.bam
## ${sample}_bqsr.bai 
#
# cmd:
## bash ./test_pipe.sh -i /gatk/resource/reference/index/bwa/hg19 -r /gatk/resource/reference/genome/hg19/hg19.fa -s /gatk/my_dir/reference/gatk_bundle/dbsnp_138.hg19.vcf -d /gatk/my_dir/reference/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -n NPC10F-T -a /gatk/my_dir/my_project/hl_call_cnv/data/tumor-NT-NPC-WES/fastqData/NPC10F-T_1.fastq.gz -b /gatk/my_dir/my_project/hl_call_cnv/data/tumor-NT-NPC-WES/fastqData/NPC10F-T_2.fastq.gz
## bash ./test_pipe.sh -i -r -s -d -n -a -b




while getopts "i:r:s:d:n:a:b:" arg #选项后面的冒号表示该选项需要参数
do
    case $arg in
        i)
            idx=$OPTARG #参数存在$OPTARG中
            ;;
        r)
            ref=$OPTARG #参数存在$OPTARG中
            ;;
        s)
            snp=$OPTARG
            ;;
        d)
            indel=$OPTARG
            ;;
        n)  
            sample=$OPTARG
            ;;
        a)
            fq1=$OPTARG
            ;;
        b)
            fq2=$OPTARG
        ;;
    esac
done

echo ''
echo 'action! mkdir ...'

mkdir ${sample}

echo ''
echo 'action! bwa ...'
bwa mem -t 8 -R '@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina' \
${idx} \
${fq1} \
${fq2} \
> ${sample}/${sample}.sam 2> ${sample}/${sample}.bwa.log
# ================================================
echo 'samtools sort ......'

samtools view -S ${sample}/${sample}.sam -b -o ./${sample}/${sample}.bam
samtools sort -@ 8 -o ${sample}/${sample}.sorted.bam ${sample}/${sample}.bam

echo ''
echo 'MarkDuplicates ......'
/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' MarkDuplicates \
-I ${sample}/${sample}.sorted.bam  \
-O ${sample}/${sample}_sorted_marked.bam  \
-M ${sample}/$sample.metrics  1> ${sample}/${sample}_log.mark 2>&1

echo 'FixMateInformation ......'
/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' FixMateInformation \
-I ${sample}/${sample}_sorted_marked.bam \
-O ${sample}/${sample}_sorted_marked_fixed.bam \
-SO coordinate  1> ${sample}/${sample}_log.fix 2>&1
echo 'samtools index ......'
samtools index ${sample}/${sample}_sorted_marked_fixed.bam

echo ''
echo 'BaseRecalibrator ......'
/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' BaseRecalibrator \
-R $ref \
-I ${sample}/${sample}_sorted_marked_fixed.bam \
--known-sites $snp \
--known-sites $indel \
-O ${sample}/${sample}_recal.table  1> ${sample}/${sample}_log.recal 2>&1

echo ''
echo 'ApplyBQSR ......'
/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' ApplyBQSR \
-R $ref \
-I ${sample}/${sample}_sorted_marked_fixed.bam -bqsr ${sample}/${sample}_recal.table \
-O ${sample}/${sample}_bqsr.bam  1> ${sample}/${sample}_log.ApplyBQSR 2>&1


