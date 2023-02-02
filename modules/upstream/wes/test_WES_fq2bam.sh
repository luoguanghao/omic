


while getopts "r:s:d:n:a:b" arg #选项后面的冒号表示该选项需要参数
do
    case $arg in
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
        ;;
    esac
done
#===========

#GATK=/gatk/gatk
#ref=/gatk/resource/reference/genome/hg19/hg19.fa
#snp=/gatk/my_dir/reference/gatk_bundle/dbsnp_138.hg19.vcf
#indel=/gatk/my_dir/reference/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
#sample=NPC10F-T

echo 'action! bwa ...'
bwa mem -t 8 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" \
${ref} \
${fq1} \
${fq2} \
> ${sample}.sam 2> ${sample}.bwa.log
# ================================================
echo 'samtools sort ......'

samtools view -S ${sample}.sam -b -o ./${sample}.bam
samtools sort -@ 8 -o ${sample}.sorted.bam ${sample}.bam

echo 'MarkDuplicates ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${sample}.sorted.bam  \
-O ${sample}_sorted_marked.bam  \
-M $sample.metrics  1>${sample}_log.mark 2>&1
echo 'FixMateInformation ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
-I ${sample}_sorted_marked.bam \
-O ${sample}_sorted_marked_fixed.bam \
-SO coordinate  1>${sample}_log.fix 2>&1

echo 'samtools index ......'
samtools index ${sample}_sorted_marked_fixed.bam


echo 'BaseRecalibrator ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam \
--known-sites $snp \
--known-sites $indel \
-O ${sample}_recal.table  1>${sample}_log.recal 2>&1

echo 'ApplyBQSR ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam -bqsr ${sample}_recal.table \
-O ${sample}_bqsr.bam  1>${sample}_log.ApplyBQSR 2>&1








################################
# only for try
################################





echo 'action! bwa ...'
echo "bwa mem -t 8 -R '@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina' \
${ref} \
${fq1} \
${fq2} \
> ${sample}.sam 2> ${sample}.bwa.log"
# ================================================
echo 'samtools sort ......'

echo "samtools view -S ${sample}.sam -b -o ./${sample}.bam"
echo "samtools sort -@ 8 -o ${sample}.sorted.bam ${sample}.bam"

echo 'MarkDuplicates ......'
echo "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' MarkDuplicates \
-I ${sample}.sorted.bam  \
-O ${sample}_sorted_marked.bam  \
-M $sample.metrics  1>${sample}_log.mark 2>&1"

echo 'FixMateInformation ......'
echo "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' FixMateInformation \
-I ${sample}_sorted_marked.bam \
-O ${sample}_sorted_marked_fixed.bam \
-SO coordinate  1>${sample}_log.fix 2>&1"
echo 'samtools index ......'
echo "samtools index ${sample}_sorted_marked_fixed.bam"


echo 'BaseRecalibrator ......'
echo "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' BaseRecalibrator \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam \
--known-sites $snp \
--known-sites $indel \
-O ${sample}_recal.table  1>${sample}_log.recal 2>&1"

echo '\nApplyBQSR ......'
echo "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' ApplyBQSR \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam -bqsr ${sample}_recal.table \
-O ${sample}_bqsr.bam  1>${sample}_log.ApplyBQSR 2>&1"







## =====
echo "/gatk/gatk --java-options '-Xmx20G -Djava.io.tmpdir=./' BaseRecalibrator \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam \
--known-sites $snp \
--known-sites $indel \
-O ${sample}_recal.table  1>${sample}_log.recal 2>&1"



############
############
GATK=/gatk/gatk
ref=/gatk/resource/reference/genome/hg19/hg19.fa
snp=/gatk/my_dir/reference/gatk_bundle/dbsnp_138.hg19.vcf
indel=/gatk/my_dir/reference/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
sample=NPC10F-T

echo 'action! bwa ...'
bwa mem -t 8 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" \
/gatk/resource/reference/index/bwa/hg19 \
/gatk/my_dir/my_project/hl_call_cnv/data/tumor-NT-NPC-WES/fastqData/NPC10F-T_1.fastq.gz \
/gatk/my_dir/my_project/hl_call_cnv/data/tumor-NT-NPC-WES/fastqData/NPC10F-T_2.fastq.gz \
> ${sample}.sam 2> ${sample}.bwa.log
# ================================================
echo 'samtools sort ......'

samtools view -S ${sample}.sam -b -o ./${sample}.bam
samtools sort -@ 8 -o ${sample}.sorted.bam ${sample}.bam

echo 'MarkDuplicates ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${sample}.sorted.bam  \
-O ${sample}_sorted_marked.bam  \
-M $sample.metrics  1>${sample}_log.mark 2>&1
echo 'FixMateInformation ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
-I ${sample}_sorted_marked.bam \
-O ${sample}_sorted_marked_fixed.bam \
-SO coordinate  1>${sample}_log.fix 2>&1
echo 'samtools index ......'
samtools index ${sample}_sorted_marked_fixed.bam


echo 'BaseRecalibrator ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam \
--known-sites $snp \
--known-sites $indel \
-O ${sample}_recal.table  1>${sample}_log.recal 2>&1

echo 'ApplyBQSR ......'
/gatk/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
-R $ref \
-I ${sample}_sorted_marked_fixed.bam -bqsr ${sample}_recal.table \
-O ${sample}_bqsr.bam  1>${sample}_log.ApplyBQSR 2>&1