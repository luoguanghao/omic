# cmd:
## bash ./hisat2_htcount.sh -i -r -g -n -a -b
#
#
#

while getopts "i:r:g:n:a:b:" arg #选项后面的冒号表示该选项需要参数
do
    case $arg in
        i)
            idx=$OPTARG #参数存在$OPTARG中
            ;;
        r)
            ref=$OPTARG #参数存在$OPTARG中
            ;;
        g)
            gtf=$OPTARG
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
echo '#####'

echo "${fq2}"

echo '#####'

echo ''
echo 'action! mkdir ...'
mkdir ${sample}

echo ''
echo 'action! hisat2 ...'
echo "hisat2 -p 8 -t -x ${idx} -1 ${fq1} -2 ${fq2} -S ${sample}/${sample}.sam 1> ${sample}/${sample}_log.hisat2 2>&1"

hisat2 -p 8 -t -x ${idx} -1 ${fq1} -2 ${fq2} -S ${sample}/${sample}.sam 1> ${sample}/${sample}_log.hisat2 2>&1

echo ''
echo 'action! samtools ...'

samtools view  -S ${sample}/${sample}.sam -b > ${sample}/${sample}.bam
samtools sort ${sample}/${sample}.bam -o ${sample}/${sample}_sorted.bam 1> ${sample}/${sample}_log.samtools_sort 2>&1
samtools index ${sample}/${sample}_sorted.bam ${sample}/${sample}_sorted.bai 1> ${sample}/${sample}_log.samtools_index 2>&1

rm ${sample}/${sample}.sam
rm ${sample}/${sample}.bam

#echo ''
#echo 'action! htseq-count ...'
#echo "htseq-count -r pos -f bam -s no -i gene_name ${sample}/${sample}_sorted.bam ${gtf} > ${sample}/${sample}.count.txt 1> ${sample}/${sample}_log.htseq 2>&1"

#htseq-count -r pos -f bam -s no -i gene_name ${sample}/${sample}_sorted.bam ${gtf} 1> ${sample}/${sample}.count.txt 2> ${sample}/${sample}_log.htseq

echo ''
echo 'done !'






