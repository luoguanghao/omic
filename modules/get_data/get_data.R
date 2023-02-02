library(tidyverse)


###### prepare for ascp download
df = read_tsv('./filereport_read_run_PRJNA723578_tsv.txt')

cat(sprintf(
    'nohup ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/%s/%s/%s ./ &',
        sapply(strsplit(df$sra_ftp,'/'),'[[',4),
        sapply(strsplit(df$sra_ftp,'/'),'[[',5),
        sapply(strsplit(df$sra_ftp,'/'),'[[',6)), sep='\n')


###### unpack




# sh脚本
#!/bin/sh
for i in *sra
do
echo $i
/y/home/lgh/software/sratoolkit/bin/fasterq-dump --split-3 $i
done



