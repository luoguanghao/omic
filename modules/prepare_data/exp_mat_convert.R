
#######################
# 计算基因转录本长度
#######################

library(GenomicRanges)
library(rtracklayer)

gtf_dir = '/home/lgh/reference/hl_reference/gtf/ensembl/Homo_sapiens.GRCh38.98.gtf.gz'

############

cout_transcript_length<-function(gtf_dir){
library(GenomicRanges)
library(rtracklayer)
suppressMessages(library(tidyverse))
  gtf <- rtracklayer::import(gtf_dir)
  gtf_df=as.data.frame(gtf) 
  gtf_df_exon<-subset(gtf_df,type=="exon")
  transcripts<-as.character(unique(as.character(gtf_df_exon$transcript_id)))
  
  out_tab<-data.frame()
  for (i in transcripts){
    transcripts=i
    rt<-gtf_df_exon %>% filter(transcript_id==i)
    if (nrow(rt)>1){
      lengths=sum(rt$width)
    }else{
      lengths=as.numeric(rt$width)
    }
    out_tab=rbind(out_tab,cbind(transcripts,lengths))
  }
  return(out_tab)
}


trans_len_df = cout_transcript_length(gtf_dir)

### output ###
write_tsv(trans_len_df,'trans_len_df.tsv')
#######################
#######################


#######################
# count2tpm
#######################








#######################
# count2rpkm
#######################




