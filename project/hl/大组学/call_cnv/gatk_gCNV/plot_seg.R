library(data.table)

# === function ===
ReadTSV = function(tsv_file) {
    # We need to filter out header lines beginning with '@';
    # however, the standard 'fread("grep ...")' causes issues with the default Docker container, so we use a temporary file.
    # See https://github.com/broadinstitute/gatk/issues/4140.
    temp_file = tempfile()
    system(sprintf('grep -v ^@ "%s" > %s', tsv_file, temp_file))
    return(suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
}

PlotCopyRatiosWithModeledSegments = function(denoised_copy_ratios_df, modeled_segments_df, contig_names, contig_starts, point_size=0.2) {
   points_start_index = 1
   for (s in 1:nrow(modeled_segments_df)) {
       #skip segments with no points
       num_points = modeled_segments_df[s, "NUM_POINTS_COPY_RATIO"]
       if (num_points == 0) {
           next
       }
       points_end_index = points_start_index + num_points

       contig = modeled_segments_df[s, "CONTIG"]
       offset = contig_starts[match(contig, contig_names)]
       segment_start = offset + modeled_segments_df[s, "START"]
       segment_end = offset + modeled_segments_df[s, "END"]
       genomic_coordinates = offset + denoised_copy_ratios_df[points_start_index:points_end_index, "MIDDLE"]

       denoised_copy_ratios = denoised_copy_ratios_df[points_start_index:points_end_index, "COPY_RATIO"]

       colors = c("coral", "dodgerblue")
       points(x=genomic_coordinates, y=denoised_copy_ratios, col=colors[s %% 2 + 1], pch=".", cex=point_size)

       copy_ratio_posterior_10 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_10"]
       copy_ratio_posterior_50 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_50"]
       copy_ratio_posterior_90 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_90"]
       segments(x0=segment_start, y0=copy_ratio_posterior_50, x1=segment_end, y1=copy_ratio_posterior_50, col="black", lwd=2, lty=1)
       rect(xleft=segment_start, ybottom=copy_ratio_posterior_10, xright=segment_end, ytop=copy_ratio_posterior_90, lwd=1, lty=1)

       points_start_index = points_start_index + num_points
   }
}

SetUpPlot = function(sample_name, y.lab, y.min, y.max, x.lab, contig_names, contig_starts, contig_ends, do_label_contigs) {
    num_contigs = length(contig_names)
    contig_centers = (contig_starts + contig_ends) / 2
    genome_length = contig_ends[num_contigs]
    suppressWarnings(par(mar=c(3.1, 3.6, 3.6, 0), mgp=c(2, -0.2, -1.1)))
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main=sample_name, xaxt="n")
    mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
    mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

    if (do_label_contigs) {
        mtext(text=contig_names[1:num_contigs], side=1, line=ifelse(c(1:num_contigs) %% 2 == 1, -0.45, 0.0),
        at = contig_centers[1:num_contigs], las=1, cex=par("cex.axis") * par("cex") * 0.7)
    }

    for (i in 1:num_contigs) {
        use.col = ifelse(i %% 2 == 1, "grey90", "white")
        rect(xleft=contig_starts[i], ybottom=y.min, xright=contig_ends[i], ytop=y.max, col=use.col, border=NA)
    }
}

# === input ===

dirPath <- '/home/lgh/my_upstream_project/hl_call_cnv/gatk_sCNV/try_again_5/denoise'
finalSegFile <- list.files(dirPath, pattern ='*.modelFinal.seg$')
finalSegdir <- file.path(dirPath, finalSegFile)

dirPath <- '/home/lgh/my_upstream_project/hl_call_cnv/gatk_sCNV/try_again_5/denoise'
denoisedCRFile <- list.files(dirPath, pattern ='*.denoised_copy_ratios.tsv$')
denoisedCRdir <- file.path(dirPath, denoisedCRFile)



i = 1

modeled_segments_file = finalSegdir[i]
denoised_copy_ratios_file = denoisedCRdir[i]
ref_dict_file = '/home/lgh/reference/genome/hg19/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.dict'
#modeled_segments_file = '/home/lgh/my_upstream_project/hl_call_cnv/gatk_sCNV/old/sandbox/FJ20190320_FDHE21H002338_T_clean.modelFinal.seg'
#denoised_copy_ratios_file = '/home/lgh/my_upstream_project/hl_call_cnv/gatk_sCNV/old/sandbox/FJ20190320_FDHE21H002338.denoisedCR.tsv'
#ref_dict_file = '/home/lgh/reference/genome/hg19/old_hg19/hg19.sub.dict'

allelic_counts_file = NULL
maximum_copy_ratio = 4

point_size_copy_ratio = 0.2

output_file='seg.png'


# === work ===

ref_dict_df = read.csv(ref_dict_file,sep='\t',skip=1,header=FALSE)

ref_dict_df = ref_dict_df%>%filter( sapply(strsplit(V2,':'),'[[',2)%in%as.character(c(1:22)) )
ref_dict_df = ref_dict_df[order(as.numeric(sapply(strsplit(ref_dict_df$V2,':'),'[[',2))),]

contig_names = sapply(strsplit(ref_dict_df[[2]],':'),'[[',2)
contig_lengths = sapply(strsplit(ref_dict_df[[3]],':'),'[[',2)
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))


num_plots = ifelse(all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file))), 2, 1) # <<<
# png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
# par(mfrow=c(num_plots, 1), cex=0.75, las=1)


denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)


#transform to linear copy ratio
denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

#determine copy-ratio midpoints
denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

#plot up to maximum_copy_ratio (or full range, if maximum_copy_ratio = Infinity)
maximum_denoised_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(denoised_copy_ratios_df[["COPY_RATIO"]])

modeled_segments_df = ReadTSV(modeled_segments_file)  # <<<

# === plot ==
if(FALSE){
    png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(num_plots, 1), cex=0.75, las=1)
    ## Contig Background for data
    SetUpPlot('sample_name', "denoised copy ratio", 0, maximum_denoised_copy_ratio, "contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatiosWithModeledSegments(denoised_copy_ratios_df, modeled_segments_df, contig_names, contig_starts, point_size_copy_ratio)
    dev.off()
}
if(TRUE){
    ## Contig Background for data
    SetUpPlot('sample_name', "denoised copy ratio", 0, maximum_denoised_copy_ratio, "contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatiosWithModeledSegments(denoised_copy_ratios_df, modeled_segments_df, contig_names, contig_starts, point_size_copy_ratio)
}





































