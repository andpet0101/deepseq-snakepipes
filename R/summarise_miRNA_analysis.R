#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(scales)
library(VennDiagram)

#############
# Functions #
#############

bp_labels = function(l){format(big.mark=",",l)}
kb_labels = function(l,r=0){format(big.mark=",",round(l/1000,r))}
Mb_labels = function(l,r=0){format(big.mark=",",round(l/1000000,r))}
Gb_labels = function(l,r=0){format(big.mark=",",round(l/1000000000,r))}

# parses the bowtie logs
parse_bowtie_log = function(bowtie_log){
	bowtie_lines = readLines(bowtie_log)
	bowtie_lines = bowtie_lines[grep("^#",bowtie_lines)]

	if(length(bowtie_lines)==0){
		return()
	}

	# file
	bowtie_summary = data.frame(library=gsub('\\.\\S+$','',basename(bowtie_log)))
	
	# reads processed
	line = bowtie_lines[grep("^# reads processed",bowtie_lines)][1]
	bowtie_summary$total_reads = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
	bowtie_summary[is.na(bowtie_summary$total_reads),"total_reads"] = 0
	
	# reads aligned
	line = bowtie_lines[grep("^# reads with at least one reported alignment",bowtie_lines)][1]
	bowtie_summary$reads_aligned = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
	bowtie_summary[is.na(bowtie_summary$reads_aligned),"reads_aligned"] = 0

	# reads failed to align
	line = bowtie_lines[grep("^# reads that failed to align",bowtie_lines)][1]
	bowtie_summary$reads_unaligned = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
	bowtie_summary[is.na(bowtie_summary$reads_unaligned),"reads_unaligned"] = 0

	# alignment multimappers (according bowtie -m)
	line = bowtie_lines[grep("^# reads with alignments suppressed due to",bowtie_lines)][1]
	bowtie_summary$reads_multimappers = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
	bowtie_summary[is.na(bowtie_summary$reads_multimappers),"reads_multimappers"] = 0

	# reads aligned incl multimappers
	bowtie_summary$reads_aligned_with_multi = bowtie_summary$reads_aligned+bowtie_summary$reads_multimappers	
	
	bowtie_summary
}


# from Matthias' DESeq methods
plotCorrelation <- function(deseqmatrix, conditions = "", cormethod = "pearson", cluster = FALSE, filename = "") {
  require(pheatmap)
  require(RColorBrewer)
  
  cormatrix <- cor(deseqmatrix, use = "pairwise.complete.obs", method = cormethod)
  cornumbers <- round(cormatrix, 3)
  # colnames(cormatrix) <- NULL
  if (cormethod == "pearson") {
    colors <- colorRampPalette(brewer.pal(8, "Blues"))(255)
    title <- "Pearson Correlation"  
  }
  if (cormethod == "spearman") {
    colors <- colorRampPalette(brewer.pal(8, "Reds"))(255)
    title <- "Spearman Correlation"
  }
  clust <- ifelse(cluster, TRUE, FALSE)
  
  if (!is.data.frame(conditions)) {
    conditions <- NA
  }
  
  cols <- ncol(cormatrix)
  if (cols <= 10) {
    fontsize <- 8
    pheatmap(cormatrix,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontsize, fontsize_row = fontsize+2, fontsize_col = fontsize+2,
             fontsize_number = fontsize+1,
             # cellwidth = 20, cellheight = 20,
             display_numbers = cornumbers, number_color = "Black",
             filename = filename
    )
  } else if((10 < cols) & (cols <= 20)) {
    fontsize <- 7
    pheatmap(cormatrix,
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontsize, fontsize_row = fontsize, fontsize_col = fontsize,
             fontsize_number = 0.9 * fontsize,
             cellwidth = 20, cellheight = 20,
             display_numbers = cornumbers, number_color = "Black",
             filename = filename
    )
  } else if(cols > 20) {
    fontsize <- 4
    pheatmap(cormatrix, 
             cluster_rows = clust, cluster_cols = clust,
             color=colors, annotation_col = conditions,
             main = title, border_color = NA,
             fontsize = fontsize, fontsize_row = fontsize, fontsize_col = fontsize,
             cellwidth = 20, cellheight = 20,
             display_numbers = F,
             filename = filename
    )
  }
}

# colour brewer palette for 24 samples
color_brewer_qual_palette = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7","#ffffb3",
"#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")

# repeat 100 times = 2400 samples max
color_brewer_qual_palette = rep(color_brewer_qual_palette,100)

########
# Main #
########

arguments = commandArgs(trailingOnly = TRUE)
if(length(arguments)!=3){
	stop("Needs the bfx id, the data directory and the pdf directory as arguments\n", call.=FALSE)
}

bfx_id = arguments[1]
data_directory = arguments[2]
pdf_directory = arguments[3]

#########################################
# 1. plot initial genomic mapping rates #
#########################################
unfiltered_mapping_logs = list.files(path=data_directory,pattern=".unfiltered_mapping.txt$",full.names=TRUE)
unfiltered_mapping_summary = ldply(unfiltered_mapping_logs,function(x){parse_bowtie_log(x)})


unfiltered_mapping_summary_m = gather(unfiltered_mapping_summary[,c("library","reads_aligned_with_multi","reads_unaligned")],"metric","value",-library)
unfiltered_mapping_summary_plot1 = ggplot(unfiltered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"unfiltered_reads_genomic_mapping1.pdf",sep="_"),unfiltered_mapping_summary_plot1,path=pdf_directory)

unfiltered_mapping_summary_m = gather(unfiltered_mapping_summary[,c("library","reads_aligned","reads_multimappers","reads_unaligned")],"metric","value",-library)
unfiltered_mapping_summary_plot2 = ggplot(unfiltered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","multi-mapper (>1x)","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"unfiltered_reads_genomic_mapping2.pdf",sep="_"),path=pdf_directory)

#######################
# 2. plot read status #
#######################
tRNA_filter_files = list.files(path=data_directory,pattern=".tRNA.txt$",full.names=TRUE)
tRNA_len_table = ldply(tRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	length_table
})
rRNA_filter_files = list.files(path=data_directory,pattern=".rRNA.txt$",full.names=TRUE)
rRNA_len_table = ldply(rRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	length_table
})
otherRNA_filter_files = list.files(path=data_directory,pattern=".otherRNA.txt$",full.names=TRUE)
otherRNA_len_table = ldply(otherRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	length_table
})
miRNA_length_dist_files = list.files(path=data_directory,pattern=".miRNA_data.txt$",full.names=TRUE)
miRNA_data_len_table = ldply(miRNA_length_dist_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	length_table
})
	

read_status_len_table = bind_rows(tRNA_len_table,rRNA_len_table,otherRNA_len_table,miRNA_data_len_table)
read_status_len_table$plottype = ifelse(read_status_len_table$type=="tRNA" | read_status_len_table$type=="Mt_tRNA","tRNA",
				ifelse(read_status_len_table$type=="rRNA" | read_status_len_table$type=="Mt_rRNA","rRNA",
				ifelse(read_status_len_table$type=="protein_coding","mRNA",
				ifelse(read_status_len_table$type=="clean","passed","ncRNA"))))
read_status_len_table$plottype = factor(read_status_len_table$plottype,levels=c("passed","tRNA","rRNA","mRNA","ncRNA"))

read_status_len_table = read_status_len_table %>% 
			group_by(library) %>% mutate(frac_count = count/sum(count)) %>%
			as.data.frame() 
read_status_summary = read_status_len_table %>%
			group_by(library,plottype) %>% summarise(count=sum(count)) %>%
			as.data.frame()

# total counts
read_status_summary_plot1 = ggplot(read_status_summary,aes(x=library,y=count,fill=plottype)) +
geom_bar(stat="identity",position="stack",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Number of reads",label=comma) +
scale_fill_brewer("Read status",type="qual",palette=2,direction=-1) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"reads_filter_status1.pdf",sep="_"),read_status_summary_plot1,path=pdf_directory)

# percent
read_status_summary_plot2 = ggplot(read_status_summary,aes(x=library,y=count,fill=plottype)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Read status",type="qual",palette=2,direction=-1) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"reads_filter_status2.pdf",sep="_"),read_status_summary_plot2,path=pdf_directory)

#################################
# 3. plot read status vs length #
#################################
read_status_length_plot = ggplot(read_status_len_table,aes(x=length,y=frac_count,fill=plottype,colour=plottype)) +
geom_bar(stat="identity",position="stack") +
theme_bw() +
scale_x_continuous("Read length in bp",limits=c(min(read_status_len_table$length),max(read_status_len_table$length))) +
scale_y_continuous("Fraction of reads") +
scale_fill_brewer("Read status",type="qual",palette=2,direction=-1) + 
scale_colour_brewer("Read status",type="qual",palette=2,direction=-1) + 
facet_wrap(~ library,scales="free_x",ncol=3)

ggsave(paste(bfx_id,"reads_filter_status_length_distribution.pdf",sep="_"),read_status_length_plot,width=10,height=7,path=pdf_directory)

##########################################
# 4. plot filtered genomic mapping rates #
##########################################
filtered_mapping_logs = list.files(path=data_directory,pattern=".filtered_mapping.txt$",full.names=TRUE)
filtered_mapping_summary = ldply(filtered_mapping_logs,function(x){parse_bowtie_log(x)})

filtered_mapping_summary_m = gather(filtered_mapping_summary[,c("library","reads_aligned_with_multi","reads_unaligned")],"metric","value",-library)
filtered_mapping_summary_plot1 = ggplot(filtered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"filtered_reads_genomic_mapping1.pdf",sep="_"),filtered_mapping_summary_plot1,path=pdf_directory)

filtered_mapping_summary_m = gather(filtered_mapping_summary[,c("library","reads_aligned","reads_multimappers","reads_unaligned")],"metric","value",-library)
filtered_mapping_summary_plot2 = ggplot(filtered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","multi-mapper (>1)","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"filtered_reads_genomic_mapping2.pdf",sep="_"),filtered_mapping_summary_plot2,path=pdf_directory)

################################################################################
# 6. read in mirdeep star tables, add library and normalised expression levels #
################################################################################
mirdeep_star_files = list.files(path=data_directory,pattern=".mirdeep_star.csv$",full.names=TRUE)
mirdeep_star_data_table = ldply(mirdeep_star_files,function(file){
	mirdeep_table = read.table(file,header=T,sep="\t")
	mirdeep_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	mirdeep_table
})

colnames(mirdeep_star_data_table) = c("miRNA","score","chr","strand","location_precursor","count","location_miRNA","precursor_seq","precursor_sec_struct","mature_seq_precursor_seq","star_loop_seq","UMI_reads","UMIs","UMI_distribution","library")
mirdeep_star_data_table = mirdeep_star_data_table[,c("library","miRNA","score","chr","strand","location_precursor","count","location_miRNA","precursor_seq","precursor_sec_struct","mature_seq_precursor_seq","star_loop_seq","UMI_reads","UMIs","UMI_distribution")]

# add normalised counts
mirdeep_star_data_table = mirdeep_star_data_table %>%
			group_by(library) %>% mutate(norm_count = count*1000000/sum(count)) %>%
			as.data.frame()

# filter by score (known miRNAs: >=0, novelmiRNAs: >=7)
mirdeep_star_table = subset(mirdeep_star_data_table,(grepl('novelMiR',mirdeep_star_data_table$miRNA) & mirdeep_star_data_table$score>=7) | (!grepl('novelMiR',mirdeep_star_data_table$miRNA) & mirdeep_star_data_table$score>=0))

# remove novel miRNAs (at least initially)
mirdeep_star_table = subset(mirdeep_star_table,!grepl('novelMiR',mirdeep_star_table$miRNA))

# set character vector
mirdeep_star_table$UMI_distribution = as.character(mirdeep_star_table$UMI_distribution)

################################
# 7. summarise profiled miRNAs #
################################
mirdeep_star_summary = mirdeep_star_table %>% 
			group_by(library) %>% summarise(All=sum(norm_count>0),count_0=sum(norm_count>0 & norm_count<=10),count_10=sum(norm_count>10 & norm_count<=100),count_100=sum(norm_count>100 & norm_count<=1000),count_1000=sum(norm_count>1000 & norm_count<=10000),count_10000=sum(norm_count>10000)) %>%
			as.data.frame()
			
colnames(mirdeep_star_summary) = c("library","All","0-10","11-100","101-1000","1001-10000",">10000")
mirdeep_star_summary_m = gather(mirdeep_star_summary,"metric","value",-library)
mirdeep_star_summary_m$metric = factor(mirdeep_star_summary_m$metric,levels=c("All","0-10","11-100","101-1000","1001-10000",">10000"))

number_of_profiled_miRNA_plot = ggplot(subset(mirdeep_star_summary_m,metric=="All"),aes(x=library,y=value,fill=library)) +
geom_bar(stat="identity",colour="black") +
scale_x_discrete("Library") +
scale_y_continuous("Number of profiled miRNAs") +
theme_bw() +
scale_fill_manual("Library",values=color_brewer_qual_palette) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"number_of_profiled_miRNAs.pdf",sep="_"),number_of_profiled_miRNA_plot,path=pdf_directory)

miRNA_read_counts_plot = ggplot(subset(mirdeep_star_summary_m,metric!="All"),aes(x=library,y=value,fill=library)) +
geom_bar(stat="identity",colour="black") +
scale_x_discrete("Library") +
scale_y_continuous("Number of profiled miRNAs") +
theme_bw() +
facet_wrap(~ metric,scales="free_y",ncol=3) +
scale_fill_manual("Library",values=color_brewer_qual_palette) +
theme(axis.text.x=element_blank())
ggsave(paste(bfx_id,"number_of_profiled_miRNAs_with_counts.pdf",sep="_"),miRNA_read_counts_plot,width=10,height=7,path=pdf_directory)

mirRNA_ranked_table = mirdeep_star_table[,c("library","count","norm_count")] %>% group_by(library) %>% 
mutate(rank=rank(-count,ties.method="random")) %>% arrange(rank) %>% mutate(frac_cum_sum=cumsum(count)/sum(count))

miRNA_ranked_counts_plot = ggplot(mirRNA_ranked_table,aes(x=rank,y=frac_cum_sum,colour=library)) + 
geom_line() +
theme_bw() +
scale_x_continuous("MiRNA ranked from high to low expressed") +
scale_y_continuous("Fraction of sequence data",limits=c(0,1)) +
scale_colour_manual("Library",values=color_brewer_qual_palette)
ggsave(paste(bfx_id,"miRNAs_ranked_counts_plot.pdf",sep="_"),miRNA_ranked_counts_plot,width=8,height=7,path=pdf_directory)

#miRNA_ranked_counts_zoom_plot = ggplot(mirRNA_ranked_table,aes(x=rank,y=frac_cum_sum,colour=library)) + 
#geom_line() +
#theme_bw(10) +
#scale_x_continuous("MiRNA ranked from high to low expressed",limits=c(0,30)) +
#scale_y_continuous("Fraction of sequence data",limits=c(0,1)) +
#scale_colour_manual("Library",values=color_brewer_qual_palette) +
#theme(legend.position="none")
#
#miRNA_ranked_counts_adv_plot = miRNA_ranked_counts_plot + 
#annotation_custom(ggplot_gtable(ggplot_build(miRNA_ranked_counts_zoom_plot)),xmin=0.2*max(mirRNA_ranked_table$rank),ymin=0.1,ymax=0.6,xmax=max(mirRNA_ranked_table$rank))
#ggsave(paste(bfx_id,"miRNAs_ranked_counts_plot_with_zoom.pdf",sep="_"),miRNA_ranked_counts_adv_plot,width=9,height=7,path=pdf_directory)

########################
# 8. expression tables #
########################

miRNA_expression_table = spread(mirdeep_star_table[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"cpm_normalised_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

miRNA_expression_table = spread(mirdeep_star_table[,c("library","miRNA","count")],library,count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"raw_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

#############################
# 9. expression correlation #
#############################

miRNA_expression_table = spread(mirdeep_star_table[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
miRNA_expression_table$miRNA = NULL
plotCorrelation(miRNA_expression_table, conditions="", cormethod="spearman", cluster=FALSE, filename=paste(pdf_directory,paste(bfx_id,"miRNA_expression_spearman_correlation.pdf",sep="_"),sep="/"))
plotCorrelation(miRNA_expression_table, conditions="", cormethod="pearson", cluster=FALSE, filename=paste(pdf_directory,paste(bfx_id,"miRNA_expression_pearson_correlation.pdf",sep="_"),sep="/"))

#########################
# 10. UMI distributions #
#########################
#miRNA_UMI_table = mirdeep_star_table[,c("library","miRNA","count","norm_count","UMI_distribution")]

#miRNA_UMI_table = miRNA_UMI_table %>% rowwise() %>% 
#	do({
#	umis = unlist(strsplit(c(.$UMI_distribution), ","))
#	umi_seqs = as.character(gsub(":\\d+","",umis))
#	umi_counts = as.numeric(gsub("\\S+:","",umis))
#	data.frame(library=rep(.$library,length(umi_counts)),miRNA=rep(.$miRNA,length(umi_counts)),count=rep(.$count,length(umi_counts)),
#	norm_count=rep(.$norm_count,length(umi_counts)),UMI_distribution=rep("",length(umi_counts)),UMI_seq=umi_seqs,UMI_count=umi_counts)
#	}) %>% as.data.frame()

UMIs_norm_count_lm = mirdeep_star_table %>% 
	group_by(library) %>% 
	do(mod = lm(UMIs ~ norm_count, data = .),x_pos=max(.$norm_count),y_pos=max(.$UMIs)) %>% 
	mutate(R2 = summary(mod)$r.squared) %>% as.data.frame()
UMIs_norm_count_lm$mod = NULL
UMIs_norm_count_lm$x_pos = unlist(UMIs_norm_count_lm$x_pos)
UMIs_norm_count_lm$y_pos = unlist(UMIs_norm_count_lm$y_pos)
UMIs_norm_count_lm$R2 = round(UMIs_norm_count_lm$R2,2)


expression_UMI_relation_plot = ggplot(mirdeep_star_table,aes(x=norm_count,y=UMIs,colour=library)) + 
geom_point() + 
geom_smooth(method="lm",se=F) +
geom_text(data=UMIs_norm_count_lm,aes(x=x_pos,y=y_pos,label=paste("R^2 == ",R2)),vjust=0,hjust=1,parse=TRUE) +
theme_bw() +
scale_x_continuous("Normalised count") +
scale_y_continuous("Number of UMIs") +
scale_colour_discrete("Library")

ggsave(paste(bfx_id,"expression_UMI_relation_plot.pdf",sep="_"),expression_UMI_relation_plot,width=9,height=7,path=pdf_directory)





