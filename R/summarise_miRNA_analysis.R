#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(scales)
library(VennDiagram)
library(gtools)
library(GenomicRanges)
library(grid)

#############
# arguments #
#############

arguments = commandArgs(trailingOnly = F)
# get path to script
script_index = grep("--file",arguments)
script_dir = dirname(sub("--file=","",arguments[script_index]))
# get rid of leading arguments
arguments_length = length(arguments)
arguments_index = grep("--args",arguments)[1]
if(is.na(arguments_index)){
	arguments_index=arguments_length
} else{
	arguments_index=arguments_index+1
}
arguments = arguments[arguments_index:arguments_length]

if(length(arguments)!=3){
	stop("Needs the bfx id, the data directory and the pdf directory as arguments\n", call.=FALSE)
}

bfx_id = arguments[1]
data_directory = arguments[2]
pdf_directory = arguments[3]

#bfx_id = "bfx700"
#data_directory = "miRNA/report/data"
#pdf_directory = "miRNA/report/pdf"

#############
# Functions #
#############

source(paste(script_dir,"functions_settings.R",sep="/"))

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

#########################################
# 1. plot initial genomic mapping rates #
#########################################
unfiltered_mapping_logs = list.files(path=data_directory,pattern="\\.unfiltered_mapping.txt$",full.names=TRUE)
unfiltered_mapping_summary = ldply(unfiltered_mapping_logs,function(x){parse_bowtie_log(x)})
library_names = fixLibraryNames(unfiltered_mapping_summary$library)
unfiltered_mapping_summary$library = factor(library_names,levels=unique(library_names))

unfiltered_mapping_summary_m = gather(unfiltered_mapping_summary[,c("library","reads_aligned_with_multi","reads_unaligned")],"metric","value",-library)
unfiltered_mapping_summary_plot1 = ggplot(unfiltered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="stack",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Number of reads",labels=comma) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Mapping numbers of input reads")
ggsave(paste(bfx_id,"unfiltered_reads_genomic_mapping1.pdf",sep="_"),unfiltered_mapping_summary_plot1,path=pdf_directory,width=8)

unfiltered_mapping_summary_plot2 = ggplot(unfiltered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Mapping percentages of input reads")
ggsave(paste(bfx_id,"unfiltered_reads_genomic_mapping2.pdf",sep="_"),path=pdf_directory,width=8)

#######################
# 2. plot read status #
#######################
tRNA_filter_files = list.files(path=data_directory,pattern="\\.tRNA.txt$",full.names=TRUE)
tRNA_len_table = ldply(tRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t",stringsAsFactors=F)
	if(nrow(length_table)>0){
		length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	}else{
		length_table$library = c()	
	}
	length_table
})
rRNA_filter_files = list.files(path=data_directory,pattern="\\.rRNA.txt$",full.names=TRUE)
rRNA_len_table = ldply(rRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t",stringsAsFactors=F)
	if(nrow(length_table)>0){
		length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	}else{
		length_table$library = c()	
	}
	length_table
})
otherRNA_filter_files = list.files(path=data_directory,pattern="\\.otherRNA.txt$",full.names=TRUE)
otherRNA_len_table = ldply(otherRNA_filter_files,function(file){
	length_table = read.table(file,header=T,sep="\t",stringsAsFactors=F)
	if(nrow(length_table)>0){
		length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	}else{
		length_table$library = c()	
	}
	length_table
})
miRNA_length_dist_files = list.files(path=data_directory,pattern="\\.miRNA_data.txt$",full.names=TRUE)
miRNA_data_len_table = ldply(miRNA_length_dist_files,function(file){
	length_table = read.table(file,header=T,sep="\t",stringsAsFactors=F)
	if(nrow(length_table)>0){
		length_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	}else{
		length_table$library = c()	
	}
	length_table
})
	
read_status_len_table = bind_rows(tRNA_len_table,rRNA_len_table,otherRNA_len_table,miRNA_data_len_table)
library_names = fixLibraryNames(read_status_len_table$library)
read_status_len_table$library = factor(library_names,levels=unique(library_names))
read_status_len_table$type = ifelse(read_status_len_table$type=="clean","passed",read_status_len_table$type)
read_status_len_table$type = ifelse(read_status_len_table$type=="protein_coding","mRNA",read_status_len_table$type)
read_status_len_table$type = factor(read_status_len_table$type)
read_status_len_table$type = relevel(read_status_len_table$type,"passed")

read_status_len_table = read_status_len_table %>% 
			group_by(library) %>% mutate(frac_count = count/sum(count)) %>%
			as.data.frame() 
read_status_summary = read_status_len_table %>%
			group_by(library,type) %>% summarise(count=sum(count)) %>%
			as.data.frame()

read_status_summary_plot1 = ggplot(read_status_summary,aes(x=library,y=count,fill=type)) +
geom_bar(stat="identity",position="stack",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Number of reads",label=comma) +
scale_fill_manual("Read status",values=color_brewer_qual_palette) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Read status numbers")
ggsave(paste(bfx_id,"reads_filter_status1.pdf",sep="_"),read_status_summary_plot1,path=pdf_directory)

read_status_summary_plot2 = ggplot(read_status_summary,aes(x=library,y=count,fill=type)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_manual("Read status",values=color_brewer_qual_palette) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Read status percentages")
ggsave(paste(bfx_id,"reads_filter_status2.pdf",sep="_"),read_status_summary_plot2,path=pdf_directory)

read_status_len_table$plottype = ifelse(read_status_len_table$type=="tRNA" | read_status_len_table$type=="Mt_tRNA","tRNA",
				ifelse(read_status_len_table$type=="rRNA" | read_status_len_table$type=="Mt_rRNA","rRNA",
				ifelse(read_status_len_table$type=="mRNA","mRNA",
				ifelse(read_status_len_table$type=="passed","passed","ncRNA"))))
read_status_len_table$plottype = factor(read_status_len_table$plottype,levels=c("passed","tRNA","rRNA","mRNA","ncRNA"))

read_status_len_table = read_status_len_table %>% 
			group_by(library) %>% mutate(frac_count = count/sum(count)) %>%
			as.data.frame() 
read_status_summary = read_status_len_table %>%
			group_by(library,plottype) %>% summarise(count=sum(count)) %>%
			as.data.frame()

read_status_summary_plot3 = ggplot(read_status_summary,aes(x=library,y=count,fill=plottype)) +
geom_bar(stat="identity",position="stack",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Number of reads",label=comma) +
scale_fill_manual("Read status",values=color_brewer_qual_palette) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Read status numbers")
ggsave(paste(bfx_id,"reads_filter_status3.pdf",sep="_"),read_status_summary_plot3,path=pdf_directory)

read_status_summary_plot4 = ggplot(read_status_summary,aes(x=library,y=count,fill=plottype)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_manual("Read status",values=color_brewer_qual_palette) + 
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Read status percentages")
ggsave(paste(bfx_id,"reads_filter_status4.pdf",sep="_"),read_status_summary_plot4,path=pdf_directory)

#################################
# 3. plot read status vs length #
#################################
num_columns = 4

read_status_length_plot = ggplot(read_status_len_table,aes(x=length,y=frac_count,fill=plottype,colour=plottype)) +
geom_bar(stat="identity",position="stack") +
theme_bw() +
scale_x_continuous("Read length in bp",limits=c(min(read_status_len_table$length),max(read_status_len_table$length))) +
scale_y_continuous("Fraction of reads") +
scale_fill_brewer("Read status",type="qual",palette=2,direction=-1) + 
scale_colour_brewer("Read status",type="qual",palette=2,direction=-1) + 
facet_wrap(~ library,scales="free_x",ncol=num_columns) +
theme(legend.position="bottom") +
ggtitle("Read status and length distribution")

# caclulate width and height based on the assumptions:
# - landscape, inch units
# - four columns, each facet plot has a width of 11/4 inches
# - any number of rows, each facet plot has a height of 6.5/3 (so that 12 plots and a legend fit onto one A4 page)
# - the number of columns (=width) is fixed so that only the number of rows (height) needs to be adjusted

required_width = 11/4*num_columns
required_height = 6.5/3*ceiling(length(levels(read_status_len_table$library))/num_columns)
ggsave(paste(bfx_id,"reads_filter_status_length_distribution.pdf",sep="_"),read_status_length_plot,width=required_width,height=required_height,path=pdf_directory)

##########################################
# 4. plot filtered genomic mapping rates #
##########################################
filtered_mapping_logs = list.files(path=data_directory,pattern="\\.filtered_mapping.txt$",full.names=TRUE)
filtered_mapping_summary = ldply(filtered_mapping_logs,function(x){parse_bowtie_log(x)})
library_names = fixLibraryNames(filtered_mapping_summary$library)
filtered_mapping_summary$library = factor(library_names,levels=unique(library_names))

filtered_mapping_summary_m = gather(filtered_mapping_summary[,c("library","reads_aligned_with_multi","reads_unaligned")],"metric","value",-library)
filtered_mapping_summary_plot1 = ggplot(filtered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="stack",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Numbers of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(bfx_id,"filtered_reads_genomic_mapping1.pdf",sep="_"),filtered_mapping_summary_plot1,path=pdf_directory,width=8)
ggtitle("Mapping numbers")

filtered_mapping_summary_plot2 = ggplot(filtered_mapping_summary_m,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill",colour="black") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Mapping status",type="qual",palette=2,labels=c("aligned","unaligned")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1)) +
ggtitle("Mapping percentages of miRNA reads")
ggsave(paste(bfx_id,"filtered_reads_genomic_mapping2.pdf",sep="_"),filtered_mapping_summary_plot2,path=pdf_directory,width=8)

################################################################################
# 6. read in mirdeep star tables, add library and normalised expression levels #
################################################################################
mirdeep_star_files = list.files(path=data_directory,pattern="\\.mirdeep_star.csv$",full.names=TRUE)
mirdeep_star_data_table = ldply(mirdeep_star_files,function(file){
	mirdeep_table = read.table(file,header=T,sep="\t")
	mirdeep_table$library = gsub('_R\\d$','',gsub('\\.\\S+$','',basename(file)))
	mirdeep_table
})

colnames(mirdeep_star_data_table) = c("miRNA","score","chr","strand","location_precursor","count","location_miRNA","precursor_seq","precursor_sec_struct","mature_seq_precursor_seq","star_loop_seq","UMI_reads","UMIs","UMI_distribution","library")
library_names = fixLibraryNames(mirdeep_star_data_table$library)
mirdeep_star_data_table$library = factor(library_names,levels=unique(library_names))
mirdeep_star_data_table$type = factor(ifelse(grepl('novelMiR',mirdeep_star_data_table$miRNA),"novel","known"),levels=c("known","novel"))
mirdeep_star_data_table = mirdeep_star_data_table[,c("library","miRNA","type","score","chr","strand","location_precursor","count","location_miRNA","precursor_seq","precursor_sec_struct","mature_seq_precursor_seq","star_loop_seq","UMI_reads","UMIs","UMI_distribution")]


# add normalised counts (cpm)
mirdeep_star_data_table = mirdeep_star_data_table %>%
			group_by(library) %>% mutate(norm_count = count*1000000/sum(count)) %>%
			as.data.frame()

# filter by score (known miRNAs: >=0, novelmiRNAs: >=7)
mirdeep_star_table = subset(mirdeep_star_data_table,(type=="novel" & score>=7) | (type=="known" & score>=0))

# remove novel miRNAs (at least initially)
# mirdeep_star_table = subset(mirdeep_star_table,!grepl('novelMiR',mirdeep_star_table$miRNA))

# set character vector
mirdeep_star_table$UMI_distribution = as.character(mirdeep_star_table$UMI_distribution)

# output filtered tables
#for(lib in unique(mirdeep_star_table$library)){
#	write.table(subset(mirdeep_star_table,library==lib),paste(data_directory,paste(bfx_id,"putative_novel_miRNAs.csv",sep="_"),
#}

################################
# 7. summarise profiled miRNAs #
################################
mirdeep_star_summary = mirdeep_star_table %>% 
			complete(library,type,fill=list(norm_count=0)) %>%
			group_by(library,type) %>% summarise(All=sum(norm_count>0),count_0=sum(norm_count>0 & norm_count<=10),count_10=sum(norm_count>10 & norm_count<=100),count_100=sum(norm_count>100 & norm_count<=1000),count_1000=sum(norm_count>1000 & norm_count<=10000),count_10000=sum(norm_count>10000)) %>%		
			as.data.frame()
			
colnames(mirdeep_star_summary) = c("library","type","All","0-10","11-100","101-1000","1001-10000",">10000")
mirdeep_star_summary_m = gather(mirdeep_star_summary,"metric","value",-c(library,type))
mirdeep_star_summary_m$metric = factor(mirdeep_star_summary_m$metric,levels=c("All","0-10","11-100","101-1000","1001-10000",">10000"))

number_of_profiled_miRNA_plot = ggplot(subset(mirdeep_star_summary_m,metric=="All"),aes(x=library,y=value,fill=library,alpha=type)) +
geom_bar(stat="identity",colour="black") +
scale_x_discrete("Library") +
scale_y_continuous("Number of identified miRNAs") +
theme_bw() +
scale_fill_manual("Library",values=color_brewer_qual_palette) +
scale_alpha_manual("Type of miRNA",values=c(1,0.2),labels=c("known","putative novel")) +
#scale_linetype_manual("Type of miRNA",values=c("solid","dashed"),labels=c("known","putative novel"),guide=guide_legend(override.aes = list(fill="white"))) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1),legend.position="bottom") +
ggtitle("Identified MiRNAs")
ggsave(paste(bfx_id,"number_of_profiled_miRNAs.pdf",sep="_"),number_of_profiled_miRNA_plot,path=pdf_directory,width=9,height=8)

miRNA_read_counts_plot = ggplot(subset(mirdeep_star_summary_m,metric!="All"),aes(x=library,y=value,fill=library,alpha=type)) +
geom_bar(stat="identity",colour="black") +
scale_x_discrete("Library") +
scale_y_continuous("Number of identified miRNAs") +
theme_bw() +
facet_wrap(~ metric,scales="free_y",ncol=3) +
scale_fill_manual("Library",values=color_brewer_qual_palette) +
scale_alpha_manual("Type of miRNA",values=c(1,0.2),labels=c("known","putative novel")) +
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") +
ggtitle("MiRNAs by expression level")
ggsave(paste(bfx_id,"number_of_profiled_miRNAs_with_counts.pdf",sep="_"),miRNA_read_counts_plot,width=9,height=8,path=pdf_directory)

mirRNA_ranked_table = mirdeep_star_table[,c("library","count","norm_count")] %>% group_by(library) %>% 
mutate(rank=rank(-count,ties.method="random")) %>% arrange(rank) %>% mutate(frac_cum_sum=cumsum(count)/sum(count))
max_rank = max(mirRNA_ranked_table$rank)

miRNA_ranked_counts_plot = ggplot(mirRNA_ranked_table,aes(x=rank,y=frac_cum_sum,colour=library)) + 
geom_line() +
theme_bw() +
scale_x_continuous("miRNAs ranked from high to low expressed",trans=log10_trans(),limits=c(1,max_rank),breaks=c(1,10,50,100,seq(200,round(max_rank,-2),100))) +
scale_y_continuous("Fraction of miRNA reads",limits=c(0,1)) +
scale_colour_manual("Library",values=color_brewer_qual_palette) +
theme(legend.position="bottom") +
ggtitle("MiRNAs ranked by expression (high to low)")
ggsave(paste(bfx_id,"miRNAs_ranked_counts_plot.pdf",sep="_"),miRNA_ranked_counts_plot,width=9,height=8,path=pdf_directory)

#########################################################
# 8. expression tables and correlation for known miRNAs #
#########################################################

mirdeep_star_table_known = subset(mirdeep_star_table,type=="known")
miRNA_expression_table = spread(mirdeep_star_table_known[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"known_miRNAs_cpm_normalised_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

miRNA_expression_table = spread(mirdeep_star_table_known[,c("library","miRNA","count")],library,count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"known_miRNAs_raw_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

miRNA_expression_table = spread(mirdeep_star_table_known[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
miRNA_expression_table$miRNA = NULL

if(ncol(miRNA_expression_table)>1){
	plotCorrelation(miRNA_expression_table, conditions="", cormethod="spearman", cluster=FALSE, title = "Spearman correlation of counts of known miRNAs",filename=paste(pdf_directory,paste(bfx_id,"known_miRNA_expression_spearman_correlation.pdf",sep="_"),sep="/"))
	plotCorrelation(miRNA_expression_table, conditions="", cormethod="pearson", cluster=FALSE, title = "Pearson correlation of counts of known miRNAs",filename=paste(pdf_directory,paste(bfx_id,"known_miRNA_expression_pearson_correlation.pdf",sep="_"),sep="/"))
} else {
	ggsave(grid.text("No known_miRNA_expression_spearman_correlation plot possible"),file=paste(pdf_directory,paste(bfx_id,"known_miRNA_expression_spearman_correlation.pdf",sep="_"),sep="/"))
	ggsave(grid.text("No known_miRNA_expression_pearson_correlation plot possible"),file=paste(pdf_directory,paste(bfx_id,"known_miRNA_expression_pearson_correlation.pdf",sep="_"),sep="/"))
}

################################################################################
# 9. deal with novel miRNAs (find and group overlapping, filter and output)   #
################################################################################

# a) create data frame with library, miRNA, chr, start, stop, strand, accession
mirdeep_star_table_novel = subset(mirdeep_star_table,type=="novel")
novel_miRNAs = mirdeep_star_table_novel[,c("library","miRNA","chr","location_precursor","strand")]
novel_miRNAs = separate(novel_miRNAs,location_precursor,c("start","end"),sep="-")
novel_miRNAs$start = as.numeric(novel_miRNAs$start)
novel_miRNAs$end = as.numeric(novel_miRNAs$end)
novel_miRNAs$accession = paste(novel_miRNAs$library,novel_miRNAs$miRNA,sep="#")

# b) create GRanges object, merge with reduce and assign to groups
novel_miRNAs_gr = GRanges(seqnames=Rle(novel_miRNAs$chr),
			ranges=IRanges(start=novel_miRNAs$start,end=novel_miRNAs$end,names=novel_miRNAs$accession),
			strand=Rle(novel_miRNAs$strand),library=novel_miRNAs$library)

novel_miRNAs_merged_gr = reduce(novel_miRNAs_gr,ignore.strand=F,min.gapwidth=0L)
novel_miRNAs_merged_gr$group = paste("novel-mir",1:length(novel_miRNAs_merged_gr),sep="-")
overlaps_between_novel_miRNAs = findOverlaps(novel_miRNAs_gr,novel_miRNAs_merged_gr,ignore.strand=F)
novel_miRNAs_gr$group = novel_miRNAs_merged_gr[subjectHits(overlaps_between_novel_miRNAs),]$group
mirdeep_star_table_novel$group = novel_miRNAs_gr$group
mirdeep_star_table_novel = merge(mirdeep_star_table_novel,as.data.frame(novel_miRNAs_merged_gr)[,c("group","start","end")],by="group")


# c) collapse original data frame by group and output a table for novel miRNAs
novel_miRNAs_table = mirdeep_star_table_novel %>% group_by(group) %>%
summarise(type=first(type),
	chr=first(chr),
	start=first(start),
	end=first(end),
	strand=first(strand),
	num_libraries=length(library),
	mean_score=mean(score),
	median_score=median(score),
	mean_norm_count=mean(norm_count),
	median_norm_count=mean(norm_count),
	libraries=paste(library,collapse=","),
	put_mature_seq_precursor_seqs=paste(unique(mature_seq_precursor_seq),collapse=","),
	put_precursor_seqs=paste(unique(precursor_seq),collapse=",")
	) %>%
as.data.frame()
novel_miRNAs_table_colnames = colnames(novel_miRNAs_table)
novel_miRNAs_table_colnames[1] = "put_miRNA"
colnames(novel_miRNAs_table) = novel_miRNAs_table_colnames

novel_miRNAs_table$put_miRNA = factor(novel_miRNAs_table$put_miRNA,levels=mixedsort1(unique(as.character(novel_miRNAs_table$put_miRNA)),decreasing=T))
novel_miRNAs_table$chr = factor(novel_miRNAs_table$chr,levels=mixedsort1(unique(as.character(novel_miRNAs_table$chr))))
novel_miRNAs_table = novel_miRNAs_table[order(novel_miRNAs_table$put_miRNA),]

write.table(novel_miRNAs_table,paste(data_directory,paste(bfx_id,"putative_novel_miRNAs.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

####################################################################
# 10. expression tables and correlation for putative novel miRNAs  #
####################################################################

# in 9. it can happen that from the same sample multiple predictions were merged so that now we have multiple count values per merged putative miRNA and sample 
#
# in this case take mean values of counts
novel_miRNAs_count_table = mirdeep_star_table_novel[,c("library","group","count","norm_count")] %>% group_by(library,group) %>% summarise(count=mean(count),norm_count=mean(norm_count)) %>% as.data.frame()
colnames(novel_miRNAs_count_table) = c("library","miRNA","count","norm_count")

miRNA_expression_table = spread(novel_miRNAs_count_table[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"novel_miRNAs_cpm_normalised_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

miRNA_expression_table = spread(novel_miRNAs_count_table[,c("library","miRNA","count")],library,count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
write.table(miRNA_expression_table,paste(data_directory,paste(bfx_id,"novel_miRNAs_raw_counts.csv",sep="_"),sep="/"),quote=F,sep="\t",row.names=F,col.names=T)

miRNA_expression_table = spread(novel_miRNAs_count_table[,c("library","miRNA","norm_count")],library,norm_count,fill=0)
row.names(miRNA_expression_table) = miRNA_expression_table$miRNA
miRNA_expression_table$miRNA = NULL

if(ncol(miRNA_expression_table)>1){
	plotCorrelation(miRNA_expression_table, conditions="", cormethod="spearman", cluster=FALSE, title = "Spearman correlation of counts of novel miRNAs", filename=paste(pdf_directory,paste(bfx_id,"novel_miRNA_expression_spearman_correlation.pdf",sep="_"),sep="/"))
	plotCorrelation(miRNA_expression_table, conditions="", cormethod="pearson", cluster=FALSE,  title = "Spearman correlation of counts of novel miRNAs", filename=paste(pdf_directory,paste(bfx_id,"novel_miRNA_expression_pearson_correlation.pdf",sep="_"),sep="/"))
} else {
	ggsave(grid.text("No novel_miRNA_expression_spearman_correlation plot possible"),file=paste(pdf_directory,paste(bfx_id,"novel_miRNA_expression_spearman_correlation.pdf",sep="_"),sep="/"))
	ggsave(grid.text("No novel_miRNA_expression_pearson_correlation plot possible"),file=paste(pdf_directory,paste(bfx_id,"novel_miRNA_expression_pearson_correlation.pdf",sep="_"),sep="/"))
}

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

# only done for known miRNA
mirdeep_star_table_known = subset(mirdeep_star_table,type=="known")

UMIs_count_lm = mirdeep_star_table_known %>% 
	group_by(library) %>% 
	do(mod = lm(count ~ UMIs, data = subset(.,count<=64000)),x_pos=max(.$count),y_pos=max(.$UMIs)) %>% 
	mutate(R2 = summary(mod)$r.squared,intercept = coef(mod)[1],slope=coef(mod)[2]) %>% as.data.frame()
UMIs_count_lm$mod = NULL
UMIs_count_lm$x_pos = unlist(UMIs_count_lm$x_pos)
UMIs_count_lm$y_pos = unlist(UMIs_count_lm$y_pos)
UMIs_count_lm$R2 = round(UMIs_count_lm$R2,2)


expression_UMI_relation_plot = ggplot(mirdeep_star_table_known,aes(x=count,y=UMIs,colour=library)) + 
geom_point() + 
geom_abline(data=UMIs_count_lm,aes(intercept=intercept,slope=slope)) +
facet_wrap(~library) +
#geom_text(data=UMIs_count_lm,aes(x=x_pos,y=y_pos,label=paste("R^2 == ",R2)),vjust=0,hjust=1,parse=TRUE) +
theme_bw() +
scale_x_continuous("Normalised count",trans="log10") +
scale_y_continuous("Number of UMIs",trans="log10") +
scale_colour_discrete("Library") +
facet_wrap(~library) +
theme(legend.position="bottom") +
ggtitle("Expression level vs UMIs")

required_width = 11/4*num_columns
required_height = 6.5/3*ceiling(length(levels(mirdeep_star_table$library))/num_columns)
ggsave(paste(bfx_id,"expression_UMI_relation_plot.pdf",sep="_"),expression_UMI_relation_plot,width=required_width,height=required_height,path=pdf_directory)







