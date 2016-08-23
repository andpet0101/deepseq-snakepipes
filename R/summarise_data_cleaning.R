#!/usr/bin/env Rscript 

library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(scales)

#############
# Arguments #
#############

arguments = commandArgs(trailingOnly = TRUE)
if(length(arguments)!=3){
	stop("Needs the bfx id, the data directory and the pdf directory as arguments\n", call.=FALSE)
}


bfx_id = arguments[1]
data_directory = arguments[2]
pdf_directory = arguments[3]

#############
# Functions #
#############

bp_labels = function(l){format(big.mark=",",l)}
kb_labels = function(l,r=0){format(big.mark=",",round(l/1000,r))}
Mb_labels = function(l,r=0){format(big.mark=",",round(l/1000000,r))}
Gb_labels = function(l,r=0){format(big.mark=",",round(l/1000000000,r))}


parse_cutadapt_summary = function(cutadapt_file){
	cutadapt_lines = readLines(cutadapt_file)
	cutadapt_lines = cutadapt_lines[grep("^\\s*(Total|Read|Quality|Pairs)",cutadapt_lines,perl=TRUE)]
	
	if(length(cutadapt_lines)==0){
		return()
	}	
	
	# file
	cutadapt_summary = data.frame(library=gsub('\\.cutadapt\\.txt$','',basename(cutadapt_file)))		
	
	# paired end data
	if(length(grep("^(Total read pairs processed)",cutadapt_lines))>0){
		# total read pairs and reads
		line = cutadapt_lines[grep("^Total read pairs processed",cutadapt_lines)][1]
		read_pairs = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
		read_pairs[is.na(read_pairs)] = 0
		cutadapt_summary$total_reads = read_pairs*2
		
		# read pairs filtered
		line = cutadapt_lines[grep("^Pairs written",cutadapt_lines)][1]
		read_pairs = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		read_pairs[is.na(read_pairs)] = 0
		cutadapt_summary$reads_passed = read_pairs*2
		
		# with adapter (R1)
		line = cutadapt_lines[grep("^\\s+Read 1 with adapter",cutadapt_lines,perl = TRUE)][1]
		with_adapter_R1 = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		print(with_adapter_R1)
		with_adapter_R1[is.na(with_adapter_R1)] = 0
		
		# with adapter (R2)
		line = cutadapt_lines[grep("^\\s+Read 2 with adapter",cutadapt_lines,perl = TRUE)][1]
		with_adapter_R2 = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		with_adapter_R2[is.na(with_adapter_R2)] = 0

		cutadapt_summary$with_adapter = with_adapter_R1 + with_adapter_R2

		# too short after adapter (= primer dimer)
		line = cutadapt_lines[grep("^Pairs that were too short",cutadapt_lines)][1]
		read_pairs = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		read_pairs[is.na(read_pairs)] = 0
		cutadapt_summary$too_short = read_pairs*2


	# single end data
	} else {	
		# total reads
		line = cutadapt_lines[grep("^Total reads processed",cutadapt_lines)][1]
		cutadapt_summary$total_reads = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
		cutadapt_summary[is.na(cutadapt_summary$total_reads),"total_reads"] = 0

		# reads filtered
		line = cutadapt_lines[grep("^Reads written",cutadapt_lines)][1]
		cutadapt_summary$reads_passed = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		cutadapt_summary[is.na(cutadapt_summary$reads_passed),"reads_passed"] = 0	

		# with adapter
		line = cutadapt_lines[grep("^Reads with adapters",cutadapt_lines,perl = TRUE)][1]
		cutadapt_summary$with_adapter_R1 = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		cutadapt_summary[is.na(cutadapt_summary$with_adapter_R1),"with_adapter_R1"] = 0
		
		# too short after adapter (= primer dimer)
		line = cutadapt_lines[grep("^Reads that were too short",cutadapt_lines)][1]
		cutadapt_summary$too_short = as.numeric(gsub(',','',gsub('.+\\s(\\S+)\\s\\(\\S+','\\1',line)))
		cutadapt_summary[is.na(cutadapt_summary$too_short),"too_short"] = 0	
	}


	# total bp
	line = cutadapt_lines[grep("^Total basepairs processed",cutadapt_lines)][1]
	cutadapt_summary$total_bp = as.numeric(gsub('bp','',gsub(',','',gsub('.+\\s(\\S+) bp$','\\1',line))))
	cutadapt_summary[is.na(cutadapt_summary$total_bp),"total_bp"] = 0	

	# bp filtered
	line = cutadapt_lines[grep("^Total written",cutadapt_lines)][1]
	cutadapt_summary$bp_passed = as.numeric(gsub('bp','',gsub(',','',gsub('.+\\s(\\S+)\\sbp\\s\\(\\S+','\\1',line))))
	cutadapt_summary[is.na(cutadapt_summary$bp_passed),"bp_passed"] = 0

	cutadapt_summary		
}

parse_umi_summary = function(umi_file){
	umi_lines = readLines(umi_file)
	umi_lines = umi_lines[grep("^(Total|Fragment)",umi_lines)]
	
	if(length(umi_lines)==0){
		return()
	}	

	# file
	umi_summary = data.frame(library=gsub('\\.umi\\.txt$','',basename(umi_file)))	
	
	# total reads
	line = umi_lines[grep("^(Total fragments)",umi_lines)][1]
	umi_summary$total_reads = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
	umi_summary[is.na(umi_summary$total_reads),"total_reads"] = 0	

	# reads with UMI
	line = umi_lines[grep("^(Fragments with UMI)",umi_lines)][1]
	umi_summary$with_umi = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
	umi_summary[is.na(umi_summary$with_umi),"with_umi"] = 0

	# reads with UMI - total length (without the actual UMI)
	line = umi_lines[grep("^(Fragment length with UMI)",umi_lines)][1]
	umi_summary$with_umi_bp = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
	umi_summary[is.na(umi_summary$with_umi_bp),"with_umi_bp"] = 0

	# reads too short after UMI
	line = umi_lines[grep("^(Fragments that were too short)",umi_lines)][1]
	umi_summary$too_short = as.numeric(gsub(',','',gsub('.+\\s(\\S+)$','\\1',line)))
	umi_summary[is.na(umi_summary$too_short),"too_short"] = 0

	umi_summary$without_umi = umi_summary$total_reads - umi_summary$with_umi - umi_summary$too_short

	umi_summary
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


########
# Main #
########

# read in cutadapt and prepare trimming data frame
cutadapt_files = list.files(path=data_directory,pattern=".cutadapt.txt$",full.names=TRUE)

# short side remark: trying to convert from plyr -> dplyr but apparently there is no ldply equivalent in dplyr
cutadapt_summary = ldply(cutadapt_files,function(x){parse_cutadapt_summary(x)})
trimming_summary = cutadapt_summary

# if available, read in umi files - only process if there is at least one nonempty file
umi_files = list.files(path=data_directory,pattern=".umi.txt$",full.names=TRUE)

have_umi_data = length(umi_files)>0 & sum(file.size(umi_files))>0

if(have_umi_data){
	umi_summary = ldply(umi_files,function(x){parse_umi_summary(x)})
	umi_summary$too_short_umi = umi_summary$too_short
	umi_summary$total_reads = NULL
	umi_summary$too_short = NULL

	trimming_summary = merge(cutadapt_summary,umi_summary,by="library",all=T)
	trimming_summary$with_umi[is.na(trimming_summary$with_umi)] = 0
	trimming_summary$with_umi_bp[is.na(trimming_summary$with_umi_bp)] = 0
	trimming_summary$too_short_umi[is.na(trimming_summary$too_short_umi)] = 0
	trimming_summary$without_umi[is.na(trimming_summary$without_umi)] = 0
	
	# recalc reads passed and bp passed
	trimming_summary$reads_passed = ifelse(trimming_summary$with_umi>0,trimming_summary$with_umi,trimming_summary$reads_passed)
	trimming_summary$bp_passed = ifelse(trimming_summary$with_umi>0,trimming_summary$with_umi_bp,trimming_summary$bp_passed)
	
	# add reads that are too short after UMI to the general too short category
	trimming_summary$too_short = trimming_summary$too_short + trimming_summary$too_short_umi
	trimming_summary$too_short_umi = NULL
}

trimming_summary$reads_filtered = trimming_summary$total_reads - trimming_summary$reads_passed
trimming_summary$bp_filtered = trimming_summary$total_bp - trimming_summary$bp_passed

trimming_summary_m = gather(trimming_summary,"metric","value",-library)

# plot read input data 
input_data_reads_plot = ggplot(subset(trimming_summary_m,metric %in% c("total_reads")),aes(x=library,y=value)) +
geom_bar(stat="identity") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Number of reads",label=comma) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(pdf_directory,paste(bfx_id,"input_reads.pdf",sep="_"),sep="/"),input_data_reads_plot)

# plot bp input data
input_data_bp_plot = ggplot(subset(trimming_summary_m,metric %in% c("total_bp")),aes(x=library,y=value)) +
geom_bar(stat="identity") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Total sequence data in Mb",label=Mb_labels) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(pdf_directory,paste(bfx_id,"input_bp.pdf",sep="_"),sep="/"),input_data_bp_plot)


# read fate plot (only before and after)
plot_data = subset(trimming_summary_m,metric %in% c("reads_filtered","reads_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("reads_passed","reads_filtered"))
read_fate_plot = ggplot(plot_data,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Read status",type="qual",palette=2,labels=c("passed","filtered")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(pdf_directory,paste(bfx_id,"reads_removed.pdf",sep="_"),sep="/"),read_fate_plot)


# bp fate plot (only before and after)
plot_data = subset(trimming_summary_m,metric %in% c("bp_filtered","bp_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("bp_passed","bp_filtered"))
bp_fate_plot = ggplot(plot_data,aes(x=library,y=value,fill=metric)) +
geom_bar(stat="identity",position="fill") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of bp",label=percent) +
scale_fill_brewer("Bp status",type="qual",palette=2,labels=c("passed","filtered")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(pdf_directory,paste(bfx_id,"bp_removed.pdf",sep="_"),sep="/"),bp_fate_plot)

# read fate plot (complex)
plot_data = subset(trimming_summary_m,metric %in% c("without_adapter","too_short","without_umi","reads_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("reads_passed","without_adapter","too_short","without_umi"))
plot_data = subset(plot_data,!(metric=="without_umi" & value==0))
read_complex_fate_plot = ggplot(plot_data,aes(x=library,y=value,fill=metric,order=metric)) +
geom_bar(stat="identity",position="fill") +
theme_bw() +
scale_x_discrete("Library") +
scale_y_continuous("Percentage of reads",label=percent) +
scale_fill_brewer("Read status",type="qual",palette=2,labels=c("passed","without adapter","too short/primer dimer","without UMI")) +
theme(axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
ggsave(paste(pdf_directory,paste(bfx_id,"reads_removed_by_cause.pdf",sep="_"),sep="/"),read_complex_fate_plot)

# plot length distribution
length_dist_files = list.files(path=data_directory,pattern=".length_distribution.txt$",full.names=TRUE)

length_dist_table = ldply(length_dist_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.length_distribution\\.txt$','',basename(file)))
	length_table
})


length_dist_table = length_dist_table[,c("library","dir","length","reads")]
length_dist_table =  length_dist_table %>% 
			group_by(library) %>% mutate(frac_reads=reads/sum(reads)) %>%
			as.data.frame()

length_dist_plot = ggplot(length_dist_table,aes(x=length,y=frac_reads)) +
geom_bar(stat="identity",colour="black",fill="black") +
theme_bw() +
scale_x_continuous("Read length in bp",limits=c(min(length_dist_table$length),max(length_dist_table$length))) +
scale_y_continuous("Fraction of reads") +
facet_wrap(~ library,scales="free_x",ncol=3)

ggsave(paste(pdf_directory,paste(bfx_id,"clean_reads_length_distribution.pdf",sep="_"),sep="/"),length_dist_plot,width=10,height=7)


# plot umi distribution
if(have_umi_data){
	umi_distribution_files = list.files(path=data_directory,pattern=".umi.csv.gz$",full.names=TRUE)
	umi_distribution_files = umi_distribution_files[file.size(umi_distribution_files)>0]
	umi_distribution_table = ldply(umi_distribution_files,function(file){
		distribution_table = read.table(gzfile(file),header=T,sep="\t")
		distribution_table$library = gsub('\\.umi\\.csv\\.gz$','',basename(file))
		distribution_table
	})
	
	umi_distribution_table = umi_distribution_table[,c("library","UMI","count")]
	umi_distribution_table = umi_distribution_table %>% 
				group_by(library) %>% mutate(norm_count=count*1000000/sum(count)) %>%
				as.data.frame()

	# summary
	umi_distribution_summary = umi_distribution_table %>% group_by(library) %>% 
				summarise(umis=length(UMI),umi10x=sum(norm_count>0 & norm_count<=10),umi100x=sum(norm_count>10 & norm_count<=100),
				umi1000x=sum(norm_count>100 & norm_count<=1000),umi10000x=sum(norm_count>1000)) %>%
				as.data.frame()
	umi_distribution_summary_m = gather(umi_distribution_summary,"class","value",-library)
	umi_distribution_summary_m$class = factor(umi_distribution_summary_m$class,levels=c("umis","umi10x","umi100x","umi1000x","umi10000x"))

	umi_distribution_summary_plot = ggplot(umi_distribution_summary_m,aes(x=class,y=value,fill=library)) +
	geom_bar(stat="identity",position="dodge",colour="black") +
	theme_bw() +
	scale_x_discrete("UMI normalised count (rpm)",labels=c("All","0-10","11-100","101-1000",">1000")) +
	scale_y_continuous("Number of UMIs",labels=comma) +
	scale_fill_discrete("Library")
	ggsave(paste(pdf_directory,paste(bfx_id,"umi_summary.pdf",sep="_"),sep="/"),umi_distribution_summary_plot)

	ggplot(umi_distribution_table,aes(x=count,colour=library)) + 
	geom_freqpoly(binwidth=10) + scale_x_continuous(limits=c(0,300))
	
	# rank plot
	num_top = 1000
	top_per_library = umi_distribution_table %>% 
	group_by(library) %>% 
	do({
		top_counts = head(.[order(-.$norm_count),],num_top)
		top_counts$rank = 1:num_top
		top_counts
	}) %>%
	as.data.frame()
	
	umi_rank_plot = ggplot(top_per_library,aes(x=rank,y=norm_count,colour=library)) +
	geom_line() +
	theme_bw() +
	scale_x_continuous("UMI Rank") +
	scale_y_continuous("UMI normalised count (rpm)") +
	scale_colour_discrete("Library")
	ggsave(paste(pdf_directory,paste(bfx_id,"umi_ranked_top_counts.pdf",sep="_"),sep="/"),umi_rank_plot)

	umi_rank_10_plot = ggplot(top_per_library,aes(x=rank,y=norm_count,colour=library)) +
	geom_line() +
	theme_bw(10) +
	scale_x_continuous("UMI Rank",limits=c(0,10),breaks=0:10) +
	scale_y_continuous("UMI normalised count (rpm)") +
	scale_colour_discrete("Library") +
	theme(legend.position="none")

	umi_adv_rank_plot = umi_rank_plot + 
	annotation_custom(ggplot_gtable(ggplot_build(umi_rank_10_plot)),xmin=0.3*num_top,ymin=0.5*max(top_per_library$norm_count),xmax=num_top,ymax=max(top_per_library$norm_count))

	# correlation plot (only spearman rank correlation)
	umi_distribution_table_w = spread(umi_distribution_table[,c("library","UMI","norm_count")],library,norm_count,fill=0)
	row.names(umi_distribution_table_w) = umi_distribution_table_w$UMI
	umi_distribution_table_w$UMI = NULL
	plotCorrelation(umi_distribution_table_w, conditions="", cormethod="spearman", cluster=FALSE, filename=paste(pdf_directory,paste(bfx_id,"umi_spearman_correlation.pdf",sep="_"),sep="/"))
}





