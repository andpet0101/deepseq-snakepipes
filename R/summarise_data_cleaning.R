#!/usr/bin/env Rscript 

# reorders the library paths so that user paths will be the last place to look for libraries
lib_paths = .libPaths()
home_last_order = order(grepl('^/home',.libPaths()))
.libPaths(lib_paths[home_last_order])


library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(scales)
library(grid)
theme_set(theme_bw(14))


#############
# Arguments #
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

#############
# Functions #
#############
source(paste(script_dir,"functions_settings.R",sep="/"))

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

########
# Main #
########

# read in cutadapt and prepare trimming data frame
cutadapt_files = list.files(path=data_directory,pattern=".cutadapt.txt$",full.names=TRUE)

# short side remark: trying to convert from plyr -> dplyr but apparently there is no ldply equivalent in dplyr
cutadapt_summary = ldply(cutadapt_files,function(x){parse_cutadapt_summary(x)})
library_names = fixLibraryNames(cutadapt_summary$library)
cutadapt_summary$library = factor(library_names,levels=unique(library_names))
trimming_summary = cutadapt_summary

# if available, read in umi files - only process if there is at least one nonempty file
umi_files = list.files(path=data_directory,pattern=".umi.txt$",full.names=TRUE)

have_umi_data = length(umi_files)>0 & sum(file.size(umi_files))>0

if(have_umi_data){
	umi_summary = ldply(umi_files,function(x){parse_umi_summary(x)})
	library_names = fixLibraryNames(umi_summary$library)
	umi_summary$library = factor(library_names,levels=unique(library_names))
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
inputdf <- subset(trimming_summary_m,metric == "total_reads")
inputdf$Libraries <- 'Libraries'
input_data_reads_plot <- ggplot(inputdf, aes(x=Libraries, y=value)) + geom_boxplot(notch = T, fill = 'darkgoldenrod1') + 
  scale_y_continuous("Number of reads",label=comma) + ggtitle("Sequenced reads") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste(pdf_directory,paste(bfx_id,"input_reads.pdf",sep="_"),sep="/"),input_data_reads_plot)

# plot bp input data
inputdf <- subset(trimming_summary_m,metric == "total_bp")
inputdf$Libraries <- 'Libraries'
input_data_bp_plot = ggplot(inputdf ,aes(x=Libraries,y=value)) + geom_boxplot(notch = T, fill = 'darkgoldenrod1') + 
  scale_y_continuous("Total sequence data in Mb",label=Mb_labels) + ggtitle("Sequenced bps") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste(pdf_directory,paste(bfx_id,"input_bp.pdf",sep="_"),sep="/"),input_data_bp_plot)


# read fate plot (only before and after)
plot_data = subset(trimming_summary_m,metric %in% c("reads_filtered","reads_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("reads_passed","reads_filtered"))
plot_data$percent <- plot_data$value/max(plot_data$value)
read_fate_plot = ggplot(plot_data,aes(x=metric,y=percent, fill=metric)) + geom_boxplot(notch = T) + 
  scale_x_discrete("Read status", labels = c('Passed', 'Filtered')) + 
  scale_y_continuous("Percentage of reads", labels = percent) + 
  scale_fill_brewer("Read status",type="qual",palette=2) +
  ggtitle("Reads removed after data cleaning") + theme(legend.position = 'none') 
ggsave(paste(pdf_directory,paste(bfx_id,"reads_removed.pdf",sep="_"),sep="/"),read_fate_plot)


# bp fate plot (only before and after)
plot_data = subset(trimming_summary_m,metric %in% c("bp_filtered","bp_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("bp_passed","bp_filtered"))
plot_data$percent <- plot_data$value/max(plot_data$value)
bp_fate_plot = ggplot(plot_data,aes(x=metric,y=percent, fill=metric)) + geom_boxplot(notch = T) + 
  scale_x_discrete("Read status", labels = c('Passed', 'Filtered')) +
  scale_y_continuous("Percentage of reads", labels = percent) + 
  scale_fill_brewer("Read status",type="qual",palette=2) +
  ggtitle("Bps removed after data cleaning") + theme(legend.position = 'none') 
ggsave(paste(pdf_directory,paste(bfx_id,"bp_removed.pdf",sep="_"),sep="/"),bp_fate_plot)

# read fate plot (complex)
plot_data = subset(trimming_summary_m,metric %in% c("without_adapter","too_short","without_umi","reads_passed"))
plot_data$metric = factor(plot_data$metric,levels=c("reads_passed","without_adapter","too_short","without_umi"))
plot_data = subset(plot_data,!(metric=="without_umi" & value==0))
plot_data$percent <- plot_data$value/max(plot_data$value)

read_complex_fate_plot = ggplot(plot_data,aes(x=metric,y=percent, fill=metric, order = metric)) + geom_boxplot(notch = T) + 
  scale_fill_brewer("Read status",type="qual",palette=2,labels=c("passed","without adapter","too short/primer dimer","without UMI")) +
  scale_y_continuous("Percentage of reads", labels = percent) + 
  ggtitle("Reads removed by reason") + theme(legend.position = 'none')
ggsave(paste(pdf_directory,paste(bfx_id,"reads_removed_by_cause.pdf",sep="_"),sep="/"),read_complex_fate_plot)

# plot length distribution
length_dist_files = list.files(path=data_directory,pattern=".length_distribution.txt$",full.names=TRUE)

length_dist_table = ldply(length_dist_files,function(file){
	length_table = read.table(file,header=T,sep="\t")
	length_table$library = gsub('_R\\d$','',gsub('\\.length_distribution\\.txt$','',basename(file)))
	length_table
})
library_names = fixLibraryNames(length_dist_table$library)
length_dist_table$library = factor(library_names,levels=unique(library_names))

length_dist_table = length_dist_table[,c("library","dir","length","reads")]
length_dist_table =  length_dist_table %>% 
			group_by(library) %>% mutate(frac_reads=reads/sum(reads)) %>%
			as.data.frame()

num_columns = 4

length_dist_table$length <- factor(length_dist_table$length)

length_dist_plot = ggplot(length_dist_table,aes(x=length,y=frac_reads)) +
  geom_boxplot(notch = F,fill="black") + scale_y_continuous("Fraction of reads") +
  xlab("Read length in bp") + ggtitle("Length distribution of clean reads")

required_width = 11/4*num_columns
required_height = 6.5/3*ceiling(length(levels(length_dist_table$library))/num_columns)


ggsave(paste(pdf_directory,paste(bfx_id,"clean_reads_length_distribution.pdf",sep="_"),sep="/"),length_dist_plot)
# ,width=required_width,height=required_height)


# plot umi distribution
if(have_umi_data){
	umi_distribution_files = list.files(path=data_directory,pattern=".umi.csv.gz$",full.names=TRUE)
	umi_distribution_files = umi_distribution_files[file.size(umi_distribution_files)>0]
	umi_distribution_table = ldply(umi_distribution_files,function(file){
		distribution_table = read.table(gzfile(file),header=T,sep="\t")
		distribution_table$library = gsub('\\.umi\\.csv\\.gz$','',basename(file))
		distribution_table
	})
	library_names = fixLibraryNames(umi_distribution_table$library)
	umi_distribution_table$library = factor(library_names,levels=unique(library_names))

	
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
	scale_fill_discrete("Library") +
	ggtitle("UMI distribution") +
	theme(legend.position="bottom")
	ggsave(paste(pdf_directory,paste(bfx_id,"umi_summary.pdf",sep="_"),sep="/"),umi_distribution_summary_plot,width=9,height=8)

	#ggplot(umi_distribution_table,aes(x=count,colour=library)) + 
	#geom_freqpoly(binwidth=10) + scale_x_continuous(limits=c(0,300))
	
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
	scale_colour_discrete("Library") +
	ggtitle("UMI ranked by occurrence (high to low)") +
	theme(legend.position="bottom")
	ggsave(paste(pdf_directory,paste(bfx_id,"umi_ranked_top_counts.pdf",sep="_"),sep="/"),umi_rank_plot,width=9,height=8)


	# correlation plot (only spearman rank correlation)
	umi_distribution_table_w = spread(umi_distribution_table[,c("library","UMI","norm_count")],library,norm_count,fill=0)
	row.names(umi_distribution_table_w) = umi_distribution_table_w$UMI
	umi_distribution_table_w$UMI = NULL

	if(ncol(umi_distribution_table)>1){
		plotCorrelation(umi_distribution_table_w, conditions="", cormethod="spearman", cluster=FALSE,  title = "Spearman correlation of UMI counts", filename=paste(pdf_directory,paste(bfx_id,"umi_spearman_correlation.pdf",sep="_"),sep="/"))
	} else {
		ggsave(grid.text("No umi_spearman_correlation plot possible"),file=paste(pdf_directory,paste(bfx_id,"umi_spearman_correlation.pdf",sep="_"),sep="/"))
	}
}





