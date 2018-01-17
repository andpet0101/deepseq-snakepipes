#!/usr/bin/env Rscript

# reorders the library paths so that user paths will be the last place to look for libraries
lib_paths = .libPaths()
home_last_order = order(grepl('^/home',.libPaths()))
.libPaths(lib_paths[home_last_order])


library(ggplot2)
library(plyr)
library(scales)
theme_set(theme_bw(12))
library(reshape2)
library(dplyr)
library(tidyr)
library(grid)

##############
# Arguments  #
##############

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

#bfx_id = "bfx735"
#data_directory = "fastqc/report/data"
#pdf_directory = "fastqc/report/pdf"

#############
# Functions #
#############
source(paste(script_dir,"functions_settings.R",sep="/"))

# fixes library names
fixLibraryNames = function(library_names){
  library_names_wo_bfx = as.character(gsub("^bfx\\d+_","",library_names))
  library_names_wo_bfx_lib = as.character(gsub("^L\\d+_","",library_names_wo_bfx))
  
  if(length(unique(library_names_wo_bfx))==length(unique(library_names_wo_bfx_lib))){
    library_names_wo_bfx_lib
  }else{
    library_names_wo_bfx
  }
}


parse_fastqc_file = function(fastqc_file){
	# if zip then extract
	is_zipped = grepl(".zip$",fastqc_file)
	if(is_zipped){
		fastqc_unzipped_directory = gsub('\\.zip$','',basename(fastqc_file))
		unzip(fastqc_file)
                fastqc_file = paste(fastqc_unzipped_directory ,"fastqc_data.txt",sep="/")
	}

	# read in fastqc file
	fastqc_results = list()
	fastqc_lines = readLines(fastqc_file)
	fastqc_lines = fastqc_lines[2:length(fastqc_lines)]
	fastqc_lines = fastqc_lines[!grepl('^>>END_MODULE',fastqc_lines)]
  	
	
	if(is_zipped){
		fastqc_basename = gsub('_fastqc$','',dirname(fastqc_file))
	} else{
		fastqc_basename = gsub('_fastqc.txt$','',basename(fastqc_file))
	}
	
	fastqc_results[['Library']] = gsub('_(R\\d)$','',fastqc_basename)
	read_match = gsub('\\S+_(R\\d)$','\\1',gsub('_fastqc.txt$','',basename(fastqc_file)))
	fastqc_results[['Read']] = ifelse(grepl('^R\\d$',read_match),read_match,'R1')
  
	# get summary lines
	summary_lines = gsub('^>>','',fastqc_lines[grepl('^>>',fastqc_lines)])
	summary_table = read.table(text=summary_lines,sep="\t",col.names=c("Metric","Status"))
	summary_table$Variable = tolower(gsub('\\s+','_',summary_table$Metric))
	fastqc_results[['summary']] = summary_table
  
	# convert the remaining lines into list of vectors
	fastqc_lines = fastqc_lines[grepl('^\\S',fastqc_lines)]
  	
	tables_list = list()
	for(line in fastqc_lines){
		# '>>' starts a new table followed by a '#' header
		if(grepl('^>>',line,perl=T)){
			tables_list = append(tables_list,c(line))
		}
		# this is the only case where a '#' starts a new table without preceding '>>'
		else if(grepl('^#Duplication Level',line,perl=T)){
			tables_list = append(tables_list,c(gsub('^#','',line)))		
		}
		else if(grepl('^#',line,perl=T)){
			tables_list[[length(tables_list)]]  = append(tables_list[[length(tables_list)]] ,c(gsub('^#','',line)))
		}
		else{
			tables_list[[length(tables_list)]] = append(tables_list[[length(tables_list)]],line)
		}
	}
	
	# sequence_duplication_levels contains two tables
	table_names = as.vector(summary_table$Variable)
	table_names_length = length(table_names)
	table_names_insert_index = grep("sequence_duplication_levels",table_names)
	table_names = c(table_names[1:table_names_insert_index-1],"total_deduplicated_levels",table_names[table_names_insert_index:table_names_length])
  
	# vectors into tables and store in fastqc_results, note: we use the summary column 'Variable' to name the data frames
	for(index in 1:length(table_names)){
		variable = table_names[index]
		table_text = tables_list[[index]]
		table_text = table_text[!grepl('^>>',table_text,perl=T)]
		if(length(table_text)>0){
			if(variable=="total_deduplicated_levels"){
				fastqc_results[[variable]] = read.table(text=table_text,sep="\t",header=F,col.names=c("Measure","Value"),stringsAsFactors=F)
			} else if(variable=="sequence_duplication_levels") {
				fastqc_results[[variable]] = read.table(text=table_text,sep="\t",header=T,stringsAsFactors=F)
			} else {
				fastqc_results[[variable]] = read.table(text=table_text,sep="\t",header=T,stringsAsFactors=F)
			}
		}
	}
	
	# remove extracted files
	if(is_zipped){
		unlink(fastqc_unzipped_directory,recursive=T)	
	}
	fastqc_results
}

# function defines a (optimal)number of columns, number of rows, width and height for a facetted plot
calc_facet_plot_size = function(num_plots,max_plots=20){
	if(num_plots==1){
		num_columns = 1
		num_rows = 1
		width = num_columns*4
		height = num_rows*4
	} else if(num_plots>1 & num_plots<=4){
		num_columns = 2
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*3
		height = num_rows*3
	} else if(num_plots>4 & num_plots<10){
		num_columns = 3
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*3
		height = num_rows*3	
	} else{
		num_columns = 4
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*2
		height = num_rows*2
	}
	
	plot_size_data = list()
	plot_size_data[['num_columns']] = num_columns
	plot_size_data[['num_rows']] = num_rows
	plot_size_data[['width']] = width
	plot_size_data[['height']] = height
	plot_size_data
}

########
# Main #
########

# read in fastqc files and prepare a large fastq_result object that contains different tables
fastqc_files = list.files(path=data_directory,pattern="_fastqc\\.(txt|zip)$",full.names=TRUE)
fastqc_result = list()

for(f in fastqc_files){
	fastqc = parse_fastqc_file(f)
	tables = names(fastqc)
	tables = tables[!grepl('Library|Read',tables)]
	for(t in tables){
		fastqc[[t]][['Library']] = fastqc$Library
		fastqc[[t]][['Read']] = fastqc$Read
		if(t %in% names(fastqc_result)){
			fastqc_result[[t]] = rbind(fastqc_result[[t]],fastqc[[t]])
		} else {
			fastqc_result[[t]] = fastqc[[t]]
	    	}
	}
}

# names(fastqc_result)
# [1] "summary"                      "basic_statistics"            
# [3] "per_base_sequence_quality"    "per_tile_sequence_quality"   
# [5] "per_sequence_quality_scores"  "per_base_sequence_content"   
# [7] "per_sequence_gc_content"      "per_base_n_content"          
# [9] "sequence_length_distribution" "sequence_duplication_levels" 
#[11] "overrepresented_sequences"    "adapter_content"             
#[13] "kmer_content" 

# set Library, Read and Library_read factor levels
all_library_names = c()
all_read_dirs = c("R1","R2")
for(t in names(fastqc_result)){
	library_names = fixLibraryNames(fastqc_result[[t]][['Library']])
	fastqc_result[[t]][['Library']] = library_names
	all_library_names = unique(c(library_names,all_library_names))
}

for(t in names(fastqc_result)){
	fastqc_result[[t]][['Library']] = factor(fastqc_result[[t]][['Library']],levels=all_library_names)
	fastqc_result[[t]][['Read']] = factor(fastqc_result[[t]][['Read']],levels=all_read_dirs)
	fastqc_result[[t]][['Library_read']] = factor(paste(fastqc_result[[t]][['Library']],fastqc_result[[t]][['Read']]),
		levels=apply(expand.grid(all_library_names,all_read_dirs),1,paste,collapse=" "))
	column_names = colnames(fastqc_result[[t]])
	column_names = c("Library","Read","Library_read",setdiff(column_names,c("Library","Read","Library_read")))
	fastqc_result[[t]] = fastqc_result[[t]][,column_names]
}

# do some minor data fixing
fastqc_result$summary$Status = factor(fastqc_result$summary$Status,levels=c("pass","warn","fail"))
fastqc_result$summary$Metric = factor(fastqc_result$summary$Metric,levels=unique(fastqc_result$summary$Metric))
fastqc_result$per_base_sequence_quality$Base = factor(fastqc_result$per_base_sequence_quality$Base,levels=unique(fastqc_result$per_base_sequence_quality$Base),ordered=T)
fastqc_result$per_base_sequence_content$Base = factor(fastqc_result$per_base_sequence_content$Base,levels=unique(fastqc_result$per_base_sequence_content$Base),ordered=T)
fastqc_result$per_base_n_content$Base = factor(fastqc_result$per_base_n_content$Base,levels=unique(fastqc_result$per_base_n_content$Base),ordered=T)
fastqc_result$adapter_content$Position = factor(fastqc_result$adapter_content$Position,levels=levels(fastqc_result$per_base_n_content$Base))
fastqc_result$kmer_content$Max.Obs.Exp.Position = factor(fastqc_result$kmer_content$Max.Obs.Exp.Position,levels=unique(fastqc_result$kmer_content$Max.Obs.Exp.Position))
fastqc_result$sequence_length_distribution$Length = factor(fastqc_result$sequence_length_distribution$Length,levels=unique(fastqc_result$sequence_length_distribution$Length))

fastqc_result$per_sequence_quality_scores = fastqc_result$per_sequence_quality_scores %>% group_by(Library_read) %>% mutate(Fraction=Count/sum(Count+1)) %>% as.data.frame()
fastqc_result$per_base_n_content = fastqc_result$per_base_n_content %>% group_by(Library_read) %>% mutate(Fraction=N.Count/sum(N.Count+1)) %>% as.data.frame()
fastqc_result$per_sequence_gc_content = fastqc_result$per_sequence_gc_content %>% group_by(Library_read) %>% mutate(Fraction=Count/sum(Count+1)) %>% as.data.frame()
fastqc_result$sequence_length_distribution = fastqc_result$sequence_length_distribution %>% group_by(Library,Read,Library_read) %>% mutate(Fraction=Count/sum(Count+1)) %>% as.data.frame()



# thats it - now do plots

# too many datasets? - do different plots
too_many_datasets = length(fastqc_files) > 20


# summary (note: fastqc_result$summary = fastqc_result[['summary']])
num_plots = length(unique(fastqc_result$summary$Read))
plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))

if(!too_many_datasets){
	fastqc_overview_plot = ggplot(fastqc_result$summary,aes(x=Metric,y=Library)) + 
	geom_tile(aes(fill=Status),colour="black") + 
	theme_bw(12) +
	scale_x_discrete("Metric") + 
	scale_y_discrete("Library") +
	scale_fill_manual("Status",values=c("#91cf60","#ffffbf","#fc8d59")) +
	facet_wrap(~ Read,ncol=plot_size$num_columns) +
	ggtitle("FastQC summary") +
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")
} else {
	fastqc_overview_plot = ggplot(fastqc_result$summary,aes(x=Metric,fill=Status)) +
	geom_bar(position="fill",colour="black") +
	theme_bw(12) +
	scale_x_discrete("Metric") + 
	scale_y_continuous("Percentage of libraries",labels=percent) +
	scale_fill_manual("Status",values=c("#91cf60","#ffffbf","#fc8d59")) +
	facet_wrap(~ Read,ncol=plot_size$num_columns) +
	ggtitle("FastQC summary") +
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")
}
ggsave(paste(pdf_directory,paste(bfx_id,"fastqc_summary.pdf",sep="_"),sep="/"),plot=fastqc_overview_plot,width=plot_size$width,height=plot_size$height)
write.table(fastqc_result$summary[,-3],paste(data_directory,paste(bfx_id,"fastqc_summary.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)

# per base sequence quality
if("per_base_sequence_quality" %in% names(fastqc_result)){
	x_axis_breaks = levels(fastqc_result$per_base_sequence_quality$Base)
	x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})

	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(fastqc_result$per_base_sequence_quality$Library_read)))
		per_base_sequence_quality_plot = ggplot(fastqc_result$per_base_sequence_quality,aes(x=Base)) + 
		geom_boxplot(aes(ymin=X10th.Percentile,lower=Lower.Quartile,middle=Median,upper=Upper.Quartile,ymax=X90th.Percentile),stat="identity",fill="yellow") +
		geom_line(aes(y=Mean,group=Library),colour="red") +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Quality",limits=c(0,max(fastqc_result$per_base_sequence_quality$X90th.Percentile))) +
		facet_wrap(~Library_read,ncol=plot_size$num_columns,scales="free_x") +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Per base sequence quality distribution")
	} else {
		per_base_sequence_quality_m = gather(fastqc_result$per_base_sequence_quality,Metric,Value,-c(Library,Read,Library_read,Base))
		per_base_sequence_quality_m$Metric_read = factor(paste(per_base_sequence_quality_m$Metric,per_base_sequence_quality_m$Read),
		                                                 levels=apply(expand.grid(c("X10th.Percentile","Lower.Quartile","Mean","Median","Upper.Quartile","X90th.Percentile"),unique(per_base_sequence_quality_m$Read)),1,paste,collapse=" "))		
		
		plot_size = calc_facet_plot_size(length(unique(per_base_sequence_quality_m$Metric)))	
		per_base_sequence_quality_plot = ggplot(per_base_sequence_quality_m,aes(x=Base,y=Value)) + 
		geom_boxplot(fill="yellow",outlier.size=1) +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Quality",limits=c(0,max(fastqc_result$per_base_sequence_quality$X90th.Percentile))) +
		facet_wrap(~Metric_read,ncol=plot_size$num_columns) +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Per base sequence quality distribution")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"per_base_sequence_quality.pdf",sep="_"),sep="/"),plot=per_base_sequence_quality_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$per_base_sequence_quality[,-3],paste(data_directory,paste(bfx_id,"per_base_sequence_quality.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"per_base_sequence_quality.csv",sep="_"),sep="/"))
	ggsave(grid.text("No per base sequence quality plot possible"),file=paste(pdf_directory,paste(bfx_id,"per_base_sequence_quality.pdf",sep="_"),sep="/"))
}

# per sequence quality scores
if("per_base_sequence_quality" %in% names(fastqc_result)){
	num_plots = length(unique(fastqc_result$per_sequence_quality_scores$Read))
	plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))

	if(!too_many_datasets){
		per_sequence_quality_scores_plot = ggplot(fastqc_result$per_sequence_quality_scores,aes(x=Quality,y=Fraction,colour=Library)) +
		geom_line() +
		theme_bw(12) +
		scale_x_continuous("Quality",limits=c(0,max(fastqc_result$per_sequence_quality_scores$Quality)+1)) +
		scale_y_continuous("Fraction of reads") +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		ggtitle("Per sequence quality scores") +
		facet_wrap(~Read,scales="free_x",ncol=plot_size$num_columns) +
		theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="bottom")
	} else {
		mean_sequence_quality_scores = fastqc_result$per_sequence_quality_scores %>% group_by(Quality) %>% summarise(Fraction=mean(Fraction)) %>% as.data.frame()
		per_sequence_quality_scores_plot = ggplot(fastqc_result$per_sequence_quality_scores,aes(x=Quality,y=Fraction,group=Quality)) +
		geom_boxplot(fill="yellow") +
		geom_point(data=mean_sequence_quality_scores,aes(x=Quality,y=Fraction,group=1),colour="red",shape=3) +	
		theme_bw(12) +
		scale_x_continuous("Quality",limits=c(0,max(fastqc_result$per_sequence_quality_scores$Quality)+1)) +
		scale_y_continuous("Fraction of reads") +
		facet_wrap(~Read,ncol=plot_size$num_columns) +
		ggtitle("Per sequence quality scores")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"per_sequence_quality_scores.pdf",sep="_"),sep="/"),plot=per_sequence_quality_scores_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$per_sequence_quality_scores[,-3],paste(data_directory,paste(bfx_id,"per_sequence_quality_scores.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"per_sequence_quality_scores.csv",sep="_"),sep="/"))
	ggsave(grid.text("No per sequence quality scores plot possible"),file=paste(pdf_directory,paste(bfx_id,"per_sequence_quality_scores.pdf",sep="_"),sep="/"))
}

# per sequence base content
if("per_base_sequence_content" %in% names(fastqc_result)){
	per_base_sequence_content_m = gather(fastqc_result$per_base_sequence_content,Nt,Percent,-c(Base,Library,Read,Library_read))
	x_axis_breaks = levels(per_base_sequence_content_m$Base)
	x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(per_base_sequence_content_m$Library_read)))
		per_base_sequence_content_plot = ggplot(per_base_sequence_content_m,aes(x=Base,y=Percent,colour=Nt,group=Nt)) +
		geom_line() +
		geom_hline(yintercept=25,colour="gray",linetype="dashed") +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Percentage of reads") +
		scale_colour_manual("Nt",values=c("#1f78b4","#e31a1c","#6a3d9a","#33a02c")) +
		ggtitle("Per base sequence content") +
		facet_wrap(~Library_read,scales="free_x",ncol=plot_size$num_columns) +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,vjust=1,hjust=1))
	} else {
		per_base_sequence_content_m = per_base_sequence_content_m %>% group_by(Read,Base,Nt) %>% summarise(Percent = mean(Percent)) %>% as.data.frame()
		per_base_sequence_content_plot = ggplot(per_base_sequence_content_m,aes(x=Base,y=Percent,colour=Nt,fill=Nt,group=paste0(Read,Base,Nt))) +
		geom_bar(position="stack",stat="identity") +
		theme_bw(10) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels,expand=c(0,0)) +
		scale_y_continuous("Percentage of reads",expand=c(0,0)) +
		scale_fill_manual("Nt",values=c("#1f78b4","#e31a1c","#6a3d9a","#33a02c")) +
		scale_colour_manual("Nt",values=c("#1f78b4","#e31a1c","#6a3d9a","#33a02c")) +
		facet_grid(~Read) +
		theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1)) + 
		ggtitle("Per base sequence content")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"per_base_sequence_content.pdf",sep="_"),sep="/"),plot=per_base_sequence_content_plot,width=8,height=11)
	write.table(fastqc_result$per_base_sequence_content[,-3],paste(data_directory,paste(bfx_id,"per_base_sequence_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"per_base_sequence_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No per base sequence content plot possible"),file=paste(pdf_directory,paste(bfx_id,"per_base_sequence_content.pdf",sep="_"),sep="/"))
}

# per sequence GC content
if("per_sequence_gc_content" %in% names(fastqc_result)){
	num_plots = length(unique(fastqc_result$per_sequence_gc_content$Read))
	plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))	
	if(!too_many_datasets){
		per_sequence_gc_content_plot = ggplot(fastqc_result$per_sequence_gc_content,aes(x=GC.Content,y=Fraction,colour=Library)) +
		geom_line() +
		theme_bw(12) +
		scale_x_continuous("GC content",limits=c(0,100)) +
		scale_y_continuous("Fraction of reads") +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		theme(legend.position="bottom") +
		ggtitle("Per sequence GC content")
	} else {
		mean_gc_content_fraction = fastqc_result$per_sequence_gc_content %>% group_by(GC.Content) %>% summarise(Fraction=mean(Fraction)) %>% as.data.frame()
		per_sequence_gc_content_plot = ggplot(fastqc_result$per_sequence_gc_content,aes(x=GC.Content,y=Fraction,group=GC.Content)) +
		geom_boxplot(fill="yellow") +
		geom_point(data=mean_gc_content_fraction,aes(x=GC.Content,y=Fraction,group=1),colour="red",shape=3) +	
		theme_bw(12) +
		scale_x_continuous("GC content",limits=c(0,100)) +
		scale_y_continuous("Fraction of reads") +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		ggtitle("Per sequence GC content")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"per_sequence_gc_content.pdf",sep="_"),sep="/"),plot=per_sequence_gc_content_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$per_sequence_gc_content[,-3],paste(data_directory,paste(bfx_id,"per_sequence_gc_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"per_sequence_gc_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No per sequence GC content plot possible"),file=paste(pdf_directory,paste(bfx_id,"per_sequence_gc_content.pdf",sep="_"),sep="/"))
}


# per base N content
if("per_base_n_content" %in% names(fastqc_result)){
	x_axis_breaks = levels(fastqc_result$per_base_n_content$Base)
	x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})
	num_plots = length(unique(fastqc_result$per_base_n_content$Read))
	plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))	
	if(!too_many_datasets){
		per_base_n_content_plot = ggplot(fastqc_result$per_base_n_content,aes(x=Base,y=Fraction,colour=Library)) +
		geom_line() +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Fraction of reads",limits=c(0,ifelse(max(fastqc_result$per_base_n_content$Fraction)==0,1,max(fastqc_result$per_base_n_content$Fraction)))) +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		facet_wrap(~Read,scales="free_x",ncol=plot_size$num_columns) +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Per sequence N content")
	} else {
		mean_per_base_content_fraction = fastqc_result$per_base_n_content %>% group_by(Base) %>% summarise(Fraction=mean(Fraction)) %>% as.data.frame()
		per_base_n_content_plot = ggplot(fastqc_result$per_base_n_content,aes(x=Base,y=Fraction,group=Base)) +
		geom_boxplot(fill="yellow") +
		geom_point(data=mean_per_base_content_fraction,aes(x=Base,y=Fraction,group=1),colour="red",shape=3) +	
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Fraction of reads",limits=c(0,ifelse(max(fastqc_result$per_base_n_content$Fraction)==0,1,max(fastqc_result$per_base_n_content$Fraction)))) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Per base N content")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"per_base_n_content.pdf",sep="_"),sep="/"),plot=per_base_n_content_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$per_base_n_content[,-3],paste(data_directory,paste(bfx_id,"per_base_n_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"per_base_n_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No per base N content plot possible"),file=paste(pdf_directory,paste(bfx_id,"per_base_n_content.pdf",sep="_"),sep="/"))
}


# sequence length distribution
if("sequence_length_distribution" %in% names(fastqc_result)){
	x_axis_breaks = levels(fastqc_result$sequence_length_distribution$Length)
	x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})
	num_plots = length(unique(fastqc_result$sequence_length_distribution$Read))
	plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))
	num_columns = plot_size$num_columns
	required_width = plot_size$width
	required_height = plot_size$height
	if(!too_many_datasets){
		sequence_length_distribution_plot = ggplot(fastqc_result$sequence_length_distribution,aes(x=Length,y=Fraction,color=Library,group=1)) +
		geom_line() +
		geom_point() +
		theme_bw(12) +
		scale_x_discrete("Length",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Fraction of reads",limits=c(0,1)) +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		theme(legend.position="bottom") +
		ggtitle("Sequence length distribution")
	} else {
		mean_sequence_length_distribution_fraction = fastqc_result$sequence_length_distribution %>% group_by(Length) %>% summarise(Fraction=mean(Fraction)) %>% as.data.frame()
		sequence_length_distribution_plot = ggplot(fastqc_result$sequence_length_distribution,aes(x=Length,y=Fraction,group=Length)) +
		geom_boxplot(fill="yellow") +
		geom_point(data=mean_sequence_length_distribution_fraction,aes(x=Length,y=Fraction,group=1),colour="red",shape=3) +	
		theme_bw(12) +
		scale_x_discrete("Length",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Fraction of reads",limits=c(0,1)) +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		ggtitle("Sequence length distribution")
	  }
	ggsave(paste(pdf_directory,paste(bfx_id,"sequence_length_distribution.pdf",sep="_"),sep="/"),plot=sequence_length_distribution_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$sequence_length_distribution[,-3],paste(data_directory,paste(bfx_id,"sequence_length_distribution.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"sequence_length_distribution.csv",sep="_"),sep="/"))
	ggsave(grid.text("No sequence length distribution plot possible"),file=paste(pdf_directory,paste(bfx_id,"sequence_length_distribution.pdf",sep="_"),sep="/"))
}

# total deduplicated levels
if("total_deduplicated_levels" %in% names(fastqc_result)){
	num_plots = length(unique(fastqc_result$total_deduplicated_levels$Read))
	plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))
	num_columns = plot_size$num_columns
	required_width = plot_size$width
	required_height = plot_size$height
	if(!too_many_datasets){
		total_deduplicated_levels_plot = ggplot(fastqc_result$total_deduplicated_levels,aes(x=Library,y=Value,fill=Library)) +
		geom_bar(stat="identity",colour="black") +
		theme_bw(12) +
		scale_x_discrete("Library") +
		scale_y_continuous("Percentage of reads",limits=c(0,100)) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		scale_fill_manual("Library",values=color_brewer_qual_palette) +
		theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
		ggtitle("Total reads after deduplication")
	} else {
		total_deduplicated_levels_plot = ggplot(fastqc_result$total_deduplicated_levels,aes(x=as.factor(1),y=Value,group=1)) +
		geom_boxplot(fill="yellow") +
		theme_bw(12) +
		scale_x_discrete("",labels="Libraries") +
		scale_y_continuous("Percentage of reads",limits=c(0,100)) +
		scale_colour_manual("Library",values=color_brewer_qual_palette) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		ggtitle("Total reads after deduplication")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"total_deduplicated_levels.pdf",sep="_"),sep="/"),plot=total_deduplicated_levels_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$total_deduplicated_levels[,-3],paste(data_directory,paste(bfx_id,"total_deduplicated_levels.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"total_deduplicated_levels.csv",sep="_"),sep="/"))
	ggsave(grid.text("No total deduplicated levels plot possible"),file=paste(pdf_directory,paste(bfx_id,"total_deduplicated_levels.pdf",sep="_"),sep="/"))
}

# sequence duplication levels
if("sequence_duplication_levels" %in% names(fastqc_result)){
	sequence_duplication_levels_m = gather(fastqc_result$sequence_duplication_levels,Metric,Value,-c(Duplication.Level,Library,Read,Library_read))
	sequence_duplication_levels_m$Duplication.Level = factor(sequence_duplication_levels_m$Duplication.Level,levels=unique(sequence_duplication_levels_m$Duplication.Level))
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(sequence_duplication_levels_m$Library_read)))
		sequence_duplication_levels_plot = ggplot(sequence_duplication_levels_m,aes(x=Duplication.Level,y=Value,colour=Metric,group=Metric)) +
		geom_line() +
		theme_bw(10) +
		scale_x_discrete("Duplication level") +
		scale_y_continuous("Percentage of reads",limits=c(0,100)) +
		facet_wrap(~Library_read,ncol=plot_size$num_columns,scales="free_x") +
		scale_colour_manual("Duplication status",values=c("blue","red"),labels=c("After deduplication","Total reads")) +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
		ggtitle("Sequence duplication levels")
	} else {
		num_plots = length(unique(fastqc_result$sequence_duplication_levels$Read))
		plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))
		mean_sequence_duplication_levels_m = sequence_duplication_levels_m %>% group_by(Duplication.Level,Metric) %>% summarise(Value=mean(Value))
		sequence_duplication_levels_plot = ggplot(sequence_duplication_levels_m,aes(x=Duplication.Level,y=Value,fill=Metric)) +
		geom_boxplot() +
		geom_line(data=mean_sequence_duplication_levels_m,aes(x=Duplication.Level,y=Value,colour=Metric,group=Metric)) +
		theme_bw(12) +
		scale_x_discrete("Duplication level") +
		scale_y_continuous("Percentage of reads",limits=c(0,100)) +
		scale_colour_manual("Duplication status",values=c("blue","red"),labels=c("After deduplication","Total reads")) +
		scale_fill_manual("Duplication status",values=c("blue","red"),labels=c("After deduplication","Total reads")) +
		facet_wrap(~Read,ncol=plot_size$num_columns,scales="free_x") +
		theme(legend.position="bottom") +
		ggtitle("Sequence duplication levels")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"sequence_duplication_levels.pdf",sep="_"),sep="/"),plot=sequence_duplication_levels_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$sequence_duplication_levels[,-3],paste(data_directory,paste(bfx_id,"sequence_duplication_levels.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"sequence_duplication_levels.csv",sep="_"),sep="/"))
	ggsave(grid.text("No sequence duplication levels plot possible"),file=paste(pdf_directory,paste(bfx_id,"sequence_duplication_levels.pdf",sep="_"),sep="/"))
}


# overrepresented sequences
if("overrepresented_sequences" %in% names(fastqc_result)){
	top10_overrepresented_sequences_by_library = fastqc_result$overrepresented_sequences %>% group_by(Library) %>% mutate(Rank = rank(Percentage)) %>% arrange(desc(Rank)) %>% as.data.frame() %>% filter(Rank<=10) %>% as.data.frame()
	top10_overrepresented_sequences_by_library$Sequence = factor(top10_overrepresented_sequences_by_library$Sequence)
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(top10_overrepresented_sequences_by_library$Library_read)))
		top10_overrepresented_sequences_plot = ggplot(top10_overrepresented_sequences_by_library,aes(x=-Rank,y=Percentage)) +
		geom_bar(stat="identity",fill="orange",colour="black") +
		geom_text(aes(x=-Rank,y=0,label=Sequence,colour=Sequence),angle=90,hjust=0,family="mono",fontface="bold",size=2) +
		theme_bw(12) +
		scale_x_discrete("Max top10 overrepresented sequences in dataset") +
		scale_y_continuous("Percentage of reads") +
		theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="none") +
		facet_wrap(~Library_read,scales="free_x",ncol=plot_size$num_columns) +
		ggtitle("Overrepresented sequences")
	} else {
		top10_overrepresented_sequences = as.data.frame(table(subset(top10_overrepresented_sequences_by_library,Read=="R1")$Sequence))
		top10_overrepresented_sequences$Read = "R1"
		colnames(top10_overrepresented_sequences) = c("Sequence","Count","Read")
		
		if("R2" %in% top10_overrepresented_sequences_by_library$Read){
			top10_overrepresented_sequences_R2 = as.data.frame(table(subset(top10_overrepresented_sequences_by_library,Read=="R2")$Sequence))
			top10_overrepresented_sequences_R2$Read = "R2"
			colnames(top10_overrepresented_sequences_R2) = c("Sequence","Count","Read")
			top10_overrepresented_sequences = rbind(top10_overrepresented_sequences,top10_overrepresented_sequences_R2)
		}
		
		top10_overrepresented_sequences = top10_overrepresented_sequences[order(-top10_overrepresented_sequences$Count),]
		top10_overrepresented_sequences$Position = factor(1:nrow(top10_overrepresented_sequences),levels=1:(nrow(top10_overrepresented_sequences)))

		num_plots = length(unique(top10_overrepresented_sequences$Read))
		plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))

		top10_overrepresented_sequences_plot = ggplot(head(top10_overrepresented_sequences,20),aes(x=Position,y=Count)) + 
		geom_bar(stat="identity",fill="orange",colour="black") + 
		geom_text(aes(x=Position,y=0,label=Sequence),angle=90,hjust=0,family="mono",fontface="bold") +
		theme_bw(12) +
		scale_x_discrete("Sequences ordered by decreasing appearence in datasets") +
		scale_y_continuous("Number of datasets with sequence in top10") +
		theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
		facet_wrap(~Read,scales="free_x",ncol=plot_size$num_columns) +
		ggtitle("Overrepresented sequences")
	}

	ggsave(paste(pdf_directory,paste(bfx_id,"overrepresented_sequences.pdf",sep="_"),sep="/"),plot=top10_overrepresented_sequences_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$overrepresented_sequences[,-3],paste(data_directory,paste(bfx_id,"overrepresented_sequences.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"overrepresented_sequences.csv",sep="_"),sep="/"))
	ggsave(grid.text("No overrepresented sequences plot possible"),file=paste(pdf_directory,paste(bfx_id,"overrepresented_sequences.pdf",sep="_"),sep="/"))
}


# adapter content
if("adapter_content" %in% names(fastqc_result)){
	adapter_content_m = tidyr::gather(fastqc_result$adapter_content,Adapter,Content,-c(Position,Library,Read,Library_read))
	x_axis_breaks = levels(fastqc_result$adapter_content$Position)
	x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(adapter_content_m$Library_read)))
		adapter_content_plot = ggplot(adapter_content_m,aes(x=Position,y=Content,colour=Adapter)) +
		geom_line() +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Percentage of reads") +
		scale_colour_manual("Adapter",values=color_brewer_qual_palette) +
		facet_wrap(~Library_read,ncol=plot_size$num_columns) +
		theme(legend.position="bottom",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Adapter content")
	} else {
		num_plots = unique(length(adapter_content_m$Adapter))*unique(length(adapter_content_m$Read))
		plot_size = list(num_columns = num_plots,num_rows=1,height=10,width=10)
		adapter_content_plot = ggplot(adapter_content_m,aes(x=Position,y=Content,fill=Adapter,group=Position)) +
		geom_boxplot() +
		theme_bw(12) +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Percentage of reads") +
		scale_fill_manual("Adapter",values=color_brewer_qual_palette) +
		facet_grid(Adapter~Read) +
		theme(legend.position="none",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("Adapter content")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"adapter_content.pdf",sep="_"),sep="/"),plot=adapter_content_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$adapter_content,paste(data_directory,paste(bfx_id,"adapter_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
	
	# note: fastqc provides the cummulative adapter content seen at each position - thus, to get a total, take the max, also calculate a total adapter content
	adapter_content_total = adapter_content_m %>% group_by(Library,Read,Library_read,Adapter) %>% summarise(Content=max(Content)) %>% 
					group_by(Library,Read,Library_read) %>% 
					do(rbind(.,data.frame(Library=first(.$Library),Read=first(.$Read),Library_read=first(.$Library_read),Adapter=c("Total"),Content=sum(.$Content)))) %>% 
					as.data.frame()
	adapter_content_total$Adapter = relevel(factor(adapter_content_total$Adapter),"Total")
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(adapter_content_total$Library_read)))
		adapter_content_plot = ggplot(adapter_content_total,aes(x=Adapter,y=Content,fill=Adapter)) +
		geom_bar(stat="identity") +
		theme_bw(12) +
		scale_x_discrete("Adapter") +
		scale_y_continuous("Percentage of reads") +
		scale_fill_manual("Adapter",values=color_brewer_qual_palette) +
		facet_wrap(~Library_read,ncol=plot_size$num_columns) +
		theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") +
		ggtitle("Total adapter content")

	} else {
		plot_size = calc_facet_plot_size(length(unique(adapter_content_total$Read)))
		adapter_content_plot = ggplot(adapter_content_total,aes(x=Adapter,y=Content,fill=Adapter,group=Adapter)) +
		geom_boxplot() +
		theme_bw(12) +
		scale_x_discrete("Adapter") +
		scale_y_continuous("Percentage of reads") +
		scale_fill_manual("Adapter",values=color_brewer_qual_palette) +
		facet_wrap(~Read,ncol=plot_size$num_columns) +
		theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") +
		ggtitle("Adapter content")
	}
	ggsave(paste(pdf_directory,paste(bfx_id,"total_adapter_content.pdf",sep="_"),sep="/"),plot=adapter_content_plot,width=plot_size$width,height=plot_size$height)
	write.table(spread(adapter_content_total[,-3],Adapter,Content,fill=0),paste(data_directory,paste(bfx_id,"total_adapter_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"adapter_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No total adapter content plot possible"),file=paste(pdf_directory,paste(bfx_id,"total_adapter_content.pdf",sep="_"),sep="/"))
	file.create(paste(data_directory,paste(bfx_id,"total_adapter_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No adapter content plot possible"),file=paste(pdf_directory,paste(bfx_id,"adapter_content.pdf",sep="_"),sep="/"))
}

# k-mer content  - select only one kmer per position
if("kmer_content" %in% names(fastqc_result) && length(fastqc_result$kmer_content$Max.Obs.Exp.Position)>0){
	kmer_content_plot_data = fastqc_result$kmer_content %>% group_by(Library,Read,Library_read,Max.Obs.Exp.Position) %>% mutate(Rank = rank(-Obs.Exp.Max)) %>% filter(Rank==1) %>%
	group_by(Library,Read,Library_read) %>% mutate(Plot_row=rank(Max.Obs.Exp.Position)) %>% as.data.frame()
	kmer_content_plot_data$Sequence = factor(kmer_content_plot_data$Sequence)
	x_axis_breaks = levels(fastqc_result$kmer_content$Max.Obs.Exp.Position)
	x_axis_labels = x_axis_breaks	
	#x_axis_labels = sapply(1:length(x_axis_breaks),function(x){ifelse(x%%5==1,x_axis_breaks[x],"")})
	if(!too_many_datasets){
		plot_size = calc_facet_plot_size(length(unique(kmer_content_plot_data$Library_read)))
		kmer_content_plot = ggplot(kmer_content_plot_data,aes(y=Plot_row,x=Max.Obs.Exp.Position,colour=Sequence,label=Sequence)) + 
		geom_text(hjust=0,family="mono") + 
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		facet_wrap(~Library_read,scales="free",ncol=plot_size$num_columns) +
		theme(legend.position="none",axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("K-mer content")

	} else {
		kmer_content_plot_data_summarised = kmer_content_plot_data %>% group_by(Read,Max.Obs.Exp.Position) %>% do(as.data.frame(table(Sequence=.$Sequence))) %>% as.data.frame()
		kmer_content_plot_data_summarised$Sequence = factor(kmer_content_plot_data_summarised$Sequence)
		num_plots = length(unique(kmer_content_plot_data_summarised$Read))
		plot_size = list(num_columns = num_plots,num_rows=1,height=8,width=ifelse(num_plots==1,8,12))

	  	kmer_content_plot = ggplot(subset(kmer_content_plot_data_summarised,Freq>0),aes(x=Max.Obs.Exp.Position,y=Freq,fill=Sequence,label=Sequence)) +
		geom_bar(stat="identity",position="stack",colour="black") +
		geom_text(position="stack",hjust=1,family="mono",angle=90,fontface="bold") +
		scale_x_discrete("Position in read (bp)",limits=x_axis_breaks,breaks=x_axis_breaks,labels=x_axis_labels) +
		scale_y_continuous("Number of data sets") +
		facet_wrap(~Read,scales="free",ncol=plot_size$num_columns) +
		theme(legend.position="none",axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
		ggtitle("K-mer content")
	 }

	ggsave(paste(pdf_directory,paste(bfx_id,"kmer_content.pdf",sep="_"),sep="/"),plot=kmer_content_plot,width=plot_size$width,height=plot_size$height)
	write.table(fastqc_result$kmer_content[,-3],paste(data_directory,paste(bfx_id,"kmer_content.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
}else{
	file.create(paste(data_directory,paste(bfx_id,"kmer_content.csv",sep="_"),sep="/"))
	ggsave(grid.text("No kmer content plot possible"),file=paste(pdf_directory,paste(bfx_id,"kmer_content.pdf",sep="_"),sep="/"))
}


