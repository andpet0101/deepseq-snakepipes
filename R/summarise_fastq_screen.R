#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(scales)
theme_set(theme_bw(12))
library(reshape2)
library(dplyr)
library(tidyr)

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

#bfx_id = "nobfx"
#data_directory = "fastq_screen/report/data/"
#pdf_directory = "fastq_screen/report/pdf/"

#############
# Functions #
#############
source(paste(script_dir,"functions_settings.R",sep="/"))
#source("~/git/deepseq-snakepipes/R/functions_settings.R")

parse_fastq_screen_file = function(fastq_screen_file){
	fastq_screen_result = list()
	# library and read
	fastq_screen_basename = gsub('_screen.csv$','',basename(fastq_screen_file))
	fastq_screen_result[['Library']] = gsub('_(R\\d)$','',fastq_screen_basename)
	read_match = gsub('\\S+_(R\\d)$','\\1',gsub('_screen.csv$','',basename(fastq_screen_file)))
	fastq_screen_result[['Read']] = ifelse(grepl('^R\\d$',read_match),read_match,'R1')
	
	# read content
	fastq_screen_lines = readLines(fastq_screen_file)
	# skip first line (comment)
	fastq_screen_lines = fastq_screen_lines[-1]
	# fix the '#' and the '%' symbols
	fastq_screen_lines[1] = gsub("#","",fastq_screen_lines[1])
	fastq_screen_lines[1] = gsub("%","perc_",fastq_screen_lines[1])
	
	# "no hit" line
	no_hit_line_index = grep("Hit_no_libraries",fastq_screen_lines)
	no_hit_line = fastq_screen_lines[no_hit_line_index]
	fastq_screen_result[["No_hit"]] = as.double(gsub("^\\S+\\s+","",no_hit_line,perl=TRUE))
	fastq_screen_lines = fastq_screen_lines[-no_hit_line_index]
	
	# the rest contains the libraries
	fastq_screen_result[["databases"]] = read.table(text=fastq_screen_lines,sep="\t",header=T)
	col_names = colnames(fastq_screen_result[["databases"]])
	col_names[1] = "Database"
	colnames(fastq_screen_result[["databases"]]) = col_names
	fastq_screen_result
}

# function defines a (optimal)number of columns, number of rows, width and height for a facetted plot
calc_facet_plot_size = function(num_plots,max_plots=20){
	if(num_plots==1){
		num_columns = 1
		num_rows = 1
		width = num_columns*5 # 5
		height = num_rows*6 # 6
	} else if(num_plots>1 & num_plots<=4){
		num_columns = 2
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*3 # 6
		height = num_rows*4 # 8
	} else if(num_plots>4 & num_plots<10){
		num_columns = 3
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*3 # 9
		height = num_rows*4 # 12	
	} else{
		num_columns = 4
		num_rows = ceiling(num_plots/num_columns)
		width = num_columns*3.5 # 8
		height = num_rows*4 # 12
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
fastq_screen_files = list.files(path=data_directory,pattern="_screen\\.csv$",full.names=TRUE)
fastq_screen_result = list()
all_library_names = c()
all_read_dirs = c("R1","R2")

for(f in fastq_screen_files){
	fastq_screen = parse_fastq_screen_file(f)
	fastq_screen[["Library"]] = fixLibraryNames(fastq_screen[["Library"]])
	fastq_screen[["databases"]][["Library"]] = fastq_screen[["Library"]]
	fastq_screen[["databases"]][["Read"]] = fastq_screen[["Read"]]
	if("databases" %in% names(fastq_screen_result)){
		fastq_screen_result[["databases"]] = rbind(fastq_screen_result[["databases"]],fastq_screen[["databases"]])
	} else {
		fastq_screen_result[["databases"]] = fastq_screen[["databases"]]
	}
	if("no_hit" %in% names(fastq_screen_result)){
		fastq_screen_result[["no_hit"]] = rbind(fastq_screen_result[["no_hit"]],data.frame("Library"=fastq_screen[["Library"]],"Read"=fastq_screen[["Read"]],"No_hit"=fastq_screen[["No_hit"]]))
	} else {
		fastq_screen_result[["no_hit"]] = data.frame("Library"=fastq_screen[["Library"]],"Read"=fastq_screen[["Read"]],"No_hit"=fastq_screen[["No_hit"]])
	}
}

all_library_names = unique(c(fixLibraryNames(fastq_screen_result[["databases"]][["Library"]]),fixLibraryNames(fastq_screen_result[["no_hit"]][["Library"]])))
all_read_dirs = unique(c(fixLibraryNames(fastq_screen_result[["databases"]][["Read"]]),fixLibraryNames(fastq_screen_result[["no_hit"]][["Read"]])))

for(t in c("databases","no_hit")){
	fastq_screen_result[[t]][["Library"]] = factor(fastq_screen_result[[t]][["Library"]],levels=all_library_names)
	fastq_screen_result[[t]][["Read"]] = factor(fastq_screen_result[[t]][["Read"]],levels=all_read_dirs)
	fastq_screen_result[[t]][["Library_read"]] = factor(paste(fastq_screen_result[[t]][["Library"]],fastq_screen_result[[t]][["Read"]]),levels=apply(expand.grid(all_library_names,all_read_dirs),1,paste,collapse=" "))
	col_names = colnames(fastq_screen_result[[t]])
	col_names = c("Library","Read","Library_read",setdiff(col_names,c("Library","Read","Library_read")))
	fastq_screen_result[[t]] = fastq_screen_result[[t]][,col_names]
}

# plot database distribution (only perc_One_hit_one_library and perc_Multiple_hits_one_library)
libraries_per_row = 20
library_row = data.frame(Library=levels(fastq_screen_result$databases$Library),Row=ceiling(1:length(levels(fastq_screen_result$databases$Library))/libraries_per_row))
fastq_screen_result$databases$Row = NULL
fastq_screen_result$databases = merge(fastq_screen_result$databases,library_row,by="Library")
plot_size = calc_facet_plot_size(length(unique(fastq_screen_result$databases$Row)))
database_specific_hits_plot = ggplot(fastq_screen_result$databases,aes(x=Library,y=perc_One_hit_one_library+perc_Multiple_hits_one_library,fill=Database)) + 
	geom_bar(stat="identity",position="stack") +
	theme_bw(11) +
	scale_x_discrete("Library") + 
	scale_y_continuous("Percentage of reads",limits=c(0,100)) +
	scale_fill_manual("Database",values=color_brewer_qual_palette) +
	facet_wrap(~ Row,scales="free",ncol=plot_size$num_columns) +
	theme(strip.background=element_blank(),strip.text=element_blank(),legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
	ggtitle("Contamination: reads exclusive to database")

ggsave(paste(pdf_directory,paste(bfx_id,"contamination_database_specific_hits.pdf",sep="_"),sep="/"),plot=database_specific_hits_plot,width=plot_size$width,height=plot_size$height)
write.table(spread(fastq_screen_result$databases[,c("Library","Read","Library_read","Database","perc_One_hit_one_library")],Database,perc_One_hit_one_library,fill=0),paste(data_directory,paste(bfx_id,"contamination_database_specific_hits.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)

# plot No_hit stats
libraries_per_row = 20
fastq_screen_result$no_hit$Row = NULL
fastq_screen_result$no_hit = merge(fastq_screen_result$no_hit,library_row,by="Library")
plot_size = calc_facet_plot_size(length(unique(fastq_screen_result$no_hit$Row)))
no_hits_plot = ggplot(fastq_screen_result$no_hit,aes(x=Library,y=No_hit)) + 
	geom_bar(stat="identity") +
	theme_bw(11) +
	scale_x_discrete("Library") + 
	scale_y_continuous("Percentage of reads",limits=c(0,100)) +
	facet_wrap(~ Row,scales="free",ncol=plot_size$num_columns) +
	theme(strip.background=element_blank(),strip.text=element_blank(),legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
	ggtitle("Contamination: reads exclusive to database")
ggsave(paste(pdf_directory,paste(bfx_id,"no_hits.pdf",sep="_"),sep="/"),plot=no_hits_plot,width=plot_size$width,height=plot_size$height)
write.table(fastq_screen_result$no_hit[,c("Library","Read","Library_read","No_hit")],paste(data_directory,paste(bfx_id,"no_hits.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)


