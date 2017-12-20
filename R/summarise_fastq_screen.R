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
library(magrittr)

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

#bfx_id = "bfx853"
#data_directory = "fastq_screen/report/data"
#pdf_directory = "fastq_screen/report/pdf"

#############
# Functions #
#############
source(paste(script_dir,"functions_settings.R",sep="/"))
#source("~/git/deepseq-snakepipes/R/functions_settings.R")

parse_fastq_screen_file = function(fastq_screen_file){
	fastq_screen_result = list()
	# library and read
	fastq_screen_basename = gsub('_screen.csv$','',basename(fastq_screen_file))
	fastq_screen_result[['Basename']] = gsub('_(R\\d)$','',fastq_screen_basename)
	read_match = gsub('\\S+_(R\\d)$','\\1',gsub('_screen.csv$','',basename(fastq_screen_file)))
	libraryid = gsub('^(L\\d+)_\\S+','\\1',fastq_screen_result[['Basename']])
	sample = gsub('^L\\d+_','',fastq_screen_result[['Basename']])

	fastq_screen_result[['Library']] = libraryid
	fastq_screen_result[['Sample']] = sample
	fastq_screen_result[['Read']] = ifelse(grepl('^R\\d$',read_match),read_match,'R1')
	
	# read content
	fastq_screen_lines = readLines(fastq_screen_file)
	# skip first line (comment)
	fastq_screen_lines = fastq_screen_lines[-1]
	# fix the '#' and the '%' symbols
	fastq_screen_lines[1] = gsub("#","",fastq_screen_lines[1])
	fastq_screen_lines[1] = gsub("%","perc_",fastq_screen_lines[1])
	
	# "no hit" line
	no_hit_line_index = grep("Hit_no_genomes",fastq_screen_lines)
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
		if(num_plots<20){num_columns = 4}
		else if(num_plots<30){num_columns = 5}
		else{num_columns = 6}
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
if(length(fastq_screen_files)==0){
  stop("No *_screen.csv files found in the data directory!\n", call.=FALSE)
}
fastq_screen_result = list()
all_library_names = c()
all_read_dirs = c("R1","R2")

for(f in fastq_screen_files){
	fastq_screen = parse_fastq_screen_file(f)
	fastq_screen[["databases"]][["Library"]] = fastq_screen[["Library"]]
	fastq_screen[["databases"]][["Sample"]] = fastq_screen[["Sample"]]
	fastq_screen[["databases"]][["Read"]] = fastq_screen[["Read"]]
	if("databases" %in% names(fastq_screen_result)){
		fastq_screen_result[["databases"]] = rbind(fastq_screen_result[["databases"]],fastq_screen[["databases"]])
	} else {
		fastq_screen_result[["databases"]] = fastq_screen[["databases"]]
	}
	if("no_hit" %in% names(fastq_screen_result)){
		fastq_screen_result[["no_hit"]] = rbind(fastq_screen_result[["no_hit"]],data.frame("Library"=fastq_screen[["Library"]],"Sample"=fastq_screen[["Sample"]],"Read"=fastq_screen[["Read"]],"Reads_processed"=fastq_screen[["databases"]][["Reads_processed"]][1],"No_hit"=fastq_screen[["No_hit"]]))
	} else {
		fastq_screen_result[["no_hit"]] = data.frame("Library"=fastq_screen[["Library"]],"Sample"=fastq_screen[["Sample"]],"Read"=fastq_screen[["Read"]],"Reads_processed"=fastq_screen[["databases"]][["Reads_processed"]][1],"No_hit"=fastq_screen[["No_hit"]])
	}
}

all_libraries = unique(c(as.character(fastq_screen_result[["databases"]][["Library"]]),as.character(fastq_screen_result[["no_hit"]][["Library"]])))
all_samples = unique(c(as.character(fastq_screen_result[["databases"]][["Sample"]]),as.character(fastq_screen_result[["no_hit"]][["Sample"]])))
all_read_dirs = unique(c(as.character(fastq_screen_result[["databases"]][["Read"]]),as.character(fastq_screen_result[["no_hit"]][["Read"]])))

for(t in c("databases","no_hit")){
	fastq_screen_result[[t]][["Library"]] = factor(fastq_screen_result[[t]][["Library"]],levels=all_libraries)
	fastq_screen_result[[t]][["Sample"]] = factor(fastq_screen_result[[t]][["Sample"]],levels=all_samples)
	fastq_screen_result[[t]][["Read"]] = factor(fastq_screen_result[[t]][["Read"]],levels=all_read_dirs)
	fastq_screen_result[[t]][["Library_read"]] = factor(paste(fastq_screen_result[[t]][["Library"]],fastq_screen_result[[t]][["Read"]]),levels=apply(expand.grid(all_libraries,all_read_dirs),1,paste,collapse=" "))
	col_names = colnames(fastq_screen_result[[t]])
	col_names = c("Library","Sample","Read","Library_read",setdiff(col_names,c("Library","Sample","Read","Library_read")))
	fastq_screen_result[[t]] = fastq_screen_result[[t]][,col_names]
}

# add to database distribution: no hits to any database
databases_plot_data = fastq_screen_result$databases[,c("Database","Sample","Library","Read","Library_read","Reads_processed","perc_One_hit_one_genome","perc_Multiple_hits_one_genome")]
databases_plot_data = rbind(databases_plot_data,data.frame(Database="No hit",
          fastq_screen_result$no_hit[,c("Sample","Library","Read","Library_read","Reads_processed")],
           perc_One_hit_one_genome=fastq_screen_result$no_hit$No_hit,
           perc_Multiple_hits_one_genome=0))
# add to database distribution: ambiguous hits (cannot decide)
ambig = databases_plot_data %>% group_by(Sample,Library,Read,Library_read,Reads_processed) %>%
  summarise(Ambiguous=100-sum(perc_One_hit_one_genome+perc_Multiple_hits_one_genome)) %>%
  as.data.frame()
databases_plot_data = rbind(databases_plot_data,data.frame(Database="Ambiguous",
          ambig[,c("Sample","Library","Read","Library_read","Reads_processed")],
          perc_One_hit_one_genome=ambig$Ambiguous,
          perc_Multiple_hits_one_genome=0))

# sort data and database levels by decreasing mean percentage
mean_perc_per_database = databases_plot_data %>% group_by(Database) %>%
  summarise(Mean_perc = mean(perc_One_hit_one_genome+perc_Multiple_hits_one_genome)) %>%
  arrange(desc(Mean_perc)) %>% as.data.frame()
databases_plot_data$Database = factor(as.character(databases_plot_data$Database),levels = as.character(mean_perc_per_database$Database))

# group 20 libraries per plot
libraries_per_row = 20
library_row = data.frame(Library=levels(databases_plot_data$Library),Row=ceiling(1:length(levels(databases_plot_data$Library))/libraries_per_row))
databases_plot_data$Row = NULL
databases_plot_data = merge(databases_plot_data,library_row,by="Library")
plot_size = calc_facet_plot_size(length(unique(databases_plot_data$Row)))

# plot database distribution
database_specific_hits_plot = ggplot(databases_plot_data,aes(x=Sample,y=perc_One_hit_one_genome+perc_Multiple_hits_one_genome,fill=Database)) + 
 	geom_bar(stat="identity",position=position_stack(reverse=T)) +
 	theme_bw(11) +
 	scale_x_discrete("Sample") + 
 	scale_y_continuous("Percentage of reads",limits=c(0,100)) +
 	scale_fill_manual("Species",values=color_brewer_qual_palette) +
 	facet_wrap(~ Row,scales="free",ncol=plot_size$num_columns) +
 	theme(strip.background=element_blank(),strip.text=element_blank(),legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
 	ggtitle("Reads assigned to species")

ggsave(paste(pdf_directory,paste(bfx_id,"contamination_species_specific_hits.pdf",sep="_"),sep="/"),plot=database_specific_hits_plot,width=plot_size$width,height=plot_size$height)

one_genome_hits_table = data.frame(Library=databases_plot_data$Library,
				Sample=databases_plot_data$Sample,
				Read=databases_plot_data$Read,
				Database=databases_plot_data$Database,
				Reads_processed=databases_plot_data$Reads_processed,
				hits=databases_plot_data$perc_One_hit_one_genome+databases_plot_data$perc_Multiple_hits_one_genome)
one_genome_hits_table$Database = factor(one_genome_hits_table$Database,levels=sort(unique(as.character(one_genome_hits_table$Database))))
one_genome_hits_table = spread(one_genome_hits_table,Database,hits,fill=0)
write.table(one_genome_hits_table,paste(data_directory,paste(bfx_id,"contamination_species_specific_hits.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)

# plot No_hit stats
libraries_per_row = 20
fastq_screen_result$no_hit$Row = NULL
fastq_screen_result$no_hit = merge(fastq_screen_result$no_hit,library_row,by="Library")
plot_size = calc_facet_plot_size(length(unique(fastq_screen_result$no_hit$Row)))
no_hits_plot = ggplot(fastq_screen_result$no_hit,aes(x=Library,y=No_hit)) + 
	geom_bar(stat="identity") +
	theme_bw(11) +
	scale_x_discrete("Sample",labels=unique(databases_plot_data[order(databases_plot_data$Library),"Sample"])) + 
	scale_y_continuous("Percentage of reads",limits=c(0,100)) +
	facet_wrap(~ Row,scales="free",ncol=plot_size$num_columns) +
	theme(strip.background=element_blank(),strip.text=element_blank(),legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
	ggtitle("Reads with no hits")
ggsave(paste(pdf_directory,paste(bfx_id,"no_hits.pdf",sep="_"),sep="/"),plot=no_hits_plot,width=plot_size$width,height=plot_size$height)
write.table(fastq_screen_result$no_hit[,c("Library","Sample","Read","No_hit")],paste(data_directory,paste(bfx_id,"no_hits.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)


