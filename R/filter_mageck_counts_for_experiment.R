#!/usr/bin/env Rscript --vanilla

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
if(length(arguments)!=5){
	stop("Needs the mageck counts table, the treatements (separated by comma), the controls (separated by comma), the minimum read count per row and the name of the filtered counts table\n", call.=FALSE)
}

total_counts_table_name = arguments[1]
treatments = arguments[2]
controls = arguments[3]
min_read_count = as.integer(arguments[4])
filtered_counts_table_name = arguments[5]

#total_counts_table_name = "mageck/count/all.count.txt"
#treatments = "Erastin"
#controls = "Control"
#min_read_count = 10
#filtered_counts_table_name = "test.csv"


# read in total counts, then select sgRNA, Gene, treatments and controls columns, filter and write to stdout
counts_table = read.table(total_counts_table_name,header=T,sep="\t",check.names=F)
columns_to_select = c("sgRNA","Gene",strsplit(treatments,",")[[1]],strsplit(controls,",")[[1]])
counts_table = counts_table[,columns_to_select]

row_sums = rowSums(counts_table[,-1:-2])
counts_table = counts_table[which(row_sums>=min_read_count),]

write.table(counts_table,filtered_counts_table_name,sep="\t",row.names=F,col.names=T,quote=F)

