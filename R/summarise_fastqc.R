#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(scales)
theme_set(theme_bw(12))
library(reshape2)

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

# read in cutadapt and prepare trimming data frame
fastqc_files = list.files(path=data_directory,pattern="_fastqc.txt$",full.names=TRUE)

fastqc_overview = ldply(fastqc_files,function(f){
	fastqc_out = read.table(pipe(paste('cat ',f,' | grep "^>>" | grep -v "END" | sed "s/>>//"')),sep="\t",header=F)
	colnames(fastqc_out) = c("metric","status")
	fastqc_out$library = gsub('_(R\\d)$','',gsub('_fastqc.txt$','',basename(f)))
	read_match = gsub('\\S+_(R\\d)$','\\1',gsub('_fastqc.txt$','',basename(f)))
	fastqc_out$read = ifelse(grepl('^R\\d$',read_match),read_match,'R1')
	fastqc_out = fastqc_out[,c("library","read","metric","status")]
	fastqc_out
})
	
fastqc_overview$library = factor(fastqc_overview$library,levels=unique(fastqc_overview$library))
fastqc_overview$metric = factor(fastqc_overview$metric,levels=unique(fastqc_overview$metric))
fastqc_overview$status = factor(fastqc_overview$status,levels=c("pass","warn","fail"))
fastqc_overview$read = factor(fastqc_overview$read,levels=unique(fastqc_overview$read))
	
fastqc_overview_plot = ggplot(fastqc_overview,aes(x=metric,y=library)) + 
geom_tile(aes(fill=status),colour="black") + 
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
scale_x_discrete("Metric") + 
scale_y_discrete("Library") +
scale_fill_manual("Status",values=c("#91cf60","#ffffbf","#fc8d59")) +
facet_wrap(~ read)
	
ggsave(paste(pdf_directory,paste(bfx_id,"fastqc_overview.pdf",sep="_"),sep="/"),plot=fastqc_overview_plot,width=10,height=8)
write.table(fastqc_overview,paste(data_directory,paste(bfx_id,"fastqc_overview.csv",sep="_"),sep="/"),col.names=T,row.names=F,sep="\t",quote=F)


