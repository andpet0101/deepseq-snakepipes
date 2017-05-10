#!/usr/bin/env Rscript --vanilla

library(plyr)
library(ggplot2)
theme_set(theme_bw(12))
library(reshape2)
library(scales)
library(data.table)


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

if(length(arguments)!=4){
	stop("Needs the bfx id, the data directory, the fastq directory and the pdf directory as arguments\n", call.=FALSE)
}

bfx_id = arguments[1]
data.dir = arguments[2]
fastq.dir = arguments[3]
plots.dir = arguments[4]

## load files 
files = gsub(
  "(^.*)\\.cls$",
  "\\1.txt",
  list.files(path=data.dir,pattern="^.*\\.cls$",full.names=TRUE)
)

counts = ldply(files, function(x){
  df <- fread(x, skip = 1)[,c(1,3), with=FALSE]
  setnames(df, gsub(" ", ".", colnames(df)))
  df$file <- rep(gsub("(.*)\\.txt", "\\1", basename(x)), nrow(df))  

  return(df)
  }
)

# ## create the factors for plotting:
counts$file = gsub("_Track.+","",counts$file)
counts$file <- factor(counts$file)

p.bar <- ggplot(counts, aes(x=file)) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), axis.title.y=element_text(vjust=0.1)) +
  geom_bar() + 
  scale_y_continuous(labels = comma) +
  ggtitle("Number of detected sequences\n(clonotypes)") + 
  xlab("Library") + 
  ylab("Number of sequences") 

ggsave(paste(plots.dir,paste(bfx_id,"number_of_sequences_per_sample.pdf",sep="_"),sep="/"), p.bar)

p.den <- ggplot(counts, aes(x=Read.count)) + 
  geom_line(stat="density") +
  facet_wrap( ~ file) +
  scale_x_log10("Number of reads",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous("Number of sequences") +
  ggtitle("Distribution of reads per sequence\n(clonotypes)")

ggsave(paste(plots.dir,paste(bfx_id,"distribution_of_reads_per_sample_density.pdf",sep="_"),sep="/"), p.den, height=7, width=9)

## histogram
p.hist <- ggplot(counts, aes(x=Read.count)) + geom_histogram(binwidth=1, position="dodge") + 
  facet_wrap(~file) + 
  ggtitle("Distribution of reads per sequence\n(clonotypes)") +
  scale_y_continuous("Number of sequences") +
  scale_x_log10("Number of reads",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))

ggsave(paste(plots.dir,paste(bfx_id,"distribution_of_reads_per_sample_histogram.pdf",sep="_"),sep="/"), p.hist, height=7, width=9)

# proper frequency plot:
p.freq <- ggplot(counts, aes(x=Read.count)) + 
  geom_freqpoly(binwidth = 1) + 
  scale_y_continuous(labels = comma) + 
  facet_wrap(~file) + 
  ggtitle("Distribution of reads per sequence\n(clonotypes)") +
  scale_x_log10("Number of reads",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous("Number of sequences",labels=comma)

ggsave(paste(plots.dir,paste(bfx_id,"distribution_of_reads_per_sample_freqpoly.pdf",sep="_"),sep="/"), p.freq, height=7, width=9)

## summarize sum of counts
counts_sum <- ddply(counts, c("file"), summarize, count.sums=sum(Read.count))
counts_sum$file = factor(counts_sum$file)
counts_sum$lib = gsub("^(L\\d+)_.+","\\1",counts_sum$file)
counts_sum$lib = factor(counts_sum$lib)

p.sum <- ggplot(counts_sum, aes(x=file, y=count.sums, fill=file)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9)) + 
  scale_y_continuous(labels = comma) +  
  ggtitle("Total number of assigned reads") + 
  xlab("Library") + 
  ylab("Number of reads") + 
  scale_fill_discrete("Library") +
  geom_text(aes(label=comma(count.sums)), vjust=0.5,hjust=1, size=3,angle=90) +
  # facet_wrap(~library) + 
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1),legend.position="none")


ggsave(paste(plots.dir,paste(bfx_id,"total_reads_per_sample.pdf",sep="_"),sep="/"), p.sum, height=7, width=10)


#########################
## reads count table:
fastqs=list.files(path=fastq.dir,pattern=".fastq.gz",full.names=TRUE)
libids=gsub("(.*)\\.fastq.gz", "\\1", basename(fastqs))


n_reads <- ldply(fastqs, function(x){
  lib=gsub("(.*)\\.fastq.gz", "\\1", basename(x))
  print(lib)
  command=paste("zcat ", x, " | echo $((`wc -l`/4))", sep="")
  total_reads=as.integer(system(command, intern=TRUE))
  print(total_reads)
  df=data.frame(
    file=lib,
    total.sequenced.reads=total_reads)
  return(df)
  }
)
n_reads$file <- factor(n_reads$file)
n_reads$lib = gsub("^(L\\d+)_.+","\\1",n_reads$file)
n_reads$lib <- factor(n_reads$lib)

final_table <- merge(n_reads, counts_sum, by="lib")
final_table <- transform(final_table, percentage=count.sums/total.sequenced.reads*100)

final_table = final_table[,c("lib","file.x","file.y","total.sequenced.reads","count.sums","percentage")]
colnames(final_table) = c("ID","InternalName","ExternalName","Reads","Counts","Percentage")

n_clones <- ddply(counts, c("file"), summarize, number_clonotypes=length(CDR3.nucleotide.sequence))
n_clones$ID = gsub("^(L\\d+)_.+","\\1",n_clones$file)
n_clones$ID <- factor(n_clones$ID)
n_clones$file = NULL

final_table <- merge(final_table, n_clones, by="ID")
colnames(final_table) = c("ID","InternalName","ExternalName","Reads","Counts","Percentage","Num_Clonotypes")

write.table(final_table, file=paste(data.dir,paste(bfx_id,"percentage_of_reads_used_in_MiTCR.csv",sep="_"),sep="/"), sep="\t", row.names = FALSE, quote=FALSE)
