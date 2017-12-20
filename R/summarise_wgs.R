data_directory = "."
insertsize_files = list.files(path=data_directory,pattern="insertsize_metrics.csv",full.names=TRUE)

insertsize_result = list()

for(f in insertsize_files){
	insertsize_metrics = read.table(f,skip=6,header=T,nrows=1,sep="\t")
	insertsize_metrics$Library = gsub("^\\S+_L\\d+_","",gsub("_insertsize_metrics\\.csv","",basename(f)))
	if("metrics" %in% names(insertsize_result)){
		insertsize_result[["metrics"]] = rbind(insertsize_result[["metrics"]],insertsize_metrics)
	} else {
		insertsize_result[["metrics"]] = insertsize_metrics
	}	
	
	insertsize_distribution = read.table(f,skip=10,header=T,sep="\t")
	insertsize_distribution$fraction = insertsize_distribution$All_Reads.fr_count/sum(insertsize_distribution$All_Reads.fr_count)
	insertsize_distribution$Library = gsub("^\\S+_L\\d+_","",gsub("_insertsize_metrics\\.csv","",basename(f)))
	if("distribution" %in% names(insertsize_result)){
		insertsize_result[["distribution"]] = rbind(insertsize_result[["distribution"]],insertsize_distribution)
	} else {
		insertsize_result[["distribution"]] = insertsize_distribution
	}
}

color_brewer_qual_palette = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7","#ffffb3",
"#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
color_brewer_qual_palette = rep(color_brewer_qual_palette,100)

ggplot(insertsize_result$distribution,aes(x=insert_size,y=fraction,colour=Library)) + 
geom_line() +
theme_bw(12) +
scale_x_continuous("Insert size [bp]") +
scale_y_continuous("Fraction of fragments") +
scale_colour_manual("Library",values=color_brewer_qual_palette)+
ggtitle("Insert size distribution") +
theme(legend.position="none")


