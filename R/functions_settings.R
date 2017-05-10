

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


# labels for plots
bp_labels = function(l){format(big.mark=",",l)}
kb_labels = function(l,r=0){format(big.mark=",",round(l/1000,r))}
Mb_labels = function(l,r=0){format(big.mark=",",round(l/1000000,r))}
Gb_labels = function(l,r=0){format(big.mark=",",round(l/1000000000,r))}


# colour brewer palette for 24 samples; repeat 100 times = 2400 samples max
color_brewer_qual_palette = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7","#ffffb3",
"#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
color_brewer_qual_palette = rep(color_brewer_qual_palette,100)


# from Matthias' DESeq methods
plotCorrelation <- function(deseqmatrix, conditions = "", cormethod = "pearson", cluster = FALSE, filename = "", title = "Correlation") {
  require(pheatmap)
  require(RColorBrewer)
  
  cormatrix <- cor(deseqmatrix, use = "pairwise.complete.obs", method = cormethod)
  cornumbers <- round(cormatrix, 3)
  # colnames(cormatrix) <- NULL
  if (cormethod == "pearson") {
    colors <- colorRampPalette(brewer.pal(8, "Blues"))(255)
  }
  if (cormethod == "spearman") {
    colors <- colorRampPalette(brewer.pal(8, "Reds"))(255)
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

# mixedorder and mixedsort needed to be redefined from gtools: [.] should be [\\.] to avoid POSIX collating errors
mixedorder1 = function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE, 
    numeric.type = c("decimal", "roman"), roman.case = c("upper", 
        "lower", "both")) 
{
    numeric.type <- match.arg(numeric.type)
    roman.case <- match.arg(roman.case)
    if (length(x) < 1) 
        return(NULL)
    else if (length(x) == 1) 
        return(1)
    if (!is.character(x)) 
        return(order(x, decreasing = decreasing, na.last = na.last))
    delim = "\\$\\@\\$"
    if (numeric.type == "decimal") {
        regex <- "((?:(?i)(?:[-+]?)(?:(?=[\\.]?[0123456789])(?:[0123456789]*)(?:(?:[\\.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
        numeric <- function(x) as.numeric(x)
    }
    else if (numeric.type == "roman") {
        regex <- switch(roman.case, both = "([IVXCLDMivxcldm]+)", 
            upper = "([IVXCLDM]+)", lower = "([ivxcldm]+)")
        numeric <- function(x) roman2int(x)
    }
    else stop("Unknown value for numeric.type: ", numeric.type)
    nonnumeric <- function(x) {
        ifelse(is.na(numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    delimited <- gsub(regex, paste(delim, "\\1", delim, sep = ""), 
        x, perl = TRUE)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) x[x > ""])
    suppressWarnings(step1.numeric <- lapply(step1, numeric))
    suppressWarnings(step1.character <- lapply(step1, nonnumeric))
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) sapply(step1.numeric, 
        function(x) x[i]))
    step1.character.t <- lapply(1:maxelem, function(i) sapply(step1.character, 
        function(x) x[i]))
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) as.numeric(factor(x)))
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric), 
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric, 
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0) 
        if (is.na(na.last)) 
            order.frame[which.nas, ] <- NA
        else if (na.last) 
            order.frame[which.nas, ] <- Inf
        else order.frame[which.nas, ] <- -Inf
    if (length(which.blanks) > 0) 
        if (is.na(blank.last)) 
            order.frame[which.blanks, ] <- NA
        else if (blank.last) 
            order.frame[which.blanks, ] <- 1e+99
        else order.frame[which.blanks, ] <- -1e+99
    order.frame <- as.list(order.frame)
    order.frame$decreasing <- decreasing
    order.frame$na.last <- NA
    retval <- do.call("order", order.frame)
    return(retval)
}

mixedsort1 = function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE, 
    numeric.type = c("decimal", "roman"), roman.case = c("upper", 
        "lower", "both")) 
{
    ord <- mixedorder1(x, decreasing = decreasing, na.last = na.last, 
        blank.last = blank.last, numeric.type = numeric.type, 
        roman.case = roman.case)
    x[ord]
}



