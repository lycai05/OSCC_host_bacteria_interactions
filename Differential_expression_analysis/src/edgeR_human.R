##########
##### identify differentially expressed genes
##########

args = commandArgs(trailingOnly=TRUE)
# meta csv file
meta = args[1]
# working directory
wk_dir = args[2]
# output prefix
prefix = args[3]
setwd(paste(wk_dir, "featureCounts", sep = "/"))

#pk_extension = as.numeric(args[4])

print(meta)
print(wk_dir)
print(prefix)
print(getwd())

library(edgeR)
library(ggplot2)
library("pheatmap")
library(plyr)


gene_names <- scan("/home/liuyang/Projects/DO_DN/results/diff_expr/featureCounts/gene_names.txt", what = "character")
map_gene_names <- read.table("/home/liuyang/Annotation/hg38/genename_to_symbol.txt", header = F)
map_gene_names[[1]] <- as.character(map_gene_names[[1]])
map_gene_names[[2]] <- as.character(map_gene_names[[2]])

#diff_expr <- function(meta, prefix, wk_dir, gene_names,  map_gene_names) {
	meta <- read.csv(meta, header = T)
	meta[[1]] <- as.character(meta[[1]])
	meta[[2]] <- as.factor(meta[[2]])
	merged_counts <- lapply(meta[[1]], function(x){print(head(x));a <- scan(list.files(path = paste(wk_dir, "featureCounts", sep = "/"), pattern = x), what = "numeric"); c(x, a)})
	merged_counts_mx <- do.call(cbind, merged_counts)
	merged_counts_col <- merged_counts_mx[1,]
	merged_counts_mx <- merged_counts_mx[-1,]

	rownames(merged_counts_mx) <- gene_names
	colnames(merged_counts_mx) <- merged_counts_col
	class(merged_counts_mx) <- "numeric"


	# mx <- read.table("DCM_NF_merged_mx_for_diffbind.txt", header = F)
	# mx_mx <- as.matrix(mx[, 2:7])
	# rownames(mx_mx) <- mx[[1]]
	# colnames(mx_mx) <- rev(c("Y20_HK", "Y123_HK", "Y234_HK", "Y591_HK", "Y600_HK", "Y647_HK"))
	#cpm_log <- cpm(mx_mx, log = TRUE)
	#median_log2_cpm <- apply(cpm_log, 1, median)
	#hist(median_log2_cpm)
	#expr_cutoff <- -1
	#abline(v = expr_cutoff, col = "red", lwd = 3)

	#sum(median_log2_cpm > expr_cutoff)
	#data_clean <- mx_mx[median_log2_cpm > expr_cutoff, ]
	data_clean <- merged_counts_mx[rowSums(merged_counts_mx) > 0, ]


	my_sample_col <- data.frame(meta$Condition)
	#my_sample_col$Condition <- factor(meta$Condition)
	colnames(my_sample_col) <- "Condition"
	my_colour = list(
	    Condition = c(Control = "#5977ff", Treatment = "#f74747")
	)
	row.names(my_sample_col) <- colnames(data_clean)

	# sample heatmap
	cpm_log <- cpm(data_clean, log = TRUE)
	pdf(paste(prefix, "rnaseq_heatmap_cor_clutsering_plot.pdf", sep = "_"), width = 8, height = 8)
	pheatmap(cor(cpm_log), annotation_col = my_sample_col, annotation_colors = my_colour)
	dev.off()



	# PCA
	pca <- prcomp(t(cpm_log))
	pdf(paste(prefix, "rnaseq_PCA_plot.pdf", sep = "_"), width = 8, height = 8)
	plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
	text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
	dev.off()
	#summary(pca)

	#group <- c(rep("DCM", 3), rep("NF", 3))
	#group <- factor(c(2,2,2,1,1,1))
	#gender <- rev(c("F", "F", "M", "F", "F", "M"))

	y <- DGEList(data_clean)
	y <- calcNormFactors(y)

	meta$Condition <- factor(meta$Condition, levels =c("Control", "Treatment"))
	design <- model.matrix(~meta$Condition)
	design

	y <- estimateDisp(y, design)
	fit <- glmFit(y, design)
	lrt <- glmLRT(fit)
	resFilt <- topTags(lrt, n = nrow(data_clean), sort.by = "none")
   

    expr_table <- resFilt@.Data[[1]]
    expr_table$gene_symbols <- mapvalues(rownames(expr_table), map_gene_names[[1]], map_gene_names[[2]], warn_missing = FALSE)
    cpm_counts_log <- cpm(data_clean, log = TRUE)
	expr_df <- merge(expr_table, cpm_counts_log, by = "row.names", all = F)
    write.table(expr_df, paste(prefix, "diff_expression_edgeR_qvalue.txt", sep = "_"), col.names = T, row.names = T, quote = F, sep = "\t")


    expr_df2 <- merge(expr_df, data_clean, by.x = "Row.names", by.y = "row.names", all = F)

	# normal down
	# up-regulated
	sigUpReg <- expr_df2[expr_df2$FDR<0.05 & expr_df2$logFC > 1,]
print(head(sigUpReg))
	# normal high
	# down regulated
	sigDownReg <- expr_df2[expr_df2$FDR<0.05 & expr_df2$logFC < -1,]


	# sigUpReg2_kp_cols <- (apply(sigUpReg[,c(39:47)], 1, function(x) sum(x>10))>0.5)
	# sigUpReg2 <- sigUpReg[names(sigUpReg2_kp_cols[sigUpReg2_kp_cols==T]),]


	# sigDownReg2_kp_cols <- (apply(sigDownReg[,c(28:38)], 1, function(x) sum(x>10))>0.5)
	# sigDownRe2 <- sigDownReg[names(sigDownReg2_kp_cols[sigDownReg2_kp_cols==T]),]


	# sigUpReg3_kp_cols <- (apply(sigUpReg[,c(28:47)], 1, function(x) sum(x>10))>16)
	# sigUpReg3 <- sigUpReg[names(sigUpReg3_kp_cols[sigUpReg3_kp_cols==T]),]


	# sigDownReg3_kp_cols <- (apply(sigDownReg[,c(28:47)], 1, function(x) sum(x>10))>16)
	# sigDownReg3 <- sigDownReg[names(sigDownReg3_kp_cols[sigDownReg3_kp_cols==T]),]




 # Keep only those columns
 	expr_df2$diff_expr <- apply(expr_df2, 1, function(x){if(x[["Row.names"]] %in% sigUpReg[["Row.names"]]){"Treatment"}
 		else if(x[["Row.names"]] %in% sigDownReg[["Row.names"]]){"Control"}
  		else{"common"}})
    print(head(expr_df2))

   #  write.table(expr_df2, paste(prefix, "diff_expression_edgeR_qvalue_w_raw_counts_updated.txt", sep = "_"), col.names = T, row.names = F, quote = F, sep = "\t")



	#expr_df <- resFilt@.Data[[1]]
	#expr_df$ENSEMBL <- gsub("\\.[0-9]*", "", rownames(expr_df))


	# MA plot for differentially expressed genes
	theme_pub <- theme_bw() + theme(panel.grid.major = element_blank(),
	                                panel.grid.minor = element_blank(),
	                                axis.text = element_text(size = 20),
	                                #axis.title.y = element_text(vjust=1.5),
	                                #axis.title.x = element_text(hjust=0),
	                                axis.title = element_text(size = 24),
	                                #axis.line = element_line(color = "black", size = 1),
	                                #axis.ticks =element_line(size = 0.8),
	                                #axis.title.x = element_blank(),
	                                #axis.title.y = element_blank(),
	                                #axis.ticks.x=element_blank(),
	                                #axis.text.x = element_blank(),
	                                strip.text = element_text(size = 20),
	                                panel.border = element_rect(color = "black", size = 1),
	                                legend.position="none")

	plot_diff_express_genes <- function(file, anno_genes, qvalue, log2fc, methods) {
	  #file$color <- apply(file, 1, function(x){if(is.na(x['FDR'])){"Common"}
	  #  else {a <- as.numeric(as.character(x['FDR']));b <- as.numeric(as.character(x['logFC']));
	  #  if((a <= qvalue) & (b >= log2fc)){"Tumor"}
	  #  else if((a <= qvalue) & (b <= -log2fc)){"Control"}
	  #  else{"Common"}}})
	  file$color <- file$diff_expr
	  file <- subset(file, color != "NA")
	  #print(head(file))
	  file$color <- factor(file$color, levels = c("Treatment", "Control", "common"))
  	  n_tumor <- nrow(subset(file, color == "Treatment"))
  	  n_an <- nrow(subset(file, color == "Control"))
      n_common <- nrow(subset(file, color == "common"))
	  #genes2 <- paste(anno_genes, collapse="|")
	  #file1 <- file[grepl(genes2, file$SYMBOL), ]
	  #file1$SYMBOL  <- as.character(file1$SYMBOL)
	  #file1 <- file1[!(file1$SYMBOL %in% c("APOBEC2", "APOBEC4", "APOBEC3B-AS1", "CDKN2B-AS1", "APOBEC3G", "APOBEC3D", "APOBEC3A", "APOBEC3B", "APOBR", "APOBEC3C", "APOBEC4", "APOBEC2", "BRAFP1", "GGT5", "GGT7", "CDKN2AIPNL", "NARG2", "PCGT1B", "UGGT2", "UGGT1", "RABGGTB", "CDKN2AIP", "RABGGTA", "GGT2", "GGTLC1", 
	  #                                     "GGTLC2", "ADGRB2", "PGGT1B", "CEBPA-AS1", "CYP1B1-AS1", "GGTA1P", "ADH5P4", "APOBEC3H", "APOBEC3F", "BCRP6", "BCRP1", "SAA2-SAA4", "NCCRP1", "BCRP3", "BCRP2", "TECRP1", "BCRP8")), ]
	  #for (j in 1:nrow(file1)) {
	  #  file1$SYMBOL[j] <- grep(genes2, unlist(strsplit(file1$SYMBOL[j], split=",")), value = T)
	  #}
	  o <- ggplot(file, aes(x = logCPM, y = logFC, color = color)) + geom_point() +
	    scale_color_manual(values = c("#BE0005", "#206095", "grey"))+ guides(color=FALSE) + xlab("Log2(Mean)") +ylab("Log2(Tumor/Control)") +
	    #geom_point(data = file1, aes(x = log2(baseMean), y = log2FoldChange), color = "orange", size = 3)  +
	    scale_y_continuous(limits = c(-8,8), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))+ theme_pub +
		annotate("text", label = paste("N = ", c(n_an, n_common,  n_tumor), sep = ""), x = 10, y = c(-8, 0, 8), size = 5) 
	    # + #scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5), labels = c(1, 10, 100, 1000, 10000, 100000)) +  #scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6), labels = c(2, 4, 8, 16, 32, 6))
	    #geom_label_repel(data = file1, aes(x = log2(baseMean), y = log2FoldChange, label =  SYMBOL), color = "black", label.size = 3)  
	  ggsave(plot = o, filename = paste(prefix, "report_MA_plot", methods ,".pdf", sep = "_"), width = 7, height = 5)
	  return(0)
	}
	plot_diff_express_genes(expr_df2, "", 0.1, 0, "edgeR")


genename_to_biotypes <- read.table("/home/liuyang/Annotation/hg38/genename_to_biotypes.txt", header = F)
expr_df2 <- merge(expr_df2, genename_to_biotypes, by.x = "Row.names", by.y = "V1", all.x = T)


    write.table(expr_df2, paste(prefix, "diff_expression_edgeR_qvalue_w_raw_counts.txt", sep = "_"), col.names = T, row.names = T, quote = F, sep = "\t")














save.image(paste(prefix, ".RData", sep = ""))

#awk 'FS="_"{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' sigDownReg.txt > sigDownReg.bed
#awk 'FS="_"{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' sigUpReg.txt > sigUpReg.bed




