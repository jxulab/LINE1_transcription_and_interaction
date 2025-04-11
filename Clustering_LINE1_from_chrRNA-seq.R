###############draw heatmap of all LINE1 and for each LINE1 family (exclude sense LINE1).
library(gplots,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
norm_method=c('TPM')
for (n in 1:length(norm_method)){
	normed_line1_exp <- read.table(paste0('Exp_matrix/line1_',norm_method[n],'_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt'),header=TRUE,sep="\t")[,-c(9,10,11,12)]
	normed_line1_exp_matrix <- as.matrix(normed_line1_exp)
	data_for_fig <- log2(normed_line1_exp_matrix+1)
	data_for_fig[data_for_fig>1]=1
	line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(4,6)]
	subfamily_anno <- merge(normed_line1_exp_matrix,line1,by.x=0,by.y='V4')
	subfamily_anno <- subfamily_anno[,c(1,dim(subfamily_anno)[2])]
	line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5')
	#subfamily_color <- c('black','#33ABC1','#8D73BA','red','grey')
	subfamily_color <- c('red','green','blue','yellow','grey')
	subfamily_anno[,3] <- 'white'
	for (i in 1:length(line1_subfamily_need_color)){
		subfamily_anno[subfamily_anno[,2]==line1_subfamily_need_color[i],3] <- subfamily_color[i]
	}
	pdf(paste0('Heatmap/All_full_length_LINE1_clustering_distance_euclidean_merged_cell_lines_add_subfamily_anno_rm_sense_LINE1_',norm_method[n],'_10celllines_darkcolor.pdf'), width = 5, height = 5,useDingbats=FALSE)
	heatmap.2(data_for_fig,
			  scale = "none",  # You can choose "row", "column", or "none" for scaling
			  main = paste0("all LINE1 (",norm_method[n],")"),
			  distfun = function(x) dist(x, method = "euclidean"),  # Use correlation as the distance metric
			  hclustfun = function(x) hclust(x, method = "complete"),  # Specify clustering method
			  margins = c(7, 7),  # Add extra space for labels
			  #col = rev(colorRampPalette(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"))(50)),
			  col = rev(colorRampPalette(c("#67001f", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
			  #col = rev(colorRampPalette(c("#de4128", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
			  key = TRUE,  # Display the color legend
			  trace = "none",cexRow = 0.1,cexCol=0.5,RowSideColors=subfamily_anno[,3])
	dev.off()
	##draw heatmap of all full-length LINE1 with subfamily annotation.

	######only draw heatmap of L1HS,L1PA2 and L1PA3.
	young_idx=(subfamily_anno[,2]=='L1HS')|(subfamily_anno[,2]=='L1PA2')|(subfamily_anno[,2]=='L1PA3')
	subfamily_anno=subfamily_anno[young_idx,]
	data_for_fig=data_for_fig[young_idx,]
	pdf(paste0('Heatmap/Young_full_length_LINE1_clustering_distance_euclidean_merged_cell_lines_add_subfamily_anno_rm_sense_LINE1_',norm_method[n],'_10celllines_darkcolor.pdf'), width = 5, height = 5,useDingbats=FALSE)
	heatmap.2(data_for_fig,
			  scale = "none",  # You can choose "row", "column", or "none" for scaling
			  main = paste0("Young LINE1 (",norm_method[n],")"),
			  distfun = function(x) dist(x, method = "euclidean"),  # Use correlation as the distance metric
			  hclustfun = function(x) hclust(x, method = "complete"),  # Specify clustering method
			  margins = c(7, 7),  # Add extra space for labels
			  #col = rev(colorRampPalette(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"))(50)),
			  col = rev(colorRampPalette(c("#67001f", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
			  #col = rev(colorRampPalette(c("#de4128", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
			  key = TRUE,  # Display the color legend
			  trace = "none",cexRow = 0.1,cexCol=0.5,RowSideColors=subfamily_anno[,3])
	dev.off()
	####only draw heatmap of L1HS,L1PA2 and L1PA3.
	
	subfamily_anno[subfamily_anno[,3]=='white',2]='Other LINE1'
	line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5','Other LINE1')
	pdf(paste0('Heatmap/line1_subfamily_clustering_distance_euclidean_merged_cell_lines_add_subfamily_anno_rm_sense_LINE1_',norm_method[n],'_10celllines.pdf'), width = 5, height = 5, onefile = TRUE,useDingbats=FALSE)
	for (i in 1:length(line1_subfamily_need_color)){
		idx <- subfamily_anno[,2]==line1_subfamily_need_color[i]
		heatmap.2(data_for_fig[idx,c(19,20,9,10,5,6,1,2,15,16,3,4,7,8,17,18,13,14,11,12)],
				  scale = "none",  # You can choose "row", "column", or "none" for scaling
				  Colv=FALSE,
				  main = paste0("all ",line1_subfamily_need_color[i],"(",norm_method[n],")"),
				  distfun = function(x) dist(x, method = "euclidean"),  # Use correlation as the distance metric
				  hclustfun = function(x) hclust(x, method = "complete"),  # Specify clustering method
				  margins = c(7, 7),  # Add extra space for labels
				  col = rev(colorRampPalette(c("#67001f", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
				  #col = rev(colorRampPalette(c("#de4128", "#d6604d", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de"))(50)),
				  key = TRUE,  # Display the color legend
				  trace = "none",cexRow = 0.1,cexCol=0.5,RowSideColors=subfamily_anno[idx,3])
	}
	dev.off()
	###draw heatmap of each LINE1 subfamily.
}

module unload bedtools ## bedtools seems to have conflict with ggplot2.
module purge && module load R/4.2.2-shlib
R

new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
.libPaths(c(new_lib_path, .libPaths()))
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))
library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")


setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
norm_method=c('TPM')
exp_cutoff=list(c(0.5))
pdf(paste0('Subfamily_analysis/line1_subfamily_fraction_of_expressed_LINE1_merged_cell_lines_rm_sense_LINE1_10celllines.pdf'), width = 7, height = 4,onefile = TRUE,useDingbats=FALSE)
for (n in 1:length(norm_method)){
	for (e in 1:length(exp_cutoff[[n]])){
		normed_line1_exp <- read.table(paste0('Exp_matrix/line1_',norm_method[n],'_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt'),header=TRUE,sep="\t")[,-c(9,10,11,12)]
		normed_line1_exp_matrix <- as.matrix(normed_line1_exp)
		data_for_fig <- log2(normed_line1_exp_matrix+1)
		#data_for_fig[data_for_fig>8]=8
		line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(4,6)]
		subfamily_anno <- merge(normed_line1_exp_matrix,line1,by.x=0,by.y='V4')
		subfamily_anno <- subfamily_anno[,c(1,dim(subfamily_anno)[2])]
		line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5')
		#subfamily_color <- c('black','#33ABC1','#8D73BA','red','grey')
		subfamily_color <- c('red','green','blue','yellow','grey')
		subfamily_anno[,3] <- 'white'
		for (i in 1:length(line1_subfamily_need_color)){
			subfamily_anno[subfamily_anno[,2]==line1_subfamily_need_color[i],3] <- subfamily_color[i]
		}
		subfamily_anno[subfamily_anno[,3]=='white',2]='Other LINE1'
		fraction_of_expressed_line1=NULL
		cellline=gsub('_rep1','',colnames(data_for_fig)[seq(1,dim(data_for_fig)[2],2)])
		line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5','Other LINE1')
		for (i in 1:length(line1_subfamily_need_color)){
			idx <- subfamily_anno[,2]==line1_subfamily_need_color[i]
			data_for_fig_temp <- data_for_fig[idx,]
			for (j in 1:length(cellline)){
				data_for_fig_temp_spe <- data_for_fig_temp[,(2*j-1):(2*j)]
				fraction <- c(line1_subfamily_need_color[i],cellline[j],round(sum(apply(data_for_fig_temp_spe,1,mean)>log2(exp_cutoff[[n]][e]+1))/dim(data_for_fig_temp)[1],4))
				fraction_of_expressed_line1 <- rbind(fraction_of_expressed_line1,fraction)
			}
		}
		colnames(fraction_of_expressed_line1)=c('LINE1_subfamily','Cell_Line','Fraction_of_expressed_LINE1')
		fraction_of_expressed_line1 <- as.data.frame(fraction_of_expressed_line1)
		fraction_of_expressed_line1[,3] <- as.numeric(fraction_of_expressed_line1[,3])
		p <- ggplot(data=fraction_of_expressed_line1, aes(x=LINE1_subfamily, y=Fraction_of_expressed_LINE1, group=Cell_Line,color=Cell_Line)) + geom_line() + ggtitle(paste0('Normalized_by_',norm_method[n],'_exp_cutoff',exp_cutoff[[n]][e])) + ylim(0, 1)
		final_plot <- p+scale_color_manual(values=c('#699ECA','#FF8C00','#F898CB','#4DAF4A','#D65190','#731A73','#FFCB5B','#EC3E31','#0076B9','#3D505A'))
		print(final_plot)
		write.table(fraction_of_expressed_line1,'Subfamily_analysis/line1_subfamily_fraction_of_expressed_LINE1_merged_cell_lines_rm_sense_LINE1_10celllines.txt',row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
	}
}
dev.off()
##############draw the fraction of expressed LINE1 for each subfamily.
###############draw heatmap of all LINE1 and for each LINE1 family (exclude sense LINE1).