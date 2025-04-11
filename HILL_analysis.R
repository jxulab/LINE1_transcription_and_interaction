##########Define highly interacted LINE1 loci (Hills) and test if they tend to have higher expression. 
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
name=c('MCF7_merged','K562_merged','T47D_081423')
K <- 5000
B <- 50
p_cutoff=0.0001
for (i in 1:length(name)){
	all_pair_intra=read.table(paste0(name[i], '_', B, 'bins_binsize_', K, '_line1_intra_pairs_p_value.txt'),header=F,sep="\t")
	all_pair_inter=read.table(paste0(name[i],'_binsize_1M_line1_inter_pairs_p_value.txt'),header=F,sep="\t")
	sig_pair_intra=all_pair_intra[all_pair_intra[,5]<p_cutoff,]
	sig_pair_inter=all_pair_inter[all_pair_inter[,5]<p_cutoff,]
	colnames(sig_pair_intra)=colnames(sig_pair_inter)
	sig_pair=rbind(sig_pair_intra,sig_pair_inter)

	line1_freq=table(sig_pair[,6])
	sig_pair_output=merge(sig_pair,line1_freq,by.x='V6',by.y=1)
	colnames(sig_pair_output)=c('LINE1_ID','chrom1','pos1','chrom2','pos2','p value','PET_num')
	write.table(sig_pair_output,paste0(name[i],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_line1_with_PET_num.txt'),row.names=F,col.names=T,sep='\t',quote=F)
}

library(ggplot2)
library("cowplot")
cellline=c('MCF7','K562','T47D')
samples=c('MCF7_merged','K562_merged','T47D_081423')
line1_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt",header=TRUE,sep="\t")
line1_counts=log2(line1_counts+1)
col_list=list(c(13,14),c(7,8),c(19,20))
hill_list=list()
for (i in 1:length(cellline)){
  setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
  sig_pair=unique(read.table(paste0(samples[i],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_line1_with_PET_num.txt'),header=T,sep="\t")[,1:6])
  line1_PET_num=sort(table(sig_pair[,1]))
  x <- 1:length(line1_PET_num)
  y <- as.vector(line1_PET_num)

  scaled_x <- (x - min(x)) / (max(x) - min(x))
  scaled_y <- (y - min(y)) / (max(y) - min(y))

  # Calculate the slope for each point
  slopes <- diff(scaled_y) / diff(scaled_x)

  # Find the point with a slope closest to 1
  closest_slope_index <- which.min(abs(slopes - 1))
  closest_point_x <- scaled_x[closest_slope_index + 1]
  closest_point_y <- scaled_y[closest_slope_index + 1]

  # Map the index from scaled data to the original data
  original_x_value <- x[closest_slope_index + 1]
  original_y_value <- y[closest_slope_index + 1]

  pdf('All_full_length_LINE1_PET_num_scaled_curve_p0point0001.pdf')
  # Plot the scaled data
  # Calculate the tangent line
  tangent_slope <- 1
  tangent_intercept <- scaled_y[closest_slope_index] - tangent_slope * scaled_x[closest_slope_index]

  # Plot the scaled data
  plot(scaled_x, scaled_y, pch = 20, col = 'red', bg = 'pink', lwd = 1, xlab = 'Scaled X', ylab = 'Scaled Y')

  # Highlight the point with slope closest to 1
  points(scaled_x[closest_slope_index], scaled_y[closest_slope_index], pch = 19, col = 'green')

  # Plot the tangent line
  abline(a = tangent_intercept, b = tangent_slope, col = 'blue', lty = 2)
  dev.off()

  hil=names(line1_PET_num)[original_x_value:length(line1_PET_num)]
  hill_list[[i]]=hil
  write.table(hil,paste0(samples[i],'HILL_line1.txt'),row.names=F,col.names=F,sep="\t",quote=F)
  ############# Define highly interacted LINE1 using super-enhancer similar method.


  #############check wether Hil LINE1 has higher LINE1 expression.
  line1_counts_spe=line1_counts[,col_list[[i]]] 

  line1_counts_spe[rownames(line1_counts_spe)%in%hil,3]='Hill'
  line1_counts_spe[!(rownames(line1_counts_spe)%in%hil),3]='non-Hill'
  line1_counts_spe[,4]=(line1_counts_spe[,1]+line1_counts_spe[,2])/2
  colnames(line1_counts_spe)=c('Rep1','Rep2','Hil_label','Mean_of_2Rep')
  t_test_result=t.test(line1_counts_spe[line1_counts_spe[,3]=='Hill',4],line1_counts_spe[line1_counts_spe[,3]=='non-Hill',4])
  p_value=t_test_result$p.value

  # Create a data frame with your scaled data
  data <- data.frame(scaled_x, scaled_y)

  # Create the plot
  pdf(paste0(cellline[i],'_HILL_LINE1_exp_p0point0001.pdf'),useDingbats=FALSE)
  a <- ggplot(data, aes(x = scaled_x, y = scaled_y)) +
    geom_point(shape = 21, color = "red", fill = "pink", size = 3) +
    geom_point(data = data[closest_slope_index, , drop = FALSE], shape = 19, color = "green", size = 3) +
    geom_abline(intercept = tangent_intercept, slope = tangent_slope, linetype = "dashed", color = "blue") +
    labs(x = "Scaled X", y = "Scaled Y") + theme_minimal() + labs(x="LINE1",y="scaled PET number",title=cellline[i])
  b=ggplot(line1_counts_spe, aes(x=Hil_label,  y=Mean_of_2Rep,color=Hil_label)) + geom_boxplot() + scale_color_manual(values=c("red","grey")) + ylim(0,2) + labs(x="",y="log2(TPM)",title=cellline[i]) + annotate("text",x=1.5,y = 2, label = paste0('p value = ',signif(p_value,digits=2)))
  #b=ggplot(line1_counts_spe, aes(x=Hil_label,  y=Mean_of_2Rep,color=Hil_label)) + geom_boxplot() + scale_color_manual(values=c("red","grey")) + ylim(0,1) + labs(x="",y="log2(TPM)",title=cellline[i]) + annotate("text",x=1.5,y = 1, label = paste0('p value = ',signif(p_value,digits=2)))
  final_plot=plot_grid(a,b,labels = c("A", "B"),ncol = 2, nrow = 2)
  print(final_plot)
  dev.off()
}
#### HILLs has higher LINE1 expression.


#############Test whether HILL-interacting genes has higher expression.  
module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0
library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
library(cowplot,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
genes=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/genes.bed',header=F,sep="\t")
genes_plus=genes[genes[,6]=="+",]
genes_plus[,2]=genes_plus[,2]-4000
genes_minus=genes[genes[,6]=="-",]
genes_minus[,3]=genes_minus[,3]+4000
genes=rbind(genes_plus,genes_minus)
genes[genes[,2]<0,2]=1

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
PETs_name=c('MCF7_merged','K562_merged','T47D_081423')
cellline=c('MCF7','K562','T47D')
gene_PET_num_list=list()
for (c in 1:length(PETs_name)){
  HILL_line1=read.table(paste0(PETs_name[c],'HILL_line1_intra_chrom.txt'),header=F,sep="\t")
  our_pair=unique(read.table(paste0(PETs_name[c],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_line1_with_PET_num.txt'),header=T,sep="\t")[,c(2:5)])
  our_pair_left=our_pair[,1:2]; our_pair_left[,3]=our_pair_left[,2]+1
  our_pair_right=our_pair[,3:4]; our_pair_right[,3]=our_pair_right[,2]+1
  our_pair_left[,4]=1:dim(our_pair)[1]
  our_pair_right[,4]=1:dim(our_pair)[1]

  HILL_line1_info=line1_info[line1_info[,4]%in%HILL_line1[,1],1:4]

  our_pair_left_olp_line1=bedtoolsr::bt.intersect(our_pair_left,HILL_line1_info,wa=T,wb=T)

  gene_bed=genes[,1:4]
  our_pair_right_olp_gene=bedtoolsr::bt.intersect(our_pair_right,gene_bed,wa=T,wb=T)
  genes_pairs_with_HILL_line1=our_pair_right_olp_gene[our_pair_right_olp_gene[,4]%in%our_pair_left_olp_line1[,4],c(1,2,3,4,8)]
  gene_PET_num_with_HILL_line1=table(genes_pairs_with_HILL_line1[,5])
  gene_PET_num_list[[c]]=gene_PET_num_with_HILL_line1
  HILL_line1_gene_pair=merge(our_pair_left_olp_line1,our_pair_right_olp_gene,by.x='V4',by.y='V4')[,c(2,3,8,9,10,15)]
  colnames(HILL_line1_gene_pair)=c('Left_chrom','Left_pos','LINE1','Right_chrom','Right_pos','Gene')
  write.table(HILL_line1_gene_pair,paste0(cellline[c],'_HILL_line1_gene_pair.txt'),row.names=F,col.names=T,sep="\t",quote=F)
}
save(list = ls(), file = "HILL_line1_interacting_genes.RData")


setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
load("HILL_line1_interacting_genes.RData")

new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
.libPaths(c(new_lib_path, .libPaths()))
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))
library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
library("cowplot")
cellline=c('MCF7','K562','T47D')
gene_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/gene_TPM_normalized_readcount_merged_cell_lines.txt",header=TRUE,sep="\t")
gene_counts=log2(gene_counts+1)
col_list=list(c(13,14),c(7,8),c(19,20))
gene_counts_spe_list=list()
p_value_list=list()
for (i in 1:length(cellline)){
  gene_counts_spe=gene_counts[,col_list[[i]]] 

  gene_counts_spe[rownames(gene_counts_spe)%in%names(gene_PET_num_list[[i]]),3]='Hill'
  gene_counts_spe[!(rownames(gene_counts_spe)%in%names(gene_PET_num_list[[i]])),3]='non-Hill'
  gene_counts_spe[,4]=(gene_counts_spe[,1]+gene_counts_spe[,2])/2
  colnames(gene_counts_spe)=c('Rep1','Rep2','Hil_label','Mean_of_2Rep')
  t_test_result=t.test(gene_counts_spe[gene_counts_spe[,3]=='Hill',4],gene_counts_spe[gene_counts_spe[,3]=='non-Hill',4])
  p_value=t_test_result$p.value

  gene_counts_spe_list[[i]]=gene_counts_spe
  p_value_list[[i]]=p_value
}

pdf('3celllines_HILL_interacting_gene_exp_p0point0001_intra_chrom_no_xlimit.pdf',useDingbats=FALSE, height = 6)
a=ggplot(gene_counts_spe_list[[1]], aes(x=Hil_label,  y=Mean_of_2Rep,color=Hil_label)) + geom_boxplot() + scale_color_manual(values=c("red","grey")) +  labs(x="",y="log2(TPM)",title=cellline[1]) + annotate("text",x=1.5,y = 15, label = paste0('p value = ',signif(p_value_list[[1]],digits=2)))
b=ggplot(gene_counts_spe_list[[2]], aes(x=Hil_label,  y=Mean_of_2Rep,color=Hil_label)) + geom_boxplot() + scale_color_manual(values=c("red","grey")) +  labs(x="",y="log2(TPM)",title=cellline[2]) + annotate("text",x=1.5,y = 15, label = paste0('p value = ',signif(p_value_list[[2]],digits=2)))
c=ggplot(gene_counts_spe_list[[3]], aes(x=Hil_label,  y=Mean_of_2Rep,color=Hil_label)) + geom_boxplot() + scale_color_manual(values=c("red","grey")) +  labs(x="",y="log2(TPM)",title=cellline[3]) + annotate("text",x=1.5,y = 15, label = paste0('p value = ',signif(p_value_list[[3]],digits=2)))
final_plot=plot_grid(a,b,c,labels = c("A","B","C"),ncol = 2, nrow = 2)
print(final_plot)
dev.off()
########### HILL-interacting genes do have higher expression.


######Check what pathways does HILL-interacting genes enriched for.
library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
load("HILL_line1_interacting_genes.RData")
name=c('MCF7','K562','T47D')
for (i in 1:length(name)){
  write.table(gene_PET_num_list[[i]],paste0(name[i],'_HILL_interacting_genes.txt'),row.names=F,col.names=T,sep="\t",quote=F)
}

gene_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/genes.bed',header=F,sep="\t")
anno=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Data_Genomics/RNA-seq/Erythroid_Metabolome/CPDB_gene2path.txt',header=T,sep="\t")
uni_path=as.vector(unique(anno[,2]))
all_genes=unique(anno[,1])
N=length(all_genes)
enrich_path=matrix(NA,length(uni_path),7)
enrich_path[,1]=uni_path
for (i in 1:length(name)){
  k=NULL; n=NULL; M=NULL;
  gene_in_cluster=intersect(all_genes,names(gene_PET_num_list[[i]]))
  for (w in 1:dim(enrich_path)[1]){
    gene_in_path=intersect(unique(as.vector(anno[anno[,2]==enrich_path[w,1],1])),all_genes)
    M[w]=length(gene_in_path)
    n[w]=length(gene_in_cluster)
    k[w]=length(intersect(gene_in_path,gene_in_cluster))
  }                                               
  p=1-phyper(k-1,M,N-M,n)
  ratio=k/M
  enrich_path[,2:7]=cbind(N,M,n,k,p,ratio)
  enrich_path=enrich_path[order(p),]
  colnames(enrich_path)=c('Path_class','#Gene_in_all_path','#Gene_in_Path','#Gene_in_Cluster','#Overlap','P_value','Percentage')
  write.table(enrich_path,paste0(name[i],'_HILL_interacting_genes_enrich_CPDB_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
}
######Check what pathways does HILL-interacting genes enriched for.





