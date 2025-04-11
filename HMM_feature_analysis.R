###prepare cellmarkfiletable.txt file.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/ChromHMM')
cellline=c('MCF7','K562','MOLM13','CaCo2','T47D','Huh7','GM12878','THP1','MDA-MB-231','H460')
mark_name=list(c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H4K20me1','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H4K20me1','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3'),
c('H3K4me1','H3K4me3','H3K9ac','H3K27ac','H3K27me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H4K20me1','H3K27ac','H3K27me3','H3K36me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3'),
c('CTCF','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3'),
c('H3K4me3','H3K27ac','H3K27me3'))
for (i in 1:length(cellline)){
	for (j in 1:length(mark_name[[i]])){
		print(file.exists(paste0('Input_bedfiles/',cellline[i],'_',mark_name[[i]][j],'_merged.bed')))
	}
}	
output=NULL
for (i in 1:length(cellline)){
	output_temp=matrix(NA,length(mark_name[[i]]),3)
	output_temp[,1]=cellline[i]
	for (j in 1:length(mark_name[[i]])){
		output_temp[j,2]=mark_name[[i]][j]
		output_temp[j,3]=paste0(cellline[i],'_',mark_name[[i]][j],'_merged.bed')
	}
	output=rbind(output,output_temp)
}
write.table(output,'cellmarkfiletable.txt',row.names=F,col.names=F,sep="\t",quote=F)

module load java
work_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/ChromHMM
chromsize_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/chromSize
chromHMM_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Software/ChromHMM
cd $work_dir
mkdir output_15states
java -mx1600M -jar $chromHMM_dir/ChromHMM.jar BinarizeBed $chromsize_dir/hg38.chromSize.clean $work_dir/Input_bedfiles $work_dir/cellmarkfiletable.txt $work_dir/binarized
java -mx1600M -jar $chromHMM_dir/ChromHMM.jar LearnModel -noimage $work_dir/binarized $work_dir/output_15states 15 hg38
### annotate the output of ChromHMM manually and classify them into seven groups.

### 7/30/2024
###############HMM features of expressed LINE1 vs non-expressed LINE1.
module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0

library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")

cellline=c('MCF7','K562','MOLM13','CaCo2','T47D','Huh7','GM12878','THP1','MDA-MB-231','H460')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/ChromHMM/output_15states')
Emission=paste0('E',1:15)
label=c('Heterochrom','Promoter','Promoter','Promoter','Enhancer','Enhancer','Repetitive_or_CNV','Heterochrom','Heterochrom','Polycomb_repressed','Insulator','Enhancer','Transcription','Transcription','Transcription')
anno=as.data.frame(cbind(Emission,label))
chromHMM_anno_list=list()
for (i in 1:length(cellline)){
	dat=read.table(paste0(cellline[i],'_15_segments.bed'),header=F,sep="\t")
	anno_dat=merge(dat,anno,by.x='V4',by.y='Emission')[,2:5]	
	chromHMM_anno_list[[i]]=anno_dat
}

line1_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt",header=TRUE,sep="\t")[,-c(9,10,11,12)]
line1_counts=log2(line1_counts+1)

line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.5UTR-w-genomic-flanks_del_L1PA2_1008.bed',header=FALSE,sep="\t")
line1_plus=line1[line1[,6]=='+',]
line1_minus=line1[line1[,6]=='-',]
line1_plus[,2]=line1_plus[,2]-1000
line1_minus[,3]=line1_minus[,3]+1000
line1_5UTR_plus_promoter=rbind(line1_plus,line1_minus)
line1_full=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=FALSE,sep="\t")

cellline=c('MCF7','K562','MOLM13','CaCo2','T47D','Huh7','GM12878','THP1','MDA','H460')
exp_cutoff=log2(0.6+1)
data_for_fig_all=NULL
for (i in 1:length(cellline)){
	line1_counts_spe=line1_counts[,grepl(cellline[i],colnames(line1_counts))]  
	exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]>exp_cutoff)&(line1_counts_spe[,2]>exp_cutoff)]
	exp_line1_promoter_info=line1_5UTR_plus_promoter[line1_5UTR_plus_promoter[,4]%in%exp_line1,]
	exp_line1_full_info=line1_full[line1_full[,4]%in%exp_line1,]
	no_exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]<exp_cutoff)&(line1_counts_spe[,2]<exp_cutoff)]
	no_exp_line1_promoter_info=line1_5UTR_plus_promoter[line1_5UTR_plus_promoter[,4]%in%no_exp_line1,]
	no_exp_line1_full_info=line1_full[line1_full[,4]%in%no_exp_line1,]

	uni_anno=unique(chromHMM_anno_list[[i]][,4])
	data_for_fig=matrix(0,length(uni_anno),6)
	for (j in 1:length(uni_anno)){
		uni_anno_spe=chromHMM_anno_list[[i]][chromHMM_anno_list[[i]][,4]==uni_anno[j],]
		if (length(grep('Promoter',uni_anno[j]))!=0){
			exp_line1_olp_HMM=bedtoolsr::bt.intersect(exp_line1_promoter_info,uni_anno_spe,wa=T,wb=T)		
			data_for_fig[j,2]=length(unique(exp_line1_promoter_info[,4]))
			if (dim(exp_line1_olp_HMM)[1]!=0){
                data_for_fig[j,1]=length(unique(exp_line1_olp_HMM[,4]))
				data_for_fig[j,3]=round(length(unique(exp_line1_olp_HMM[,4]))/length(unique(exp_line1_promoter_info[,4])),4)
			}
			no_exp_line1_olp_HMM=bedtoolsr::bt.intersect(no_exp_line1_promoter_info,uni_anno_spe,wa=T,wb=T)
            data_for_fig[j,4:5]=c(length(unique(no_exp_line1_olp_HMM[,4])),length(unique(no_exp_line1_promoter_info[,4])))
			data_for_fig[j,6]=round(length(unique(no_exp_line1_olp_HMM[,4]))/length(unique(no_exp_line1_promoter_info[,4])),4)
		}
		else if (length(grep('Heterochrom',uni_anno[j]))!=0) {
			exp_line1_olp_HMM=bedtoolsr::bt.intersect(exp_line1_full_info,uni_anno_spe,wa=T,wb=T,f=1)
			data_for_fig[j,2]=length(unique(exp_line1_full_info[,4]))
			if (dim(exp_line1_olp_HMM)[1]!=0){
				data_for_fig[j,1]=length(unique(exp_line1_olp_HMM[,4]))
				data_for_fig[j,3]=round(length(unique(exp_line1_olp_HMM[,4]))/length(unique(exp_line1_full_info[,4])),4)
			}
			no_exp_line1_olp_HMM=bedtoolsr::bt.intersect(no_exp_line1_full_info,uni_anno_spe,wa=T,wb=T,f=1)
            data_for_fig[j,4:5]=c(length(unique(no_exp_line1_olp_HMM[,4])),length(unique(no_exp_line1_full_info[,4])))
			data_for_fig[j,6]=round(length(unique(no_exp_line1_olp_HMM[,4]))/length(unique(no_exp_line1_full_info[,4])),4)
		} else {
			exp_line1_olp_HMM=bedtoolsr::bt.intersect(exp_line1_full_info,uni_anno_spe,wa=T,wb=T)
			data_for_fig[j,2]=length(unique(exp_line1_full_info[,4]))
			if (dim(exp_line1_olp_HMM)[1]!=0){
                data_for_fig[j,1]=length(unique(exp_line1_olp_HMM[,4]))
				data_for_fig[j,3]=round(length(unique(exp_line1_olp_HMM[,4]))/length(unique(exp_line1_full_info[,4])),4)
			}
			no_exp_line1_olp_HMM=bedtoolsr::bt.intersect(no_exp_line1_full_info,uni_anno_spe,wa=T,wb=T)
            data_for_fig[j,4:5]=c(length(unique(no_exp_line1_olp_HMM[,4])),length(unique(no_exp_line1_full_info[,4])))
			data_for_fig[j,6]=round(length(unique(no_exp_line1_olp_HMM[,4]))/length(unique(no_exp_line1_full_info[,4])),4)
		}
	}
	data_for_fig=as.data.frame(data_for_fig)
	data_for_fig[,7]=data_for_fig[,3]/data_for_fig[,6]
	data_for_fig[,8]=uni_anno; 
	data_for_fig_all=rbind(data_for_fig_all,cbind(data_for_fig,cellline[i]))
}
L=dim(data_for_fig_all)[2]
data_for_fig_all[,L+1]=data_for_fig_all[,2]-data_for_fig_all[,1]
data_for_fig_all[,L+2]=data_for_fig_all[,5]-data_for_fig_all[,4]
colnames(data_for_fig_all)=c('expressed_olp_HMM','expressed_all','expressed_ratio','non-expressed_olp_HMM','non-expressed_all','non_expressed_ratio','FC','HMM_features','Cell_Line','expressed_no_olp_HMM','non-expressed_no_olp_HMM')
data_for_fisher=data_for_fig_all[,c(1,10,4,11)]
result <- apply(data_for_fisher, 1, function(x) fisher.test(matrix(x, 2, 2)))
L=dim(data_for_fig_all)[2]
for (i in 1:length(result)){
data_for_fig_all[i,L+1]=result[[i]]$p
}
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/HMM_feature_analysis')
write.table(data_for_fig_all,'expressed_LINE1_olp_HMM_features_10celllines.txt',row.names=F,col.names=T,sep="\t",quote=F)

colnames(data_for_fig_all)=c(colnames(data_for_fig_all)[1:11],'P_value')
#data_for_fig_all$log2_FC <- log2(data_for_fig_all$FC)

# Ensure 'Cell_Line' and 'HMM_features' are factors
data_for_fig_all[data_for_fig_all[,12]<1e-5,12]=1e-5
data_for_fig_all[,12]=(-log10(data_for_fig_all[,12]))
data_for_fig_all[data_for_fig_all[,7]<1,12]= (-data_for_fig_all[data_for_fig_all[,7]<1,12])
data_for_fig_all$Cell_Line <- as.factor(data_for_fig_all$Cell_Line)
data_for_fig_all$HMM_features <- as.factor(data_for_fig_all$HMM_features)

# Create a matrix of P-values with 'Cell_Line' as rows and 'HMM_features' as columns
pvalue_matrix <- with(data_for_fig_all, tapply(P_value, list(Cell_Line, HMM_features), mean))
pvalue_matrix <- pvalue_matrix[,-8]

# Create a heatmap
library(gplots,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2')
pdf('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/HMM_feature_analysis/expressed_LINE1_HMM_features_heatmap_10celllines_TPM_cutoff_0point6.pdf',width=10,height=7)
heatmap.2(pvalue_matrix, 
        scale = "none",  # You can choose "row", "column", or "none" for scaling
        main = "Heatmap of P-values",
        margins = c(12, 7),  # Add extra space for labels
        col = colorRampPalette(c("blue", "white", "red"))(50),   # You can customize the color scale
		key = TRUE,  # Display the color legend
        key.title = "-log10(P-value)",
        trace="none")  # Set the legend title
dev.off()
####The heatmap looks good.  The expressed LINE1 are enriched with promoter mark and depleted in Heterochrmatin mark.
###############HMM features of expressed LINE1 vs non-expressed LINE1.


###8/27/2024
###############Test whether the interacting regions of expressed LINE1 are enriched with enhancer mark etc.
module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0

library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")

cellline=c('MCF7','K562','T47D')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/ChromHMM/output_15states')
Emission=paste0('E',1:15)
label=c('Heterochrom','Promoter','Promoter','Promoter','Enhancer','Enhancer','Repetitive_or_CNV','Heterochrom','Heterochrom','Polycomb_repressed','Insulator','Enhancer','Transcription','Transcription','Transcription')
anno=as.data.frame(cbind(Emission,label))
chromHMM_anno_list=list()
for (i in 1:length(cellline)){
	dat=read.table(paste0(cellline[i],'_15_segments.bed'),header=F,sep="\t")
	anno_dat=merge(dat,anno,by.x='V4',by.y='Emission')[,2:5]	
	chromHMM_anno_list[[i]]=anno_dat
}

line1_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt",header=TRUE,sep="\t")
line1_counts=log2(line1_counts+1)

cellline=c('MCF7','K562','T47D')
PETs_name=c('MCF7_merged','K562_merged','T47D_081423')
col_list=list(c(13,14),c(7,8),c(19,20))
exp_cutoff=0.5
data_for_fig_all=NULL

library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
library(cowplot,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
uni_anno=unique(chromHMM_anno_list[[1]][,4])
p_value=matrix(1,length(PETs_name),length(uni_anno),dimnames=list(cellline,uni_anno))
for (i in 1:length(PETs_name)){
	line1_counts_spe=line1_counts[,col_list[[i]]]  
	exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]>=exp_cutoff)&(line1_counts_spe[,2]>=exp_cutoff)]
	no_exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]<exp_cutoff)|(line1_counts_spe[,2]<exp_cutoff)]
	PETs=read.table(paste0('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance/',PETs_name[i],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_sig_other_end_intra.bed'),header=F,sep="\t")
	exp_line1_with_PETs=intersect(exp_line1,as.vector(PETs[,4]))  #746/975
	no_exp_line1_with_PETs=intersect(no_exp_line1,as.vector(PETs[,4])) #3456/8030

	uni_anno=unique(chromHMM_anno_list[[i]][,4])
	exp_line1_PET_num=matrix(0,length(exp_line1_with_PETs),length(uni_anno))
	no_exp_line1_PET_num=matrix(0,length(no_exp_line1_with_PETs),length(uni_anno))
	for (j in 1:length(uni_anno)){
		uni_anno_spe=chromHMM_anno_list[[i]][chromHMM_anno_list[[i]][,4]==uni_anno[j],]
		region_size=sum(uni_anno_spe[,3]-uni_anno_spe[,2])
		PETs_olp_HMM=bedtoolsr::bt.intersect(PETs,uni_anno_spe,wa=T,wb=T)
		for (k in 1:length(exp_line1_with_PETs)){
			if (sum(PETs_olp_HMM[,4]==exp_line1_with_PETs[k])>0){
				PET_num_spe=dim(unique(PETs_olp_HMM[PETs_olp_HMM[,4]==exp_line1_with_PETs[k],1:3]))[1]
				all_PET_num_spe=dim(unique(PETs[PETs[,4]==exp_line1_with_PETs[k],1:3]))[1]
				exp_line1_PET_num[k,j]=(PET_num_spe/region_size)/(all_PET_num_spe/3031042417)
			}
		}
		
		for (k in 1:length(no_exp_line1_with_PETs)){
			if (sum(PETs_olp_HMM[,4]==no_exp_line1_with_PETs[k])>0){
				PET_num_spe=dim(unique(PETs_olp_HMM[PETs_olp_HMM[,4]==no_exp_line1_with_PETs[k],1:3]))[1]
				all_PET_num_spe=dim(unique(PETs[PETs[,4]==no_exp_line1_with_PETs[k],1:3]))[1]
				no_exp_line1_PET_num[k,j]=(PET_num_spe/region_size)/(all_PET_num_spe/3031042417)
			}
		}
	
	}
	
	colnames(exp_line1_PET_num)=uni_anno
	colnames(no_exp_line1_PET_num)=uni_anno
	rownames(exp_line1_PET_num)=exp_line1_with_PETs
	rownames(no_exp_line1_PET_num)=no_exp_line1_with_PETs
	#colSums(exp_line1_PET_num>0)/dim(exp_line1_PET_num)[1]
	#colSums(no_exp_line1_PET_num>0)/dim(no_exp_line1_PET_num)[1]

	
	for (j in 1:length(uni_anno)){
		t_test_result=t.test(exp_line1_PET_num[,j],no_exp_line1_PET_num[,j])
		if (t_test_result$statistic>0){
			p_value[i,j]=(-log10(t_test_result$p.value))
		}else{
			p_value[i,j]=(log10(t_test_result$p.value))
		}
	}

	data_for_fig_list=list()
	for (j in 1:length(uni_anno)){
		temp=cbind(exp_line1_PET_num[,j],1)
		colnames(temp)=c('PET_num','Exp_label')
		olp=merge(temp,line1_counts_spe,by.x=0,by.y=0)
		olp[dim(olp)[2]+1]=(olp[dim(olp)[2]-1]+olp[dim(olp)[2]])/2
		rownames(olp)=olp[,1]
		olp=olp[,-1]

		sorted_score=sort(olp[,dim(olp)[2]],decreasing=T)
		L=length(sorted_score)
		olp[(olp[,dim(olp)[2]]>=sorted_score[round(L*0.50)]),2]=2
		olp[(olp[,dim(olp)[2]]>=sorted_score[round(L*0.25)]),2]=3
		olp[(olp[,dim(olp)[2]]>=sorted_score[round(L*0.1)]),2]=4
		olp[(olp[,dim(olp)[2]]>=sorted_score[round(L*0.05)]),2]=5
		
		temp=cbind(no_exp_line1_PET_num[,j],0); colnames(temp)=c('PET_num','Exp_label')
		data_for_fig_list[[j]]=rbind(temp,olp[,1:2])
	}
	
}

save(list = ls(), file = "expressed_LINE1_interacting_PETs_HMM_features_all_variables.RData")
library(gplots,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/')
p_value[p_value>5]=5; p_value[p_value<(-5)]=(-5)
pdf('expressed_LINE1_interacting_PETs_HMM_features_heatmap_sig_PETs_chromatin_RNA_seq_10celllines.pdf',width=10,height=5)
heatmap.2(p_value, 
        scale = "none",  # You can choose "row", "column", or "none" for scaling
        main = "Heatmap of P-values",
        margins = c(12, 7),  # Add extra space for labels
        col = colorRampPalette(c("blue", "white", "red"))(50),   # You can customize the color scale
		key = TRUE,  # Display the color legend
        key.title = "-log10(P-value)",
        trace="none",
		cexRow = 2)  # Set the legend title
dev.off()
###################8/29/2024
###############Test whether the interacting regions of expressed LINE1 are enriched with enhancer mark etc.
