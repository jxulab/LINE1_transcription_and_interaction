##7/11/2024
###do motif analysis for expressed LINE1 versus non-expressed LINE1.
line1_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt",header=TRUE,sep="\t")
line1_counts=log2(line1_counts+1)

line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.5UTR-w-genomic-flanks_del_L1PA2_1008.bed',header=F,sep="\t")
cellline=c('MCF7','K562','MOLM13','CaCo2','T47D','Huh7','GM12878','THP1','MDA','H460')
exp_cutoff=log2(0.6+1)
data_for_fig_all=NULL
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/motif_analysis')
for (i in 1:length(cellline)){
	line1_counts_spe=line1_counts[,grepl(cellline[i],colnames(line1_counts))]  #
	exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]>exp_cutoff)&(line1_counts_spe[,2]>exp_cutoff)]
	exp_line1_info=line1_info[line1_info[,4]%in%exp_line1,]

	no_exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]<exp_cutoff)&(line1_counts_spe[,2]<exp_cutoff)]
	no_exp_line1_info=line1_info[line1_info[,4]%in%no_exp_line1,]

	line1_counts_spe[,3]=(line1_counts_spe[,1]+line1_counts_spe[,2])/2
	# sorted_line1=rownames(line1_counts_spe)[order(line1_counts_spe[,3])]
	# control_line1=sorted_line1[1:length(exp_line1)]
	# control_line1_info=line1_info[line1_info[,4]%in%control_line1,]
	control_line1=rownames(line1_counts_spe)[line1_counts_spe[,3]==0]
	control_line1_info=line1_info[line1_info[,4]%in%control_line1,]
	
	write.table(exp_line1_info,paste0(cellline[i],'_expressed_line1.bed'),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(no_exp_line1_info,paste0(cellline[i],'_no_expressed_line1.bed'),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(control_line1_info,paste0(cellline[i],'_control_line1.bed'),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

software_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/BisiMiao/Motif/bin
work_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/motif_analysis/expressed_vs_control
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/BisiMiao/Motif/bin
cellline=(MCF7 K562 MOLM13 CaCo2 T47D LNCaP MCF10A Huh7 GM12878 THP1 MDA H460)
for ((i=0;i<${#cellline[*]};i++))
do
mkdir ${cellline[i]}
$software_dir/findMotifsGenome.pl $work_dir/${cellline[i]}_expressed_line1.bed hg38 $work_dir/${cellline[i]} -size -1000,1000 -bg $work_dir/${cellline[i]}_control_line1.bed
done


#######################draw the motif enrichment result of expressed versus zero-expressed LINE1.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/motif_analysis/expressed_vs_control/size_plus1k_minus1k')
cellline=c('MCF7','K562','MOLM13','CaCo2','T47D','Huh7','GM12878','THP1','MDA','H460')
sig_motif_all=NULL
p_cutoff=0.001;  percent_diff_cutoff=5;  percent_cutoff=95; 
for (i in 1:length(cellline)){
	motif_result=read.table(paste0(cellline[i],'/knownResults.txt'),header=T,sep="\t",comment.char = "")
	motif_result[,7]=gsub("%","",motif_result[,7]); mode(motif_result[,7])="numeric"
	motif_result[,9]=gsub("%","",motif_result[,9]); mode(motif_result[,9])="numeric"
	motif_result[,10]=motif_result[,7]-motif_result[,9]
	
	sig_motif=motif_result[(motif_result[,5]<p_cutoff)&(motif_result[,10]>percent_diff_cutoff)&(motif_result[,7]>percent_cutoff),1]
	sig_motif_all=union(sig_motif_all,sig_motif)
	if (i==1){
		p_mtx=motif_result[,c(1,3)]
		percent_mtx=motif_result[,c(1,7)]
		percent_diff_mtx=motif_result[,c(1,10)]
	} 
	if (i>1){
		p_mtx=merge(p_mtx,motif_result[,c(1,3)],by.x='Motif.Name',by.y='Motif.Name')
		percent_mtx=merge(percent_mtx,motif_result[,c(1,7)],by.x='Motif.Name',by.y='Motif.Name')
		percent_diff_mtx=merge(percent_mtx,motif_result[,c(1,10)],by.x='Motif.Name',by.y='Motif.Name')
	}
}
sig_motif_p_mtx=p_mtx[p_mtx[,1]%in%sig_motif_all,]
sig_motif_percent_mtx=percent_mtx[percent_mtx[,1]%in%sig_motif_all,]
sig_motif_diff_percent_mtx=percent_diff_mtx[percent_diff_mtx[,1]%in%sig_motif_all,]
rownames(sig_motif_p_mtx)=gsub("/Homer","",sig_motif_p_mtx[,1]); sig_motif_p_mtx=sig_motif_p_mtx[,-1]
rownames(sig_motif_percent_mtx)=gsub("/Homer","",sig_motif_percent_mtx[,1]); sig_motif_percent_mtx=sig_motif_percent_mtx[,-1]
rownames(sig_motif_diff_percent_mtx)=gsub("/Homer","",sig_motif_diff_percent_mtx[,1]); sig_motif_diff_percent_mtx=sig_motif_diff_percent_mtx[,-1]
colnames(sig_motif_p_mtx)=cellline;  colnames(sig_motif_diff_percent_mtx)=cellline; colnames(sig_motif_percent_mtx)=cellline

# new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
# .libPaths(c(new_lib_path, .libPaths()))
# new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
# .libPaths(c(new_lib_path, .libPaths()))
# library(corrplot)
# library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
# library(Hmisc)
# library(RColorBrewer)

sig_motif_p_mtx=as.matrix(sig_motif_p_mtx)
mode(sig_motif_p_mtx)='numeric'
sig_motif_diff_percent_mtx=as.matrix(sig_motif_diff_percent_mtx)
mode(sig_motif_diff_percent_mtx)='numeric'
sig_motif_p_mtx[sig_motif_p_mtx<1e-22]=1e-22
order_flag=order(rowSums(-log10(sig_motif_p_mtx)),decreasing=T)
sig_motif_p_mtx=sig_motif_p_mtx[order_flag,c(10,4,6,5,1,9,2,7,3,8)]
sig_motif_diff_percent_mtx=sig_motif_diff_percent_mtx[order_flag,c(10,4,6,5,1,9,2,7,3,8)]
if ((dim(sig_motif_p_mtx)[1]>30)&(dim(sig_motif_p_mtx)[1]<100)){font_size=0.3}else if (dim(sig_motif_p_mtx)[1]>100){font_size=0.2}else{font_size=1}
pdf(paste0('expressed_vs_control_percent_cutoff',percent_cutoff,'_percent_diff_cutoff',percent_diff_cutoff,'_p_cutoff',gsub('\\.','point',p_cutoff),'_show_p_value_10celllines.pdf'), useDingbats=FALSE,width=10,height=7)
corrplot(-log10(sig_motif_p_mtx),type="full",is.corr = FALSE, tl.cex=font_size,tl.col="black", tl.srt=45,col=brewer.pal(n=9, name="OrRd")[3:8],p.mat=sig_motif_p_mtx, insig = "blank", sig.level = 0.01)
dev.off()
#use useDingbats=FALSE to avoid circle become square (crossed box) in illustrator.
#######################draw the motif enrichment result of expressed versus zero-expressed LINE1.


######1/13/2025
########I calculate the p value for the read distribution of FOXA1 binding for expressed and non-expression LINE1.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
normed_line1_exp <- read.table(paste0('Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt'),header=T,sep="\t")
normed_line1_exp_matrix <- as.matrix(normed_line1_exp)
normed_line1_exp_matrix <- log2(normed_line1_exp_matrix+1)
line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(4,6)]
subfamily_anno <- merge(normed_line1_exp_matrix,line1,by.x=0,by.y='V4')
subfamily_anno <- subfamily_anno[,c(1,dim(subfamily_anno)[2])]
line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5')
subfamily_color <- c('red','green','blue','yellow','grey')
subfamily_anno[,3] <- 'white'
for (i in 1:length(line1_subfamily_need_color)){
        subfamily_anno[subfamily_anno[,2]==line1_subfamily_need_color[i],3] <- subfamily_color[i]
}
subfamily_anno[subfamily_anno[,3]=='white',2]='Other_LINE1'


#####In order to compare the peak height of all individual regions using t test, I need to recalculate the RPKM of all regions.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
normed_line1_exp <- read.table(paste0('Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt'),header=T,sep="\t")
normed_line1_exp_matrix <- as.matrix(normed_line1_exp)
normed_line1_exp_matrix <- log2(normed_line1_exp_matrix+1)
line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(4,6)]
subfamily_anno <- merge(normed_line1_exp_matrix,line1,by.x=0,by.y='V4')
subfamily_anno <- subfamily_anno[,c(1,dim(subfamily_anno)[2])]
line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5')
subfamily_color <- c('red','green','blue','yellow','grey')
subfamily_anno[,3] <- 'white'
for (i in 1:length(line1_subfamily_need_color)){
        subfamily_anno[subfamily_anno[,2]==line1_subfamily_need_color[i],3] <- subfamily_color[i]
}
subfamily_anno[subfamily_anno[,3]=='white',2]='Other_LINE1'

bedfile=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/BrandonChen/Data_BrandonChen/Datasets_for_Michael/ChIP-seq/MCF7/FOXA1_rep1_MCF7/bedfiles/FOXA1_rep1_MCF7_k10000Aligned.sortedByCoord.out.primary.bed",header=F,sep="\t")
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(1:4,6,5)]
type=unique(subfamily_anno[,2])
work_dir="/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq"
p_value_ttest=NULL
for (i in 1:length(type)){
        exp_spe_type=normed_line1_exp_matrix[rownames(normed_line1_exp_matrix)%in%subfamily_anno[subfamily_anno[,2]==type[i],1],13:14]  ##MCF7
        exp_spe_type=as.data.frame(exp_spe_type)
        exp_spe_type[,3]=(exp_spe_type[,1]+exp_spe_type[,2])/2
        exp_spe_type=exp_spe_type[order(exp_spe_type[,3],decreasing=T),]
        regions=list()
        regions[[1]]=rownames(exp_spe_type)[exp_spe_type[,3]>log2(0.5+1)]
        regions[[2]]=rownames(exp_spe_type)[exp_spe_type[,3]<log2(0.5+1)]
        #regions[[2]]=rownames(exp_spe_type)[exp_spe_type[,3]==0]

        for (j in 1:length(regions)){
                regions[[j]]=line1_info[line1_info[,4]%in%regions[[j]],]
        }
        regions_list=regions

        data_for_fig=NULL
        RPKM_all_regions=NULL
        bin_num_region=40; bin_num_up=5; bin_num_down=5  # divide the region between TSS and TES into bin_num bins (default=2000).
        up=5000  # include the upstream region in the distribution (default=5000).
        down=5000 # include the downstream region in the distribution (default=5000).
        region_length_min=1000
        source('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/scripts/R_scripts/read_distribution_around_a_group_of_genes.r')
        for (j in 1:length(regions_list)){
                regions_plus=regions_list[[j]][regions_list[[j]][,6]=='+',1:3]
                regions_minus=regions_list[[j]][regions_list[[j]][,6]=='-',1:3]
                output_plus=read_distribution_around_a_group_of_genes(regions_plus,bedfile,bin_num_region,bin_num_up,bin_num_down,up,down,region_length_min)
                RPKM_plus=output_plus[[1]]; RPKM_all_regions_plus=output_plus[[2]]
                output_minus=read_distribution_around_a_group_of_genes(regions_minus,bedfile,bin_num_region,bin_num_up,bin_num_down,up,down,region_length_min)
                RPKM_minus=output_minus[[1]]; RPKM_all_regions_minus=output_minus[[2]]
                RPKM=RPKM_plus;  RPKM[,2]=(RPKM_plus[,2]+rev(RPKM_minus[,2]))/2
                RPKM_all_regions=cbind(RPKM_all_regions,rbind(cbind(RPKM_all_regions_plus,RPKM_all_regions_minus),j))
                data_for_fig=rbind(data_for_fig,cbind(rep(j,dim(RPKM)[1]),RPKM))
        }
        peak_strength_positive=RPKM_all_regions[which.max(data_for_fig[data_for_fig[,1]==1,3]),RPKM_all_regions[51,]==1]
        peak_strength_negative=RPKM_all_regions[which.max(data_for_fig[data_for_fig[,1]==1,3]),RPKM_all_regions[51,]==2]
        ttest_result=t.test(peak_strength_positive,peak_strength_negative)
        p_value_ttest[i]=ttest_result$p.value
        #write.table(data_for_fig,paste(work_dir,'/read_distribution/data_for_fig_',type[i],'_reverse_minus_strand_full_length_LINE1_two_groups_rep1_plus_minus_mean_TPM0point5_control.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
}
write.table(cbind(p_value_ttest,p_value_wilcox),paste(work_dir,'/read_distribution/reverse_minus_strand_full_length_LINE1_two_groups_rep1_plus_minus_mean_p_TPM0point5_control.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
########## use full length of LINE1 instead of 5UTR of LINE1 to draw the lineplot of MCF7 FOXA1 ChIP-seq binding.
########I calculate the p value for the read distribution of expressed and non-expression LINE1.
p_value for rep1 comparing with TPM<0.5 control: 
[1] 3.073665e-02 5.115824e-02 3.273184e-09 5.139687e-05 1.933851e-01 1.618120e-01

# Define labels
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
normed_line1_exp <- read.table(paste0('Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt'),header=T,sep="\t")
normed_line1_exp_matrix <- as.matrix(normed_line1_exp)
normed_line1_exp_matrix <- log2(normed_line1_exp_matrix+1)
line1=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")[,c(4,6)]
subfamily_anno <- merge(normed_line1_exp_matrix,line1,by.x=0,by.y='V4')
subfamily_anno <- subfamily_anno[,c(1,dim(subfamily_anno)[2])]
line1_subfamily_need_color <- c('L1HS','L1PA2','L1PA3','L1PA4','L1PA5')
subfamily_color <- c('red','green','blue','yellow','grey')
subfamily_anno[,3] <- 'white'
for (i in 1:length(line1_subfamily_need_color)){
        subfamily_anno[subfamily_anno[,2]==line1_subfamily_need_color[i],3] <- subfamily_color[i]
}
subfamily_anno[subfamily_anno[,3]=='white',2]='Other_LINE1'
type=unique(subfamily_anno[,2])

p_value=c(0.031,0.051,3.27e-09,5.14e-05,0.193,0.162)  ##TPM<0.5 control
work_dir="/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq"
label <- c('expressed', 'non-expressed')
# Create PDF for each type
pdf(paste(work_dir, '/read_distribution/FOXA1_ChIP-seq_read_distribution_full_length_line1_subfamilies_two_groups_rep1_plus_minus_mean.pdf', sep=""),onefile = TRUE,useDingbats=FALSE)
# Loop over each type
# Set up margins (bottom, left, top, right)
par(mar = c(10, 10, 10, 10))  # Increase margins as needed
for (g in 1:length(type)) {
    data_for_fig <- read.table(paste(work_dir, '/read_distribution/data_for_fig_', type[g], '_reverse_minus_strand_full_length_LINE1_two_groups_rep1_plus_minus_mean_TPM0point5_control.txt', sep=""), header=FALSE, sep="\t")

    ntrees <- max(data_for_fig[,1])  # Number of trees

    xrange <- range(data_for_fig[,2])  # X-axis range
    yrange <- range(data_for_fig[,3])  # Y-axis range

    # Set up the plot
    plot(xrange, yrange, type="n", xlab="LINE1", ylab="FOXA1 ChIP-seq RPKM", ylim=c(0, 7))

    colors <- c('red','grey')
	
    linetype <- c(1:ntrees)    # Line types
    plotchar <- seq(18, 18 + ntrees, 1)  # Plot characters

    # Add smoothed lines
    for (i in 1:ntrees) {
        tree <- data_for_fig[data_for_fig[,1] == i,]
        smoothed_values <- filter(tree[,3], rep(1, 1), sides=2)  # Applying a simple moving average with window size 5
        lines(tree[,2], smoothed_values, type="l", lwd=1.5, col=colors[i])  # Plot smoothed line
    }

    # Add legend
    legend(15, 6, label, cex=0.8, col=colors, lty=1)

    # Add title
    title(type[g])

    text(x = 25, y = 6.5, labels = paste("p =", p_value[g]), cex = 1.2)
}
dev.off()  # Close PDF device
