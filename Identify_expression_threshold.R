################Justify why we choose TPM > 0.5 as cutoff to define expressed LINE1.
gene_gtf_file=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/from_Michael/hg38_refGene.gtf',header=F,sep="\t")
gene_gtf_file[,10]=gsub("gene_id ","",sapply(strsplit(as.vector(as.matrix(gene_gtf_file[,9])),';'),function(x) x[1]))
chrom=paste0('chr',c(1:22,'X','Y','M'))
gene_bed_file=unique(gene_gtf_file[gene_gtf_file[,1]%in%chrom,c(1,4,5,10,8,7)])
write.table(gene_bed_file,"/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/from_Michael/hg38_refGene_clean.bed",row.names=F,col.names=F,sep="\t",quote=F)

module load bedtools
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38
micheal_ref_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/from_Michael
cat intron_hg38.bed cds_hg38.bed utr3_hg38.bed utr5_hg38.bed $micheal_ref_dir/hg38_refGene_clean.bed > intragenic_hg38.bed
bedtools subtract -a Whole_genome_hg38.bed -b intragenic_hg38.bed > intergenic_hg38.bed

ref_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38
ref_dir2=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/RepeatMasker/hg38
line1_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table

bedtools subtract -a $ref_dir/intergenic_hg38.bed -b $ref_dir2/RepeatMasker_LINE1.bed  > $ref_dir/intergenic_hg38_rm_all_LINE1.bed

###################
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Bedfiles
samples=(STAR_CaCo2_chrRNA_rep1 STAR_CaCo2_chrRNA_rep2 STAR_H460_S2R_chrRNA_rep1 STAR_H460_S2R_chrRNA_rep2 STAR_Huh7_chrRNA_rep1 STAR_Huh7_chrRNA_rep2 STAR_LNCaP_chrRNA_rep1 STAR_LNCaP_chrRNA_rep2 STAR_MCF10A_chrRNA_rep1 STAR_MCF10A_chrRNA_rep2 STAR_MCF7_chrRNA_rep1 STAR_MCF7_chrRNA_rep2 STAR_MDA-MB-231_chrRNA_rep1 STAR_MDA-MB-231_chrRNA_rep2 STAR_MOLM13_chrRNA_rep1 STAR_MOLM13_chrRNA_rep2 STAR_T47D_chrRNA_rep1 STAR_T47D_chrRNA_rep2 STAR_THP1_chrRNA_rep1 STAR_THP1_chrRNA_rep2 STAR_K562_chrRNA_rep2 STAR_K562_chrRNA_rep1 STAR_GM12878_chrRNA_rep1 STAR_GM12878_chrRNA_rep2)
for ((i=0;i<${#samples[*]};i++))
do
bedtools intersect -wa -a ${samples[i]}_clean.bed -b $ref_dir/intergenic_hg38_rm_all_LINE1.bed | uniq > ${samples[i]}_olp_intergenic_region_rm_all_LINE1.bed
done
##############

module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project')
line1_exp_mtx=read.table('chromatin_RNA_seq/Exp_matrix/line1_raw_readcount_merged_cell_lines.txt',header=T,sep="\t")[,-c(9,10,11,12)]
gene_exp_mtx=read.table('chromatin_RNA_seq/Exp_matrix/gene_raw_readcount_merged_cell_lines.txt',header=T,sep="\t")[,-c(9,10,11,12)]
combine_exp_mtx=rbind(gene_exp_mtx,line1_exp_mtx)
line1=read.table('annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
line1[,dim(line1)[2]+1]=line1[,3]-line1[,2]
line1_size=line1[,c(4,7)]
gene_size_old=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/gene_exonic_sizes.txt',header=F,sep="\t")
gene_size=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/from_Michael/gene_exonic_sizes.txt',header=F,sep="\t")
###the gene exon size of these two versions of gtf files are different (median difference are 304bp).  Michael's version has smaller size.
####I decide to use Michael's version as he used it to calculate the read count.
colnames(line1_size)=colnames(gene_size)
combine_size=rbind(gene_size,line1_size)

#TPM_all=combine_size[,1,drop=FALSE]
T=NULL
for (i in 1:dim(combine_exp_mtx)[2]){
	readcount_spe=combine_exp_mtx[,i,drop=FALSE]
	readcount_spe=merge(readcount_spe,combine_size,by.x=0,by.y='V1')
	t=readcount_spe[,2]/readcount_spe[,3]
	T[i]=sum(t)
	#readcount_spe[,4]=round(t*10^6/T,4)
	#TPM_all=merge(TPM_all,readcount_spe[,c(1,4)],by.x='V1',by.y='Row.names')
}
###get the T value for all samples.

library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Bedfiles')
samples=c('CaCo2_chrRNA_rep1','CaCo2_chrRNA_rep2','GM12878_chrRNA_rep1','GM12878_chrRNA_rep2','Huh7_chrRNA_rep1','Huh7_chrRNA_rep2','K562_chrRNA_rep1','K562_chrRNA_rep2','MCF7_chrRNA_rep1','MCF7_chrRNA_rep2','MDA-MB-231_chrRNA_rep1','MDA-MB-231_chrRNA_rep2','MOLM13_chrRNA_rep1','MOLM13_chrRNA_rep2','T47D_chrRNA_rep1','T47D_chrRNA_rep2','THP1_chrRNA_rep1','THP1_chrRNA_rep2','H460_S2R_chrRNA_rep1','H460_S2R_chrRNA_rep2')
#samples=c('CaCo2_chrRNA_rep1','CaCo2_chrRNA_rep2','Huh7_chrRNA_rep1','Huh7_chrRNA_rep2','LNCaP_chrRNA_rep1','LNCaP_chrRNA_rep2','MCF10A_chrRNA_rep1','MCF10A_chrRNA_rep2','MCF7_chrRNA_rep1','MCF7_chrRNA_rep2','MDA-MB-231_chrRNA_rep1','MDA-MB-231_chrRNA_rep2','MOLM13_chrRNA_rep1','MOLM13_chrRNA_rep2','T47D_chrRNA_rep1','T47D_chrRNA_rep2','THP1_chrRNA_rep1','THP1_chrRNA_rep2')
#total_reads=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Quality_checking/total_reads_mapping_to_line1_or_gene_merged_cell_lines.txt',header=T,sep="\t")
#total_reads=total_reads[-c(3,4,7,8,23,24),]
intergenic_region=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/intergenic_hg38_rm_all_LINE1.bed',header=F,sep="\t")
intergenic_region[,4]=intergenic_region[,3]-intergenic_region[,2]
filtered_intergenic_region=intergenic_region[(intergenic_region[,4]>5000)&(intergenic_region[,4]<9000),]
background_TPM=matrix(NA,dim(filtered_intergenic_region)[1],length(samples),dimnames=list(c(),samples))
t_all=matrix(NA,dim(filtered_intergenic_region)[1],length(samples),dimnames=list(c(),samples))
for (i in 1:length(samples)){
	dat=read.table(paste0('STAR_',samples[i],'_olp_intergenic_region_rm_all_LINE1.bed'),header=F,sep="\t")
	filtered_intergenic_region_anno=bedtoolsr::bt.intersect(filtered_intergenic_region,dat,wa=TRUE,c=TRUE,F=1)
    t=filtered_intergenic_region_anno[,5]/filtered_intergenic_region_anno[,4]
	background_TPM[,i]=round(t*10^6/T[i],4)
	t_all[,i]=t
}
#blacklist_flag=rowSums(t_all*100>5)>0
blacklist_flag=rowSums(t_all*100)>10
####it means rowSums(t_all)>0.1
###> sum(apply(t_all, 2, function(x) quantile(x, 0.90)))
##[1] 0.08736557
filtered_background_TPM=background_TPM[!blacklist_flag,]
write.table(filtered_background_TPM,'intergenic_regions_rm_all_LINE1_over5k_lessthan9k_TPM_rm_depth_over10_10celllines.txt',row.names=FALSE,col.names=TRUE,sep="\t",quote=F)

setwd('/Volumes/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Bedfiles')
background_TPM=read.table('intergenic_regions_rm_all_LINE1_over5k_lessthan9k_TPM_rm_depth_over10_10celllines.txt',header=TRUE,sep="\t")

background_TPM[,dim(background_TPM)[2]+1]=rowMeans(background_TPM)
data_for_fig=NULL
for (i in 1:dim(background_TPM)[2]){
    data_for_fig=rbind(data_for_fig,cbind(background_TPM[,i],colnames(background_TPM)[i]))
}
data_for_fig=as.data.frame(data_for_fig)
data_for_fig[,2]=gsub('chrRNA_','',data_for_fig[,2])
data_for_fig[,2]=gsub('V21','Mean',data_for_fig[,2])
data_for_fig[,1]=as.numeric(data_for_fig[,1])
colnames(data_for_fig)=c('TPM','Cell_lines')

data_for_fig=data_for_fig[data_for_fig[,2]=='Mean',]
quantile_95 <- round(quantile(data_for_fig$TPM, 0.95),1)
setwd('/Volumes//jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq')
pdf('background_TPM_rm_all_LINE1_over5k_lessthan9k_rm_depth_over10_10cellines.pdf', width = 5, height = 4,useDingbats=FALSE)
ggplot(data_for_fig, aes(x=TPM)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white",bins = 100)+
 geom_vline(xintercept = quantile_95, linetype = "dashed", color = "red", size = 0.5) +  
 annotate("text", x = quantile_95, y = 0.1, label = "p < 0.05", color = "red", vjust = -0.5, hjust = 0.4) + 
 geom_density(alpha=.2, fill="#FF6666") + xlim(0,0.8) + ylim(0,25) + ggtitle('rm_all_LINE1_rm_depth_over10_10celllines')
dev.off()
################Justify why we choose TPM > 0.5 as cutoff to define expressed LINE1.