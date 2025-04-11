#########use chromunity to call multi-way interactions for LINE1 project.
############sliding window  chromunity
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))

library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
library(GenomicRanges)
library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(magrittr)
library(gUtils)

#IMPORTANT STEP, move all chunks of .pore_c.parquet files into a new folder "porec_parquet_only"
## This is the path to directory with all chunks
parquet_files = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/MCF7_nextflow_output/chromunity/null_chromunity_parquet"
#parquet_files = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/K562_nextflow_output/chromunity/null_chromunity_parquet"
#parquet_files = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/T47D_nextflow_output/chromunity/null_chromunity_parquet"
#parquet_files = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/5UTR902_nextflow_output"
#parquet_files = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/SCR306_nextflow_output"

## Generating GRanges
name='MCF7_merged'
name='K562_merged'
name='T47D_081423'
name='5UTR902'
name='SCR306'
grange_file = parquet2gr(parquet_files, mc.cores = 10)
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity')
#setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity')
saveRDS(grange_file, paste0("grange_",name,".rds"))

#Loading an the grange file that I already converted from poreC parquet files 
this_gr <-readRDS(paste0("grange_",name,".rds"))

## Removing low quality monomers, using gUtils functions 
this_gr = gUtils::`%Q%`(this_gr, this_gr$pass_filter == TRUE)
output=as.data.frame(this_gr)
write.table(output[,c(1:3,6)],paste0("grange_",name,"_clean.txt"),row.names=F,col.names=T,sep='\t',quote=F)

set.seed(125)
all_concatemers = unique(this_gr$read_idx)
subsample_concatemers = sample(all_concatemers, length(all_concatemers)/2)
this_gr_subsample = this_gr %Q% (read_idx %in% subsample_concatemers)
this_gr_subsample$cid = this_gr_subsample$read_idx

load('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/mutiway_interactions/chromunity/gc_frag_cov.Rdata')
#This covariate object will be used to assign the respective values to each bin in bin-sets.

#"Sliding window" Chromunity.
chrom=paste0('chr',c(1:22,'X','Y'))
window_size=2e7
for (i in 1:length(chrom)){
	this_sliding_chrom = sliding_window_chromunity(concatemers = this_gr_subsample, resolution = 5e4, window.size = window_size, take_sub_sample = TRUE, chr  = chrom[i], subsample.frac = 0.5, mc.cores = 10)
	#saveRDS(this_sliding_chrom, paste0(name,"_",chrom[i],"_chromunity.rds"))
	write.table(this_sliding_chrom$binsets,paste0(name,"_",chrom[i],"_binsets_result_20M_window.txt"),row.names=T,col.names=T,sep="\t",quote=F)

	# Testing called communities for synergy (significant enrichment of high order interactions)
	#Step 1
	#Creating covariate object
	#set covariates
	#In this demonstration we will use GC content and number of fragments as covariates

	#Step2
	#Annotation of bin-sets with concatemer counts and covariates
	
	## concatemers used for testing
	this.chrom = gr2dt(this_sliding_chrom$concatemers)

	## Concatemers for testing
	this_gr_testing = dt2gr(gr2dt(this_gr_subsample)[!read_idx %in% unique(this.chrom$read_idx)]) 
	annotated_chrom = chromunity::annotate(binsets = this_sliding_chrom$binsets,
							   k = 10,
							   concatemers = this_gr_testing,
							   covariates = gc_frag_cov, resolution = 5e4,
							   mc.cores = 10) 
	set.seed(198)
	back_gr = sliding_window_background(chromosome= chrom[i], binsets = this_sliding_chrom$binsets, n = 1000,resolution = 5e4)
	back_gr[, V1 := NULL]
	back_gr = na.omit(back_gr) 
	## Adding few filters to remove outlier simulation, can be customized
	## Removing bins less than resolution and lying out of bounds

	## Getting seqlengths of each chromosome
	upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T) 
	setkeyv(back_gr, c("seqnames", "start"))
	back_gr = back_gr[!bid %in% back_gr[width < (5e4-1)]$bid]
	back_gr = gr2dt(gr.reduce(dt2gr(back_gr), by = "bid"))
	back_gr$bid <- as.factor(back_gr$bid)
	back_gr = merge(back_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
	back_gr = back_gr[end < V2][start < V2]
	back_gr[, overall.cardinality := .N, by = bid]
	back_gr = back_gr[overall.cardinality > 1]
	back_gr = dt2gr(back_gr)
	head(back_gr)
	set.seed(198)
	annotated_back = chromunity::annotate(binsets = back_gr,
							  k = 5,
							  concatemers = this_gr_testing,
							  covariates = gc_frag_cov, resolution = 5e4,
							  mc.cores = 10) 
	annotated_back = annotated_back[!bid %in% annotated_back[, .(sum(count)), by = bid][V1 == 0]$bid]

	set.seed(198)
	back_model = fit(annotated_back)

	annotated_chrom = sscore(annotated_chrom, model = back_model) 

	set.seed(198)
	this_synergy = na.omit(synergy(binsets = this_sliding_chrom$binsets,
								   annotated.binsets = annotated_chrom, model = back_model))
	head(this_synergy)
	## Synergy 2022-08-17 20:25:44: Scoring binsets

	this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)

	library("readr")
	write.table(this_synergy, paste0(name,"_",chrom[i],"_synergy_result_20M_window.txt"),row.names=F,col.names=T,sep="\t")
	write.table(this.chrom, paste0(name,"_",chrom[i],"_chromosome_concatemers_for_testing_20M_window.txt"),row.names=F,col.names=T,sep="\t")
	#saveRDS(this_gr_testing, paste0(name,"_",chrom[i],"_this_gr_testing.rds"))
	write.table(annotated_chrom, paste0(name,"_",chrom[i],"_bin-sets_chrom_annotations_20M_window.txt"),row.names=F,col.names=T,sep="\t")
	write.table(annotated_back, paste0(name,"_",chrom[i],"_background_bin-sets_chrom_annotations_20M_window.txt"),row.names=F,col.names=T,sep="\t")

	synergy=read.table(paste0(name,"_",chrom[i],"_synergy_result_20M_window.txt"),header=T,sep="\t")
	concatemer=read.table(paste0(name,"_",chrom[i],"_chromosome_concatemers_for_testing_20M_window.txt"),header=T,sep="\t")
	binsets=read.table(paste0(name,"_",chrom[i],"_binsets_result_20M_window.txt"),header=T,sep="\t")
	synergy_info=merge(binsets,synergy,by.x='bid',by.y='bid')
	write.table(synergy_info, paste0(name,"_",chrom[i],"_synergy_result_add_chrom_coordinates_20M_window.txt"),row.names=F,col.names=T,sep="\t")

}
############sliding window  chromunity

############regulatory element (RE) chromunity 
######Use three criterions to pre-defined interested regions to do RE chromunity analysis.
##(1) expressed LINE1; (2) genes that interact with expressed LINE1; (3)enhancers that interact with expressed LINE1.
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
gene_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/gene_TPM_normalized_readcount_merged_cell_lines.txt",header=TRUE,sep="\t")
gene_counts=log2(gene_counts+1)
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
genes=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/genes.bed',header=F,sep="\t")
genes_plus=genes[genes[,6]=="+",]
genes_plus[,2]=genes_plus[,2]-4000
genes_minus=genes[genes[,6]=="-",]
genes_minus[,3]=genes_minus[,3]+4000
genes=rbind(genes_plus,genes_minus)
genes[genes[,2]<0,2]=1

col_list=list(c(13,14),c(7,8),c(19,20))
PETs_name=c('MCF7_merged','K562_merged','T47D_081423')
exp_cutoff=0.5
gene_exp_cutoff=5
PET_cutoff=3
for (c in 1:length(PETs_name)){
	line1_counts_spe=line1_counts[,col_list[[c]]]  
	exp_line1=rownames(line1_counts_spe)[(line1_counts_spe[,1]>=exp_cutoff)&(line1_counts_spe[,2]>=exp_cutoff)]
	if (c==3){
		MCF7_interest_regions=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity/MCF7_expressed_line1_and_interact_genes_enhancers.bed',header=F,sep="\t")
		MCF7_interest_regions_info=merge(MCF7_interest_regions,line1_info,by.x=c('V1','V2','V3'),by.y=c('V1','V2','V3'))
		MCF7_exp_line1=as.vector(MCF7_interest_regions_info[,5])
		exp_line1=intersect(exp_line1,MCF7_exp_line1)  ## We do this because T47D has two many interest regions and chromunity RE can not works. So we used MCF7 to further filter expressed LINE1.
	}
	our_pair=unique(read.table(paste0('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance/',PETs_name[c],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_line1_with_PET_num.txt'),header=T,sep="\t")[,c(2:5)])
	our_pair_left=our_pair[,1:2]; our_pair_left[,3]=our_pair_left[,2]+1
	our_pair_right=our_pair[,3:4]; our_pair_right[,3]=our_pair_right[,2]+1
	our_pair_left[,4]=1:dim(our_pair)[1]
	our_pair_right[,4]=1:dim(our_pair)[1]

	exp_line1_info=line1_info[line1_info[,4]%in%exp_line1,1:4]

	#our_pair_left_olp_line1=bedtoolsr::bt.intersect(our_pair_left,exp_line1_info,wa=T)
	#line1_other_end=our_pair_right[our_pair_right[,4]%in%our_pair_left_olp_line1[,4],]
	#####didn't keep the names of LINE1 for the overlapped PETs.

	our_pair_left_olp_line1=bedtoolsr::bt.intersect(our_pair_left,exp_line1_info,wa=T,wb=T)
	line1_other_end=unique(merge(our_pair_right,our_pair_left_olp_line1,by.x='V4',by.y='V4')[,c(2:4,11)])
	#####keep the names of LINE1 for the overlapped PETs.

	gene_counts_spe=gene_counts[,col_list[[c]]]  
	exp_gene=rownames(gene_counts_spe)[(gene_counts_spe[,1]>=gene_exp_cutoff)&(gene_counts_spe[,2]>=gene_exp_cutoff)]
	#if (c==1){MCF7_exp_gene=exp_gene}
	# if (c==3){exp_gene=intersect(exp_gene,MCF7_exp_gene)}  #this way is still too many.
	if (c==3){
		MCF7_interest_regions_info=merge(MCF7_interest_regions,genes,by.x=c('V1','V2','V3'),by.y=c('V1','V2','V3'))
		MCF7_exp_gene=as.vector(MCF7_interest_regions_info[,5])
		exp_gene=intersect(exp_gene,MCF7_exp_gene)
	}
	
	expressed_gene=genes[genes[,4]%in%exp_gene,1:4]

	expressed_gene_olp_line1_PETs=bedtoolsr::bt.intersect(expressed_gene,line1_other_end,wa=T,wb=T)
	expressed_gene_with_line1_PETs=NULL
	uni_gene=unique(expressed_gene_olp_line1_PETs[,4])
	for (i in 1:length(uni_gene)){
		temp_idx=expressed_gene_olp_line1_PETs[,4]==uni_gene[i]
		max_PETs_num=max(table(expressed_gene_olp_line1_PETs[temp_idx,8]))
		expressed_gene_with_line1_PETs=rbind(expressed_gene_with_line1_PETs,cbind(unique(expressed_gene_olp_line1_PETs[temp_idx,1:4]),max_PETs_num))
	}
	expressed_gene_with_line1_PETs=expressed_gene_with_line1_PETs[expressed_gene_with_line1_PETs[,5]>PET_cutoff,]
	######require each gene has >3 PETs with one of the LINE1.

	chromHMM=chromHMM_anno_list[[c]]
	enhancer=chromHMM[chromHMM[,4]=='Enhancer',]
	promoter=chromHMM[chromHMM[,4]=='Promoter',]
	enhancer[,4]=paste0(enhancer[,4],'_',c(1:dim(enhancer)[1]))
	promoter[,4]=paste0(promoter[,4],'_',c(1:dim(promoter)[1]))

	enhancer_olp_line1_PETs=bedtoolsr::bt.intersect(enhancer,line1_other_end,wa=T,wb=T)
	enhancer_with_line1_PETs=NULL
	uni_enhancer=unique(enhancer_olp_line1_PETs[,4])
	for (i in 1:length(uni_enhancer)){
		temp_idx=enhancer_olp_line1_PETs[,4]==uni_enhancer[i]
		max_PETs_num=max(table(enhancer_olp_line1_PETs[temp_idx,8]))
		enhancer_with_line1_PETs=rbind(enhancer_with_line1_PETs,cbind(unique(enhancer_olp_line1_PETs[temp_idx,1:4]),max_PETs_num))
	}
	if (c==3){PET_cutoff=5}
	enhancer_with_line1_PETs=enhancer_with_line1_PETs[enhancer_with_line1_PETs[,5]>PET_cutoff,]

	exp_line1_info[,4]='line1'; expressed_gene_with_line1_PETs[,4]='gene'
	enhancer_with_line1_PETs[,4]='Enhancer'
	interest_line1_gene=rbind(exp_line1_info[,1:4],expressed_gene_with_line1_PETs[,1:4])
	interest_enhancer=enhancer_with_line1_PETs[,1:4]

	anno_dat=unique(rbind(interest_line1_gene,interest_enhancer))
	setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity')
	if (c!=3){
		write.table(anno_dat,paste0(cellline[c],'_expressed_line1_and_interact_genes_enhancers.bed'),row.names=F,col.names=F,sep="\t",quote=F)
	}else{
		write.table(anno_dat,paste0(cellline[c],'_expressed_line1_and_interact_genes_filtered_by_MCF7_enhancers_PET_over5.bed'),row.names=F,col.names=F,sep="\t",quote=F)
	}

	interest_promoter=unique(bedtoolsr::bt.intersect(promoter,interest_line1_gene,wa=T))
	interest_promoter[,4]='Promoter'
	anno_dat=unique(rbind(interest_promoter,interest_enhancer))
	write.table(anno_dat,paste0(cellline[c],'_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers.bed'),row.names=F,col.names=F,sep="\t",quote=F)
}

##use the filtered line1, gene and enhancer to run RE chromunity.
module load bedtools
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity
cellline=(MCF7 K562 T47D)
for ((i=0;i<${#cellline[*]};i++))
do
	bedtools sort -i ${cellline[i]}_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers.bed > ${cellline[i]}_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers_sorted.bed
	bedtools merge -i ${cellline[i]}_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers_sorted.bed > ${cellline[i]}_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers_merged.bed
done
############################################################################################


####do RE chromunity for three cell lines and ASO KD data.
############################################################################################
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity_SCR369_5UTR902_rep0/grange_SCR369.rds grange_SCR369.rds
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity_SCR369_5UTR902_rep0/grange_5UTR902.rds grange_5UTR902.rds

module purge && module load R/4.2.2-shlib

new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))

library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
library(GenomicRanges)
library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(magrittr)
library(gUtils)

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity')
#setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity')
cellline='MCF7'
#name=c('MCF7_merged','K562_merged','T47D_081423')
#name=c('SCR369')
#name=c('5UTR902')
name=c('MCF7_merged')
## getting E-P annotations
#anno_dat=read.table(paste0(cellline,'_promoter_enhancer_olp_expressed_line1_and_interact_genes_enhancers_merged'),header=F,sep="\t")
#anno_dat[,4]='interest_region'
anno_dat=read.table(paste0('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity/',cellline,'_expressed_line1_and_interact_genes_enhancers.bed'),header=F,sep="\t")
colnames(anno_dat)=c('seqnames','start','end','name')
chrom=paste0('chr',c(1:22,'X','Y'))
anno_dat=anno_dat[anno_dat[,1]%in%chrom,]
chmm <- anno_dat %>% dt2gr()
chmm_ep = dt2gr(gr2dt(chmm))

set.seed(198)
this_gr <- readRDS(paste0("grange_",name,".rds"))
this_gr = gUtils::`%Q%`(this_gr, this_gr$pass_filter == TRUE)  ## Removing low quality monomers, using gUtils functions 
reads_frag_count=table(this_gr$read_name)
reads_with_small_frag_count=names(reads_frag_count)[reads_frag_count<50]
this_gr = this_gr %Q% (read_name %in% reads_with_small_frag_count)

## resolution of the pad
resolution = 2.5e4  #this can not change. Otherwise it gave me error saying 'Shave' is not available.
tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)

## Targets
this_EP = chmm_ep
targets_EP = gr.reduce(this_EP+resolution)

## Subsample
this_gr = this_gr %&% targets_EP

## Subsampling concatemers (use 1/4 as training and another 1/4 as testing)
#out_label='2over4'
out_label='2over6'  
out_label='2over12'  
set.seed(125)
all_concatemers = unique(this_gr$read_idx)
training_concatemers = sample(all_concatemers, length(all_concatemers)/4)
#training_concatemers = sample(all_concatemers, length(all_concatemers)/6)
#training_concatemers = sample(all_concatemers, length(all_concatemers)/12)
remaining_concatemers = setdiff(all_concatemers, training_concatemers)
testing_concatemers = sample(remaining_concatemers, length(all_concatemers)/4)
#testing_concatemers = sample(remaining_concatemers, length(all_concatemers)/6)
#testing_concatemers = sample(remaining_concatemers, length(all_concatemers)/12)

this_gr_training = this_gr %Q% (read_idx %in% training_concatemers)
this_gr_testing = this_gr %Q% (read_idx %in% testing_concatemers)

##
this_gr_training$cid = this_gr_training$read_idx
this_gr_testing$cid = this_gr_testing$read_idx

## Running Chromunity
reso=3*5e4
this_re_chrom = re_chromunity(concatemers = this_gr_training, windows = targets_EP, piecewise = FALSE, shave = TRUE, resolution = reso, k.knn = 50,mc.cores = 10,bthresh = 3, cthresh = 3)
write.table(this_re_chrom$binsets,paste0(name,"_LINE1_RE_bin-sets_subsample_",out_label,".txt"),row.names=T,col.names=T,sep="\t")
#In this demonstration we will use GC content and number of fragments as covariates

load(file = "/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/mutiway_interactions/chromunity/gc_frag_cov.Rdata")
#This covariate object will be used to assign the respective values to each bin in bin-sets.

##suggestion from https://github.com/mskilab-org/chromunity/issues/3: 
##Make sure the k is set to 3. Higher k's will explode the problem space.
#If there are some very large chromunities, try excluding them before running synergy
#Try smaller bin sizes such as 25kb/10kb

set.seed(198)
annotated_re_chrom = chromunity::annotate(binsets = this_re_chrom$binsets,
           k = 3,
           concatemers = this_gr_testing,
           covariates = gc_frag_cov, resolution = reso,
           mc.cores = 10) 
write.table(this_re_chrom$binsets, paste0(name,"_LINE1_RE_bin-sets_chrom_annotations2_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")

set.seed(198)
back_re_gr = gr2dt(dt2gr(re_background(binsets = this_re_chrom$binsets, n=100000,resolution = reso)))
## Adding few filters to remove outlier simulation, can be customized
## Removing bins less than resolution and lying out of bounds

## Getting seqlengths of each chromosome
upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T) 
setkeyv(back_re_gr, c("seqnames", "start"))
back_re_gr = back_re_gr[!bid %in% back_re_gr[width < (reso-1)]$bid]
back_re_gr = gr2dt(gr.reduce(dt2gr(back_re_gr), by = "bid"))
back_re_gr$bid <- as.factor(back_re_gr$bid)
back_re_gr = merge(back_re_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
back_re_gr = back_re_gr[end < V2][start < V2]
back_re_gr[, overall.cardinality := .N, by = bid]
back_re_gr = back_re_gr[overall.cardinality > 1]
back_re_gr = dt2gr(back_re_gr)		
   
set.seed(198)
annotated_re_back = chromunity::annotate(binsets = back_re_gr,
           k = 3,
           concatemers = this_gr_testing,
           covariates = gc_frag_cov, resolution = reso,
           mc.cores = 5) 

## Removing edge cases with no counts
annotated_re_back = annotated_re_back[!bid %in% annotated_re_back[, .(sum(count)), by = bid][V1 == 0]$bid]
##The outputs of this_sliding_chrom function is a Chromunity object with a concatemer GRanges where each row is a monomer indexed by read_id as well as community id or chid. The bin-sets are created and stored in bin-sets GRanges object that maps community id to bin-id or bid.

set.seed(198)
back_re_model = fit(annotated_re_back)

annotated_re_chrom = sscore(annotated_re_chrom, model = back_re_model) 
head(annotated_re_chrom)

set.seed(198)
this_synergy = na.omit(synergy(binsets = this_re_chrom$binsets,
                annotated.binsets = annotated_re_chrom, model = back_re_model))
## Synergy 2022-08-17 20:30:29: Scoring binsets

this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)

save.image(file = paste0(name,"_chromunity.RData"))

write.table(this_synergy, paste0(name,"_LINE1_RE_synergy_result_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")
write.table(annotated_re_chrom, paste0(name,"_LINE1_RE_bin-sets_chrom_annotations_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")
write.table(this_re_chrom$binsets, paste0(name,"_LINE1_RE_bin-sets_chrom_annotations2_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")
write.table(annotated_re_back, paste0(name,"_LINE1_RE_background_bin-sets_chrom_annotations_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")

synergy=read.table(paste0(name,"_LINE1_RE_synergy_result_subsample_",out_label,".txt"),header=T,sep="\t")
binsets=read.table(paste0(name,"_LINE1_RE_bin-sets_chrom_annotations2_subsample_",out_label,".txt"),header=T,sep="\t")
synergy_info=merge(binsets,synergy,by.x='bid',by.y='bid')
write.table(synergy_info, paste0(name,"_LINE1_RE_synergy_result_add_chrom_coordinates2_subsample_",out_label,".txt"),row.names=F,col.names=T,sep="\t")
############################################################################################
############regulatory element (RE) chromunity 


####10/14/2024
######################################compare HILL to synergy events.
module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0

library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
chrom_size <- read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/chromSize/hg38.chromSize.clean', header = FALSE, sep = "\t")
whole_genome_size <- sum(chrom_size[,2])
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity')
cellline=c('MCF7','K562','T47D')
name=c('run04b_run06','K562_merged','T47D_081423')
name2=c('MCF7_run04b_run06','K562_merged','T47D_081423')
chrom=c(paste0('chr',c(1:22,'X')))
syn_regions_list=list()
cover_size_percent=NULL
for (i in 1:length(cellline)){
	syn_regions=NULL
	for (c in 1:length(chrom)){
		if (i==1){
			result=read.table(paste0('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/mutiway_interactions/chromunity/',name[i],'_',chrom[c],'_synergy_result_add_chrom_coordinates_20M_window.txt'),header=T,sep="\t")
		}else{
			result=read.table(paste0(name[i],'_',chrom[c],'_synergy_result_add_chrom_coordinates_20M_window.txt'),header=T,sep="\t")
		}
		#sig_result=result[result$p<0.05,]
		sig_result=result[result$fdr<0.01,]
		syn_regions=rbind(syn_regions,sig_result[,2:4])
	}
	
	if (i==1){
		RE_syn_regions=read.table(paste0(name2[i],'_LINE1_RE_synergy_result_add_chrom_coordinates2_subsample_2over4.txt'),header=T,sep="\t")
	}else{
		RE_syn_regions=read.table(paste0(name2[i],'_LINE1_RE_synergy_result_add_chrom_coordinates2_subsample_2over12.txt'),header=T,sep="\t")
	}
	sig_RE_syn_regions=RE_syn_regions[RE_syn_regions$fdr<0.01,2:4]
	syn_regions=rbind(syn_regions,sig_RE_syn_regions)
	merged_syn_regions=bedtoolsr::bt.merge(bedtoolsr::bt.sort(syn_regions))
	syn_regions_list[[i]]=merged_syn_regions
}

########highlight synergistic genes in the scatterplot.
name=c('MCF7_merged','K562_merged','T47D_081423')
line1_counts=read.table("/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/chromatin_RNA_seq/Exp_matrix/line1_TPM_normalized_readcount_merged_cell_lines_rm_sense_LINE1.txt",header=TRUE,sep="\t")
line1_counts=log2(line1_counts+1)
col_list=list(c(13,14),c(7,8),c(19,20))
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
for (i in 1:length(cellline)){
	setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
	sig_pair=unique(read.table(paste0(name[i],'_intra_p_cutoff0point0001_inter_p_cutoff0point0001_line1_with_PET_num.txt'),header=T,sep="\t")[,1:6])
	sig_pair=sig_pair[sig_pair[,2]==sig_pair[,4],]  ####only use intra-chrom interactions.
	line1_PET_num=sort(table(sig_pair[,1]))
	x <- 1:length(line1_PET_num)
	y <- as.vector(line1_PET_num)
	
	all_line1=names(line1_PET_num)
	all_line1_info=line1_info[line1_info[,4]%in%all_line1,1:4]
	all_line1_olp_chromunity=bedtoolsr::bt.intersect(all_line1_info,syn_regions_list[[i]],wa=T,f=1)

	scaled_x <- (x - min(x)) / (max(x) - min(x))
	scaled_y <- (y - min(y)) / (max(y) - min(y))
	slopes <- diff(scaled_y) / diff(scaled_x)
	closest_slope_index <- which.min(abs(slopes - 1))
	#tangent_slope <- 1
	tangent_slope <- slopes[closest_slope_index]
	tangent_intercept <- scaled_y[closest_slope_index] - tangent_slope * scaled_x[closest_slope_index]
	save.image(file = paste0(cellline[i],"_data_for_fig_intra_chrom_fdr0point01_K562_T47D_2over12.RData"))
}
	
##do this to avoid R package conflict.
#module unload R
#module unload bedtools
#module purge && module load R/4.2.2-shlib

cellline=c('MCF7','K562','T47D')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/interaction_significance')
for (i in 1:length(cellline)){
	load(file = paste0(cellline[i],"_data_for_fig_intra_chrom_fdr0point01_K562_T47D_2over12.RData"))
	
	data <- data.frame(scaled_x, scaled_y)
	data$chromunity <- 0
	data[names(line1_PET_num) %in% all_line1_olp_chromunity[, 4], "chromunity"] <- 1
	data$chromunity=as.factor(data$chromunity)

	new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
	.libPaths(c(new_lib_path, .libPaths()))
	new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
	.libPaths(c(new_lib_path, .libPaths()))
	library(ggplot2,lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2")
	pdf(paste0(cellline[i], '_highlight_both_slinding_window_and_RE_synergistic_LINE1_p0point0001_intra_chrom_fdr0point01_K562_T47D_2over12.pdf'), useDingbats = FALSE,width=5,height=5)
	ggplot(data, aes(x = scaled_x, y = scaled_y, color =chromunity )) +
		geom_point(shape = 21, size = 3) +  # Default color mapping will be based on 'chromunity'
		geom_point(data = data[closest_slope_index, , drop = FALSE], shape = 19, color = "green", size = 3) +
		geom_abline(intercept = tangent_intercept, slope = tangent_slope, linetype = "dashed", color = "blue") +
		scale_color_manual(values = c("0" = "grey", "1" = "violetred")) +  # Map 0 to pink and 1 to red
		labs(x = "Scaled X", y = "Scaled Y", title = cellline[i]) +
		theme_minimal()
	dev.off()
}
########highlight synergistic genes in the scatterplot.

########Calculate p value for the overlapping of HILL with synergy.
cellline=c('MCF7','K562','T47D')
name=c('MCF7_merged','K562_merged','T47D_081423')
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
N=dim(line1_info)[1]
for (i in 1:length(cellline)){
	hill_line1=read.table(paste0('HILL_analysis/',name[i],'_HILL_line1_intra_chrom.txt'),header=F,sep="\t")
	hill_line1_info=line1_info[line1_info[,4]%in%hill_line1[,1],1:4]
	hill_line1_olp_chromunity=bedtoolsr::bt.intersect(hill_line1_info,syn_regions_list[[i]],wa=T,f=1)

	all_line1_olp_chromunity=bedtoolsr::bt.intersect(line1_info[,1:4],syn_regions_list[[i]],wa=T,f=1)
	n=dim(unique(hill_line1_info))[1]
	k=dim(unique(hill_line1_olp_chromunity))[1]
	M=dim(unique(all_line1_olp_chromunity))[1]
	p=1-phyper(k-1,M,N-M,n)
}
###
MCF7: p=0
K562: p=1.482956e-10
T47D: p=9.185536e-08
########Calculate p value for the overlapping of HILL with synergy.
######################################compare HILL to synergy events.

#################################draw circosplot for MCF7,K562 and T47D as well as SCR369 and 5UTR902. Use fdr<0.01 to define significant synergy.
#######only show genes with the largest number of supporting concatemers (So the figure is not too crowd). 
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/mutiway_interactions/chromunity/grange_run04b_run06_clean.txt grange_MCF7_run04b_run06_clean.txt

cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity_SCR369_5UTR902_rep0/grange_SCR369_clean.txt grange_SCR369_clean.txt
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity_SCR369_5UTR902_rep0/grange_5UTR902_clean.txt grange_5UTR902_clean.txt


module purge && module load R/4.2.2-shlib
module load bedtools/2.30.0

library('bedtoolsr',lib.loc="/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2")
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
.libPaths(c(new_lib_path, .libPaths()))
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))
library(circlize,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old')

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C/chromunity')
#setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity_SCR369_5UTR902_rep0')
#setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/ASO_KD_poreC/chromunity')
chrom=paste0('chr',c(1:22,'X','Y'))
name='MCF7_run04b_run06'
#name='K562_merged'
name='T47D_081423'
#name='SCR369'
#name='5UTR902'
out_label='2over4'
#out_label='2over6'
out_label='2over12'
ref_gene=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/genes.bed',header=F,sep="\t")
top_gene_num=50
#sig_synergy_num=matrix(0,length(chrom),length(name),dimnames=list(chrom,name))
for (n in 1:length(name)){
    concatemer=read.table(paste0("grange_",name[n],"_clean.txt"),header=T,sep="\t")
    concatemer_bed=concatemer
    synergy=read.table(paste0(name[n],"_LINE1_RE_synergy_result_add_chrom_coordinates2_subsample_",out_label,".txt"),header=T,sep="\t")
    #sig_synergy=synergy[synergy$fdr<0.00001,]
	  sig_synergy=synergy[synergy$fdr<0.01,]
    #sig_synergy_num[i,n]=dim(sig_synergy)[1]
    uni_bid=unique(sig_synergy[,1])
    bin_bed1=NULL; bin_bed2=NULL
    for (b in 1:length(uni_bid)){
        sig_synergy_temp=sig_synergy[sig_synergy[,1]==uni_bid[b],]
        L=dim(sig_synergy_temp)[1]
        bin_bed1_temp=NULL; bin_bed2_temp=NULL
        for (l in 1:(L-1)){
            for (ll in (l+1):L){
                bin_bed1_temp=rbind(bin_bed1_temp,sig_synergy_temp[l,c(2:4,1)])
                bin_bed2_temp=rbind(bin_bed2_temp,sig_synergy_temp[ll,c(2:4,1)])
            }

        }
        bin_bed1=rbind(bin_bed1,bin_bed1_temp)
        bin_bed2=rbind(bin_bed2,bin_bed2_temp)
    }

    bin_bed1_olp_concatemer=bedtoolsr::bt.intersect(bin_bed1,concatemer_bed,wa=T,wb=T)
    bin_bed2_olp_concatemer=bedtoolsr::bt.intersect(bin_bed2,concatemer_bed,wa=T,wb=T)
    concatemer_in_bin=unique(rbind(bin_bed1_olp_concatemer[,5:7],bin_bed2_olp_concatemer[,5:7]))
    support_concatemer_num=NULL
    for (b in 1:dim(bin_bed1)[1]){
        olp_contatemer1=merge(bin_bed1[b,,drop=FALSE],bin_bed1_olp_concatemer,by.x=c('seqnames','start','end'),by.y=c('V1','V2','V3'))
        olp_contatemer2=merge(bin_bed2[b,,drop=FALSE],bin_bed2_olp_concatemer,by.x=c('seqnames','start','end'),by.y=c('V1','V2','V3'))
        olp_contatemer=unique(intersect(olp_contatemer1[,9],olp_contatemer2[,9]))
        support_concatemer_num[b]=length(olp_contatemer)
    }
    
    save.image(file = paste0(name[n],"_RE_synergy_result_circosplot_fdr0point01.RData"))
    load(paste0(name[n],"_RE_synergy_result_circosplot_fdr0point01.RData"))

    bin_bed1=bin_bed1[support_concatemer_num>1,]
    bin_bed2=bin_bed2[support_concatemer_num>1,]
    support_concatemer_num=support_concatemer_num[support_concatemer_num>1]

    #CTCF=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/BrandonChen/Data_BrandonChen/Datasets_for_Michael/ChIP-seq/MCF7/CTCF_MCF7_Stamatoyannopoulos_SE36_rep2/peaks/CTCF_MCF7_Stamatoyannopoulos_SE36_rep2_k10000_peaks.bed',header=F,sep="\t")
    #CTCF=CTCF[,1:3]
    L1HS=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/MichaelLeeJr/Data_MichaelLeeJr/LINE1_annotations/RepeatMasker_all_L1Hs_hg38.bed',header=F,sep="\t")[,c(1:3,5)]
    L1HS[,4]=1

    L1HS_full_length=L1HS[(L1HS[,3]-L1HS[,2])>5000,]

    interest_gene_bed=unique(bedtoolsr::bt.intersect(ref_gene,concatemer_in_bin,wa=T,c=T))
    ranked_concatemer_num=sort(interest_gene_bed[,7],decreasing=T)
    interest_gene_bed=interest_gene_bed[interest_gene_bed[,7]>ranked_concatemer_num[top_gene_num],1:6]

    col_flag=rep('darkorange',dim(bin_bed1)[1])
    lwd_flag=as.numeric(log10(support_concatemer_num))

    pdf(paste0(name[n],"_",out_label,"_RE_synergy_result_circosplot_fdr0point01_rm_CTCF_all_LINE1_ring_fdr0point01.pdf"))
    set.seed(123)
    circos.initializeWithIdeogram(plotType = NULL,chromosome.index = chrom)
    if (dim(interest_gene_bed)[1]>0){
        circos.genomicLabels(interest_gene_bed,labels.column=4, side = "outside")
    }
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.rect(xlim[1], 0, xlim[2], 1, col = 'lightgrey',border = "transparent")
        circos.text(mean(xlim), mean(ylim), gsub('chr','',chr), cex = 1.4, col = "dodgerblue1",
            facing = "inside", niceFacing = TRUE)
    }, track.height = 0.1, bg.border = NA)
    #circos.genomicDensity(CTCF, window.size = 1e6,col='red',track.height = 0.1, bg.border = 'lightgrey')

    #circos.genomicTrackPlotRegion(L1HS, ylim = c(0, 1),
    #    panel.fun = function(region, value, ...) {
    #        circos.genomicLines(region, value,type='h',col='seagreen', ...)
    #}, track.height = 0.07, bg.border = 'lightgrey')

    circos.genomicTrackPlotRegion(L1HS_full_length, ylim = c(0, 1),
        panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value,type='h',col='green', ...)
    }, track.height = 0.07, bg.border = 'lightgrey')

    circos.genomicLink(bin_bed1[,1:3], bin_bed2[,1:3], col = col_flag, lwd = lwd_flag, border = col_flag)
    dev.off()
}
#####################draw circosplot for MCF7,K562 and T47D as well as SCR369 and 5UTR902. Use fdr<0.01 to define significant synergy. 






