######calcualte the significant interactions from long-reads L1-capture-3C data.

###############The chromosome coordinates from pair files are the middle point of DNPII fragments
###############Extend the middle point to the start and end of DNPII fragments.
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/pair_files
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_061423_MCF7_DpnII-3C_PRO114_hg38basic/pairs/DpnII_run11_hg38_basic_unphased.unsorted.pairs MCF7_061423.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_080722_vs_041523_MCF7_DpnII-3C_hg38basic/pairs/DpnII_run04b_hg38_basic_unphased.unsorted.pairs MCF7_080722.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_080722_vs_041523_MCF7_DpnII-3C_hg38basic/pairs/DpnII_run06_hg38_basic_unphased.unsorted.pairs MCF7_041523.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_042623_K562_DpnII-3C_hg38basic/pairs/DpnII_run07_hg38_basic_unphased.unsorted.pairs K562_042623.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_020423_K562-DMSO-DpnII-3C_MIN114_hg38basic/pairs/DpnII_run16_hg38_basic_unphased.unsorted.pairs K562_020423.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_070523_K562_DpnII-3C_MIN114_hg38basic_PRO114reseq/pairs/DpnII_run12_hg38_basic_unphased.unsorted.pairs K562_070523.pairs
ln -s /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/DNAnexus_backup/results_081423_T47D-LNCaP-DpnII-3C_PRO114_hg38basic/pairs/DpnII_run13_hg38_basic_unphased.unsorted.pairs T47D_081423.pairs
cat *K562* > K562_merged.pairs
cat *MCF7* > MCF7_merged.pairs
name=(MCF7_merged K562_merged T47D_081423)
for ((i=0;i<${#name[*]};i++))
do
cut -f2,3,6,1 ${name[i]}.pairs > ${name[i]}_left.pairs
cut -f4,5,7,1 ${name[i]}.pairs > ${name[i]}_right.pairs
done

###R scripts.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/For_HiC')
dpnii=read.table('hg38_DPNII_frag_clean.bed',header=F,sep="\t")
dpnii[,dim(dpnii)[2]+1]=dpnii[,2]+floor((dpnii[,3]-dpnii[,2])/2)
dpnii[,dim(dpnii)[2]+1]=dpnii[,2]+round((dpnii[,3]-dpnii[,2])/2)
dpnii=dpnii[,c(1:3,7,8)]
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/pair_files')
name=c('MCF7_merged','K562_merged','T47D_081423')
for (i in 1:length(name)){
  pairs_left=read.table(paste0(name[i],'_left.pairs'),header=F,sep="\t")
  pairs_left_olp_dpnii1=merge(pairs_left,dpnii,by.x=c('V2','V3'),by.y=c('V1','V7'))
  pairs_left_olp_dpnii2=merge(pairs_left,dpnii,by.x=c('V2','V3'),by.y=c('V1','V8'))
  pairs_left_olp_dpnii=unique(rbind(pairs_left_olp_dpnii1[,c(1,5,6,3,4)],pairs_left_olp_dpnii2[,c(1,5,6,3,4)]))
  write.table(pairs_left_olp_dpnii,paste0(name[i],'_pairs_left.bed'),row.names=F,col.names=F,sep="\t",quote=F)

  pairs_right=read.table(paste0(name[i],'_right.pairs'),header=F,sep="\t")
  pairs_right_olp_dpnii1=merge(pairs_right,dpnii,by.x=c('V2','V3'),by.y=c('V1','V7'))
  pairs_right_olp_dpnii2=merge(pairs_right,dpnii,by.x=c('V2','V3'),by.y=c('V1','V8'))
  pairs_right_olp_dpnii=unique(rbind(pairs_right_olp_dpnii1[,c(1,5,6,3,4)],pairs_right_olp_dpnii2[,c(1,5,6,3,4)]))
  write.table(pairs_right_olp_dpnii,paste0(name[i],'_pairs_right.bed'),row.names=F,col.names=F,sep="\t",quote=F)
}
###############Extend the middle point to the start and end of DNPII fragments.


################overlap with LINE1.
module load bedtools
LINE1_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table
cd /research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/pair_files
out_dir=/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C
name=(MCF7_merged K562_merged T47D_081423)
name=(MCF7_merged)
for ((i=0;i<${#name[*]};i++))
do
bedtools intersect -wa -a ${name[i]}_pairs_left.bed -b $LINE1_dir/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed -v > $out_dir/${name[i]}_pairs_left_not_olp_LINE1.bed
bedtools intersect -wa -a ${name[i]}_pairs_left.bed -b $LINE1_dir/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed > $out_dir/${name[i]}_pairs_left_olp_LINE1.bed

bedtools intersect -wa -a ${name[i]}_pairs_right.bed -b $LINE1_dir/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed -v > $out_dir/${name[i]}_pairs_right_not_olp_LINE1.bed
bedtools intersect -wa -a ${name[i]}_pairs_right.bed -b $LINE1_dir/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed > $out_dir/${name[i]}_pairs_right_olp_LINE1.bed
done
################overlap with LINE1.

##R scripts
################output read pairs overlapped or not overlapped with LINE1
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/')
#name=c('MCF7_merged','K562_merged','T47D_081423')
name=c('K562_merged')
pairs_left=read.table(paste0('pair_files/',name,'_left.pairs'),header=F,sep="\t")
pairs_right=read.table(paste0('pair_files/',name,'_right.pairs'),header=F,sep="\t")
pairs_left_olp_line1=read.table(paste0('L1_capture_3C/',name,'_pairs_left_olp_LINE1.bed'),header=F,sep="\t")
pairs_right_olp_line1=read.table(paste0('L1_capture_3C/',name,'_pairs_right_olp_LINE1.bed'),header=F,sep="\t")
reads_olp_line1=union(pairs_left_olp_line1[,4],pairs_right_olp_line1[,4])

pairs_left_not_olp_line1=pairs_left[!(pairs_left[,1]%in%reads_olp_line1),]
pairs_right_not_olp_line1=pairs_right[!(pairs_right[,1]%in%reads_olp_line1),]
pairs_not_olp_line1=merge(pairs_left_not_olp_line1,pairs_right_not_olp_line1,by.x='V1',by.y='V1')
pairs_not_olp_line1_intra=pairs_not_olp_line1[pairs_not_olp_line1[,2]==pairs_not_olp_line1[,5],]
pairs_not_olp_line1_inter=pairs_not_olp_line1[pairs_not_olp_line1[,2]!=pairs_not_olp_line1[,5],]
write.table(pairs_not_olp_line1_intra,paste0('L1_capture_3C/',name,'_pairs_not_olp_line1_intra.txt'),row.names=F,col.names=F,sep="\t",quote=F)
write.table(pairs_not_olp_line1_inter,paste0('L1_capture_3C/',name,'_pairs_not_olp_line1_inter.txt'),row.names=F,col.names=F,sep="\t",quote=F)
######output read pairs not overlapped with LINE1.

pairs_left_olp_line1=pairs_left[(pairs_left[,1]%in%reads_olp_line1),]
pairs_right_olp_line1=pairs_right[(pairs_right[,1]%in%reads_olp_line1),]
pairs_olp_line1=merge(pairs_left_olp_line1,pairs_right_olp_line1,by.x='V1',by.y='V1')
pairs_olp_line1_intra=pairs_olp_line1[pairs_olp_line1[,2]==pairs_olp_line1[,5],]
pairs_olp_line1_inter=pairs_olp_line1[pairs_olp_line1[,2]!=pairs_olp_line1[,5],]
write.table(pairs_olp_line1_intra,paste0('L1_capture_3C/',name,'_pairs_olp_line1_intra.txt'),row.names=F,col.names=F,sep="\t",quote=F)
write.table(pairs_olp_line1_inter,paste0('L1_capture_3C/',name,'_pairs_olp_line1_inter.txt'),row.names=F,col.names=F,sep="\t",quote=F)
######output read pairs overlapped with LINE1.

#########calcaulte the distance between background intra-pairs and draw histgram.
dist_intra=abs(pairs_not_olp_line1_intra[,3]-pairs_not_olp_line1_intra[,6])
dist_intra_for_fig=as.data.frame(dist_intra)
write.table(dist_intra_for_fig,paste0('L1_capture_3C/',name,'_background_intra_pairs_distance.txt'),row.names=F,col.names=F,sep="\t",quote=F)

pdf(paste0('L1_capture_3C/',name,'_background_intra_pairs_distance_distribution.pdf'), width = 5, height = 4,useDingbats=FALSE)
ggplot(dist_intra_for_fig, aes(x=dist_intra)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white",bins = 200)+
 geom_density(alpha=.2, fill="#FF6666")  + xlim(0,5.0e+7) + ylim(0,3e-07) + ggtitle(name)
dev.off()
#########calculate the distance between background intra-pairs and draw histgram.

#########calculate the PET num in the nearby 50 bins around intergenic regions (use intergenic regions as background).
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
name=c('MCF7_merged','K562_merged','T47D_081423')
pairs_not_olp_line1_intra=read.table(paste0(name,'_pairs_not_olp_line1_intra.txt'),header=F,sep="\t")
K=5000; B=50
temp1=pairs_not_olp_line1_intra[,c(2,3,5,6)]; colnames(temp1)=c('V1','V2','V3','V4')
temp2=pairs_not_olp_line1_intra[,c(5,6,2,3)]; colnames(temp2)=c('V1','V2','V3','V4')
background_pairs=rbind(temp1,temp2)  ##make the pair reciprocal
intergenic_region=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/intergenic_hg38_rm_all_LINE1.bed',header=F,sep="\t")
intergenic_region[,4]=intergenic_region[,3]-intergenic_region[,2]
filtered_intergenic_region=intergenic_region[(intergenic_region[,4]>5000),]
filtered_intergenic_region[,4]=filtered_intergenic_region[,2]+floor((filtered_intergenic_region[,3]-filtered_intergenic_region[,2])/2)
filtered_intergenic_region[,2]=filtered_intergenic_region[,4]-K/2
filtered_intergenic_region[,3]=filtered_intergenic_region[,4]+K/2
filtered_intergenic_region[,4]=filtered_intergenic_region[,3]-filtered_intergenic_region[,2]
chrom=paste0('chr',c(1:22,'X'))

dist_breaks=K*c(1:B)
for (i in 1:length(chrom)){
	filtered_intergenic_region_spe=filtered_intergenic_region[filtered_intergenic_region[,1]==chrom[i],]
	background_pairs_spe=background_pairs[background_pairs[,1]==chrom[i],]
	dist_mtx=matrix(0,dim(filtered_intergenic_region_spe)[1],B)
	for (j in 1:dim(filtered_intergenic_region_spe)[1]){
		idx=(background_pairs_spe[,2]>filtered_intergenic_region_spe[j,2])&(background_pairs_spe[,2]<filtered_intergenic_region_spe[j,3])
		background_pairs_spe2=background_pairs_spe[idx,]
		dist_temp=abs(background_pairs_spe2[,4]-background_pairs_spe2[,2])
		bin_flag=floor(dist_temp/K)
		if (sum(bin_flag<B)>0){
			for (b in 1:B){
				dist_mtx[j,b]=sum((bin_flag>=(b-1))&(bin_flag<b))
			}
		}
	}
	write.table(dist_mtx,paste0('background_estimation/',name,'_intergenic_background_',B,'_',K,'bin_PETs_num_',chrom[i],'.txt'),row.names=F,col.names=F,sep="\t",quote=F)
}
#########calculate the PET num in the nearby 50 bins around intergenic regions (use intergenic regions as background).


############calculate intra-chromosomal p value.
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
.libPaths(c(new_lib_path, .libPaths()))
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))
library(data.table,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old')
library(dplyr)

setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
name=c('MCF7_merged','K562_merged','T47D_081423')
K <- 5000
B <- 50

# Read input files using fread for efficiency
pairs_olp_line1_intra <- fread(paste0(name, '_pairs_olp_line1_intra.txt'), header=FALSE, sep="\t")
line1_info <- fread('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed', header=FALSE, sep="\t")

# Process pairs to create reciprocal pairs
temp1 <- pairs_olp_line1_intra[, .(V1 = V2, V2 = V3, V3 = V5, V4 = V6)]
temp2 <- pairs_olp_line1_intra[, .(V1 = V5, V2 = V6, V3 = V2, V4 = V3)]
line1_pairs <- rbindlist(list(temp1, temp2))


# Define chromosomes
chromosomes <- paste0('chr', c(1:22, 'X'))

# Function to process each chromosome
process_chromosome <- function(chromosome) {
  # Filter for the current chromosome
  line1_info_spe <- line1_info[V1 == chromosome]
  line1_pairs_spe <- line1_pairs[V1 == chromosome]

  # Read background PET numbers for the current chromosome
  background_PET_num <- fread(paste0('background_estimation/',name, '_intergenic_background_', B, '_', K, 'bin_PETs_num_', chromosome, '.txt'), header=FALSE, sep="\t")
  
  # Prepare a list to collect results
  line1_pairs_add_p <- list()
  
  # Iterate over each LINE1 element
  for (i in seq_len(nrow(line1_info_spe))) {
    line1_ID <- line1_info_spe[i, V4]
    idx <- line1_pairs_spe$V2 > line1_info_spe[i, V2] & line1_pairs_spe$V2 < line1_info_spe[i, V3]
    
    if (sum(idx) > 0) {
      line1_pairs_spe2 <- line1_pairs_spe[idx]
      dist_temp <- abs(line1_pairs_spe2$V4 - line1_pairs_spe2$V2)
      bin_flag <- floor(dist_temp / K)
      
      bin_PET_num <- table(bin_flag)
      uni_bin_flag <- as.numeric(names(bin_PET_num))
      
      # Calculate p-values
      p_value_temp <- sapply(uni_bin_flag, function(b) {
        if (b == 0) {
          return(1)
        }
        if (b <= 50) {
          return(sum(bin_PET_num[names(bin_PET_num)==b] <= background_PET_num[[b]]) / nrow(background_PET_num))
        }
        return(sum(bin_PET_num[names(bin_PET_num)==b] <= background_PET_num[[B]]) / nrow(background_PET_num))
      })
      
      # Update p-values in the pairs data
      #line1_pairs_spe2[, V5 := p_value_temp[bin_flag + 1]]

      for (j in 1:length(uni_bin_flag)){
        line1_pairs_spe2[bin_flag==uni_bin_flag[j],V5 := p_value_temp[j]]
      }

      line1_pairs_add_p[[i]] <- cbind(line1_pairs_spe2, line1_ID)
    }
  
    # Combine all results for the current chromosome
    line1_pairs_add_p_all <- rbindlist(line1_pairs_add_p, fill=TRUE)
  }
  return(line1_pairs_add_p_all)
}

# Process all chromosomes
line1_pairs_add_p_all <- rbindlist(lapply(chromosomes, process_chromosome), fill=TRUE)

# Write output
fwrite(line1_pairs_add_p_all, paste0(name, '_', B, 'bins_binsize_', K, '_line1_intra_pairs_p_value.txt'), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
############calculate intra-chromosomal p value.


###calculate the PET num between intergenic regions and other chromosomes (use intergenic regions as background).
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old"
.libPaths(c(new_lib_path, .libPaths()))
new_lib_path <- "/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old_2"
.libPaths(c(new_lib_path, .libPaths()))
library(data.table,lib.loc='/home/yzhang57/R/x86_64-pc-linux-gnu-library/4.2-old')
library(dplyr)

# Set working directory
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')

# Read input files
name=c('MCF7_merged','K562_merged','T47D_081423')
pairs_not_olp_line1_inter <- fread(paste0(name, '_pairs_not_olp_line1_inter.txt'), header=FALSE, sep="\t")
colnames(pairs_not_olp_line1_inter) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6','V7')

# Create reciprocal pairs
temp1 <- pairs_not_olp_line1_inter %>% select(V2, V3, V5, V6) %>% setNames(c('V1', 'V2', 'V3', 'V4'))
temp2 <- pairs_not_olp_line1_inter %>% select(V5, V6, V2, V3) %>% setNames(c('V1', 'V2', 'V3', 'V4'))
background_pairs <- rbind(temp1, temp2)

# Read and process intergenic regions
intergenic_region <- fread('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/intergenic_hg38_rm_all_LINE1.bed', header=FALSE, sep="\t")
intergenic_region <- intergenic_region %>%
  mutate(V4 = V3 - V2) %>%
  filter(V4 > 5000) %>%
  mutate(V4 = V2 + floor((V3 - V2) / 2),
         V2 = V4 - 5000 / 2,
         V3 = V4 + 5000 / 2,
         V4 = V3 - V2)

# Read chromosome sizes
chrom_size <- fread('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/chromSize/hg38.chromSize.clean', header=FALSE, sep="\t")
setnames(chrom_size, c('Chrom', 'Size'))

# Define variables
chrom <- paste0('chr', c(1:22, 'X'))
K <- 10^6

# Loop through chromosomes
for (i in seq_along(chrom)) {
  chrom_current <- chrom[i]
  filtered_intergenic_region_spe <- intergenic_region %>% filter(V1 == chrom_current)
  background_pairs_spe <- background_pairs %>% filter(V1 == chrom_current)
  rest_chrom <- setdiff(chrom, chrom_current)
  
  for (rest_chrom_current in rest_chrom) {
    B <- ceiling(chrom_size[chrom_size$Chrom == rest_chrom_current, Size] / K)
    dist_mtx <- matrix(0, nrow = nrow(filtered_intergenic_region_spe), ncol = B)
    
    for (j in seq_len(nrow(filtered_intergenic_region_spe))) {
      region <- filtered_intergenic_region_spe[j,]
      idx <- background_pairs_spe$V2 > region$V2 & background_pairs_spe$V2 < region$V3
      background_pairs_spe2 <- background_pairs_spe[idx,]
      background_pairs_spe3 <- background_pairs_spe2 %>% filter(V3 == rest_chrom_current)
      
      if (nrow(background_pairs_spe3) > 0) {
        bin_flag <- ceiling(background_pairs_spe3$V4 / K)
        bin_PET_num <- table(bin_flag)
        uni_bin_flag <- as.numeric(names(bin_PET_num))
        
        dist_mtx[j, uni_bin_flag] <- bin_PET_num
      }
    }
    
    fwrite(as.data.table(dist_mtx), paste0('background_estimation/',name, '_', chrom_current, '_intergenic_background_to_', rest_chrom_current, '_1M_bins_PETs_num.txt'), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}
###calculate the PET num between intergenic regions and other chromosomes (use intergenic regions as background).



#########calculate the PET num between line1 and other chromosomes and compare with the intergenic background, calculate the p value.
setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/L1_capture_3C')
name=c('MCF7_merged','K562_merged','T47D_081423')
pairs_olp_line1_inter=read.table(paste0(name,'_pairs_olp_line1_inter.txt'),header=F,sep="\t")
temp1=pairs_olp_line1_inter[,c(2,3,5,6)]; colnames(temp1)=c('V1','V2','V3','V4')
temp2=pairs_olp_line1_inter[,c(5,6,2,3)]; colnames(temp2)=c('V1','V2','V3','V4')
line1_pairs=rbind(temp1,temp2)  ##make the pair reciprocal
line1_info=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Projects/LINE1_project/annotation_table/hg38_rmsk_all_LINE1_over5000bp.sorted.subfID.bed',header=F,sep="\t")
chrom=paste0('chr',c(1:22,'X'))
K=10^6
chrom_size=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/chromSize/hg38.chromSize.clean',header=F,sep="\t")
line1_pairs_add_p_all=NULL
for (i in 1:length(chrom)){
  line1_info_spe=line1_info[line1_info[,1]==chrom[i],]
  line1_pairs_spe=line1_pairs[line1_pairs[,1]==chrom[i],]
  rest_chrom=setdiff(chrom,chrom[i])
  for (c in 1:length(rest_chrom)){
    line1_pairs_spe2=line1_pairs_spe[line1_pairs_spe[,3]==rest_chrom[c],]
    if (dim(line1_pairs_spe2)[1]>0){
      B=ceiling(chrom_size[chrom_size[,1]==rest_chrom[c],2]/K)
      background_PET_num=read.table(paste0('background_estimation/',name,'_',chrom[i],'_intergenic_background_to_',rest_chrom[c],'_1M_bins_PETs_num.txt'),header=F,sep="\t")
      line1_pairs_add_p=NULL
      for (j in 1:dim(line1_info_spe)[1]){
        line1_ID=line1_info_spe[j,4]
        idx=(line1_pairs_spe2[,2]>line1_info_spe[j,2])&(line1_pairs_spe2[,2]<line1_info_spe[j,3])
        line1_pairs_spe3=line1_pairs_spe2[idx,]
        if (dim(line1_pairs_spe3)[1]>0){
          bin_flag=ceiling(line1_pairs_spe3[,4]/K)
          bin_PET_num=table(bin_flag)
          uni_bin_flag=as.numeric(names(bin_PET_num))
          for (b in 1:length(uni_bin_flag)){
            p_value_temp=sum(bin_PET_num[b]<=background_PET_num[,uni_bin_flag[b]])/dim(background_PET_num)[1]
            line1_pairs_spe3[bin_flag==uni_bin_flag[b],5]=p_value_temp
          }
          line1_pairs_add_p=rbind(line1_pairs_add_p,cbind(line1_pairs_spe3,line1_ID))
        }
      }
    }
    line1_pairs_add_p_all=rbind(line1_pairs_add_p_all,line1_pairs_add_p)
  }
}
write.table(line1_pairs_add_p_all,paste0(name,'_binsize_1M_line1_inter_pairs_p_value.txt'),row.names=F,col.names=F,sep="\t",quote=F)
#########calculate the PET num between line1 and other chromosomes and compare with the intergenic background, calculate the p value.


