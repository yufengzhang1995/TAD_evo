######################################################
# In this file we used three methods: 
# 1) calculated correlations between genes within a TAD
# 2) randomly choose two genes and then calculated the correlations between them
# 3) Calculate the correlation between genes within TAD and genes outside TAD(the distance between the pairwise genes are sampled from the distance between genes that are both from the same TAD)
# summary statistics
# plot violin plot
# plot normal density plot
# use Kolmogorov-Smirnov Tests
##################### SET UP ##########################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicAlignments")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
genome <- BSgenome.Hsapiens.UCSC.hg38
seqlengths(genome)[["chr1"]]
BiocManager::install("annotate")
BiocManager::install("biomaRt")
install.packages("ggpubr")
install.packages("rdist")
install.packages("DescTools")
install.packages("sets")

##################### SET UP ##########################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(GenomicAlignments)
library(org.Hs.eg.db)
library(annotate)
library(biomaRt)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(DescTools)
library(sets)

##################### SET UP ##########################

# source functions
source("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex/functions.R")

# set working directory
setwd("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex")

# reference genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# get genes
hg38_genes <- genes(txdb)

# external TAD files
TAD_original <- read.csv("all_10000_blocks.bedpe",sep = "\t",stringsAsFactors=FALSE,header = TRUE)



# preprocess
TAD <- TAD_original[,c(1:3)]
#TAD <- rescale_resolution(TAD,10000)
colnames(TAD) <- c("chrom","tad_start","tad_end")
TAD$chrom <- paste("chr",TAD$chrom,sep = "")
TAD <- TAD[- grep("#",TAD$chrom),]
View(TAD)
nrow(TAD)
# transform the data into GRange objects
tad_grange <- makeGRangesFromDataFrame(TAD,ignore.strand = T,seqnames.field = "chrom",start.field = "tad_start",end.field = "tad_end")
head(tad_grange)
######################################################
######################################################
#### The step has been processed so we could skip this ####

# This step used awk; should not be processed in R
awk 'BEGIN{OFS="\t";FS="\t"}{
if(NR==1){
    for(i=1;i<NF;i++){
        if($i=="SMTS")target_id=i;
    }
    print $0;
}
if(NR>1 && $target_id=="Colon")print $0;
}' GTEx_v7_Annotations_SampleAttributesDS.txt> GTEx_Data_V6_Anno_colon.tsv

# extract gene expressions in different tissues and read them out into tsv files#

# whole seq data
GTExFile_path<-"GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct" # whole seq data

# get tissue name
tissue_name_whole <- c("Adipose_Tissue","Blood_Vessel","Bone_Marrow","Blood","Brain","Breast",
                 "Cervix_Uteri","Colon", "Esophagus","Fallopian_Tube","Heart","Kidney",
                 "Muscle","Nerve","Pituitary","Prostate","Salivary_Gland","Skin","Stomach","Testis",
                 "Thyroid","Uterus","Vagina")
               
# gene expression for 20 different tissues
readdata <- function(tn){
  fn = paste("GTEx_Data_V6_Anno_",tn,".tsv",sep = "")
  dt_temp <- fread(fn,sep = "\t",header = T)
  tissue_ID <- dt_temp[,SAMPID]
  tissue_ID <- append(c("Name","Description"),tissue_ID)
  GTEx_reads <- fread(GTExFile_path,sep = "\t",header =T,skip = 2,select = tissue_ID)
  fwrite(GTEx_reads,paste("GTEx_Analysis_v6p_RNA_",tn,".tsv",sep = ""),sep = "\t")
}
list_file <- lapply(tissue_name,readdata)

######################################################
######################################################

# set working directory
setwd("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex")
setwd("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex/match")

# read in the files we export in the last step
temp_whole = list.files(pattern="*.tsv")
temp = list.files(pattern="*.tsv")

# get tissue name
tissue_name_whole <- c("Adipose_Tissue","Blood_Vessel","Bone_Marrow","Blood","Brain","Breast",
                       "Cervix_Uteri","Colon", "Esophagus","Fallopian_Tube","Heart","Kidney",
                       "Muscle","Nerve","Pituitary","Prostate","Salivary_Gland","Skin","Stomach","Testis",
                       "Thyroid","Uterus","Vagina")
tissue_name <- c("Adrenal_Gland","Bladder","Liver","Lung","Ovary","Pancreas",
                 "Small_Intestine","Spleen")



data_name_whole <- paste(tissue_name_whole,"_exp",sep = "")
data_name <- paste(tissue_name,"_exp",sep = "")
for (i in 1:length(temp_whole)) assign(data_name_whole[i], read.csv(temp_whole[i],sep = "\t",header = TRUE))
for (i in 1:length(temp)) assign(data_name[i], read.csv(temp[i],sep = "\t",header = TRUE))

# list all the gene expression files  
data_list_whole <- list(Adipose_Tissue_exp,Blood_Vessel_exp,Bone_Marrow_exp,Blood_exp,Brain_exp,Breast_exp,
                        Cervix_Uteri_exp,Colon_exp,Esophagus_exp,Fallopian_Tube_exp,Heart_exp,Kidney_exp,
                        Muscle_exp,Nerve_exp,Pituitary_exp,Prostate_exp,Salivary_Gland_exp,Skin_exp,Stomach_exp,
                        Testis_exp,Thyroid_exp,Uterus_exp,Vagina_exp,Adrenal_Gland_exp,Bladder_exp,Liver_exp,
                        Lung_exp,Ovary_exp,Pancreas_exp,Small_Intestine_exp,Spleen_exp)
# match GTEx data to Hi-c data in terms of tissue type(to do) #########
data_list <- list(Adrenal_Gland_exp,Bladder_exp,Liver_exp,Lung_exp,Ovary_exp,Pancreas_exp,Small_Intestine_exp,Spleen_exp)
length(data_list)
######################################################

# find genes within TAD
overlap_df <- overlap_human(tad_grange,TAD)
View(overlap_df)

# add gene symbol to overlap_df
overlap_df_2 <- get_overlap_df_2(overlap_df)

# count the number of genes within every TAD
overlap_tbl <- get_overlap_tbl(overlap_df_2)
View(overlap_tbl)

######################################################

# get correlation list
TAD_correlation <- correlation_list(overlap_tbl,overlap_df_2,data_list_whole)
# save
write.csv(TAD_correlation, file = "TAD_correlation.csv")
# convert to data frame
TAD_correlation_f <- vec_to_frame(TAD_correlation)

# summary statistics of correlation value list
summary(TAD_correlation)

# plot
boxplot(TAD_correlation)
ggplot(TAD_correlation_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "TAD",y = "correlation value") 
ggdensity(TAD_correlation_f$vec)
View(Bladder_exp)
######################################################

# random method 1
random_1_cor <- random_2_74477_cor(data_list)
random_1_cor_f <- vec_to_frame(random_1_cor)
# save
write.csv(random_1_cor, file = "random_1_cor.csv")
# summary statistics
summary(random_1_cor)

# plot
ggplot(random_1_cor_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "random method 1",y = "correlation value") 
ggdensity(random_1_cor_f$vec)

####################################################

# random method 2
random_2_cor = const_distance_correlation(overlap_tbl,overlap_df_2,data_list_whole)
random_2_cor_f <- vec_to_frame(random_2_cor)
# save
write.csv(random_2_cor, file = "random_2_cor.csv")
# summary statistics
summary(random_2_cor)

# plot
ggplot(random_2_cor_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "random method 2",y = "correlation value") 
ggdensity(random_2_cor_f$vec)

####################################################

# non-parametric test
ks.test(TAD_correlation ,random_1_cor)
ks.test(TAD_correlation,random_2_cor)


####################################################
####################################################
####################################################
# randomly partition the genome into 1931 blocks
# absolute length
# absolute_loc <- list()
# sum = 0 
# for(i in my_chr){
#  sum = sum + my_chr_size[[i]]
#  absolute_loc[[i]] = sum
# }
# absolute_loc_vec <- data.frame(as.matrix(unlist(absolute_loc)))
# colnames(absolute_loc_vec) <- c("aboslute_length")
# chr_size <- data.frame(as.matrix(unlist(my_chr_size)))
# chr_size$chr_name <- row.names(chr_size)
# colnames(chr_size) <- c("size")
# absolute_loc_vec$start <- unname(absolute_loc_vec$aboslute_length) - unname(chr_size$size) + 1
# absolute_loc_vec$end <- absolute_loc_vec$aboslute_length
# absolute_loc_vec$chr_name <- row.names(chr_size)
# write.table(absolute_loc_vec,"chrome_length.csv",sep = ",")

# read in csv processed before
# chrome_length <- read.csv("chrome_length.csv",header = TRUE, sep = ",")
# head(chrome_length)
# get the TAD size
# TAD_size <- TAD$tad_end - TAD$tad_start # for future random shuffling
# ggdensity(TAD_size)
########################################
# get the chromosome sizes for every chromosome
# my_chr <- c(1:22,'X','Y')
# my_chr <- gsub(pattern="^", replacement='chr', my_chr)
# my_chr_size <- list()
# for (i in my_chr){
#   my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg38[[i]])
# }
########################################
# whole_length <- sum(unname(unlist(my_chr_size)))
# set.seed(1984)
# random_start <- runif(1931, min = 0, max = 3088269832 - max(TAD_size))
# random_end <- random_start + sample(TAD_size)
# random_tbl <- data.frame(random_start,random_end)
# chr_name <- c()
# pseudo_start <- c()
# pseudo_end <- c() 
# for(i in 1:1931){
#   for(j in 1:24){
#     if(random_start[i] > absolute_loc_vec$start[j] & random_end[i]  < absolute_loc_vec$end[j]){
#       chr_name <- append(chr_name,absolute_loc_vec$chr_name[j])
#       pseudo_start <- append(pseudo_start,random_start[i])
#       pseudo_end <- append(pseudo_end, random_end[i])
#     }
#   }
# }
# length(pseudo_start)
# random_block <- data.frame(chr_name,pseudo_start,pseudo_end)
# random_block_2 <- merge(random_block,chrome_length, by = "chr_name")
# random_block_2$tad_start <- random_block_2$pseudo_start - random_block_2$start + 1
# random_block_2$tad_end <- random_block_2$pseudo_end - random_block_2$start + 1
# random_round <- apply(random_block_2[,c(2,3,7,8)],2,round)
# chromosome_name <- as.character(random_block_2[,c("chr_name")])
# random_block_3 <- cbind(random_round,chromosome_name)
# random_block_3 <- data.frame(random_block_3[,3:5])
# random_grange <- makeGRangesFromDataFrame(random_block_3,ignore.strand = T,
#                                           seqnames.field = "chromosome_name",
#                                           start.field = "tad_start",end.field = "tad_end")
# overlap_random <- overlap_human(random_grange,random_block_3)
# overlap_random_2 <- get_overlap_df_2(overlap_random)
# overlap_random_tbl <- get_overlap_tbl(overlap_random_2)
# colnames(overlap_random_tbl) <- c("seqnames","block_start","block_end")


# random_cor_li <- correlation_list(overlap_random_tbl,overlap_random_2)
# random_cor_fram <- vec_to_frame(random_cor_li)
# boxplot(random_cor_li)
# ggplot(random_cor_fram ,aes(x = id,y = random_cor_li)) +
#   geom_violin(fill="gray") + geom_boxplot(width=0.1)
# median(random_cor_li)
# ggdensity(random_cor_fram$random_cor_li)



# ======================================================
# tissue-specific TAD boundaries

# set working directory
setwd("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex")

# read in the files we export in the last step
tad_files = list.files(pattern="*.bed")

# 

TAD_name <- c("Adrenal_Gland_TAD","Bladder_TAD", "Lung1_TAD", "Lung2_TAD" , "Lung3_TAD",
              "Liver_TAD",  "Ovary_TAD"  ,  "Pancreas1_TAD", "Pancreas2_TAD" ,"Pancreas3_TAD",
              "Small_Intestine_TAD","Spleen1_TAD","Spleen2_TAD","Spleen3_TAD")
for (i in 1:length(tad_files)) assign(TAD_name[i], read.csv(tad_files[i],sep = "\t",header = FALSE)) 

TAD_list <- list(Adrenal_Gland_TAD,Bladder_TAD,Lung1_TAD,Lung2_TAD,Lung3_TAD,
              Liver_TAD,Ovary_TAD,Pancreas1_TAD,Pancreas2_TAD,Pancreas3_TAD,
              Small_Intestine_TAD,Spleen1_TAD,Spleen2_TAD,Spleen3_TAD)
TAD_transform_list <- lapply(TAD_list,TAD_transform)
grange_list <- lapply(TAD_transform_list, TAD_to_grange)

lapply(TAD_list,nrow)

#==================================================
# Adreanal Gland
# find genes within TAD
Adrenal_Gland_transform  <- TAD_transform(Adrenal_Gland_TAD)
head(Adrenal_Gland_transform)

Adrenal_Gland_grange <- TAD_to_grange(Adrenal_Gland_transform)
head(Adrenal_Gland_grange)

Adrenal_Gland_overlap_df <- overlap_human(Adrenal_Gland_grange,Adrenal_Gland_transform)
dim(Adrenal_Gland_overlap_df)
View(Adrenal_Gland_TAD)
Adrenal_Gland_TAD %>% 
  filter(V1 == "chr1") %>%
  arrange(desc(V2)) %>%
  head()
TAD %>% 
  filter(chrom == "chr1") %>%
  arrange(desc(tad_start)) %>%
  head()
View(TAD)

Adrenal_correlation_list <- correlation_tissue_specific(Adrenal_Gland_grange,Adrenal_Gland_transform,Adrenal_Gland_exp)





Adrenal_Gland_result <- pre_procedure(Adrenal_Gland_grange,Adrenal_Gland_transform)
summary(unlist(Adrenal_Gland_correlation_list[3]))
length(unlist(Adrenal_Gland_correlation_list[3]))
Adrenal_Gland_correlation_list_f <- vec_to_frame(unlist(Adrenal_Gland_correlation_list[3]))
# save
write.csv(Adrenal_Gland_correlation_list[3], file = "Adrenal_Gland_correlation_list.csv")
# plot
ggplot(Adrenal_Gland_correlation_list_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "AG",y = "correlation value") 
ggdensity(Adrenal_Gland_correlation_list_f$vec)



# use hECS TAD
AG_hECS_correlation_list <-correlation_tissue_specific(tad_grange,TAD,Adrenal_Gland_exp)
AG_hECS_correlation_list_f <- vec_to_frame(AG_hECS_correlation_list)
summary(AG_hECS_correlation_list)
length(AG_hECS_correlation_list)
ggplot(AG_hECS_correlation_list_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "AG-hESC",y = "correlation value") 
ggdensity(AG_hECS_correlation_list_f$vec)
# save
write.csv(AG_hECS_correlation_list, file = "AG_hESC_correlation_list.csv")




## random block
Adrenal_Gland_correlation_block <- const_distance_correlation_ts(
  data.frame(Adrenal_Gland_result[2]),
  data.frame(Adrenal_Gland_result[1]),
  Adrenal_Gland_exp)
Adrenal_correlation_block <- Adrenal_Gland_correlation_block
summary(Adrenal_Gland_correlation_block)
length(Adrenal_Gland_correlation_block)
Adrenal_Gland_correlation_block_f <- vec_to_frame(Adrenal_Gland_correlation_block)
# save
write.csv(Adrenal_Gland_correlation_block, file = "Adrenal_Gland_correlation_block.csv")
# use hESC
AG_hESC_correlation_block <- const_distance_correlation_ts(
  data.frame(overlap_tbl),
  data.frame(overlap_df_2),
  Adrenal_Gland_exp)
summary(AG_hESC_correlation_block)
length(AG_hESC_correlation_block)
AG_hESC_correlation_block_f <- vec_to_frame(AG_hESC_correlation_block)
# save
write.csv(AG_hESC_correlation_block, file = "AG_hESC_correlation_block.csv")
# plot
ggplot(Adrenal_Gland_correlation_block_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "AG_block",y = "correlation value") 
ggdensity(Adrenal_Gland_correlation_block_f$vec)
# plot
ggplot(AG_hESC_correlation_block_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "AG_hESC_block",y = "correlation value") 
ggdensity(Adrenal_Gland_correlation_block_f$vec)
#==================================================
#==================================================
# Bladder
# find genes within TAD
Bladder_transform  <- TAD_transform(Bladder_TAD)
Bladder_grange <- TAD_to_grange(Bladder_transform)
Bladder_correlation_list <-correlation_tissue_specific(Bladder_grange,Bladder_transform,Bladder_exp)
write.csv(Bladder_correlation_list, file = "Bladder_correlation_list.csv")
summary(Bladder_correlation_list)
length(Bladder_correlation_list)
Bladder_correlation_list_f <- vec_to_frame(Bladder_correlation_list)
# plot
ggplot(Bladder_correlation_list_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "BL",y = "correlation value") 
ggdensity(Bladder_correlation_list_f$vec)

# random block
Bladder_result <- pre_procedure(Bladder_grange,Bladder_transform)
Bladder_correlation_block <- const_distance_correlation_ts(
  data.frame(Bladder_result[2]),
  data.frame(Bladder_result[1]),
  Bladder_exp)
write.csv(Bladder_correlation_block, file = "Bladder_correlation_block.csv")
length(Bladder_correlation_block)
Bladder_correlation_block_f <- vec_to_frame(Bladder_correlation_block)

#==================================================
#==================================================
# Lung1
# find genes within TAD
Lung1_transform  <- TAD_transform(Lung1_TAD)
Lung1_grange <- TAD_to_grange(Lung1_transform)
Lung1_correlation_list <-correlation_tissue_specific(Lung1_grange,Lung1_transform,Lung_exp)
write.csv(Lung1_correlation_list, file = "Lung1_correlation_list.csv")
length(Lung1_correlation_list)
Lung1_correlation_list_f <- vec_to_frame(Lung1_correlation_list)


# random block
Lung1_result <- pre_procedure(Lung1_grange,Lung1_transform)
Lung1_correlation_block <- const_distance_correlation_ts(
  data.frame(Lung1_result[2]),
  data.frame(Lung1_result[1]),
  Lung_exp)
write.csv(Lung1_correlation_block, file = "Lung1_correlation_block.csv")
length(Lung1_correlation_block)
Lung1_correlation_block_f <- vec_to_frame(Lung1_correlation_block)

#==================================================
#==================================================
# Lung2
# find genes within TAD
Lung2_transform  <- TAD_transform(Lung2_TAD)
Lung2_grange <- TAD_to_grange(Lung2_transform)
Lung2_correlation_list <-correlation_tissue_specific(Lung2_grange,Lung2_transform,Lung_exp)
write.csv(Lung2_correlation_list, file = "Lung2_correlation_list.csv")
length(Lung2_correlation_list)
Lung2_correlation_list_f <- vec_to_frame(Lung2_correlation_list)

# random block
Lung2_result <- pre_procedure(Lung2_grange,Lung2_transform)
Lung2_correlation_block <- const_distance_correlation_ts(
  data.frame(Lung2_result[2]),
  data.frame(Lung2_result[1]),
  Lung_exp)
write.csv(Lung2_correlation_block, file = "Lung2_correlation_block.csv")
length(Lung2_correlation_block)
Lung2_correlation_block_f <- vec_to_frame(Lung2_correlation_block)

#==================================================
#==================================================
# Lung3
# find genes within TAD
Lung3_transform  <- TAD_transform(Lung3_TAD)
Lung3_grange <- TAD_to_grange(Lung3_transform)
Lung3_correlation_list <-correlation_tissue_specific(Lung3_grange,Lung3_transform,Lung_exp)
write.csv(Lung3_correlation_list, file = "Lung3_correlation_list.csv")
summary(Lung3_correlation_list)
Lung3_correlation_list_f <- vec_to_frame(Lung3_correlation_list)


# random block
Lung3_result <- pre_procedure(Lung3_grange,Lung3_transform)
Lung3_correlation_block <- const_distance_correlation_ts(
  data.frame(Lung3_result[2]),
  data.frame(Lung3_result[1]),
  Lung_exp)
write.csv(Lung3_correlation_block, file = "Lung3_correlation_block.csv")
summary(Lung3_correlation_block)
Lung3_correlation_block_f <- vec_to_frame(Lung3_correlation_block)

#==================================================
#==================================================
# Liver
# find genes within TAD
Liver_transform  <- TAD_transform(Liver_TAD)
Liver_grange <- TAD_to_grange(Liver_transform)
Liver_correlation_list <-correlation_tissue_specific(Liver_grange,Liver_transform,Liver_exp)
write.csv(Liver_correlation_list, file = "Liver_correlation_list.csv")
summary(Liver_correlation_list)
Liver_correlation_list_f <- vec_to_frame(Liver_correlation_list)


# random block
Liver_result <- pre_procedure(Liver_grange,Liver_transform)
Liver_correlation_block <- const_distance_correlation_ts(
  data.frame(Liver_result[2]),
  data.frame(Liver_result[1]),
  Liver_exp)
write.csv(Liver_correlation_block, file = "Liver_correlation_block.csv")
summary(Liver_correlation_block)
Liver_correlation_block_f <- vec_to_frame(Liver_correlation_block)

#==================================================
#==================================================
# Liver
# find genes within TAD
Ovary_transform  <- TAD_transform(Ovary_TAD)
Ovary_grange <- TAD_to_grange(Ovary_transform)
Ovary_correlation_list <-correlation_tissue_specific(Ovary_grange,Ovary_transform,Ovary_exp)
write.csv(Ovary_correlation_list, file = "Ovary_correlation_list.csv")
summary(Ovary_correlation_list)
Ovary_correlation_list_f <- vec_to_frame(Ovary_correlation_list)


# random block
Ovary_result <- pre_procedure(Ovary_grange,Ovary_transform)
Ovary_correlation_block <- const_distance_correlation_ts(
  data.frame(Ovary_result[2]),
  data.frame(Ovary_result[1]),
  Ovary_exp)
write.csv(Ovary_correlation_block, file = "Ovary_correlation_block.csv")
summary(Ovary_correlation_block)
Ovary_correlation_block_f <- vec_to_frame(Ovary_correlation_block)

#==================================================
#==================================================
# Pancreas1
# find genes within TAD
Pancreas1_transform  <- TAD_transform(Pancreas1_TAD)
Pancreas1_grange <- TAD_to_grange(Pancreas1_transform)
Pancreas1_correlation_list <-correlation_tissue_specific(Pancreas1_grange,Pancreas1_transform,Pancreas_exp)
write.csv(Pancreas1_correlation_list, file = "Pancreas1_correlation_list.csv")
summary(Pancreas1_correlation_list)
Pancreas1_correlation_list_f <- vec_to_frame(Pancreas1_correlation_list)


# random block
Pancreas1_result <- pre_procedure(Pancreas1_grange,Pancreas1_transform)
Pancreas1_correlation_block <- const_distance_correlation_ts(
  data.frame(Pancreas1_result[2]),
  data.frame(Pancreas1_result[1]),
  Pancreas_exp)
write.csv(Pancreas1_correlation_block, file = "Pancreas1_correlation_block.csv")
summary(Pancreas1_correlation_block)
Pancreas1_correlation_block_f <- vec_to_frame(Pancreas1_correlation_block)

#==================================================
# Pancreas2
# find genes within TAD
Pancreas2_transform  <- TAD_transform(Pancreas2_TAD)
Pancreas2_grange <- TAD_to_grange(Pancreas2_transform)
Pancreas2_correlation_list <-correlation_tissue_specific(Pancreas2_grange,Pancreas2_transform,Pancreas_exp)
write.csv(Pancreas2_correlation_list, file = "Pancreas2_correlation_list.csv")
length(Pancreas2_correlation_list)
Pancreas2_correlation_list_f <- vec_to_frame(Pancreas2_correlation_list)


# random block
Pancreas2_result <- pre_procedure(Pancreas2_grange,Pancreas2_transform)
Pancreas2_correlation_block <- const_distance_correlation_ts(
  data.frame(Pancreas2_result[2]),
  data.frame(Pancreas2_result[1]),
  Pancreas_exp)
write.csv(Pancreas2_correlation_block, file = "Pancreas2_correlation_block.csv")
summary(Pancreas2_correlation_block)
Pancreas2_correlation_block_f <- vec_to_frame(Pancreas2_correlation_block)

#==================================================
# Pancreas3
# find genes within TAD
Pancreas3_transform  <- TAD_transform(Pancreas3_TAD)
Pancreas3_grange <- TAD_to_grange(Pancreas3_transform)
Pancreas3_correlation_list <-correlation_tissue_specific(Pancreas3_grange,Pancreas3_transform,Pancreas_exp)
write.csv(Pancreas3_correlation_list, file = "Pancreas3_correlation_list.csv")
summary(Pancreas3_correlation_list)
Pancreas3_correlation_list_f <- vec_to_frame(Pancreas3_correlation_list)


# random block
Pancreas3_result <- pre_procedure(Pancreas3_grange,Pancreas3_transform)
Pancreas3_correlation_block <- const_distance_correlation_ts(
  data.frame(Pancreas3_result[2]),
  data.frame(Pancreas3_result[1]),
  Pancreas_exp)
write.csv(Pancreas3_correlation_block, file = "Pancreas3_correlation_block.csv")
summary(Pancreas3_correlation_block)
Pancreas3_correlation_block_f <- vec_to_frame(Pancreas3_correlation_block)

#==================================================
# Small_Intestine
# find genes within TAD
Small_Intestine_transform  <- TAD_transform(Small_Intestine_TAD)
Small_Intestine_grange <- TAD_to_grange(Small_Intestine_transform)
Small_correlation_list <-correlation_tissue_specific(Small_Intestine_grange,
                                                               Small_Intestine_transform,
                                                               Small_Intestine_exp)
write.csv(Small_correlation_list, file = "Small_correlation_list.csv")
summary(Small_correlation_list)
Small_Intestine_correlation_list_f <- vec_to_frame(Small_Intestine_correlation_list)


# random block
Small_Intestine_result <- pre_procedure(Small_Intestine_grange,Small_Intestine_transform)
Small_Intestine_correlation_block <- const_distance_correlation_ts(
  data.frame(Small_Intestine_result[2]),
  data.frame(Small_Intestine_result[1]),
  Small_Intestine_exp)
Small_correlation_block <- Small_Intestine_correlation_block
write.csv(Small_Intestine_correlation_block, file = "Small_Intestine_correlation_block.csv")
summary(Small_Intestine_correlation_block)
Small_Intestine_correlation_block_f <- vec_to_frame(Small_Intestine_correlation_block)

#==================================================
# Spleen1
# find genes within TAD
Spleen1_transform  <- TAD_transform(Spleen1_TAD)
Spleen1_grange <- TAD_to_grange(Spleen1_transform)
Spleen1_correlation_list <-correlation_tissue_specific(Spleen1_grange,Spleen1_transform,Spleen_exp)
write.csv(Spleen1_correlation_list, file = "Spleen1_correlation_list.csv")
write.csv(Spleen2_correlation_list, file = "Spleen2_correlation_list.csv")
write.csv(Spleen3_correlation_list, file = "Spleen3_correlation_list.csv")
length(Spleen1_correlation_list)
Pancreas1_correlation_list_f <- vec_to_frame(Spleen1_correlation_list)


# random block
Spleen1_result <- pre_procedure(Spleen1_grange,Spleen1_transform)
Spleen1_correlation_block <- const_distance_correlation_ts(
  data.frame(Spleen1_result[2]),
  data.frame(Spleen1_result[1]),
  Spleen_exp)
write.csv(Spleen1_correlation_block, file = "Spleen1_correlation_block.csv")
write.csv(Spleen2_correlation_block, file = "Spleen2_correlation_block.csv")
write.csv(Spleen3_correlation_block, file = "Spleen3_correlation_block.csv")
summary(Spleen1_correlation_block)
Spleen1_correlation_block_f <- vec_to_frame(Spleen1_correlation_block)

#==================================================
# Spleen2
# find genes within TAD
Spleen2_transform  <- TAD_transform(Spleen2_TAD)
Spleen2_grange <- TAD_to_grange(Spleen2_transform)
Spleen2_correlation_list <-correlation_tissue_specific(Spleen2_grange,Spleen2_transform,Spleen_exp)
summary(Spleen2_correlation_list)
Spleen2_correlation_list_f <- vec_to_frame(Spleen2_correlation_list)

# random block
Spleen2_result <- pre_procedure(Spleen2_grange,Spleen2_transform)
Spleen2_correlation_block <- const_distance_correlation_ts(
  data.frame(Spleen2_result[2]),
  data.frame(Spleen2_result[1]),
  Spleen_exp)
summary(Spleen2_correlation_block)
Spleen2_correlation_block_f <- vec_to_frame(Spleen2_correlation_block)

#==================================================
# Spleen3
# find genes within TAD
Spleen3_transform  <- TAD_transform(Spleen3_TAD)
Spleen3_grange <- TAD_to_grange(Spleen3_transform)
Spleen3_correlation_list <-correlation_tissue_specific(Spleen3_grange,Spleen3_transform,Spleen_exp)
summary(Spleen3_correlation_list)
Spleen3_correlation_list_f <- vec_to_frame(Spleen3_correlation_list)

# random block
Spleen3_result <- pre_procedure(Spleen3_grange,Spleen3_transform)
Spleen3_correlation_block <- const_distance_correlation_ts(
  data.frame(Spleen3_result[2]),
  data.frame(Spleen3_result[1]),
  Spleen_exp)
length(Spleen3_correlation_block)
Spleen3_correlation_block_f <- vec_to_frame(Spleen3_correlation_block)

#==================================================
correlation_tissue_sum <- get(paste(str_trim(lapply(str_split(TAD_name,"_"),function(y) y[1])),"_correlation_list",
                                sep = ''))
correlation_tissue_block_sum <- get(paste(str_trim(lapply(str_split(TAD_name,"_"),function(y) y[1])),"_correlation_block",
                                    sep = ''))
ks.test(correlation_tissue_sum,correlation_tissue_block_sum)
correlation_tissue_sum_f$group = c("tad")
correlation_tissue_block_sum_f$group = c("random_block")
corr_whole <- rbind(correlation_tissue_sum_f,correlation_tissue_block_sum_f)
ggplot(corr_whole,aes(x = vec,color = group)) + geom_density


correlation_tissue_sum_f <- vec_to_frame(correlation_tissue_sum)
correlation_tissue_block_sum_f <- vec_to_frame(correlation_tissue_block_sum)
ggplot(correlation_tissue_sum_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "whole_TAD",y = "correlation value") 
ggplot(correlation_tissue_block_sum_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  labs(x = "whole_random_block",y = "correlation value") 




ggdensity(Adrenal_Gland_correlation_list_f$vec)

summary(correlation_tissue_sum)
# 

TAD_name <- c("Adrenal_Gland_TAD","Bladder_TAD", "Lung1_TAD", "Lung2_TAD" , "Lung3_TAD",
              "Liver_TAD",  "Ovary_TAD"  ,  "Pancreas1_TAD", "Pancreas2_TAD" ,"Pancreas3_TAD",
              "Small_Intestine_TAD","Spleen1_TAD","Spleen2_TAD","Spleen3_TAD")
df_ts_cor <- data.frame('AG' = c(median(Adrenal_correlation_list),length(Adrenal_correlation_list)),
                        'BL' = c(median(Bladder_correlation_list),length(Bladder_correlation_list)),
                        'LG1' = c(median(Lung1_correlation_list),length(Lung1_correlation_list)),
                        'LG2' = c(median(Lung2_correlation_list),length(Lung2_correlation_list)),
                        'LG3' = c(median(Lung3_correlation_list),length(Lung3_correlation_list)),
                        'LI' = c(median(Liver_correlation_list),length(Liver_correlation_list)),
                        'OV' = c(median(Ovary_correlation_list),length(Ovary_correlation_list)),
                        'PA1' = c(median(Pancreas1_correlation_list),length(Pancreas1_correlation_list)),
                        'PA2' = c(median(Pancreas2_correlation_list),length(Pancreas2_correlation_list)),
                        'PA3' = c(median(Pancreas3_correlation_list),length(Pancreas3_correlation_list)),
                        'SI' = c(median(Small_correlation_list),length(Small_correlation_list)),
                        'SX1' = c(median(Spleen1_correlation_list),length(Spleen1_correlation_list)),
                        'SX2' = c(median(Spleen2_correlation_list),length(Spleen2_correlation_list)),
                        'SX3' = c(median(Spleen3_correlation_list),length(Spleen3_correlation_list)))
df_ts_cor_block <- data.frame('AG' = c(median(Adrenal_Gland_correlation_block),length(Adrenal_Gland_correlation_block)),
                        'BL' = c(median(Bladder_correlation_block),length(Bladder_correlation_block)),
                        'LG1' = c(median(Lung1_correlation_block),length(Lung1_correlation_block)),
                        'LG2' = c(median(Lung2_correlation_block),length(Lung2_correlation_block)),
                        'LG3' = c(median(Lung3_correlation_block),length(Lung3_correlation_block)),
                        'LI' = c(median(Liver_correlation_block),length(Liver_correlation_block)),
                        'OV' = c(median(Ovary_correlation_block),length(Ovary_correlation_block)),
                        'PA1' = c(median(Pancreas1_correlation_block),length(Pancreas1_correlation_block)),
                        'PA2' = c(median(Pancreas2_correlation_block),length(Pancreas2_correlation_block)),
                        'PA3' = c(median(Pancreas3_correlation_block),length(Pancreas3_correlation_block)),
                        'SI' = c(median(Small_Intestine_correlation_block),length(Small_Intestine_correlation_block)),
                        'SX1' = c(median(Spleen1_correlation_block),length(Spleen1_correlation_block)),
                        'SX2' = c(median(Spleen2_correlation_block),length(Spleen2_correlation_block)),
                        'SX3' = c(median(Spleen3_correlation_block),length(Spleen3_correlation_block)))
df_ts_cor_2 <- t(as.matrix(df_ts_cor))
colnames(df_ts_cor_2) <- c("median","pairs")
df_ts_cor_2 <- data.frame(df_ts_cor_2)
setDT(df_ts_cor_2, keep.rownames = TRUE)[]

df_ts_cor_block_2 <- t(as.matrix(df_ts_cor_block))
colnames(df_ts_cor_block_2) <- c("median","pairs")
df_ts_cor_block_2 <- data.frame(df_ts_cor_block_2)
setDT(df_ts_cor_block_2, keep.rownames = TRUE)[]

df_ts_cor_2$group <- rep("TAD",nrow(df_ts_cor_2))
df_ts_cor_block_2$group <- rep("random_block",nrow(df_ts_cor_block_2))
df_ts_whole <- rbind(df_ts_cor_2,df_ts_cor_block_2)



ggplot(df_ts_whole ,aes(x = rn, y = pairs,col = group)) +
  geom_bar(stat="identity",position=position_dodge()) +
  labs(x = "tissue",y = "number of gene pairs") 

ggplot(df_ts_cor_block_2 ,aes(x = rn, y = median)) +
  geom_bar(stat="identity") +
  labs(x = "tissue",y = "median correlation value") 


ggplot(df_ts_whole,aes(x = pairs,
                     y = median)) + geom_point() +
  geom_smooth(method = "lm",se = FALSE)
model1 <- lm(data = df_ts_whole,median ~ pairs)
summary(model1)

x = as.numeric(TAD[1,c(2,3)])
y = as.numeric(TAD[2,c(2,3)])
head(TAD)


Overlap(x,y)

#===========================================
# calculate distance between two genes by gene ID
# load('/Volumes/GoogleDrive/My\ Drive/summer/data/GTex/organize_correlation.RData')
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene_list <- unlist(lapply(str_split(Adrenal_Gland_exp$Name,"\\."), function(y) y[1]))

not_same_chr = 0
for(i in 1:1000){
  sample_ensembl_id <- sample(gene_list,2,replace = FALSE)
  gene_info <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('ensembl_gene_id'),
      values=c(sample_ensembl_id),
      mart=ensembl)
  if(nrow(gene_info) == 2){
    if(gene_info[1,1] != gene_info[1,2]){
      not_same_chr = not_same_chr +1
    }
  } 
}

not_same_chr/1000
all_boundaries <- TAD_list
all_boundaries[15] <- list(TAD)
all_tad_start <- lapply(seq(1,14,1),extract_start,chr_n = "chr1")
# inspect the TAD boundaries on chromosome 1
plot(density(extract_start(1,"chr6")))
lines(density(extract_start(2,"chr6")))
lines(density(extract_start(3,"chr6")))
lines(density(extract_start(4,"chr6")))
lines(density(extract_start(5,"chr6")))
lines(density(extract_start(6,"chr6")))
lines(density(extract_start(7,"chr6")))
lines(density(extract_start(8,"chr6")))
lines(density(extract_start(9,"chr6")))
lines(density(extract_start(10,"chr6")))
lines(density(extract_start(11,"chr6")))
lines(density(extract_start(12,"chr6")))
lines(density(extract_start(13,"chr6")))
lines(density(extract_start(14,"chr6")))
lines(density(as.numeric(data.frame(all_boundaries[15]) %>%
                filter(chrom == "chr6") %>% "$"(tad_start))))

plot(density(extract_end(1,"chr1")))
lines(density(extract_end(2,"chr1")))
lines(density(extract_end(3,"chr1")))
lines(density(extract_end(4,"chr1")))
lines(density(extract_end(5,"chr1")))
lines(density(extract_end(6,"chr1")))
lines(density(extract_end(7,"chr1")))
lines(density(extract_end(8,"chr1")))
lines(density(extract_end(9,"chr1")))
lines(density(extract_end(10,"chr1")))
lines(density(extract_end(11,"chr1")))
lines(density(extract_end(12,"chr1")))
lines(density(extract_end(13,"chr1")))
lines(density(extract_end(14,"chr1")))

          
                           


# correlation between distance
# conserved acorss tissue: tad mean anything(TAD shared among tissues would be important)
# average gene expressione check the correlation between different tissue(average the gene expression )
# find shared tad tissue (how to define the shared tad boundareis: how many percent could we define it is shared)
### # conserved acorss tissue: tad mean anything(TAD shared among tissues would be important)
# change definition of distance
# Rao 2014: different tissue Hi-C
######################################################
######################################################

# correlation between genes within TAD (do average over samples)
tad_cor_mean <- correlation_list_mean(overlap_tbl,overlap_df_2,data_list_whole)

#=======================================================
# Plot histogram for differen tad size
hist_TAD <- data.frame(as.matrix(as.numeric(TAD$tad_end)-as.numeric(TAD$tad_start)))
colnames(hist_TAD) <- c("size")
hist_TAD$group <- rep(c("TAD"),nrow(hist_TAD))
hist_Bladder <- different_TAD_size(Bladder_TAD,"Bladder_TAD")
hist_Liver <-  different_TAD_size(Liver_TAD,"Liver_TAD")
hist_Lung1 <-  different_TAD_size(Lung1_TAD,"Lung1_TAD")
hist_Lung2 <-  different_TAD_size(Lung2_TAD,"Lung2_TAD")
hist_Lung3 <-  different_TAD_size(Lung3_TAD,"Lung3_TAD")
hist_Ovary <-  different_TAD_size(Ovary_TAD,"Ovary_TAD")
hist_Pancreas1 <-  different_TAD_size(Pancreas1_TAD,"Pancreas1_TAD")
hist_Pancreas2 <-  different_TAD_size(Pancreas2_TAD,"Pancreas2_TAD")
hist_Pancreas3 <-  different_TAD_size(Pancreas3_TAD,"Pancreas3_TAD")
hist_SX <-  different_TAD_size(Small_Intestine_TAD,"SX_TAD")
hist_Spleen1 <- different_TAD_size(Spleen1_TAD,"Spleen1_TAD")
hist_Spleen2 <- different_TAD_size(Spleen2_TAD,"Spleen2_TAD")
hist_Spleen3 <- different_TAD_size(Spleen3_TAD,"Spleen3_TAD")
hist_AG <- different_TAD_size(Adrenal_Gland_TAD,"AG_TAD")
hist_list <- list(hist_Bladder,hist_Liver,hist_Lung1,hist_Lung2,hist_Lung3,hist_Ovary,
                  hist_Pancreas1,hist_Pancreas2,hist_Pancreas3,hist_SX,hist_Spleen1,hist_Spleen2,
                  hist_Spleen3,hist_AG)
hist_all <- bind_rows(hist_list)

ggplot(hist_all ,aes(x = size,color = group)) +
  geom_histogram(fill = "white",alpha = 0.5,position = "identity",binwidth = 10000) +
  labs(x = "size",y = "freq") 
ggplot(hist_TAD ,aes(x = size)) +
  geom_histogram(position = "identity") +
  labs(x = "size",y = "freq") 

#===========================================================
# define convserved if the TAD interval is overlapped by 90%
Overlap(x,y)
nrow(Spleen1_TAD)
#===========================================================
# calculate the gene co-expression correlation varied by distance

gene_distance_correlation <- function(){
  for(ch in chromosome_name){
    df <- data.frame(ch,c(1),chromosome_length[[ch]])
    colnames(df) <- c("seqnames","chr_start","chr_end")
    df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end")
    df_1 <- overlap_human(df_grange,TAD)
    df_2 <- get_overlap_df_2(df_1)
    df_tbl <- get_overlap_tbl(df_2)
  }
}

###################################
# only for reference #
correlation_list <- function(overlap_tbl,overlap_df_2,data_list){
  # initialization
  cor_li <- c()
  for(i in 1:nrow(overlap_tbl)){
    # specify a TAD range and find the genes within the TAD
    d <- overlap_df_2 %>%
      filter(seqnames == overlap_tbl[i,1] & tad_start == overlap_tbl[i,2] & tad_end == overlap_tbl[i,3])
    # get the the position of genes of the TAD in gene_expression data frame
    sample_id <- as.numeric(match(d$gene_symbl,Liver_exp$Description))
    # remove NA
    sample_id <- sample_id[!is.na(sample_id)]
    # apply to the all gene expression data frame
    res <- sapply(data_list,extract_generow,sample_id = sample_id)
    # bind these gene expression by column
    mean_all <- bind_cols(res)
    # throw away the first column and transpose it
    trans_mat <- data.matrix(t(mean_all[,c(-1)]))
    # set the rownames to null
    rownames(trans_mat) <- c()
    # convert the data type of every element to numeric and construct a data frame
    trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
    # get correlation 
    cor_mat <- cor(trans_mat)
    # make the correlation value in the matrix a long vector
    cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
    cor_mat[cor_mat_tri] <- NA
    corr_matrix <- melt(cor_mat,na.rm = TRUE)
    y = corr_matrix$value
    # add it to the list
    cor_li <- append(cor_li,y) 
  }
  return(cor_li)
}




