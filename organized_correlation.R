######################################################
# In this file we used three methods: 
# 1) calculated correlations between genes within a TAD
# 2) randomly choose two genes and then calculated the correlations between them
# 3) randomly break the sequences into bins and calculted the correlation between genes within a bin
# plot normal density plot
# use Kolmogorov-Smirnov Tests
######################################################

#Set up
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicAlignments")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#install.packages("ggpubr")
######################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(GenomicAlignments)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(rhdf5)
library(ggpubr)


#####################################################
######################################################
# working directory
setwd("/Volumes/GoogleDrive/My\ Drive/summer/data/GTex")
# reference genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# external TAD files
TAD <- read.csv("5000_blocks.bedpe",sep = "\t",header = TRUE)
TAD <- TAD[,c(1:3)]
colnames(TAD) <- c("chrom","tad_start","tad_end")
TAD$chrom <- paste("chr",TAD$chrom,sep = "")
# transform the data into GRange objects
tad_grange <- makeGRangesFromDataFrame(TAD,ignore.strand = T,seqnames.field = "chrom",start.field = "tad_start",end.field = "tad_end")
######################################################
# extract gene expressions in different tissues and read them out into tsv files#
# whole seq data
GTExFile_path<-"GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct" # whole seq data
# get tissue name
ANNOTATION <- read.csv("GTEx_v7_Annotations_SampleAttributesDS.txt",
                       sep = "\t",header = TRUE)
tissue_name <- c("Blood","Adipose_Tissue", "Muscle", "Blood_Vessel","Heart","Ovary","Uterus",
                 "Vagina","Breast","Skin","Salivary_Gland","Brain","Adrenal_Gland", "Thyroid",
                 "Lung","Spleen", "Pancreas","Esophagus","Stomach", "Colon", "Small_Intestine", 
                 "Prostate", "Testis", "Nerve","Pituitary","Liver","Kidney", "Fallopian_Tube", 
                 "Bladder","Cervix_Uteri","Bone_Marrow")
               
# gene expression for 31 different tissues
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
# read in the files we export in the last step
temp = list.files(pattern="*.tsv")
data_name <- paste(tissue_name,"_exp",sep = "")
for (i in 1:length(temp)) assign(data_name[i], read.csv(temp[i],sep = "\t",header = TRUE))
# several patients so we take the mean expression for every gene
mean_exp <-function(x){
  d <- data.frame(ID = rownames(x[sample_id,-c(1,2)]),
                  mean = rowMeans(x[sample_id,-c(1,2)]))
  colnames(d) <- c("ID","tissue")
  return(d)
}
# list all the gene expression files
data_list <- list(Blood_exp,Adipose_Tissue_exp,Muscle_exp,Blood_Vessel_exp,Heart_exp,Ovary_exp,Uterus_exp,Vagina_exp,Breast_exp,Skin_exp,
Salivary_Gland_exp,Brain_exp,Adrenal_Gland_exp,Thyroid_exp,Lung_exp,Spleen_exp,Pancreas_exp,Esophagus_exp,Stomach_exp ,Colon_exp,
Small_Intestine_exp ,Prostate_exp,Testis_exp,Nerve_exp,Pituitary_exp,Liver_exp,Kidney_exp,Fallopian_Tube_exp,Bladder_exp,Cervix_Uteri_exp,Bone_Marrow_exp)
######################################################
# create a function to find overlaps
overlape_human <- function(tad_grange){
  overlap_tad_gene <- findOverlaps(hg38_genes,tad_grange,select = "all",type = "within")
  overlap_tad_gene <- as.data.frame(overlap_tad_gene)
  overlap_tad_df <- TAD[overlap_tad_gene$subjectHits,]
  overlap_gene_df <-  as.data.frame(hg38_genes)
  overlap_gene_df <- overlap_gene_df[overlap_tad_gene$queryHits,]
  overlap_df <- cbind(overlap_gene_df,overlap_tad_df)
  return (overlap_df)
}
# add gene symbol
get_overlap_df_2 <- function(overlap_df){
	gene_symbl <- getSYMBOL(overlap_df$gene_id,data = 'org.Hs.eg')
    overlap_df_2 <- cbind(overlap_df,gene_symbl)
    return (overlap_df_2)
}
# count the genes within each TAD
get_overlap_tbl <- function(overlap_df_2){
	overlap_tbl <- overlap_df_2 %>% 
	    group_by(seqnames,tad_start,tad_end) %>%
        count() %>%
        filter(n >= 2)
    overlap_tbl <- data.frame(overlap_tbl)
    return (overlap_tbl)
}
overlap_df <- overlape_human(tad_grange)
overlap_df_2 <- get_overlap_df_2(overlap_df)
overlap_tbl <- get_overlap_tbl(overlap_df_2)
View(overlap_tbl)
######################################################
# construct a function to get correlation list
# overlap_tbl: count 
# overlap_df_2: original overlap dataframe with gene symbol
get_correlation_list <- function(overlap_tbl,overlap_df_2){
  cor_li <- c()
  for(i in 1:nrow(overlap_tbl)){
    d <- overlap_df_2 %>%
    filter(seqnames == overlap_tbl[i,1] & tad_start == overlap_tbl[i,2] & tad_end == overlap_tbl[i,3])
    sample_id <- as.numeric(match(d$gene_symbl,liver_exp$Description))
    sample_id <- sample_id[!is.na(sample_id)]
    res <- lapply(data_list,mean_exp)
    mean_all <- bind_cols(res)[seq(2,62,2)]
    trans_mat <- unname(t(mean_all[,c(-6)]))
    cor_mat <- cor(trans_mat)
    cor_mat_tri <-lower.tri(cor_mat,diag = TRUE)
    cor_mat[cor_mat_tri] <- NA
    corr_matrix <- melt(cor_mat,na.rm = TRUE)
    y = corr_matrix$value
    cor_li <- append(cor_li,y) 
  }
  return(cor_li)
}
# transform to datatype
vec_to_frame <- function(vec){
	return(data.frame(id = rep(1,length(vec)),vec))
}
# use 2 functions above
TAD_correlation <- get_correlation_list(overlap_tbl,overlap_df_2)
TAD_correlation_f <- vec_to_frame(TAD_correlation)
# summary statistics of correlation value list
summary(TAD_correlation)
# plot
boxplot(TAD_correlation)
ggplot(TAD_correlation_f ,aes(x = id,y = vec)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1)
ggdensity(TAD_correlation_f$vec)
######################################################
# randomly choose two genes and then calculated the correlations between them
random_correlation <- c()
mean_func_random <-function(x){
  sample_id <- sample(seq(1,nrow(liver_exp)),5000,replace = FALSE)
  d <- data.frame(ID = rownames(x[sample_id,-c(1,2)]),
                  mean = rowMeans(x[sample_id,-c(1,2)]))
  colnames(d) <- c("ID","tissue")
  return(d)
}

res <- lapply(data_list,mean_func_random)
mean_all <- bind_cols(res)[seq(2,62,2)]
View(mean_all)
trans_mat <- unname(t(mean_all[,c(-6)]))
cor_mat <- cor(trans_mat)
View(cor_mat)
cor_mat <- cor(t(mean_all[,c(-1)]))
cor_mat_tri <-lower.tri(cor_mat,diag = TRUE)
cor_mat[cor_mat_tri] <- NA
corr_matrix <- melt(cor_mat,na.rm = TRUE)
random_correlation = corr_matrix$value
summary(random_correlation)
boxplot(random_correlation)
cor_fram_random <- data.frame(id = rep(1,length(random_correlation)),random_correlation)
ggplot(cor_fram_random ,aes(x = id,y = cor_li_random)) +
  geom_violin(fill="gray") + geom_boxplot(width=0.1)
ggdensity(cor_fram_random$cor_li_random)

