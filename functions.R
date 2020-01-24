###################FUNCTION###########################
# functions for organized_correlation(1).R files
###################FUNCTION###########################

# Function:rescale the position with different resolutions
# Input: 1) dataframe with original start and end position
#        2) resolution
# Output: dataframe with rescaled start and end position
rescale_resolution <- function(df,resolution){
  fun <- function(x) as.numeric(x)/as.numeric(resolution)
  df_out <- df %>%
    mutate_at(c(2,3),fun)
  return(df_out)
}

######################################################

# Function:extract rows with specified gene names
# Input: 1) original gene expression dataframe extracted from GTEx
#        2) sample_id
# Output: gene expression dataframe with specified genes
extract_generow <-function(df,sample_id){
  d <- data.frame(ID = rownames(df[sample_id,-c(1,2)]),
                  df[sample_id,-c(1,2)])
  return(d)
}

######################################################

# Function:extract rows with specified gene names and then calculated means by row
# Input: 1) original gene expression dataframe extracted from GTEx
#        2) sample_id
# Output: gene expression dataframe with specified genes
extract_generow_mean <-function(df,sample_id){
  d <- data.frame(ID = rownames(df[sample_id,-c(1,2)]),
                  mean = rowMeans(df[sample_id,-c(1,2)]))
  colnames(d) <- c("ID","tissue")
  return(d)
}

######################################################


# Function:apply grange to all TAD boundaries 
# Input: TAD boundaries
# Output: gene expression dataframe with specified genes
TAD_transform <- function(TAD){
  TAD <- TAD[,c(1:3)]
  colnames(TAD) <- c("chrom","tad_start","tad_end")
  return(TAD)
}

TAD_to_grange <- function(TAD){
	tad_grange <- makeGRangesFromDataFrame(TAD,ignore.strand = T,seqnames.field = "chrom",start.field = "tad_start",end.field = "tad_end")
  return (tad_grange)
}

######################################################

# Function: find genes within TAD
# Input: 1) grange from data  
#        2) data frame with col(chromosome,start,end); 
#        1) and 2) are the same data with different forms 
# Output: data frame with col(seqnames,start,end,width,strand,gene_id,chrom,tad_start,tad_end)
overlap_human <- function(tad_grange,TAD){
  overlap_gene_df <-  as.data.frame(hg38_genes)
  overlap_tad_gene <- findOverlaps(hg38_genes,tad_grange,select = "all",type = "within")
  overlap_tad_gene <- as.data.frame(overlap_tad_gene)
  overlap_tad_df <- TAD[overlap_tad_gene$subjectHits,]
  overlap_gene_df <- overlap_gene_df[overlap_tad_gene$queryHits,]
  overlap_df <- cbind(overlap_gene_df,overlap_tad_df)
  return (overlap_df)
}

######################################################

# Function: find genes within BLOCK
# Input: 1) grange from data  
#        2) data frame with col(chromosome,start,end); 
#        1) and 2) are the same data with different forms 
# Output: data frame with col(seqnames,start,end,width,strand,gene_id,chrom,tad_start,tad_end)
overlap_human_block <- function(tad_grange,TAD){
  overlap_gene_df <-  as.data.frame(hg38_genes)
  overlap_tad_gene <- findOverlaps(hg38_genes,tad_grange,select = "all",type = "within")
  overlap_tad_gene <- as.data.frame(overlap_tad_gene)
  overlap_tad_df <- TAD[overlap_tad_gene$subjectHits,]
  overlap_gene_df <- overlap_gene_df[overlap_tad_gene$queryHits,]
  overlap_df <- cbind(overlap_gene_df,overlap_tad_df)
  return (overlap_df)
}

######################################################

# Function: add gene symbol to overlap_df
# Input: dataframe overlap_df
# Output: dataframe added with gene symbol
get_overlap_df_2 <- function(overlap_df){
  gene_symbl <- getSYMBOL(overlap_df$gene_id,data = 'org.Hs.eg')
  overlap_df_2 <- cbind(overlap_df,gene_symbl)
  return (overlap_df_2)
}

######################################################

# Function: count the number of genes within every TAD
# Input: dataframe overlap_df with gene symbol
# Output: dataframe grouped by TAD with col(seqnames,tad_start,tad_end,n)
get_overlap_tbl <- function(overlap_df_2){
  overlap_tbl <- overlap_df_2 %>% 
    group_by(seqnames,tad_start,tad_end) %>%
    count() %>%
    filter(n >= 2)
  overlap_tbl <- data.frame(overlap_tbl)
  return (overlap_tbl)
}

######################################################

# Function: combine the procedure overlap_human, get_overlap_df_2, get_overlap_tbl together
# Input: 1) grange from data  
#        2) data frame with col(chromosome,start,end); 
#        1) and 2) are the same data with different forms 
# Output: a list with [1]: overlap_df_2;[2]: overlap_tbl
pre_procedure <- function(tad_grange,TAD){
  overlap_df <- overlap_human(tad_grange,TAD)
  overlap_df_2 <- get_overlap_df_2(overlap_df)
  overlap_tbl <- get_overlap_tbl(overlap_df_2)
  results <- list(overlap_df_2,overlap_tbl)
  return(results)
}

######################################################

# Function: Get correlation list between genes that are in the same TAD
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol 
#        3) data list
# Output: correlation list

correlation_list_mean <- function(overlap_tbl,overlap_df_2,data_list){
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
    # apply to the all gene expression data frame and calculate the mean by row
    res <- sapply(data_list,extract_generow_mean,sample_id = sample_id)
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

######################################################

# Function: Get correlation list between genes that are in the same TAD (But average the gene expression)
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol 
#        3) data list
# Output: correlation list

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

######################################################

# Function: Get correlation list between genes that are in the same random block
# Input: 1) data frame grouped by ran
#        2) data frame with gene symbol 
#        3) data list
# Output: correlation list

# chromosome name
chromosome_name <- c("chr1","chr2","chr3","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                     "chr20","chr21","chr22","chrX") 
# chromosome length
chromosome_length <- lapply(chromosome_name, function(x) seqlengths(genome)[[x]])
chromosome_length <- setNames(chromosome_length,chromosome_name)

# size of TADs on different chromosomes
TAD_size <- get_TAD_size_list(TAD)

# number of TADs on different chromosomes
TAD_number <- get_TAD_number_list(TAD)

random_block_correlation <- function(chromosome_name,chromosome_length,TAD_number,TAD_size,data_list){
  cor_li <- c()
  for(ch in chromosome_name){
    tad_size <- TAD_size[[ch]]
    chrom_length <- chromosome_length[[ch]]
    random_block_start <- sample(seq(1,chrom_length,1),TAD_number[[ch]])
    random_block_end <- random_block_start + sample(tad_size)
    df <- data.frame(ch,random_block_start,random_block_end)
    colnames(df) <- c("seqnames","start","end")
    df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end")
    df_1 <- overlap_human(df_grange,TAD)
    df_2 <- get_overlap_df_2(df_1)
    df_tbl <- get_overlap_tbl(df_2)
    y <- correlation_list(df_tbl,df_2,data_list_whole)
    cor_li <- append(cor_li,y) 
  }
  return(cor_li)
}

######################################################

# Function: Get correlation list between genes that are in the same random block(average over sample)
# Input: 1) data frame grouped by random block
#        2) data frame with gene symbol 
#        3) data list
# Output: correlation list

random_block_correlation <- function(chromosome_name,chromosome_length,TAD_number,TAD_size,data_list){
  cor_li <- c()
  for(ch in chromosome_name){
    tad_size <- TAD_size[[ch]]
    chrom_length <- chromosome_length[[ch]]
    random_block_start <- sample(seq(1,chrom_length,1),TAD_number[[ch]])
    random_block_end <- random_block_start + sample(tad_size)
    df <- data.frame(ch,random_block_start,random_block_end)
    colnames(df) <- c("seqnames","start","end")
    df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end")
    df_1 <- overlap_human(df_grange,TAD)
    df_2 <- get_overlap_df_2(df_1)
    df_tbl <- get_overlap_tbl(df_2)
    y <- correlation_list_mean(df_tbl,df_2,data_list_whole)
    cor_li <- append(cor_li,y) 
  }
  return(cor_li)
}

######################################################
######################################################

# Function: Get correlation list between genes that are in the same TAD (tissue-specific)
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol 
#        3) tissue expression data
# Output: correlation list

correlation_list_tissue <- function(overlap_tbl,overlap_df_2,tissue_exp){
    
    cor_li <- c()
    for(i in 1:nrow(overlap_tbl)){
        d <- overlap_df_2 %>%
        filter(seqnames == overlap_tbl[i,1] & tad_start == overlap_tbl[i,2] & tad_end == overlap_tbl[i,3])
        sample_id <- as.numeric(match(d$gene_symbl,Liver_exp$Description))
        sample_id <- sample_id[!is.na(sample_id)]
        res <- extract_generow(tissue_exp,sample_id)
        mean_all <- res
        trans_mat <- data.matrix(t(mean_all[,c(-1)]))
        rownames(trans_mat) <- c()
        trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
        cor_mat <- cor(trans_mat)
        cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
        cor_mat[cor_mat_tri] <- NA
        corr_matrix <- melt(cor_mat,na.rm = TRUE)
        y = corr_matrix$value
        cor_li <- append(cor_li,y)
    }
    return(cor_li)
	
}

######################################################

# Function: combine pre-processing with correlation together
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol 
#        3) tissue expression data
# Output: correlation list

correlation_tissue_specific <- function(tad_grange,TAD,tissue_exp){
    res_1 <- pre_procedure(tad_grange,TAD)
    res_2 <- correlation_list_tissue(data.frame(res_1[2]),
                                     data.frame(res_1[1]),
                                        tissue_exp)
    return (res_2)
}

######################################################

# Function: transform from list to data frame in order to plot violin plot
# Input: correlation list
# Output: dataframe of correlation list
vec_to_frame <- function(vec){
  return(data.frame(id = rep(1,length(vec)),vec))
}

######################################################

# Function: Randomly choose 2 genes 150 times and calculate correlations
# Input: data list of 20 gene expression data frame
# Output: correlation list
random_2_74477_cor <-function(data_list){
  cor_li <- c()
  gene_exp_all <- bind_cols(data_list)
  #gene_exp_all_no <- gene_exp_all %>%
      #filter(GTEX.1117F.0226.SM.5GZZ7 != 0)

  for(i in 1:74477){
    sample_id <- sample(seq(1,nrow(gene_exp_all)),2,replace = FALSE)
    gene_exp_some <- data.frame(ID = rownames(gene_exp_all[sample_id,-c(1,2)]),
                                gene_exp_all[sample_id,-c(1,2)])
    trans_mat <- data.matrix(t(gene_exp_some[,c(-1)]))
    rownames(trans_mat) <- c()
    trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
    trans_mat[is.na(trans_mat)] <- 0
    cor_mat <- cor(trans_mat)
    cor_mat_tri <-lower.tri(cor_mat,diag = TRUE)
    cor_mat[cor_mat_tri] <- NA
    corr_matrix <- melt(cor_mat,na.rm = TRUE)
    y = corr_matrix$value
    cor_li <- append(cor_li,y)
  }
  return(cor_li)
}

######################################################

# Function: calculate the gene co-expression correlation varied by distance
# Input: 
# Output: 

##preparation
# get datasets including distance
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
# get gene symbols
gene_list <- unlist(lapply(str_split(Adrenal_Gland_exp$Name,"\\."), function(y) y[1]))


gene_distance_correlation <- function(){
  for(ch in chromosome_name){
    dc_df <- data.frame(0,0)
    colnames(dc_df) <- c("distance","correlation")
    # create table count
    ch = "chr1"
    df <- data.frame(ch,c(1),chromosome_length[[ch]])
    colnames(df) <- c("seqnames","chr_start","chr_end")
    df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end")
    df_1 <- overlap_human(df_grange,TAD)
    df_2 <- get_overlap_df_2(df_1)
    head(df_2)
    sample_id <- as.numeric(match(df_2$gene_symbl,Liver_exp$Description))
    sample_id <- sample_id[!is.na(sample_id)]
    sample_gene_id <- Liver_exp[sample_id,1]
    sample_gene_id_new <- unlist(lapply(str_split(sample_gene_id,"\\."), function(y) y[1]))
    res <- sapply(data_list_whole,extract_generow,sample_id = sample_id)
    mean_all <- bind_cols(res)
    for(i in 1:nrow(mean_all)){
      for(j in 1:nrow(mean_all)){
        correlation <- cor(as.numeric(mean_all[i,c(-1)]),as.numeric(mean_all[j,c(-1)]))
        gene_info <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
                                         filters=c('ensembl_gene_id'),
                                         values=c(sample_gene_id_new[c(i,j)]),
                                         mart=ensembl)
        if(nrow(gene_info) ==2){
          distance <- abs(gene_info[1,2] - gene_info[2,2])
          dc_df <- rbind(dc_df,data.frame(distance,correlation))
        }
        
    }
    }
  }
  return(dc_df)
}












######################################################

# Function: calculate distance between two genes within a TAD(all)
# Input ：
# Output: a vector of gene distance list
gene_dist_list <- function(overlap_tbl,overlap_df_2){
  dist_li <- c()
  for(i in 1:nrow(overlap_tbl)){
    d <- overlap_df_2 %>%
      filter(seqnames == overlap_tbl[i,1] & tad_start == overlap_tbl[i,2] & tad_end == overlap_tbl[i,3])
    dist_li <- append(dist_li,as.vector(dist(d$start)))
  }
  return(dist_li)
}

######################################################

# Function: calculate distance between two genes within a TAD (for every chromosome)
# Input ：
# Output: a vector of gene distance list
gene_dist_chr_list <- function(overlap_tbl,overlap_df_2){
  chr_name <- as.character(unique(overlap_tbl$seqnames))
  for (chr in chr_name){
    d <- overlap_df_2 %>%
      filter(seqnames == chr)
    dist <- as.vector(dist(d$start))
    assign(chr,dist)
  }
  chr_list <- list(chr1,chr2,chr3,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,
               chr20,chr21,chr22,chrX)
  return(chr_list)
}

######################################################

# Function: calculate TAD size (for every chromosome)
# Input ：TAD dataframe
# Output: a vector of gene distance list
get_TAD_size_list <- function(TAD){
  chr_list <- list()
  chr_name <- as.character(unique(TAD$chrom))
  for (chr in chr_name){
    d <- TAD %>%
      filter(chrom == chr)
    dist <- as.vector(as.numeric(d$tad_end) - as.numeric(d$tad_start))
    chr_list[[chr]] = dist
  }
  return(chr_list)
}

######################################################

# Function: get TAD numbers (for every chromosome)
# Input ：TAD dataframe
# Output: a vector of gene distance list
get_TAD_number_list <- function(TAD){
  tad_num <- list()
  chr_name <- as.character(unique(TAD$chrom))
  for (chr in chr_name){
    x <- TAD %>%
      filter(chrom == chr) %>% nrow()
    tad_num[[chr]] = x
  }
  return(tad_num)
}

######################################################

######################################################

# Function: convert original bed file to tad size for histogram
# Input ：TAD bed file
# Output: dataframe containing group and tad size
different_TAD_size <- function(TAD,chr_TAD){
  hist_TAD <- data.frame(as.matrix(as.numeric(TAD[,3])-as.numeric(TAD[,2])))
  colnames(hist_TAD) <- c("size")
  hist_TAD$group <- rep(c(chr_TAD),nrow(hist_TAD))
  return(hist_TAD)
}

######################################################
######################################################

# Function: Calculate the correlation in random blocks
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol
#        3) data list
# Output: correlation list
const_distance_correlation <- function(overlap_tbl,overlap_df_2,data_list){
  gene_dist <- gene_dist_list(overlap_tbl,overlap_df_2)
  cor_li <- c()
  for(j in 1:nrow(overlap_tbl)){
    d <- overlap_df_2 %>%
      filter(seqnames == overlap_tbl[j,1] & tad_start == overlap_tbl[j,2] & tad_end == overlap_tbl[j,3])
    for(i in 1:nrow(d)){
      x <- d$end[i] + sample(gene_dist,1) # new_start candidate
      y <- d$start[i] - sample(gene_dist,1) # new_end candidate
      if(x > d[i,"tad_end"]){
        df <- data.frame(d[1,1],c(d$start[i],x),c(d$end[i],x + 20000))
        colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
        df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                              seqnames.field = "seqnames",
                                              start.field = "gene_start_reg",
                                              end.field = "gene_end_reg")
        p_gene <- overlap_human_block(df_grange,df)$gene_id
        p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
        p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
        if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
          res <- sapply(data_list,extract_generow,sample_id = p_gene_id)
          mean_all <- bind_cols(res)
          trans_mat <- data.matrix(t(mean_all[,c(-1)]))
          trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
          cor_mat <- cor(trans_mat)
          cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
          cor_mat[cor_mat_tri] <- NA
          corr_matrix <- melt(cor_mat,na.rm = TRUE)
          cor_li <- append(cor_li,corr_matrix$value)
        }
        
      }
      if(y < d[i,"tad_start"]){
        df <- data.frame(d[1,1],c(d$start[i],y - 20000),c(d$end[i],y))
        colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
        df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                              seqnames.field = "seqnames",
                                              start.field = "gene_start_reg",
                                              end.field = "gene_end_reg")
        p_gene <- overlap_human_block(df_grange,df)$gene_id
        p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
        p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
        if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
          res <- sapply(data_list,extract_generow,sample_id = p_gene_id)
          mean_all <- bind_cols(res)
          trans_mat <- data.matrix(t(mean_all[,c(-1)]))
          trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
          cor_mat <- cor(trans_mat)
          cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
          cor_mat[cor_mat_tri] <- NA
          corr_matrix <- melt(cor_mat,na.rm = TRUE)
          cor_li <- append(cor_li,corr_matrix$value)
        }
      }
    }
  }
  return(cor_li)
  
}

######################################################

######################################################

# Function: Calculate the correlation between genes within TAD and genes outside TAD
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol
#        3) data list
# Output: correlation list
const_distance_correlation <- function(overlap_tbl,overlap_df_2,data_list){
    gene_dist <- gene_dist_list(overlap_tbl,overlap_df_2)
    cor_li <- c()
    for(j in 1:nrow(overlap_tbl)){
        d <- overlap_df_2 %>%
        filter(seqnames == overlap_tbl[j,1] & tad_start == overlap_tbl[j,2] & tad_end == overlap_tbl[j,3])
        for(i in 1:nrow(d)){
            x <- d$end[i] + sample(gene_dist,1) # new_start candidate
            y <- d$start[i] - sample(gene_dist,1) # new_end candidate
            if(x > d[i,"tad_end"]){
                df <- data.frame(d[1,1],c(d$start[i],x),c(d$end[i],x + 20000))
                colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
                df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                seqnames.field = "seqnames",
                start.field = "gene_start_reg",
                end.field = "gene_end_reg")
                p_gene <- overlap_human_block(df_grange,df)$gene_id
                p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
                p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
                if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
                    res <- sapply(data_list,extract_generow,sample_id = p_gene_id)
                    mean_all <- bind_cols(res)
                    trans_mat <- data.matrix(t(mean_all[,c(-1)]))
                    trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
                    cor_mat <- cor(trans_mat)
                    cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
                    cor_mat[cor_mat_tri] <- NA
                    corr_matrix <- melt(cor_mat,na.rm = TRUE)
                    cor_li <- append(cor_li,corr_matrix$value)
                }
                
            }
            if(y < d[i,"tad_start"]){
                df <- data.frame(d[1,1],c(d$start[i],y - 20000),c(d$end[i],y))
                colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
                df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                seqnames.field = "seqnames",
                start.field = "gene_start_reg",
                end.field = "gene_end_reg")
                p_gene <- overlap_human_block(df_grange,df)$gene_id
                p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
                p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
                if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
                    res <- sapply(data_list,extract_generow,sample_id = p_gene_id)
                    mean_all <- bind_cols(res)
                    trans_mat <- data.matrix(t(mean_all[,c(-1)]))
                    trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
                    cor_mat <- cor(trans_mat)
                    cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
                    cor_mat[cor_mat_tri] <- NA
                    corr_matrix <- melt(cor_mat,na.rm = TRUE)
                    cor_li <- append(cor_li,corr_matrix$value)
                }
            }
        }
    }
    return(cor_li)
    
}

######################################################

# Function: Calculate the correlation between genes within TAD and genes outside TAD(tissue-specific)
# Input: 1) data frame grouped by TAD
#        2) data frame with gene symbol 
# Output: correlation list
const_distance_correlation_ts <- function(overlap_tbl,overlap_df_2,tissue_exp){
  gene_dist <- gene_dist_list(overlap_tbl,overlap_df_2)
  cor_li <- c()
  for(j in 1:nrow(overlap_tbl)){
    d <- overlap_df_2 %>%
      filter(seqnames == overlap_tbl[j,1] & tad_start == overlap_tbl[j,2] & tad_end == overlap_tbl[j,3])  
    for(i in 1:nrow(d)){
        x <- d$end[i] + sample(gene_dist,1) # new_start candidate
        y <- d$start[i] - sample(gene_dist,1) # new_end candidate
        if(x > d[i,"tad_end"]){
          df <- data.frame(d[1,1],c(d$start[i],x),c(d$end[i],x + 20000))
          colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
          df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                                    seqnames.field = "seqnames",
                                                    start.field = "gene_start_reg",
                                                    end.field = "gene_end_reg")
          p_gene <- overlap_human_block(df_grange,df)$gene_id
          p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
          p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
          if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
            res <- extract_generow(tissue_exp,p_gene_id)
            mean_all <- res
            trans_mat <- data.matrix(t(mean_all[,c(-1)]))
            trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
            cor_mat <- cor(trans_mat)
            cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
            cor_mat[cor_mat_tri] <- NA
            corr_matrix <- melt(cor_mat,na.rm = TRUE)
            cor_li <- append(cor_li,corr_matrix$value)
          }

        }
        if(y < d[i,"tad_start"]){
          df <- data.frame(d[1,1],c(d$start[i],y - 20000),c(d$end[i],y))
          colnames(df) <- c("seqnames","gene_start_reg","gene_end_reg")
          df_grange <- makeGRangesFromDataFrame(df,ignore.strand = T,
                                                seqnames.field = "seqnames",
                                                start.field = "gene_start_reg",
                                                end.field = "gene_end_reg")
          p_gene <- overlap_human_block(df_grange,df)$gene_id
          p_gene_symbl <- getSYMBOL(p_gene,data = 'org.Hs.eg')
          p_gene_id <- as.numeric(match(p_gene_symbl,Liver_exp$Description))
          if(length(p_gene_id[!is.na(p_gene_id)]) > 1){
            res <- extract_generow(tissue_exp,p_gene_id)
            mean_all <- res
            trans_mat <- data.matrix(t(mean_all[,c(-1)]))
            trans_mat <- as.data.frame(apply(trans_mat,2, as.numeric))
            cor_mat <- cor(trans_mat)
            cor_mat_tri <- lower.tri(cor_mat,diag = TRUE)
            cor_mat[cor_mat_tri] <- NA
            corr_matrix <- melt(cor_mat,na.rm = TRUE)
            cor_li <- append(cor_li,corr_matrix$value)
          }
        }
    }
 }
  return(cor_li)
}

######################################################

# Function: combine all the process together
# Input: 1) tissue specific TAD bounadaires
#        2) tissue specific expression 
# Output:1) TAD dataframe
#        2) TAD grange
#        3) tissue specific correlation list
corr_tissue_all <- function(TAD,tissue_exp){
  TAD_to_official <- TAD_transform(TAD)
  official_to_grange <- TAD_to_grange(TAD_to_official)
  corr_li <- correlation_tissue_specific(official_to_grange,TAD_to_official,tissue_exp)
  results <- list(TAD_to_official,official_to_grange,corr_li)
  return(results)
}



######################################################

# Function: when coordinate x and y equals m and n we get the contact value and avoid null at the same time
# Input: coordinates m and n
# Output: contact value 
extract_contact <- function(m,n){
  x =  mouse_chr1_contact %>% filter(V1 == m & V2 == n) %>% "$"(V3)
  if(length(x) ==0){
    return(0)
  } else{
    return(x)
  }
}

######################################################

# Function: extract tad_start or tad_end from every tissue-specific 
# Input: index from 1~15,chromosome number
# Output: 

extract_start <- function(index,chr_n){
  col_vl <- data.frame(all_boundaries[index]) %>%
    filter(V1 == chr_n) %>% "$"(V2)
  return(col_vl)
}

extract_end <- function(index,chr_n){
  col_vl <- data.frame(all_boundaries[index]) %>%
    filter(V1 == chr_n) %>% "$"(V3)
  return(col_vl)
}

######################################################

# Function:  
# Input: dataframe 1 and dataframe 2
# Output: 
chr_name <- paste("chr",c(seq(1,22),"X","Y"),sep = "")
testify_overlap <- function(df1,df2){
  for (name in chr_name){
    df_process_1 <- df1[df1$V1 == name,][,c(2,3)] 
    df_process_2 <- df2[df2$V1 == name,][,c(2,3)] 
    for( i in 1:5){
      x <- c(df_process_1[i,1],df_process_1[i,2])
      for( j in 1:5){
        y <- c(df_process_2[j,1],df_process_2[j,2])
        pert <- Overlap(x,y)/40000
        if(pert >= 0.9){
          z = conserved_tad
          temp_list <- c(name,x,y)
          conserved_tad <- list(z,temp_list)
        }
      }
    }
  }
  return(conserved_tad)
}
  
  
# 1:nrow(df_process_1
# 1:nrow(df_process_2


#=====================================================
######################################################

# Function: replace NA with zero
# Input: dataframe containing NA
# Output: dataframe with zero
replaceNAwithZero <- function(x){
    x[is.na(x)] <- 0
    return(x)
}
#=====================================================
#ensembl <- useMart("ensembl")
#ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
#      filters=c('hgnc_symbol', 'ensembl_gene_id'),
#      values=list('TTN', 'ENSG00000155657'),
#      mart=ensembl)
