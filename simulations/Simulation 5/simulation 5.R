########################
setwd('D:/paper')
library(corncob)
# Remove rows that are all zeros
run_corncob <- function(countdata, metadata) {
  library(phyloseq)
  set.seed(1) 
  flag <- rowSums(countdata[]) > 0
  countdata <- countdata[flag, ]
  countdata <- as.matrix(countdata)
  identical(colnames(countdata), rownames(metadata))
  rows_to_keep <- intersect(colnames(countdata), rownames(metadata))
  metadata <- metadata[rows_to_keep,, drop = F]
  countdata <- countdata[, rows_to_keep]
  identical(colnames(countdata), rownames(metadata))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(countdata, taxa_are_rows = TRUE), 
                           phyloseq::sample_data(metadata))
  # Differential abundance
  corn_da <- differentialTest(formula = ~ comparison,
                              phi.formula = ~ comparison,
                              formula_null = ~ 1,
                              phi.formula_null = ~ comparison,
                              data = ps, 
                              test = "Wald", boot = FALSE,
                              fdr_cutoff = 0.05)
  
  fdr_corncob <- corn_da$significant_taxa
  dim(data.frame(fdr_corncob))
  fdr_corncob <- as.data.frame(fdr_corncob) 
  
  return(fdr_corncob)
}
corncob_feng <- run_corncob(otu, meta)
####################### Load data#####
### Feng's original data
load("D:/paper/realdata/DAA-main/Data/realdata/FengQ_2015.Rdata")
otu <- phy@otu_table
otu <- as.data.frame(otu)
taxa <- as.data.frame(phy@tax_table)
meta <- as.data.frame(phy@sam_data)
colnames(meta)[colnames(meta) == "study_condition"] <- "comparison"
run_corncob(otu, meta)
## Save Feng's original data as csv for imputation in GEMimp
write.csv(otu, 'D:/paper/feng_taxa_analysis/countdata.csv')
### Feng's imputed data
otuimpute <- read.csv("D:/paper/feng_taxa_analysis/imputed_file_real_data.csv", row.names = 1)
corncob_feng_imputed <- run_corncob(otuimpute, meta)
colnames(corncob_feng_imputed) <- 'fdr_corncob_imputed'


run_comparision <- function(rawdf, imputeddf) {
  imputeddf <- cbind("id" = rownames(imputeddf), imputeddf)
  rawdf <- cbind("id" = rownames(rawdf), rawdf)
  merged_df <- merge(imputeddf, rawdf, all = TRUE)
  result_df <- data.frame(fdr_corncob_imputed = NA, fdr_corncob = NA)
  for (i in 1:length(merged_df$fdr_corncob_imputed)) {
    if (merged_df$fdr_corncob_imputed[i] %in% merged_df$fdr_corncob) {
      match_index <- which(merged_df$fdr_corncob == merged_df$fdr_corncob_imputed[i])
      result_df <- rbind(result_df, data.frame(fdr_corncob_imputed = merged_df$fdr_corncob_imputed[i], fdr_corncob = merged_df$fdr_corncob[match_index]))
    } else {
      result_df <- rbind(result_df, data.frame(fdr_corncob_imputed = merged_df$fdr_corncob_imputed[i], fdr_corncob = NA))
    }
  }
  for (i in 1:length(merged_df$fdr_corncob)) {
    if (!(merged_df$fdr_corncob[i] %in% merged_df$fdr_corncob_imputed)) {
      result_df <- rbind(result_df, data.frame(fdr_corncob_imputed = NA, fdr_corncob = merged_df$fdr_corncob[i]))
    }
  }
  result_df <- result_df[-1, ]
  return(result_df)
}

############# Place Feng's original data results corncob_feng and imputed data results corncob_feng_imputed into this function's dataframe
feng_corncob <- run_comparision(corncob_feng, corncob_feng_imputed)
write.csv(feng_corncob, 'D:/paper/feng_taxa_analysis/feng_corncob_analysis.csv')

################################# 
run_Maaslin2 <- function(countdata, metadata) {
  library(Maaslin2)
  library(phyloseq)
  ps <- phyloseq::phyloseq(phyloseq::otu_table(countdata, taxa_are_rows = TRUE), 
                           phyloseq::sample_data(metadata))
  mas <- Maaslin2(
    input_data = data.frame(otu_table(ps)),
    input_metadata = data.frame(sample_data(ps)),
    output = "./Maaslin2_default_zeller_output",
    min_abundance = 0.0,
    min_prevalence = 0.0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.05,
    fixed_effects = "comparison",
    correction = "BH",
    standardize = FALSE,
    cores = 1)
  
  mas_res_df <- mas$results
  
  fdr_mas <- mas_res_df %>% dplyr::filter(qval < 0.05)
  
  dim(fdr_mas)
  return(fdr_mas)
}
colnames(meta)[colnames(meta) == "study_condition"] <- "comparison"
Maaslin2_feng <- run_Maaslin2(otu, meta)   # Feng's data taxa
Maaslin2_feng_imputed <- run_Maaslin2(otuimpute, meta) # Feng's imputed taxa
save(Maaslin2_feng, file = "D:/paper/feng_taxa_analysis/taxa_feng_Qin.RData")
save(Maaslin2_feng_imputed, file = "D:/paper/feng_taxa_analysis/taxa_Maaslin2_feng_imputed.RData")

###################################################
# Try with generated data
simudata <- read.csv('D:/paper/simulated_data/SyntheticMicrobiome-Counts.txt', sep = '\t', row.names = 1)
metadata <- simudata[1, ]
simudata <- simudata[-1, ]
metadata <- t(metadata)
colnames(metadata) <- 'comparison'
metadata <- as.data.frame(metadata)
Maaslin2_simu <- run_Maaslin2(simudata, metadata)

###################################### KarlssonFH data #######################
#################################### corncob imputation method ####################
otu_KarlssonFH <- phy@otu_table
otu_KarlssonFH <- as.data.frame(otu_KarlssonFH)  # count data
# Save the count data as csv for imputation in GEMimp
write.csv(otu_KarlssonFH, 'D:/paper/Karlsson_taxa_analysis/Karlsson_count.csv')
taxa <- as.data.frame(phy@tax_table)
meta <- as.data.frame(phy@sam_data)  # metadata
colnames(meta)[colnames(meta) == "study_condition"] <- "comparison"
# Import imputed KarlssonFH data
otu_KarlssonFH_imputed <- read.csv('D:/paper/Karlsson_taxa_analysis/imputed_file_real_data.csv', row.names = 1)
################### corncob method
corncob_KarlssonFH <- run_corncob(otu_KarlssonFH, meta)
corncob_KarlssonFH_imputed <- run_corncob(otu_KarlssonFH_imputed, meta)
colnames(corncob_KarlssonFH_imputed) <- 'fdr_corncob_imputed'
KarlssonFH_corncob <- run_comparision(corncob_KarlssonFH, corncob_KarlssonFH_imputed)
write.csv(KarlssonFH_corncob, 'D:/paper/Karlsson_taxa_analysis/Karlsson_corncob_analysis.csv')

################################## Maaslin2 imputation method #######################
Maaslin2_KarlssonFH <- run_Maaslin2(otu_KarlssonFH, meta) # Karlsson's data taxa
Maaslin2_KarlssonFH_imputed <- run_Maaslin2(otu_KarlssonFH_imputed, meta) # Karlsson's imputed taxa
####################### Temporarily save KarlssonFH's data before using run_comparison function to merge
save(Maaslin2_KarlssonFH, file = "D:/paper/Karlsson_taxa_analysis/taxa_Maaslin2_KarlssonFH.RData")
save(Maaslin2_KarlssonFH_imputed, file = "D:/paper/Karlsson_taxa_analysis/taxa_Maaslin2_KarlssonFH_imputed.RData")

####################################### Qin data #######################
#################################### corncob imputation method ####################
load("D:/paper/realdata/DAA-main/Data/realdata/QinJ_2012.Rdata")
otu_Qin <- phy@otu_table
otu_Qin <- as.data.frame(otu_Qin)  # count data
# Save the count data as csv for imputation in GEMimp
write.csv(otu_Qin, 'D:/paper/Qin_taxa_analysis/Qin_count.csv')  # Qin's count data
taxa <- as.data.frame(phy@tax_table)
meta <- as.data.frame(phy@sam_data)  # metadata
colnames(meta)[colnames(meta) == "study_condition"] <- "comparison"
# Import imputed Qin data
otu_Qin_imputed <- read.csv('D:/paper/Qin_taxa_analysis/imputed_file_real_data.csv', row.names = 1)
################### corncob method
corncob_Qin <- run_corncob(otu_Qin, meta)
corncob_Qin_imputed <- run_corncob(otu_Qin_imputed, meta)
colnames(corncob_Qin_imputed) <- 'fdr_corncob_imputed'
Qin_corncob <- run_comparision(corncob_Qin,corncob_Qin_imputed)
write.csv(Qin_corncob,'D:/paper/Qin_taxa_analysis/Qin_corncob_analysis.csv')

################################# Maaslin2 Imputation Method #######################
Maaslin2_Qin <- run_Maaslin2(otu_Qin, meta) # Taxa for Karlsson's data
Maaslin2_Qin_imputed <- run_Maaslin2(otu_Qin_imputed, meta) # Taxa for Karlsson's imputed data
####################### 
save(Maaslin2_Qin, file = "D:/paper/Qin_taxa_analysis/taxa_Maaslin2_Qin.RData")
save(Maaslin2_Qin_imputed, file = "D:/paper/Qin_taxa_analysis/taxa_Maaslin2_Qin_imputed.RData")

################################ zeller Data's Maaslin2 Method ################
otu_zeller <- read.csv('D:/paper/crc_zeller/crc_zeller_genus_table.tsv', sep='\t', row.names = 1) # Original data
imputeddata <- read.csv('D:/paper/imputed_file_data_imp_zeller_realdata.csv', row.names = 1)
metadata <- read.csv('D:/paper/crc_zeller/crc_zeller_metadata.tsv', sep='\t', row.names = 1)
colnames(otu_zeller) <- gsub("\\.", "-", colnames(otu_zeller))
otu_zeller <- as.matrix(otu_zeller)
identical(colnames(otu_zeller), rownames(metadata))
Maaslin2_zeller <- run_Maaslin2(otu_zeller, metadata) # Taxa for zeller's data

############################################### DEseq2 Method

run_deseq <- function(ASV_table, groupings) {
  dim(ASV_table)
  dim(groupings)
  identical(colnames(ASV_table), rownames(groupings))
  rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
  groupings <- groupings[rows_to_keep, , drop = F]
  ASV_table <- ASV_table[, rows_to_keep]
  identical(colnames(ASV_table), rownames(groupings))
  colnames(groupings)[1] <- "Groupings"
  library(DESeq2)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                        colData = groupings,
                                        design = ~ Groupings)
  dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")
  
  res <- results(dds_res,
                 tidy = T,
                 format = "DataFrame",
                 contrast = c("Groupings", "control", "T2D"))
  head(res)
  
  DEG <- res
  logFC_cutoff <- 2
  DEG$change <- as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                 ifelse(DEG$log2FoldChange > logFC_cutoff, "UP", "DOWN"),
                                 "NOT"))
  this_title <- paste0('Cutoff for logFC is ', round(logFC_cutoff, 3),
                       '\nThe number of up taxa is ', nrow(DEG[DEG$change == 'UP',]),
                       '\nThe number of down taxa is ', nrow(DEG[DEG$change == 'DOWN',]))
  DEG <- na.omit(DEG)
  library(ggplot2)
  ggplot(data = DEG, aes(x = log2FoldChange,
                         y = -log10(pvalue),
                         color = change)) +
    geom_point(alpha = 0.8, size = 3) +
    labs(x = "log2 fold change") + ylab("-log10 pvalue") +
    ggtitle(this_title) + theme_bw(base_size = 20) +
    theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    scale_color_manual(values = c('#a121f0', '#bebebe', '#ffad21')) -> p1
  p1 + xlim(NA, 10) + ylim(NA, 30) -> p2
  
  library(patchwork)
  p1 + p2
  DA <- DEG[DEG$change == 'DOWN' | DEG$change == 'UP',] # Results for zeller's crc data set
  return(DA)
}

# Import zeller's data
otu_zeller <- read.csv('D:/paper/crc_zeller/crc_zeller_genus_table.tsv', sep = '\t', row.names = 1) # Original data
imputeddata <- read.csv('D:/paper/imputed_file_data_imp_zeller_realdata.csv', row.names = 1)
metadata <- read.csv('D:/paper/crc_zeller/crc_zeller_metadata.tsv', sep = '\t', row.names = 1)
colnames(otu_zeller) <- gsub("\\.", "-", colnames(otu_zeller))
colnames(imputeddata) <- gsub("\\.", "-", colnames(imputeddata))
colnames(metadata) <- 'comparison'
zeller_deseq <- run_deseq(otu_zeller, metadata)
zeller_deseq_imputed <- run_deseq(imputeddata, metadata)
save(zeller_deseq, file = 'D:/paper/zeller_taxa_analysis/zeller_deseq.Rdata')
save(zeller_deseq_imputed, file = 'D:/paper/zeller_taxa_analysis/zeller_deseq_imputed.Rdata')

# Import feng's data
# Load data
# Feng's original data
load("D:/paper/realdata/DAA-main/Data/realdata/FengQ_2015.Rdata")
otu_feng <- phy@otu_table
otu_feng <- as.data.frame(otu_feng)
taxa_feng <- as.data.frame(phy@tax_table)
meta_feng <- as.data.frame(phy@sam_data)
colnames(meta_feng)[colnames(meta_feng) == "study_condition"] <- "comparison"
# Only extract the comparison column
meta_feng <- meta_feng['comparison']
# Save Feng's original data as csv for imputation in GEMimp
############ Feng's imputed data
otuimpute_feng <- read.csv("D:/paper/feng_taxa_analysis/imputed_file_real_data.csv", row.names = 1)
feng_deseq <- run_deseq(otu_feng, meta_feng)
feng_deseq_imputed <- run_deseq(otuimpute_feng, meta_feng)
save(feng_deseq, file = 'D:/paper/feng_taxa_analysis/feng_deseq.Rdata')
save(feng_deseq_imputed, file = 'D:/paper/feng_taxa_analysis/feng_deseq_imputed.Rdata')

# Import KarlssonFH data
otu_KarlssonFH <- phy@otu_table
otu_KarlssonFH <- as.data.frame(otu_KarlssonFH)  # count data
# Save the count data as csv for imputation in GEMimp
# write.csv(otu_KarlssonFH, 'D:/paper/Karlsson_taxa_analysis/Karlsson_count.csv')
taxa_KarlssonFH <- as.data.frame(phy@tax_table)
meta_KarlssonFH <- as.data.frame(phy@sam_data)  # metadata
colnames(meta_KarlssonFH)[colnames(meta_KarlssonFH) == "study_condition"] <- "comparison"
# Only extract the comparison column
meta_KarlssonFH <- meta_KarlssonFH['comparison']
# Import imputed KarlssonFH data
otu_KarlssonFH_imputed <- read.csv('D:/paper/Karlsson_taxa_analysis/imputed_file_real_data.csv',row.names = 1)
KarlssonFH_deseq <- run_deseq(otu_KarlssonFH,meta_KarlssonFH)
KarlssonFH_deseq_imputed <- run_deseq(otu_KarlssonFH_imputed,meta_KarlssonFH)
save(KarlssonFH_deseq,file = 'D:/paper/Karlsson_taxa_analysis/KarlssonFH_deseq.Rdata')
save(KarlssonFH_deseq_imputed,file = 'D:/paper/Karlsson_taxa_analysis/KarlssonFH_deseq_imputed.Rdata')

###################Import Qin's data
load("D:/paper/realdata/DAA-main/Data/realdata/QinJ_2012.Rdata")
otu_Qin <- phy@otu_table
otu_Qin <- as.data.frame(otu_Qin)  
#write.csv(otu_Qin,'D:/paper/Qin_taxa_analysis/Qin_count.csv')
taxa_Qin <- as.data.frame(phy@tax_table)
meta_Qin <- as.data.frame(phy@sam_data)#metadata
colnames(meta_Qin)[colnames(meta_Qin) == "study_condition"] <- "comparison"
meta_Qin <- meta_Qin[,'comparison']
otu_Qin_imputed <- read.csv('D:/paper/Qin_taxa_analysis/imputed_file_real_data.csv',row.names = 1)
Qin_deseq <- run_deseq(otu_Qin,meta_Qin)
Qin_deseq_imputed <- run_deseq(otu_Qin_imputed,meta_Qin)
save(Qin_deseq,file = 'D:/paper/Qin_taxa_analysis/Qin_deseq.Rdata')
save(Qin_deseq_imputed,file = 'D:/paper/Qin_taxa_analysis/Qin_deseq_imputed.Rdata')


####################Draw  Wayne diagram
library(ggVennDiagram)
####before imputation
length(zeller_deseq$row)
length(taxa_zeller_all_result$fdr_corncob)
x <- list(DESeq2 = na.omit(zeller_deseq$row),
          corncob = na.omit(taxa_zeller_all_result$fdr_corncob),
          Maaslin2 = na.omit(taxa_zeller_all_Maaslin2_result$fdr_Maaslin2))

overlap_data <- process_region_data(Venn(x))
overlap_data[overlap_data$id == 123,]$item
ggVennDiagram(x,label_alpha=0) +
  scale_fill_gradient(low="gray100",high = "gray60",guide="none")


#####after imputation
x2 <- list(DESeq2 = na.omit(zeller_deseq_imputed$row),
           corncob = na.omit(taxa_zeller_all_result$fdr_corncob_imputed),
           Maaslin2 = na.omit(taxa_zeller_all_Maaslin2_result$fdr_Maaslin2_imputed))

overlap_data <- process_region_data(Venn(x2))
overlap_data[overlap_data$id == 123,]$item
ggVennDiagram(x2,label_alpha=0) +
  scale_fill_gradient(low="gray100",high = "gray60",guide="none")


run_Maaslin2_pval <- function(countdata,metadata){
  library(Maaslin2)
  library(phyloseq)
  ps <- phyloseq::phyloseq(phyloseq::otu_table(countdata, taxa_are_rows = TRUE), 
                           phyloseq::sample_data(metadata))
  mas <- Maaslin2(
    input_data = data.frame(otu_table(ps)),
    input_metadata = data.frame(sample_data(ps)),
    output = "./Maaslin2_default_KarlssonFH_pval",
    min_abundance = 0.0,
    min_prevalence = 0.0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.05,
    fixed_effects = "comparison",
    correction = "BH",
    standardize = FALSE,
    cores = 1)
  
  mas_res_df <- mas$results
  
  return(mas)
  return(mas_res_df)
}
#before imputation
KarlssonFH_pval <- run_Maaslin2_pval(otu_KarlssonFH,meta_KarlssonFH)
result_KarlssonFH_pval <- KarlssonFH_pval$results


KarlssonFH_pval <- result_KarlssonFH_pval[,c('feature','pval')]

KarlssonFH_pval$interval <- cut(KarlssonFH_pval$pval, breaks = c(0.0,0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7,0.8, 0.9,1.0))

interval_counts <- table(KarlssonFH_pval$interval)

# Draw a histogram
barplot(interval_counts, 
        main = "Taxa Counts by p-value Interval raw", 
        xlab = "p-value Interval", 
        ylab = "Counts", 
        col = "transparent",  
        border = "black", 
        ylim = c(0, max(interval_counts) + 10),
        xaxt = "n",  
        space = 0)   


axis(side = 1, at = c(0,2,4,6,8,10), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), col.axis = "black")

###after imputation

KarlssonFH_pval_imputed <- run_Maaslin2_pval(otu_KarlssonFH_imputed,meta_KarlssonFH)
result_KarlssonFH_pval_imputed <- KarlssonFH_pval_imputed$results


KarlssonFH_pval_imputed <- result_KarlssonFH_pval_imputed[,c('feature','pval')]
KarlssonFH_pval_imputed$interval <- cut(KarlssonFH_pval_imputed$pval, breaks = seq(0, 1, by = 0.1))
interval_counts <- table(KarlssonFH_pval_imputed$interval)

#
barplot(interval_counts, 
        main = "Taxa Counts by p-value Interval imputed", 
        xlab = "p-value Interval", 
        ylab = "Counts", 
        col = "transparent",  
        border = "black", 
        ylim = c(0, max(interval_counts) + 10),
        xaxt = "n",  
        space = 0)   

axis(side = 1, at = c(0,2,4,6,8,10), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), col.axis = "black")

#########################################
run_pval <- function(otu,meta){
  taxa_pval <- run_Maaslin2_pval(otu,meta)
  result_taxa_pval <- taxa_pval$results
  taxa_pval <- result_taxa_pval[,c('feature','pval')]
  taxa_pval$interval <- cut(taxa_pval$pval, breaks = c(0.0,0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7,0.8, 0.9,1.0))
  interval_counts <- table(taxa_pval$interval)
  p <- barplot(interval_counts, 
               xlab = "p-value Interval", 
               ylab = "Counts", 
               col = "transparent",  
               border = "black", 
               ylim = c(0, max(interval_counts) + 10),
               xaxt = "n",  
               space = 0)   

  axis(side = 1, at = c(0,2,4,6,8,10), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), col.axis = "black")
  return(p)
}

####Import Zeller data
meta_zeller <- read.csv('D:/paper/crc_zeller/crc_zeller_metadata.tsv',sep='\t',row.names = 1)
otu_zeller_imputed <- read.csv('D:/paper/imputed_file_data_imp_zeller_realdata.csv',row.names = 1)
colnames(otu_zeller_imputed) <- gsub("\\.","-",colnames(otu_zeller_imputed))
zeller_pval <- run_pval(otu_zeller,meta_zeller) 
zeller_imputed_pval <- run_pval(otu_zeller_imputed,meta_zeller) 


#####Import Feng data
feng_pval <- run_pval(otu_feng,meta_feng) 
feng_imputed_pval <- run_pval(otuimpute_feng,meta_feng) 


#####Import Qin data
Qin_pval <- run_pval(otu_Qin,meta_Qin) 
Qin_imputed_pval <- run_pval(otu_Qin_imputed,meta_Qin) 



##############################Random Forest Cross Validation#########
####before imputation 
microdata <- t(otu_zeller)
taxa_list <- as.character(na.omit(taxa_zeller_all_Maaslin2_result$fdr_Maaslin2))
zellerdata <- microdata[,taxa_list]
zellerdata <- as.data.frame(zellerdata)
combined_data <- cbind(zellerdata, meta_zeller)

########################After imputation

microdata_imputed <- t(otu_zeller_imputed)
taxa_list_imputed <- as.character(na.omit(taxa_zeller_all_Maaslin2_result$fdr_Maaslin2_imputed))
zellerdata_imputed <- microdata_imputed[,taxa_list]
zellerdata_imputed <- as.data.frame(zellerdata_imputed)
combined_data_imputed <- cbind(zellerdata_imputed, meta_zeller)


#####################Random Forest Cross Validation
library(dplyr) 
library(data.table)
library(randomForest) 
library(caTools)
library(caret) 
library(pROC) 
library(ggplot2) 
library(ggpubr) 
library(ggprism) 


set.seed(1234) 
#before imputation
#split <- sample.split(combined_data_imputed$comparison, SplitRatio = 0.8) 
#after imputation
#split <- sample.split(combined_data_imputed$comparison, SplitRatio = 0.8) 
train_data <- subset(combined_data_imputed, split == TRUE)  
test_data <- subset(combined_data_imputed, split == FALSE)  


X_train <- train_data[, -11]
y_train <- as.factor(train_data[, 11])


model <- randomForest(x = X_train, y = y_train, ntree = 100)

print(model)


ctrl <- trainControl(method = "cv", number = 5)
grid <- expand.grid(mtry = c(2, 4, 6))  

rf_model <- train(x = X_train, y = y_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = grid)


print(rf_model)



grid <- expand.grid(mtry = c(2))  
modellist <- list()


for (ntree in c(100,200, 300)) {
  set.seed(123)
  fit <- train(x = X_train, y = y_train, method="rf", 
               metric="Accuracy", tuneGrid=grid, 
               trControl=ctrl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}


# compare results
results <- resamples(modellist)

summary(results)


final_model <- randomForest(x = X_train, y = y_train,mtry = 2,ntree = 200)
print(final_model)



X_test <- test_data[, -11]
y_test <- as.factor(test_data[, 11])
test_predictions <- predict(final_model, newdata = test_data)




confusion_matrix <- confusionMatrix(test_predictions, y_test)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]


print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))




confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("H","CRC")
rownames(confusion_matrix_df) <- c("H","CRC")
draw_data <- round(confusion_matrix_df / rowSums(confusion_matrix_df),2)
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)


ggplot(draw_data, aes(real,variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(value))) +
  scale_fill_gradient(low = "#F0F0F0", high = "#3575b5") +
  labs(x = "True", y = "Guess", title = "Confusion matrix") +
  theme_prism(border = T)+
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")




test_predictions <- predict(final_model, newdata = test_data,type = "prob")

roc_obj <- roc(response = y_test, predictor = test_predictions[, 2])
roc_auc <- auc(roc_obj)


roc_data <- data.frame(1 - roc_obj$specificities, roc_obj$sensitivities)


ggplot(roc_data, aes(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities)) +
  geom_line(color = "#0073C2FF", size = 1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "gray") +
  geom_text(aes(x = 0.8, y = 0.2, label = paste("AUC =", round(roc_auc, 2))), size = 4, color = "black") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_pubr() +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  ggtitle("ROC Curve") +
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme_prism(border = T)



y_test_numeric <- ifelse(y_test == levels(y_test)[2], 1, 0)
pr_curve <- pr.curve(scores.class0 = test_predictions[, 2], weights.class0 = y_test_numeric, curve = TRUE)
plot(pr_curve, col = "blue", xlab = "Recall", ylab = "Precision", main = "Precision-Recall Curve(after imputed)")


pr_auc <- auc(pr_curve)
print(paste("Precision-Recall AUC:", pr_auc))
