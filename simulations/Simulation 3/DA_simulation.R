#################simulated data########################
library(SparseDOSSA2)
set.seed(2024)
metadata <- matrix(rbinom(n = 100, size = 1, prob = 0.5), nrow = 1, ncol = 100)
metadata <- t(metadata)
colnames(metadata) <- 'comparison'
sim <- SparseDOSSA2(template = "Stool", n_sample = 100, n_feature= 150,new_features = TRUE,metadata_matrix = metadata, perc_feature_spiked_metadata = 0.3,spike_metadata = "abundance")
#params <- sim[["params"]][["feature_param"]]
spike <- sim[["spike_metadata"]][["feature_metadata_spike_df"]]
simudata <- sim[["simulated_data"]]
save(simudata,file = 'D:/paper/simulated_data/simudata.RData')
save(spike,file = 'D:/paper/simulated_data/predefined_taxa.RData')
save(metadata,file = 'D:/paper/simulated_data/simu_metadata.RData')
write.csv(simudata,file = 'D:/paper/simulated_data/simudata.csv')


#Several methods for differential abundance

###################Maaslin2#######################
otu <- t(simudata)
simudata <- as.data.frame(simudata)
metadata <- as.data.frame(metadata)
# Set the row name of metadata to match the column name of simudata
rownames(metadata) <- paste0("Sample", 1:nrow(metadata))
set.seed(1234)
ps <- phyloseq::phyloseq(phyloseq::otu_table(simudata, taxa_are_rows = TRUE), 
                         phyloseq::sample_data(metadata))
colnames(metadata) <- 'comparison'
#Both simudata and metadata here must be dataframes, and the row name of metadata must match the column name of simudata
library(phyloseq)
library(Maaslin2)
mas_1 <- Maaslin2(
  input_data = data.frame(otu_table(ps)),
  input_metadata = data.frame(sample_data(ps)),
  output = "./Maaslin2_simudata_output",
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
mas_res_df <- mas_1$results

fdr_mas <- mas_res_df %>%
  dplyr::filter(qval < 0.05)

dim(fdr_mas)


######################corncob#################
set.seed(1) 

#Delete rows with all zeros
flag <- rowSums(simudata[])>0
simudata <- simudata[flag,]
#colnames(countdata) <- gsub("\\.","-",colnames(countdata))
simudata <- as.matrix(simudata)
identical(colnames(simudata), rownames(metadata))
rows_to_keep <- intersect(colnames(simudata), rownames(metadata))
metadata <- metadata[rows_to_keep,,drop=F]
countdata <- countdata[,rows_to_keep]
identical(colnames(countdata), rownames(metadata))
library(corncob)

ps <- phyloseq::phyloseq(phyloseq::otu_table(simudata, taxa_are_rows = TRUE), 
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


####################deseq#####################
run_deseq <- function(ASV_table,groupings){
  dim(ASV_table)
  dim(groupings)
  identical(colnames(ASV_table), rownames(groupings))
  rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
  groupings <- groupings[rows_to_keep,,drop=F]
  ASV_table <- ASV_table[,rows_to_keep]
  identical(colnames(ASV_table), rownames(groupings))
  colnames(groupings)[1] <- "Groupings"
  library(DESeq2)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                        colData=groupings,
                                        design = ~ Groupings)
  dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")
  
  res <- results(dds_res, 
                 tidy=T, 
                 format="DataFrame",
                 contrast = c("Groupings","0","1"))
  head(res)
  
  DEG<-res
  logFC_cutoff<-2
  DEG$change<-as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                               ifelse(DEG$log2FoldChange>logFC_cutoff,"UP","DOWN"),
                               "NOT"))
  this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                       '\nThe number of up taxa is ',nrow(DEG[DEG$change =='UP',]) ,
                       '\nThe number of down taxa is ',nrow(DEG[DEG$change =='DOWN',]))
  DEG<-na.omit(DEG)
  library(ggplot2)
  ggplot(data=DEG,aes(x=log2FoldChange,
                      y=-log10(pvalue),
                      color=change))+
    geom_point(alpha=0.8,size=3)+
    labs(x="log2 fold change")+ ylab("-log10 pvalue")+
    ggtitle(this_title)+theme_bw(base_size = 20)+
    theme(plot.title = element_text(size=15,hjust=0.5),)+
    scale_color_manual(values=c('#a121f0','#bebebe','#ffad21')) -> p1
  p1+xlim(NA,10)+ylim(NA,30) -> p2
  
  library(patchwork)
  p1+p2
  DA <- DEG[DEG$change =='DOWN' | DEG$change =='UP',]
  return(DA)
  
}

colnames(metadata) <- 'Groupings'
metadata$Groupings<-factor(metadata$Groupings)
simu_deseq <- run_deseq(simudata,metadata)


#########################before imputation##############
#Identify taxa using three methods
taxa_Maaslin2 <- fdr_mas$feature
taxa_corncob <- fdr_corncob$fdr_corncob
taxa_deseq <- simu_deseq$row
#Pre-defined taxa
taxa_predefined <- spike$feature_spiked
#Calculate recall rate, accuracy rate, etc

#Maaslin2
commontaxa_Maaslin2 <- intersect(taxa_predefined, taxa_Maaslin2)
countaxa_Maaslin2 <- length(commontaxa_Maaslin2)
#corncob
commontaxa_corncob <- intersect(taxa_predefined, taxa_corncob)
countaxa_corncob <- length(commontaxa_corncob)
#deseq
commontaxa_deseq <- intersect(taxa_predefined, taxa_deseq)
countaxa_deseq <- length(commontaxa_deseq)

#maaslin2
precision_Maaslin2=round(countaxa_Maaslin2/length(taxa_predefined),5)
recall_Maaslin2=round(countaxa_Maaslin2/length(taxa_Maaslin2),5)
F1_Maaslin2=round(2*precision_Maaslin2*recall_Maaslin2/(precision_Maaslin2+recall_Maaslin2),5)
#corncob
precision_corncob=round(countaxa_corncob/length(taxa_predefined),5)
recall_corncob=round(countaxa_corncob/length(taxa_corncob),5)
F1_corncob=round(2*precision_corncob*recall_corncob/(precision_corncob+recall_corncob),5)
#deseq
precision_deseq=round(countaxa_deseq/length(taxa_predefined),5)
recall_deseq=round(countaxa_deseq/length(taxa_deseq),5)
F1_deseq=round(2*precision_deseq*recall_deseq/(precision_deseq+recall_deseq),5)



######################after imputation###########
#The next step is to use GEMimp for imputation
simudata_imputed <- read.csv('D:/paper/simulated_data/simudata_imputed.csv',row.names = 1)
simudata <- as.data.frame(simudata)
metadata <- as.data.frame(metadata)
#Using the same operation as above, significant taxa were obtained using three differential analysis methods(omitted)
#Identify taxa using three methods
taxa_Maaslin2_imputed <- fdr_mas_imputed$feature
taxa_corncob_imputed <- fdr_corncob_imputed$fdr_corncob
taxa_deseq_imputed <- simu_deseq_imputed$row