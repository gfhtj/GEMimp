simudata2000_zero <- read.csv('D:/paper/simulated_data/sequencing depth/simudata2000_zero.csv', row.names = 1)  


#### Relative Abundance
simrel2000 <- sim2000[["simulated_matrices"]][["rel"]]
simrel2000 <- as.matrix(simrel2000)
write.csv(simrel2000, file = 'D:/paper/simulated_data/sequencing depth/simrel2000.csv')

simurel2000_zero <- read.csv('D:/paper/simulated_data/sequencing depth/simurel2000_zero.csv', row.names = 1)  ## Geimpute corresponding imputed data
simurel2000_imputed <- read.csv('D:/paper/simulated_data/sequencing depth/simurel2000_imputed.csv', row.names = 1)  ## Geimpute corresponding imputed data
simurel2000 <- as.matrix(simurel2000)
simurel2000_imputed <- as.matrix(simurel2000_imputed)
simurel2000_zero <- as.matrix(simurel2000_zero)
# Calculate MSE
GE_rmse <- sqrt(mean((simurel2000_imputed - simurel2000)^2, na.rm = TRUE))
noimp_rmse <- sqrt(mean((simurel2000_zero - simurel2000)^2, na.rm = TRUE))

############ DrImpute
set.seed(13)
library(DrImpute)
simurel2000_dr <- DrImpute(as.matrix(simurel2000_zero))
# Got the data imputed by drimpute after 10% imputation, i.e., simurel_dr
write.csv(simurel_dr, file = 'D:/paper/simulated_data/sequencing depth/simurel_dr2000.csv')
DrImpute_rmse <- sqrt(mean((simurel2000_dr - simurel2000)^2, na.rm = TRUE))  # 0.01098564

###  SAVER
library(SAVER)
# simurel2000_zero <- as.matrix(simurel2000_zero)
# simurel2000_zero <- as.data.frame(simurel2000_zero)
simudata2000_zero <- as.data.frame(simudata2000_zero) 
cortex.saver2000 <- saver(simudata2000_zero, ncores = 1, estimates.only = TRUE)  
## Convert to relative abundance
#### Relative Abundance
otu <- t(cortex.saver2000)
# Calculate the OTU relative abundance for each sample
otu_rel <- otu / rowSums(otu)
# Normalize the relative abundance for each sample to ensure the sum is 1
otu_rel_normalized <- sweep(otu_rel, 1, rowSums(otu_rel), "/")
simurel2000_saver <- t(otu_rel_normalized) # Relative abundance data
# write.csv(simurelative, file = 'D:/paper/simulated_data/sequencing depth/simurel2000.csv')

write.csv(simurel2000_saver, file = 'D:/paper/simulated_data/sequencing depth/simurel_saver2000.csv')
saver_rmse <- sqrt(mean((simurel2000_saver - simrel2000)^2, na.rm = TRUE)) # 0.0112938

############ MAGIC
## Import magic imputed data
simurel2000_magic <- read.csv('D:/paper/simulated_data/sequencing depth/simurel_mgic2000.csv', row.names = 1)  ## Magic corresponding imputed data
simurel2000_magic <- as.matrix(simurel2000_magic)
magic_rmse <- sqrt(mean((simurel2000_magic - simrel2000)^2, na.rm = TRUE)) # 0.02851441

############### Additional mbimpute
### mbimpute
### For mbimpute imputation, a count matrix is needed
library(mbImpute)
library(glmnet)
library(Matrix)


simurel2000 <- read.csv("D:/paper/simulated_data/sequencing depth/simurel2000_zero.csv", row.names = 1)
# rawfile <- read.csv("D:/Reading Literature Sync Folder/mbImpute code/0825new/Yzi_0.5.csv", row.names = 1)
# rawfile <- floor(10^rawfile-1.01) # Already converted to count matrix, imputed below with mbimpute
imputed_count_matrix_list <- mbImpute(otu_tab = t(simurel2000), k=3)
# Keep the imputed log matrix
mbimp_file <- imputed_count_matrix_list$imp_count_mat_norm 
mbimp_file <- mbimp_file/1000000
mbimp_file <- as.matrix(mbimp_file)
simrel2000 <- as.matrix(simrel2000)

### Save it as a csv format
# write.csv(simurelative, file = 'D:/paper/simulated_data/sequencing depth/simurel2000.csv')
mbimpute_rmse <- sqrt(mean((simurel2000_mbimpute - simrel2000)^2, na.rm = TRUE)) # 0.02346817
mbimpute_rmse

simurel5000 <- read.csv("D:/paper/simulated_data/sequencing depth/simurel5000_zero.csv", row.names = 1)
imputed_count_matrix_list2 <- mbImpute(otu_tab = t(simurel5000), k=3)
# Keep the imputed log matrix
mbimp_file5000 <- imputed_count_matrix_list2$imp_count_mat_norm 
mbimp_file5000 <- mbimp_file5000/1000000
mbimp_file5000 <- as.matrix(mbimp_file5000)

simrel5000 <- as.matrix(simrel5000)
mbimpute_rmse <- sqrt(mean((t(mbimp_file5000) - simrel5000)^2, na.rm = TRUE)) # 0.01704186

simurel10000 <- read.csv("D:/paper/simulated_data/sequencing depth/simurel10000_zero.csv", row.names = 1)
imputed_count_matrix_list3 <- mbImpute(otu_tab = t(simurel10000), k=3)
# Keep the imputed log matrix
mbimp_file10000 <- imputed_count_matrix_list3$imp_count_mat_norm 
mbimp_file10000 <- mbimp_file10000/1000000
mbimp_file10000 <- as.matrix(mbimp_file10000)
simurel10000 <- as.matrix(simurel10000)
mbimpute_rmse10000 <- sqrt(mean((t(mbimp_file10000) - simurel10000)^2, na.rm = TRUE)) # 0.01705645
mbimpute_rmse10000