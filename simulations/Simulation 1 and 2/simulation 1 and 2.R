###########zero-inflated data###################
# Calculate the non-zero average abundance
non_zero_means_2000_2 <- apply(logdata2000_2, 2, function(x) {
  # Exclude log10(1.01) non-zero values
  non_zero_values <- x[x != log10(1.01)]
  # If there are non-zero values, calculate the mean; otherwise return 0
  if (length(non_zero_values) > 0) {
    mean(non_zero_values, na.rm = TRUE)
  } else {
    0
  }
})

# Check for NA values
na_check <- any(is.na(non_zero_means_2000_2))
print(na_check)

# If necessary, replace NA values with 0
non_zero_means_2000_2[is.na(non_zero_means_2000_2)] <- 0
# Calculate the adjusted probability of imputing 0, zj
zj_adjusted_2000_2 <- 1 - non_zero_proportions_2000_2 / (non_zero_means_2000_2 + 1)
any(is.na(zj_adjusted_2000_2))
# Ensure zj is between 0 and 1
zj_adjusted_2000_2 <- pmin(pmax(zj_adjusted_2000_2, 0), 1)
any(is.na(zj_adjusted_2000_2))
# Generate a binomial distribution Iij
set.seed(123)
Iij_adjusted_2000_2 <- matrix(rbinom(n = nrow(simudata2000_2) * ncol(simudata2000_2), size = 1, prob = zj_adjusted_2000_2), 
                              nrow = nrow(simudata2000_2), ncol = ncol(simudata2000_2))
# Check for NA values
if (any(is.na(Iij_adjusted_2000_2))) {
  warning("Iij_adjusted_2000_2 contains NA values")
  #zj_adjusted_2000_2[is.na(zj_adjusted_2000_2)] <- 0 
}
# Multiply Iij by the log-transformed matrix to impute 0s
filled_logdata_adjusted_2000_2 <- logdata2000_2 * Iij_adjusted_2000_2





###########simulation 1 and 2#############
####function#####
process_data <- function(raw_data, imputed_data) {
  # Extract the numeric columns from the data frame
  raw_values <- as.matrix(raw_data)
  imputed_values <- as.matrix(imputed_data)
  
  # Calculate the Pearson correlation coefficient
  correlation_matrix <- cor(imputed_values, raw_values)
  
  microbial_names <- colnames(raw_data)
  s_vector <- microbial_names

  bacteria_names <- c()

  for (s in s_vector) {
    split_str <- strsplit(s, "__")[[1]]
    split_str <- strsplit(split_str[2], "_")[[1]]
    # Extract the second part, i.e., the bacterial name
    bacteria_name <- split_str[1]
    bacteria_names <- c(bacteria_names, bacteria_name)
  }
  
  name_counts <- table(bacteria_names)
  sorted_counts <- sort(name_counts, decreasing = TRUE)
  
  row_name <- rownames(correlation_matrix)
  last_counts <- as.matrix(sorted_counts[sorted_counts >= 14])
  target_columns <- rownames(last_counts)
 
  matched_elements <- c()
  
  for (element in row_name) {
    if (any(grepl(paste(target_columns, collapse = "|"), element))) {
      matched_elements <- c(matched_elements, element)
    }
  }
  matched_values <- numeric(length(matched_elements))
  
  for (i in 1:length(matched_elements)) {
    col_name <- matched_elements[i]
    if (col_name %in% colnames(correlation_matrix)) {
      matched_values[i] <- correlation_matrix[col_name, col_name]
    }
  }
  
  matched_df <- data.frame(matched_elements, matched_values)
  rownames(matched_df) <- matched_df$matched_elements
  matched_df$matched_elements <- NULL
  matched_df <- na.omit(matched_df)
  matched_df$group <- gsub("^s__([^_]+).*", "\\1", rownames(matched_df))
  result <- aggregate(matched_values ~ group, data = matched_df, FUN = mean)
  return(result)
}

process_data <- function(raw_data, imputed_data) {
  raw_values <- as.matrix(raw_data)
  imputed_values <- as.matrix(imputed_data)
  
  # Calculate the Spearman correlation coefficient
  correlation_matrix <- cor(imputed_values, raw_values, method = "spearman")
  
  microbial_names <- colnames(raw_data)
  s_vector <- microbial_names
  bacteria_names <- c()

  for (s in s_vector) {
    split_str <- strsplit(s, "__")[[1]]
    split_str <- strsplit(split_str[2], "_")[[1]]
    bacteria_name <- split_str[1]
    bacteria_names <- c(bacteria_names, bacteria_name)
  }
  
  name_counts <- table(bacteria_names)
  sorted_counts <- sort(name_counts, decreasing = TRUE)
  
  row_name <- rownames(correlation_matrix)
  #to keep the same with PCC
  target_columns <- as.matrix(c('Bacteroides', 'Clostridium', 'Eubacterium', 'Lachnospiraceae', 'Lactobacillus', 'Prevotella', 'Streptococcus'))
  

  matched_elements <- c()

  for (element in row_name) {
    if (any(grepl(paste(target_columns, collapse = "|"), element))) {
      matched_elements <- c(matched_elements, element)
    }
  }
  
  matched_values <- numeric(length(matched_elements))
  
  for (i in 1:length(matched_elements)) {
    col_name <- matched_elements[i]
    if (col_name %in% colnames(correlation_matrix)) {
      matched_values[i] <- correlation_matrix[col_name, col_name]
    }
  }
  matched_df <- data.frame(matched_elements, matched_values)
  
  rownames(matched_df) <- matched_df$matched_elements
  matched_df$matched_elements <- NULL
  matched_df <- na.omit(matched_df)
  matched_df$group <- gsub("^s__([^_]+).*", "\\1", rownames(matched_df))
  result <- aggregate(matched_values ~ group, data = matched_df, FUN = mean)
  return(result)
}

###############################
#simulation 1
logdata0915_2 <- read.csv( 'D:/github/simulations/Simulation 1/raw_file.csv',row.names = 1) #170*469
rawdata0915_2 <- floor(10^logdata0915_2-1.01) 
rawdata0915_2[rawdata0915_2 == -1] = 0

# Multiply Iij by the log-transformed matrix to impute zeros
#filled_logdata_adjusted_2 <- logdata0915_2 * Iij_adjusted_2   #from zero-inflated data
# Convert the results back to the original data format
filled_rawdata_adjusted_2 <- floor(10^filled_logdata_adjusted_2 - 1.01)
filled_rawdata_adjusted_2[filled_rawdata_adjusted_2 == -1] <- 0

# Calculate the number of zeros after imputation
zero_count_after_imputation_2 <- sum(filled_logdata_adjusted_2 == 0, na.rm = TRUE)
total_elements_2 <- nrow(filled_logdata_adjusted_2) * ncol(filled_logdata_adjusted_2)
zero_proportion_after_imputation_2 <- zero_count_after_imputation_2 / total_elements_2

# Print the proportion of zeros
zero_proportion_after_imputation_2 # 0.7797818


##simulation 2
logdata0915_2 <- read.csv('D:/github/simulations/Simulation 1/raw_file.csv', row.names = 1) 
y_sim = logdata0915_2
filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
  y = y_sim[,col_i]
  n = length(y)
  nz <- sum(y <= (log10(1.01) + 1e-6))
  pz = 1 - nz/n
  test = pz - 1.96 * sqrt(pz * (1-pz)/n)
  if(nz == n || test <= 0){
    return(0)
  }else{
    return(1)
  }
})
filter_vec <- which(unlist(filter) == 1)
y_sim = y_sim[, filter_vec]
write.csv(y_sim, file = './high_data.csv') #  170*294
################################mbimpute
library(mbImpute)
library(glmnet)
library(Matrix)



imputed_count_matrix_list0915_2 <- mbImpute(otu_tab = filled_rawdata_adjusted_2,k=3, unnormalized = T)

mbimp_file0915_a <- imputed_count_matrix_list0915_2$imp_count_mat_norm

scale <- rowSums(rawdata0915_2) / 10^6
rawdata0915_2 <- rawdata0915_2 / scale
result0915_a <- process_data2(rawdata0915_2, mbimp_file0915_a)
print(result0915_a)

result0915_a <- process_data_spearman(rawdata0915_2, mbimp_file0915_a)
print(result0915_a)




######################drimpute
library(DrImpute)
drimpute_data0915_2 <- DrImpute(as.matrix(filled_logdata_adjusted_2))

drimpute_data0915_2 <- floor(10^drimpute_data0915_2-1.01) 
scale <- rowSums(drimpute_data0915_2) / 10^6
drimpute_data0915_2 <- drimpute_data0915_2 / scale
result0915_dr_a <- process_data2(rawdata0915_2, drimpute_data0915_2)
print(result0915_dr_a)

result0915_dr_a <- process_data_spearman(rawdata0915_2, drimpute_data0915_2)
print(result0915_dr_a)


##########################saver
library(SAVER)
imputed_data_bysaver0915_2 <- saver(filled_logdata_adjusted_2, ncores = 12, estimates.only = TRUE)

imputed_data_bysaver0915_2 <- floor(10^imputed_data_bysaver0915_2-1.01) 
scale <- rowSums(imputed_data_bysaver0915_2) / 10^6
imputed_data_bysaver0915_2 <- imputed_data_bysaver0915_2 / scale
result0915_saver_a <- process_data2(rawdata0915_2, imputed_data_bysaver0915_2)
print(result0915_saver_a)

result0915_saver_a <- process_data_spearman(rawdata0915_2, imputed_data_bysaver0915_2)
print(result0915_saver_a)


################################GEMimp
imputed_data_byGEMimp0915_2 <- read.csv('D:/GEMimp/0915/outptfile_bygemimp_2.csv',row.names = 1) 

imputed_data_byGEMimp0915_2 <- floor(10^imputed_data_byGEMimp0915_2-1.01) 
imputed_data_byGEMimp0915_2[imputed_data_byGEMimp0915_2 == imputed_data_byGEMimp0915_2[1,1]]=0
scale <- rowSums(imputed_data_byGEMimp0915_2) / 10^6
imputed_data_byGEMimp0915_2 <- imputed_data_byGEMimp0915_2 / scale
result0915_GEMimp_a <- process_data2(rawdata0915_2, imputed_data_byGEMimp0915_2)
print(result0915_GEMimp_a)

result0915_GEMimp_a <- process_data_spearman(rawdata0915_2, imputed_data_byGEMimp0915_2)
print(result0915_GEMimp_a)


################################magic
imputed_data_bymagic0915_2 <- read.csv('D:/GEMimp/0915/outptfile_bymagic_2.csv',row.names = 1)

imputed_data_bymagic0915_2 <- floor(10^imputed_data_bymagic0915_2-1.01) 
scale <- rowSums(imputed_data_bymagic0915_2) / 10^6
imputed_data_bymagic0915_2 <- imputed_data_bymagic0915_2 / scale
result0915_magic_a <- process_data2(rawdata0915_2, imputed_data_bymagic0915_2)
print(result0915_magic_a)

result0915_magic_a <- process_data_spearman(rawdata0915_2, imputed_data_bymagic0915_2)
print(result0915_magic_a)

