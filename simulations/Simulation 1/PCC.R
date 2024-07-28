setwd('C:/Users/szw/Desktop/GE-Impute-main/src')

raw_data <- read.csv("./raw_file.csv",row.names = 1) #raw complete data
imputed_data <- read.csv("./imputed_file.csv",row.names = 1)#imputed data using GEMimp


raw_values <- as.matrix(raw_data)
imputed_values <- as.matrix(imputed_data)

correlation_matrix <- cor(imputed_values, raw_values)

print(correlation_matrix)


microbial_names <- colnames(raw_data)

s_vector <- microbial_names

bacteria_names <- c()


for (s in s_vector) {
  # #Use the strsplit function to split a string, using two underscores as separators
  split_str <- strsplit(s, "__")[[1]]
  split_str <- strsplit(split_str[2], "_")[[1]]
  bacteria_name <- split_str[1]
  bacteria_names <- c(bacteria_names, bacteria_name)
}

print(bacteria_names)
name_counts <- table(bacteria_names)
sorted_counts <- sort(name_counts, decreasing = TRUE)
print(sorted_counts)

row_name <- rownames(correlation_matrix)
last_counts <- as.matrix(sorted_counts[sorted_counts>=14])
target_columns <-rownames(last_counts)

print(matched_elements)


############################
matched_elements <- c()

for (element in row_name) {
  # #Use the grepl function to check if the current element is included in any of the target_columns elements
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

print(matched_values)
print(matched_elements)


matched_df <- data.frame(matched_elements, matched_values)
rownames(matched_df) <- matched_df$matched_elements
matched_df$matched_elements <- NULL

print(matched_df)
#Remove certain rows that are NA
matched_df <- na.omit(matched_df)


# Extract category information from row names (for example, extract Lactobacillus)
matched_df$group <- gsub("^s__([^_]+).*", "\\1", rownames(matched_df))

# Calculate the average value by category using the aggregate function
result <- aggregate(matched_values ~ group, data = matched_df, FUN = mean)

print(result)