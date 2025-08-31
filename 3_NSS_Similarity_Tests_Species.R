# NSS: similarity tests, taxa
library(tidyverse)
library(readxl)
library(corrplot)
library(RColorBrewer)
theme_set(theme_bw())

#Read in trait data output from script 1
straits <- read.csv("TraitData_Species_NSS_18Feb25.csv")

# which species are present in all time bins?
p5_occ <- read.csv("./Feb25 species no CEE/Bin_500b/Species_PA_500b.csv", row.names = 1)
f5_occ <- read.csv("./Feb25 species no CEE/Bin_500_700/Species_PA_500_700.csv", row.names = 1)
f7_occ <- read.csv("./Feb25 species no CEE/Bin_700_900/Species_PA_700_900.csv", row.names = 1)
f9_occ <- read.csv("./Feb25 species no CEE/Bin_900_1100/Species_PA_900_1100.csv", row.names = 1)
f11_occ <- read.csv("./Feb25 species no CEE/Bin_1100_1300/Species_PA_1100_1300.csv", row.names = 1)
f13_occ <- read.csv("./Feb25 species no CEE/Bin_1300_1500/Species_PA_1300_1500.csv", row.names = 1)
f15_occ <- read.csv("./Feb25 species no CEE/Bin_1500_1700/Species_PA_1500_1700.csv", row.names = 1)
f17_occ <- read.csv("./Feb25 species no CEE/Bin_1700a/Species_PA_1700a.csv", row.names = 1)

taxmats <- list(p5_occ, f5_occ, f7_occ,f9_occ, f11_occ, f13_occ,  f15_occ,
                f17_occ)

matrices <- lapply(taxmats, function(df) {
  df[is.na(df)] <- 0
  df[] <- lapply(df, as.numeric)
  return(df)
})

#Identifying empty columns by matrix
column_sums <- lapply(matrices, colSums)

# Print all taxon column sums, by time bin
cat("Column sums for each matrix:\n")
for (i in seq_along(column_sums)) {
  cat(sprintf("Matrix %d:\n", i))
  print(column_sums[[i]])
  cat("\n")
}

#Species occurring at at least 10 sites
sum_matrix <- do.call(rbind, column_sums)

# Create a single vector with the sum of each column across all matrices
total_column_sums <- colSums(sum_matrix)
summary(total_column_sums) #remove bottom quantile
low_counts <- names(total_column_sums[total_column_sums < 8])

empty_columns <- lapply(matrices, function(mat) {
  # Find columns where the sum is zero
  zero_cols <- colSums(mat) == 0
  # Return names of empty columns
  colnames(mat)[zero_cols]
})

# To get a list of all unique empty column names across all matrices
all_empty_columns <- unique(unlist(empty_columns))

# To count how many times each column is empty across matrices
sp_empty_count <- table(unlist(empty_columns))

# Convert to a named vector for easier manipulation
empty_column_counts <- as.vector(sp_empty_count)
names(empty_column_counts) <- names(table(unlist(empty_columns)))

# Sort by number of time bins empty
empty_column_counts <- sort(empty_column_counts, decreasing = TRUE)

# Print the results
print(empty_column_counts)

#Remove species missing in more than 3 time bins
low_counts2 <- names(empty_column_counts[empty_column_counts>4]) #remove these 19 species
all_remove <- c(low_counts, low_counts2)
all_remove <- unique(all_remove)  #28 species less than 8 sites total, or missing from more than 4 time bins

p5_occ <- p5_occ[, !colnames(p5_occ) %in% all_remove]  #leaves 85 species
f5_occ <- f5_occ[, !colnames(f5_occ) %in% all_remove]
f7_occ <- f7_occ[, !colnames(f7_occ) %in% all_remove]
f9_occ <- f9_occ[, !colnames(f9_occ) %in% all_remove]
f11_occ <- f11_occ[, !colnames(f11_occ) %in% all_remove]
f13_occ <- f13_occ[, !colnames(f13_occ) %in% all_remove]
f15_occ <- f15_occ[, !colnames(f15_occ) %in% all_remove]
f17_occ <- f17_occ[, !colnames(f17_occ) %in% all_remove]

# How many sites per bin?
taxmats <- list(p5_occ, f5_occ, f7_occ,f9_occ, f11_occ, f13_occ,  f15_occ,
                f17_occ)
sitecount <- sapply(taxmats, nrow)
sitecount 

# Function to transform occurrence matrix to correlation matrix
transform_and_calculate_corr <- function(mat) {
  data_occ <- as.matrix(mat)
  data_occ[is.na(data_occ)] <- 0  #0s for numeric
  Matrix <- as.matrix(t(data_occ))  # Transposed
  sim2 <- Matrix / sqrt(rowSums(Matrix * Matrix)) #multiply
  sim2 <- sim2 %*% t(sim2)
  return(sim2)
}

#list of corr matrices for all time bins
sim2_matrices <- lapply(taxmats, transform_and_calculate_corr)

sim2_matrices[[3]][,1]

#Long format table
# Extract column and row names from the first matrix
col_names <- colnames(sim2_matrices[[1]])
row_names <- rownames(sim2_matrices[[1]])

# Initialize the dataframe with column names and row names
df <- data.frame(matrix(ncol = length(sim2_matrices) + 2, nrow = length(sim2_matrices[[6]])))
colnames(df) <- c("Taxon1", "Taxon2", paste0("Bin_", 1:length(sim2_matrices)))

#List of row, column taxon pairs
row_index_list <- list()

row_index <- 1
for (row in 1:nrow(sim2_matrices[[1]])) {
  for (col in 1:ncol(sim2_matrices[[1]])) {
    # Save the row and column names in the dataframe
    df[row_index, "Taxon1"] <- rownames(sim2_matrices[[1]])[row]
    df[row_index, "Taxon2"] <- colnames(sim2_matrices[[1]])[col]
    
    # Save the row index for this (row, col) pair
    row_index_list[[paste(row, col, sep = "_AND_")]] <- row_index
    row_index <- row_index + 1
  }
}

# Fill in the observed values from the matrices
for (i in 1:length(sim2_matrices)) {
  matuse <- sim2_matrices[[i]]
  for (row in 1:nrow(matuse)) {
    for (col in 1:ncol(matuse)) {
      row_index <- row_index_list[[paste(row, col, sep = "_AND_")]]
      df[row_index, i + 2] <- matuse[row, col]
    }
  }
}


# Identify Oceanodromous-pair, Potamodromous-pair, Ocean-Poto pairs 
colnames(species)
species <- straits %>% filter(GBIF_species %in% colnames(f17_occ))

df$FreshMarine <- NA
for (i in 1:nrow(df)) {
  taxon1 <- df$Taxon1[i]
  taxon2 <- df$Taxon2[i]
  # Get the OAS values for Taxon1 and Taxon2 from the family trait data
  oas_taxon1 <- straits[straits$GBIF_species == taxon1,10]
  oas_taxon2 <- straits[straits$GBIF_species == taxon2,10]
  
  # Check the conditions and assign the appropriate value to fresh_marine
  if (oas_taxon1 == "Oceanodromous" && oas_taxon2 == "Oceanodromous") {  #both fully marine
    df$FreshMarine[i] <- "Marine"
  } else if (oas_taxon1 == "Potamodromous" && oas_taxon2 == "Potamodromous") {  #both fully freshwater
    df$FreshMarine[i] <- "Fresh"
   } else if (oas_taxon1 == "Potamodromous" && oas_taxon2 == "Oceanodromous") {  #one of each
      df$FreshMarine[i] <- "FreshMarine"
   } else if (oas_taxon1 == "Oceanodromous" && oas_taxon2 == "Potamodromous") {  #one of each
     df$FreshMarine[i] <- "FreshMarine"
  } else {
    df$FreshMarine[i] <- "Other" 
  }
}

df %>% dplyr::count(FreshMarine)

#Shuffle matrices and repeat
# Function to shuffle the entries in each column of a matrix
shuffle_matrix <- function(mat) {
  apply(mat, 2, sample)
}

# Create a list to store the shuffled matrices
shuffled_matrices <- list()

# Generate 500 shuffled versions of each matrix
for (i in 1:length(taxmats)) {
  mat <- taxmats[[i]]
  rownames(mat) <- rownames(taxmats[[i]])
  for (j in 1:500) {
    shuffled_matrices[[paste0("Bin_", i, "_Shuffled_", j)]] <- shuffle_matrix(mat)
  }
}

#function for shuffled matrices' correlation
transform_and_calculate_corrnull <- function(mat) {
  #rownames(mat) <- raw$X
  #data_occ <- as.matrix(mat[, 2:ncol(mat)])  # Only taxa columns
  mat[is.na(mat)] <- 0  #0s for numeric
  Matrix <- as.matrix(t(mat))  # Transposed
  sim2 <- Matrix / sqrt(rowSums(Matrix * Matrix)) #multiply
  sim2 <- sim2 %*% t(sim2)
  return(sim2)
}

#Find correlation values for each shuffled matrix
sim2_shuffled_matrices <- lapply(shuffled_matrices, transform_and_calculate_corrnull)

# Create a list to store the corr values for each time bin
values_lists <- vector("list", length(taxmats))
names(values_lists) <- paste0("Bin_", 1:length(taxmats))

# Initialize the lists within values_lists
for (i in 1:length(values_lists)) {
  values_lists[[i]] <- list()
}

# Save into list the values from the shuffled matrices
for (i in 1:length(sim2_shuffled_matrices)) {
  matuse <- sim2_shuffled_matrices[[i]]
  bin_index <- ceiling(i / 500)  # Determine which time bins (first 500 are bin 1, etc.)
  for (row in 1:nrow(matuse)) {
    for (col in 1:ncol(matuse)) {
      key <- paste(rownames(matuse)[row], colnames(matuse)[col], sep = "_AND_") #save which two taxa go with value
      if (!is.null(values_lists[[bin_index]][[key]])) {
        values_lists[[bin_index]][[key]] <- c(values_lists[[bin_index]][[key]], matuse[row, col])
      } else {
        values_lists[[bin_index]][[key]] <- matuse[row, col]  #save value
      }
    }
  }
  if (i %in% c(1000,2000,3000)){
    print(i)
  }
}

# Add columns to store mean and sd values of null matrices for each bin
for (i in 1:length(taxmats)) {
  df[[paste0("Mean_Bin_", i)]] <- NA
  df[[paste0("SD_Bin_", i)]] <- NA
  df[[paste0("Significance_Bin", i)]] <- NA
}

# Calculate the mean and standard deviation for each (Taxon1, Taxon2) pair in each time bin
for (bin_index in 1:length(values_lists)) {
  for (key in names(values_lists[[bin_index]])) {
    values <- values_lists[[bin_index]][[key]]
    taxon1 <- strsplit(key, "_AND_")[[1]][1]
    taxon2 <- strsplit(key, "_AND_")[[1]][2]
    mean_value <- mean(values)
    sd_value <- sd(values)
    
    # Save mean and SD of null matrices into original dataframe
    row_index <- which(df$Taxon1 == taxon1 & df$Taxon2 == taxon2)
    df[row_index, paste0("Mean_Bin_", bin_index)] <- mean_value
    df[row_index, paste0("SD_Bin_", bin_index)] <- sd_value
    
    # Calculate the z-score of observed vs. null (in terms of SD's from mean)
    observed_value <- df[row_index, paste0("Bin_", bin_index)]
    if (!is.na(observed_value) && !is.na(mean_value) && !is.na(sd_value)) {
      z_score <- (observed_value - mean_value) / sd_value
      df[row_index, paste0("Significance_Bin", bin_index)] <- z_score
    }
  }
}

# remove pairs where two sites are the same
nrow(df)
colnames(df)
df <-  df %>% filter(!(Taxon1==Taxon2))
nrow(df)

# Consider results-- how many pairs of species have significant z score, relative to null?
head(df)
df %>% count(Significance_Bin1 > 1.96)
df %>% count(Significance_Bin2 > 1.96)
df %>% count(Significance_Bin3 > 1.96)
df %>% count(Significance_Bin4 > 1.96)
df %>% count(Significance_Bin5 > 1.96)
df %>% count(Significance_Bin6 > 1.96)
df %>% count(Significance_Bin7 > 1.96)
df %>% count(Significance_Bin8 > 1.96)

summary_df <- df %>% group_by(FreshMarine) %>% 
  summarize(
    Bin_1_over_1.96 = sum(Significance_Bin1 >= 1.96, na.rm = TRUE),
    Bin_1_under_1.96 = sum(Significance_Bin1 < 1.96, na.rm = TRUE),
    Bin_2_over_1.96 = sum(Significance_Bin2 >= 1.96, na.rm = TRUE),
    Bin_2_under_1.96 = sum(Significance_Bin2 < 1.96, na.rm = TRUE),
    Bin_3_over_1.96 = sum(Significance_Bin3 >= 1.96, na.rm = TRUE),
    Bin_3_under_1.96 = sum(Significance_Bin3 < 1.96, na.rm = TRUE),
    Bin_4_over_1.96 = sum(Significance_Bin4 >= 1.96, na.rm = TRUE),
    Bin_4_under_1.96 = sum(Significance_Bin4 < 1.96, na.rm = TRUE),
    Bin_5_over_1.96 = sum(Significance_Bin5 >= 1.96, na.rm = TRUE),
    Bin_5_under_1.96 = sum(Significance_Bin5 < 1.96, na.rm = TRUE),
    Bin_6_over_1.96 = sum(Significance_Bin6 >= 1.96, na.rm = TRUE),
    Bin_6_under_1.96 = sum(Significance_Bin6 < 1.96, na.rm = TRUE),
    Bin_7_over_1.96 = sum(Significance_Bin7 >= 1.96, na.rm = TRUE),
    Bin_7_under_1.96 = sum(Significance_Bin7 < 1.96, na.rm = TRUE),
    Bin_8_over_1.96 = sum(Significance_Bin8 >= 1.96, na.rm = TRUE),
    Bin_8_under_1.96 = sum(Significance_Bin8 < 1.96, na.rm = TRUE)
  ) %>%
  mutate(
    Ratio_Bin_1 = Bin_1_over_1.96 / Bin_1_under_1.96,
    Ratio_Bin_2 = Bin_2_over_1.96 / Bin_2_under_1.96,
    Ratio_Bin_3 = Bin_3_over_1.96 / Bin_3_under_1.96,
    Ratio_Bin_4 = Bin_4_over_1.96 / Bin_4_under_1.96,
    Ratio_Bin_5 = Bin_5_over_1.96 / Bin_5_under_1.96,
    Ratio_Bin_6 = Bin_6_over_1.96 / Bin_6_under_1.96,
    Ratio_Bin_7 = Bin_7_over_1.96 / Bin_7_under_1.96,
    Ratio_Bin_8 = Bin_8_over_1.96 / Bin_8_under_1.96
  )

#save results
write.csv(df, "Results_PairwiseCorr_Species_500null_8timebins.csv")
write.csv(summary_df, "SummaryResults_Significance_SpeciesPairs_500null_8timebins.csv")

# Plotting ####
df <- read.csv("Results_PairwiseCorr_Species_500null_8timebins.csv", row.names = 1)
summary_df <- read.csv("SummaryResults_Significance_SpeciesPairs_500null_8timebins.csv", row.names = 1)
summary_df <-  summary_df %>% mutate(
  Perc_Bin_1 = Bin_1_over_1.96 / (Bin_1_over_1.96+Bin_1_under_1.96),
  Perc_Bin_2 = Bin_2_over_1.96 / (Bin_2_over_1.96 + Bin_2_under_1.96),
  Perc_Bin_3 = Bin_3_over_1.96 / (Bin_3_over_1.96 + Bin_3_under_1.96),
  Perc_Bin_4 = Bin_4_over_1.96 / (Bin_4_over_1.96 + Bin_4_under_1.96),
  Perc_Bin_5 = Bin_5_over_1.96 / (Bin_5_over_1.96 + Bin_5_under_1.96),
  Perc_Bin_6 = Bin_6_over_1.96 / (Bin_6_over_1.96 + Bin_6_under_1.96),
  Perc_Bin_7 = Bin_7_over_1.96 / (Bin_7_over_1.96 + Bin_7_under_1.96),
  Perc_Bin_8 = Bin_8_over_1.96 / (Bin_8_over_1.96 + Bin_8_under_1.96)
) %>% mutate(
  Tot_Bin_1 = (Bin_1_over_1.96+Bin_1_under_1.96),
  Tot_Bin_2 =  (Bin_2_over_1.96 + Bin_2_under_1.96),
  Tot_Bin_3 =  (Bin_3_over_1.96 + Bin_3_under_1.96),
  Tot_Bin_4 =  (Bin_4_over_1.96 + Bin_4_under_1.96),
  Tot_Bin_5 = (Bin_5_over_1.96 + Bin_5_under_1.96),
  Tot_Bin_6 =  (Bin_6_over_1.96 + Bin_6_under_1.96),
  Tot_Bin_7 =  (Bin_7_over_1.96 + Bin_7_under_1.96),
  Tot_Bin_8 =  (Bin_8_over_1.96 + Bin_8_under_1.96)
)
look <- summary_df[,c(1, 34:41)] #how many pairs in each bin?

long_summary_df <- summary_df %>%
  pivot_longer(cols = starts_with("Ratio_Bin_"), names_to = "Bin", values_to = "Ratio") %>%
  mutate(Bin = as.numeric(gsub("Ratio_Bin_", "", Bin)))

long_summary_df_perc <- summary_df %>%
  pivot_longer(cols = starts_with("Perc_Bin_"), names_to = "Bin", values_to = "Percentage_Significant") %>%
  mutate(Bin = as.numeric(gsub("Perc_Bin_", "", Bin)))

# Plot across bins, broken down by pairwise ocean access score
ggplot(long_summary_df, aes(x = Bin, y = Ratio, color = FreshMarine, group = FreshMarine)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Marine" = "blue2", "Fresh" = "forestgreen", 
                                "FreshMarine" = "magenta2", "Other"= "grey2"))+
  labs(title = "Ratio of Significant to Nonsignificant vs. Null Correlation",
       x = "Bin Number",
       y = "Ratio",
       color = "Pairwise Shared Habitat?") +
  theme_bw()

ggplot(long_summary_df_perc, aes(x = Bin, y = Percentage_Significant, color = FreshMarine, group = FreshMarine)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Marine" = "blue2", "Fresh" = "forestgreen", 
                                "FreshMarine" = "magenta2", "Other"= "grey2"))+
  labs(
       x = "Bin Number",
       y = "Percentage Significant",
       color = "Pairwise Shared Habitat?") +
  theme_bw()

write.csv(long_summary_df_perc, "Species_FreshMarine_BinSummary_Percent.csv")

#Example plot of null distribution vs. observed
taxon1 <- "Osmerus_eperlanus"
taxon2 <- "Scomber_scombrus"
bin_index <- 3

key <- paste(taxon1, taxon2, sep = "_AND_")
values <- values_lists[[bin_index]][[key]]
plot_df <- data.frame(Value = values)

#actual corr
observed_value <- df[df$Taxon1 == taxon1 & df$Taxon2 == taxon2, paste0("Bin_", bin_index)]
# mean null corr
mean_value <- df[df$Taxon1 == taxon1 & df$Taxon2 == taxon2, paste0("Mean_Bin_", bin_index)]

ggplot(plot_df, aes(x = Value)) +
  geom_density(color = "blue", fill = "blue", alpha = 0.3) +  # Smooth curve
  geom_vline(aes(xintercept = observed_value), color = "red", linetype = "dashed", size = 1) +  # Red vertical line
  labs(title = paste("Distribution of Correlation Scores for", taxon1, "and", taxon2),
       x = "Value",
       y = "Density") +
  annotate("text", x = 0.2, y = 15, label = paste( "Mean:", round(mean_value, 2)), 
           color = "black", vjust = -0.5)+
  annotate("text", x = 0.2, y = 16, label = paste("Observed:", round(observed_value, 2)), 
           color = "black", vjust = -0.9) +
  theme_bw()


#  Example plot of single pair over time
null_res <- read.csv("Results_PairwiseCorr_Species_500null_8timebins.csv", row.names = 1)
colnames(null_res)
# Find observed values for each bin
observed_pair <- null_res[null_res$Taxon1 == taxon1 & null_res$Taxon2 == taxon2, 3:10]
# Find Significance of these observed values relative to null
sig_pair <- null_res[null_res$Taxon1 == taxon1 & null_res$Taxon2 == taxon2, c(11,14,17,20,23,26,29,32)]
# Find mean similarity across shuffled matrices
binn_pair <- null_res[null_res$Taxon1 == taxon1 & null_res$Taxon2 == taxon2, c(11,15,18,21,24,27,30,33)]

pdat <- cbind(c(1:8),as.numeric(observed_pair[1,]), as.numeric(binn_pair[1,]))
pdat <- as.data.frame(pdat)

#Plot: pink is observed, grey is null
ggplot(pdat)+
  geom_line(aes(x=V1, y=V2), color="magenta3", size=2)+
  geom_line(aes(x=V1, y=V3), color="grey", size=2)+ xlab("Time")+ ylab("Similarity")+
  labs(title = paste("Similarity Scores for", taxon1, "and", taxon2))+
  theme_bw()
