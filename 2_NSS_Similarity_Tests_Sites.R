# NSS: similarity tests, sites
# Goal is to compare urban and rural sites from same and different regions

library(tidyverse)
library(readxl)
library(corrplot)
library(RColorBrewer)
theme_set(theme_bw())

# set working directory where outputs from script 1 are
sitedata <- read.csv("Site_Metadata_SpeciesLevel_noCEE_180225.csv", row.names = 1)
# Input occurrence similarity
raw <- read.csv("Family_PA_noCEE_18Feb25.csv")

colnames(sitedata)
bins <- unique(sitedata$time_bins2) #find all time bins
bins
# select sites in each time bin
thisbin1 <- sitedata[sitedata$time_bins2 %in% bins[8],]  
p5_occ <- raw %>% filter(X %in% thisbin1$DB_Assemblage_ID)
thisbin2 <- sitedata[sitedata$time_bins2 %in% bins[3],]  
f5_occ <- raw %>% filter(X %in% thisbin2$DB_Assemblage_ID)
thisbin3 <- sitedata[sitedata$time_bins2 %in% bins[4],]  
f7_occ <- raw %>% filter(X %in% thisbin3$DB_Assemblage_ID)
thisbin4 <- sitedata[sitedata$time_bins2 %in% bins[7],]  
f9_occ <- raw %>% filter(X %in% thisbin4$DB_Assemblage_ID)
thisbin5 <- sitedata[sitedata$time_bins2 %in% bins[6],]  
f11_occ <- raw %>% filter(X %in% thisbin5$DB_Assemblage_ID)
thisbin6 <- sitedata[sitedata$time_bins2 %in% bins[1],]  
f13_occ <- raw %>% filter(X %in% thisbin6$DB_Assemblage_ID)
thisbin7 <- sitedata[sitedata$time_bins2 %in% bins[5],]  
f15_occ <- raw %>% filter(X %in% thisbin7$DB_Assemblage_ID)
thisbin8 <- sitedata[sitedata$time_bins2 %in% bins[2],]  
f17_occ <- raw %>% filter(X %in% thisbin8$DB_Assemblage_ID)

# List of individual matrices for each time bin
sitmats <- list(p5_occ, f5_occ, f7_occ,f9_occ, f11_occ, f13_occ,  f15_occ,
                f17_occ)

matrices <- lapply(sitmats, function(df) {
  rownames(df) <- df$X
  df[is.na(df)] <- 0
  df <- df[,2:ncol(df)]
  df[] <- lapply(df, as.numeric)
  return(df)
})

#Identifying count of species per site
row_sums <- lapply(matrices, rowSums)

#remove sites with <3 species
p5_remove <- rownames(matrices[[1]])[row_sums[[1]]<3]
f5_remove <- rownames(matrices[[2]])[row_sums[[2]]<3]
f7_remove <- rownames(matrices[[3]])[row_sums[[3]]<3]
f9_remove <- rownames(matrices[[4]])[row_sums[[4]]<3]
f11_remove <- rownames(matrices[[5]])[row_sums[[5]]<3]
f13_remove <- rownames(matrices[[6]])[row_sums[[6]]<3]
f15_remove <- rownames(matrices[[7]])[row_sums[[7]]<3]
f17_remove <- rownames(matrices[[8]])[row_sums[[8]]<3]

p5_occ <- p5_occ %>% filter(!(X %in% p5_remove)) 
f5_occ <- f5_occ %>% filter(!(X %in% f5_remove))
f7_occ <- f7_occ %>% filter(!(X %in% f7_remove))
f9_occ <- f9_occ %>% filter(!(X %in% f9_remove))
f11_occ <- f11_occ %>% filter(!(X %in% f11_remove))
f13_occ <- f13_occ %>% filter(!(X %in% f13_remove))
f15_occ <- f15_occ %>% filter(!(X %in% f15_remove))
f17_occ <- f17_occ %>% filter(!(X %in% f17_remove))

# List for all time bins without species with 1 or 2 occurrences in bin
simats <- list(p5_occ, f5_occ, f7_occ,f9_occ, f11_occ, f13_occ,  f15_occ,
                f17_occ)

# Function to transform occurrence matrix to correlation matrix
transform_and_calculate_corr <- function(mat) {
  rownames(mat) <- mat$X
  data_occ <- as.matrix(mat[,2:ncol(mat)])
  data_occ[is.na(data_occ)] <- 0  #0s for numeric
  Matrix <- as.matrix(data_occ)  # For sites, not transposed
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix)) #multiply
  sim <- sim %*% t(sim)
  return(sim)
}

#list of corr matrices for all time bins' sites, pairwise
sim_matrices <- lapply(simats, transform_and_calculate_corr)

head(sim_matrices[[2]])

#Transform to long format table
df_list <- list()  #create list of sim data frames
row_index_lists <- list()   #create list of row indices for each pair of sites

# Extract column and row names, different for each bin's sites
for (p in 1:8){
  col_names <- colnames(sim_matrices[[p]])
  row_names <- rownames(sim_matrices[[p]])
  # Initialize the dataframe with column names and row names
  df_list[[p]] <- data.frame(matrix(ncol = 3, nrow = length(sim_matrices[[p]])))
  colnames(df_list[[p]]) <- c("Site1", "Site2", paste0("Bin_", p,"_obs"))
  #List of row, column site pairs
  row_index_lists[[p]] <- list()
  row_index <- 1
  for (row in 1:nrow(sim_matrices[[p]])) {  #number of site in row
    for (col in 1:ncol(sim_matrices[[p]])) {  #number of site in column
      # Save the row and column names in the dataframe
      df_list[[p]][row_index, "Site1"] <- rownames(sim_matrices[[p]])[row]
      df_list[[p]][row_index, "Site2"] <- colnames(sim_matrices[[p]])[col]
      # Save the row index for this (row, col) pair
      row_index_lists[[p]][[paste(row, col, sep = "_AND_")]] <- row_index
      row_index <- row_index + 1
    }
  }
  
}

# Fill in the observed sim values from the input matrix
for (l in 1:8){
  for (row in 1:nrow(sim_matrices[[l]])) {
    for (col in 1:ncol(sim_matrices[[l]])) {
      row_index <- row_index_lists[[l]][[paste(row, col, sep = "_AND_")]]
      df_list[[l]][row_index, 3] <- sim_matrices[[l]][row, col] #save value from sim matrix into df in list
    }
  }
}

#identify pairs of sites in same and different countries that are rural and urban
for (o in 1:8){
df_list[[o]]$UrbanRural <- NA
df_list[[o]]$UrbanRural_noreg <- NA
for (i in 1:nrow(df_list[[o]])) {
  site1 <- df_list[[o]]$Site1[i]
  site2 <- df_list[[o]]$Site2[i]
  # Get the OAS values for Taxon1 and Taxon2 from the family trait data
  reg_site1 <- sitedata[sitedata$DB_Assemblage_ID == site1, 6]  #what region is site in
  reg_site2 <- sitedata[sitedata$DB_Assemblage_ID == site2, 6]
  urb_site1 <- sitedata[sitedata$DB_Assemblage_ID == site1,9] #Urban or not
  urb_site2 <- sitedata[sitedata$DB_Assemblage_ID == site2,9]
  count_site1 <- sitedata[sitedata$DB_Assemblage_ID == site1,5] #country name
  count_site2 <- sitedata[sitedata$DB_Assemblage_ID == site2,5]
  
  # Check the conditions and assign the appropriate values to pairs
  if(!(is.na(urb_site1)) && !(is.na(urb_site2))){
 if (urb_site1 == T && urb_site2== T){
      if (reg_site1 == "Western Europe" && reg_site2 == "Western Europe"){
        df_list[[o]][i,4] <- "BothUrban_WE"
      }
      else if (reg_site1 == "Britian & Ireland" && reg_site2 == "Britian & Ireland"){
        df_list[[o]][i,4] <- "BothUrban_BI"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Scandinavia"){
        df_list[[o]][i,4] <- "BothUrban_SC"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Western Europe" |
               reg_site1 == "Western Europe" && reg_site2 == "Scandinavia") {
        df_list[[o]][i,4] <- "BothUrban_SC_WE"
      }
      else if (reg_site1 == "Britian & Ireland" && reg_site2 == "Western Europe" |
               reg_site1 == "Western Europe" && reg_site2 == "Britian & Ireland") {
        df_list[[o]][i,4] <- "BothUrban_BI_WE"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Britian & Ireland" |
               reg_site1 == "Britian & Ireland" && reg_site2 == "Scandinavia") {
        df_list[[o]][i,4] <- "BothUrban_SC_BI"
      }
    }
    if (urb_site1 == F && urb_site2== F){
      if (reg_site1 == "Western Europe" && reg_site2 == "Western Europe"){
        df_list[[o]][i,4] <- "BothRural_WE"
      }
      else if (reg_site1 == "Britian & Ireland" && reg_site2 == "Britian & Ireland"){
        df_list[[o]][i,4] <- "BothRural_BI"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Scandinavia"){
        df_list[[o]][i,4] <- "BothRural_SC"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Western Europe" |
               reg_site1 == "Western Europe" && reg_site2 == "Scandinavia") {
        df_list[[o]][i,4] <- "BothRural_SC_WE"
      }
      else if (reg_site1 == "Britian & Ireland" && reg_site2 == "Western Europe" |
               reg_site1 == "Western Europe" && reg_site2 == "Britian & Ireland") {
        df_list[[o]][i,4] <- "BothRural_BI_WE"
      }
      else if (reg_site1 == "Scandinavia" && reg_site2 == "Britian & Ireland" |
               reg_site1 == "Britian & Ireland" && reg_site2 == "Scandinavia") {
        df_list[[o]][i,4] <- "BothRural_SC_BI"
      }
    }
    if (urb_site1 == T && urb_site2== F){
      if (reg_site1== "Britian & Ireland"){
        if (reg_site2 == "Britian & Ireland"){
          df_list[[o]][i,4] <- "UrbanRural_BI"
        }
        else if (reg_site2 == "Scandinavia"){
          df_list[[o]][i,4] <- "UrbanBI_RuralSC"
        }
        else if (reg_site2 == "Western Europe"){
          df_list[[o]][i,4] <- "UrbanBI_RuralWE"
        }
      }
      if (reg_site1== "Western Europe"){
          if (reg_site2 == "Western Europe"){
            df_list[[o]][i,4] <- "UrbanRural_WE"
          }
          else if (reg_site2 == "Scandinavia"){
            df_list[[o]][i,4] <- "UrbanWE_RuralSC"
          }
          else if (reg_site2 == "Britian & Ireland"){
            df_list[[o]][i,4] <- "UrbanWE_RuralBI"
          }
      }
       if (reg_site1== "Scandinavia"){
            if (reg_site2 == "Scandinavia"){
              df_list[[o]][i,4] <- "UrbanRural_SC"
            }
            else if (reg_site2 == "Western Europe"){
              df_list[[o]][i,4] <- "UrbanSC_RuralWE"
            }
            else if (reg_site2 == "Britian & Ireland"){
              df_list[[o]][i,4] <- "UrbanSC_RuralBI"
            }
       }
    }
    if (urb_site1 == F && urb_site2== T){
      if (reg_site2== "Britian & Ireland"){
        if (reg_site1 == "Britian & Ireland"){
          df_list[[o]][i,4] <- "UrbanRural_BI"
        }
        else if (reg_site1 == "Scandinavia"){
          df_list[[o]][i,4] <- "UrbanBI_RuralSC"
        }
        else if (reg_site1 == "Western Europe"){
          df_list[[o]][i,4] <- "UrbanBI_RuralWE"
        }
      }
      if (reg_site2== "Western Europe"){
        if (reg_site1 == "Western Europe"){
          df_list[[o]][i,4] <- "UrbanRural_WE"
        }
        else if (reg_site1 == "Scandinavia"){
          df_list[[o]][i,4] <- "UrbanWE_RuralSC"
        }
        else if (reg_site1 == "Britian & Ireland"){
          df_list[[o]][i,4] <- "UrbanWE_RuralBI"
        }
      }
      if (reg_site2== "Scandinavia"){
        if (reg_site1 == "Scandinavia"){
          df_list[[o]][i,4] <- "UrbanRural_SC"
        }
        else if (reg_site1 == "Western Europe"){
          df_list[[o]][i,4] <- "UrbanSC_RuralWE"
        }
        else if (reg_site1 == "Britian & Ireland"){
          df_list[[o]][i,4] <- "UrbanSC_RuralBI"
        }
      }
    }  #close all conditions
  if (df_list[[o]][i,4] %in% c("BothUrban_SC", "BothUrban_WE", "BothUrban_BI")){
      df_list[[o]][i,5] <- "BothUrban_SameReg"
  }
  else if (df_list[[o]][i,4] %in% c("BothUrban_BI_WE", "BothUrban_SC_BI", "BothUrban_SC_WE")){
    df_list[[o]][i,5] <- "BothUrban_DiffReg"
  }
  else if (df_list[[o]][i,4] %in% c("BothRural_BI_WE", "BothRural_SC_BI", "BothRural_SC_WE")){
    df_list[[o]][i,5] <- "BothRural_DiffReg"
  }
  else if (df_list[[o]][i,4] %in% c("BothRural_BI", "BothRural_SC", "BothRural_WE")){
    df_list[[o]][i,5] <- "BothRural_SameReg"
  }
  else if (df_list[[o]][i,4] %in% c("UrbanRural_BI", "UrbanRural_SC", "UrbanRural_WE")){
    df_list[[o]][i,5] <- "UrbanRural_SameReg"
  }
  else if (df_list[[o]][i,4] %in% c("UrbanBI_RuralSC", "UrbanBI_RuralWE","UrbanSC_RuralBI", 
                               "UrbanSC_RuralWE", "UrbanWE_RuralBI", "UrbanWE_RuralSC")){
    df_list[[o]][i,5] <- "UrbanRural_DiffReg"
  }
  }
  else {
    df_list[[o]][i,4] <- "Other"  #includes NAs
    df_list[[o]][i,5] <- "Other"
  }
}
print(o)
#print(df_list[[o]] %>% count(UrbanRural_noreg))
} #close loop across time bins

print(df_list[[1]] %>% count(UrbanRural))  
print(df_list[[1]] %>% count(UrbanRural_noreg))

save(df_list, file = file.path("C:/Users/abbap/OneDrive - University of Helsinki/NSS Fish Collaboration/Similarity 25/SiteBin_df_list.Rdata"))
load("C:/Users/abbap/OneDrive - University of Helsinki/NSS Fish Collaboration/Similarity 25/SiteBin_df_list.Rdata")

# Null model of occurrence matrix-- calculate pairwise sim and input into matrix 

# Function to shuffle the entries in each row (site) of input matrix
shuffle_matrix_row <- function(mat) {
  t(apply(mat, 1, sample))
}

# Generate 500 shuffled versions of matrix, for each time bin
shuffled_matrices <- list()

for (r in 1:8){
  shuffled_matrices[[r]] <- list()
  mat <- simats[[r]][2:ncol(simats[[r]])]
  rownames(mat) <- simats[[r]]$X  # to fill with shuffled matrices
  for (j in 1:500) {
    shuffled_matrices[[r]][[j]] <- shuffle_matrix_row(mat)  #
  }
}

# function for similarity of shuffled matrices
transform_and_calculate_corrnull <- function(mat) {
  mat[is.na(mat)] <- 0  #0s for numeric
  Matrix <- as.matrix((mat))  #for sites, not transposed
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix)) 
  sim <- sim %*% t(sim)
  return(sim) #resulting similarity matrix saved
}

#Find correlation values for each shuffled matrix
sim_shuffled_matrices <- list()

for (t in 1:8){
sim_shuffled_matrices[[t]] <- lapply(shuffled_matrices[[t]], transform_and_calculate_corrnull)
print(t)
}

# List of sim scores in null matrices for each pair of sites

# Create lists to store the corr values for each time bin
values_lists <- list()

# Save into list the values from the shuffled matrices
for (w in 1:8){
  values_lists[[w]] <- list()
for (i in 1:length(sim_shuffled_matrices[[w]])) {
  matuse <- sim_shuffled_matrices[[w]][[i]]
  rownames(matuse) <- simats[[w]]$X
  colnames(matuse) <- simats[[w]]$X 
 for (row in 1:nrow(matuse)) {
    for (col in 1:ncol(matuse)) {
      key <- paste(rownames(matuse)[row], colnames(matuse)[col], sep = "_AND_") #save which two sites go with value
      values_lists[[w]][[key]] <- c(values_lists[[w]][[key]], matuse[row, col])
    }
  }
  if (i %in% c(100,400)){
    print(i)
  }
}
  print(w)
}
  
length(values_lists[[4]]) #This is number of site pairs

save(values_lists, file = file.path("ShuffledSitesValues_Lists.Rdata")) #input beginning of path for your wd
load("ShuffledSitesValues_Lists.Rdata")


# Calculate the mean and standard deviation for each (Taxon1, Taxon2) pair from null and compare
for (l in 1:8){
  df_list[[l]]$Mean_Null <- NA
  df_list[[l]]$SD_Null <-  NA
  df_list[[l]]$Significance_InputcNull <- NA
  for (key in names(values_lists[[l]])){
    values <- values_lists[[l]][[key]]
    site1 <- strsplit(key, "_AND_")[[1]][1]
    site2 <- strsplit(key, "_AND_")[[1]][2]
    mean_value <- mean(values)
    sd_value <- sd(values)
    
    # Save mean and SD of null matrices into original dataframe
    row_index <- which(df_list[[l]]$Site1 == site1 & df_list[[l]]$Site2 == site2)
    df_list[[l]][row_index, 6] <- mean_value
    df_list[[l]][row_index, 7] <- sd_value
    
    # Calculate the z-score of observed vs. null (in terms of SD's from mean)
    observed_value <- df_list[[l]][row_index,3]
    if (!is.na(observed_value) && !is.na(mean_value) && !is.na(sd_value)) {
      z_score <- (observed_value - mean_value) / sd_value
      df_list[[l]][row_index, 8] <- z_score
    }
    }
  # remove pairs where two sites are the same
 # nrow(df_list[[l]])
  df_list[[l]] <-  df_list[[l]] %>% filter(!(Site1==Site2))
  print(nrow(df_list[[l]]))
}

save(df_list, file = file.path("SiteBin_Family_df_list.Rdata"))
load("SiteBin_Family_df_list.Rdata")


#Plot similarities by categories vs. null
colnames(df_list[[1]])
ggplot(df_list[[8]])+ geom_violin(aes(x= UrbanRural_noreg, y=Bin_8_obs), fill="lightblue")+
  geom_violin(aes(x= UrbanRural_noreg, y=Mean_Null), fill="pink3") + ylab("Similarity Score")
 
ggplot(df_list[[4]])+ geom_violin(aes(x= UrbanRural, y=Bin_4_obs), fill="lightblue")+
  geom_violin(aes(x= UrbanRural, y=Mean_Null), fill="pink3") + 
  ylab("Similarity Score")

# Plot significance ratio over time

# Consider results-- how many pairs of sites have significant z score
head(df_list[[3]])
Bin_1_over_1.96 <- df_list[[1]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96) %>% mutate(TimeBin= 1)
Bin_2_over_1.96 <- df_list[[2]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 2)
Bin_3_over_1.96 <- df_list[[3]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 3)
Bin_4_over_1.96 <- df_list[[4]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 4)
Bin_5_over_1.96 <- df_list[[5]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 5)
Bin_6_over_1.96 <- df_list[[6]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 6)
Bin_7_over_1.96 <- df_list[[7]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 7)
Bin_8_over_1.96 <- df_list[[8]] %>% group_by(UrbanRural) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 8)

df_sig <- rbind(Bin_1_over_1.96,Bin_2_over_1.96,Bin_3_over_1.96,Bin_4_over_1.96,
                Bin_5_over_1.96,Bin_6_over_1.96,Bin_7_over_1.96,Bin_8_over_1.96)
colnames(df_sig)[2] <- "Over1.96SDs"

# Find ratios between sig. and not, and % of significant cases in each category
Significance_Grouped <- df_sig %>%
  group_by(TimeBin, UrbanRural) %>%
  summarise(
    Ratio = (sum(n[Over1.96SDs== TRUE]) /
      sum(n[Over1.96SDs == FALSE])),
    Percentage = (sum(n[Over1.96SDs == TRUE]) /
      (sum(n[Over1.96SDs == FALSE] + sum(n[Over1.96SDs == TRUE]))))
  )
  
#save results
write.csv(df_sig, "Results_PairwiseCorr_Sites_500null_8timebins.csv")
write.csv(Significance_Grouped, "SummaryResults_Significance_SitePairs_Family_500null_8timebins.csv")
write.csv(Significance_Grouped_nr, "SummaryResults_Significance_NoReg_SitePairs_Family_500null_8timebins.csv")


# Plotting ####
# Plot across bins, broken down by pair urban/rural comparison
ggplot(Significance_Grouped, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  #scale_color_manual(values = c("Marine" = "blue2", "Fresh" = "forestgreen", 
   #                             "FreshMarine" = "magenta2", "Other"= "grey2"))+
  labs(title = "Percentage Significant",
       x = "Bin Number",
       y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

# Without region-specific comparisons
Bin_1_over_1.96 <- df_list[[1]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96) %>% mutate(TimeBin= 1)
Bin_2_over_1.96 <- df_list[[2]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 2)
Bin_3_over_1.96 <- df_list[[3]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 3)
Bin_4_over_1.96 <- df_list[[4]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 4)
Bin_5_over_1.96 <- df_list[[5]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 5)
Bin_6_over_1.96 <- df_list[[6]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 6)
Bin_7_over_1.96 <- df_list[[7]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 7)
Bin_8_over_1.96 <- df_list[[8]] %>% group_by(UrbanRural_noreg) %>% count(Significance_InputcNull >= 1.96)%>% mutate(TimeBin= 8)

df_sig_nr <- rbind(Bin_1_over_1.96,Bin_2_over_1.96,Bin_3_over_1.96,Bin_4_over_1.96,
                Bin_5_over_1.96,Bin_6_over_1.96,Bin_7_over_1.96,Bin_8_over_1.96)
colnames(df_sig_nr)[2] <- "Over1.96SDs"

# Find ratios between sig. and not, and % of significant cases in each category
Significance_Grouped_nr <- df_sig_nr %>%
  group_by(TimeBin, UrbanRural_noreg) %>%
  summarise(
    Ratio = (sum(n[Over1.96SDs== TRUE]) /
               sum(n[Over1.96SDs == FALSE])),
    Percentage = (sum(n[Over1.96SDs == TRUE]) /
                    (sum(n[Over1.96SDs == FALSE] + sum(n[Over1.96SDs == TRUE]))))
  )


ggplot(Significance_Grouped_nr, aes(x = TimeBin, y = Percentage, color = UrbanRural_noreg, 
                                    group = UrbanRural_noreg)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
 # scale_color_manual(values = c("Marine" = "blue2", "Fresh" = "forestgreen", 
  #                              "FreshMarine" = "magenta2", "Other"= "grey2"))+
  labs(
    x = "Bin Number",
    y = "Percentage Significant",
    color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

# remove comparisons of specifc regions
Significance_Grouped_nr$UrbRur <- sub("_.*", "", Significance_Grouped_nr$UrbanRural_noreg) 
Significance_Grouped_nr_group <- Significance_Grouped_nr %>% 
  group_by(UrbRur, TimeBin) %>%
  summarise(Average_Percentage = mean(Percentage))

ggplot(Significance_Grouped_nr_group, aes(x = TimeBin, y = Average_Percentage, color = UrbRur, 
                                    group = UrbRur)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  # scale_color_manual(values = c("Marine" = "blue2", "Fresh" = "forestgreen", 
  #                              "FreshMarine" = "magenta2", "Other"= "grey2"))+
  labs(
    x = "Bin Number",
    y = "Percentage Significant",
    color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

# Comparisons between groups, plots
Significance_Grouped_bothurb <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "BothUrban")) 

ggplot(Significance_Grouped_bothurb, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothUrban_BI_WE" = "#8000FF",
                                "BothUrban_SC_BI" = "forestgreen",
                                "BothUrban_SC_WE" = "#FF8000",
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanBI_RuralSC" = "#99FF66", "UrbanRural_BI"= "grey2"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

Significance_Grouped_bothrur <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "BothRural")) 

ggplot(Significance_Grouped_bothrur, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothRural_BI" = "dodgerblue",
                                "BothRural_WE" = "pink3",
                                "BothRural_SC"= "goldenrod2",
                                "BothRural_BI_WE" = "#8000FF",
                                "BothRural_SC_BI" = "forestgreen",
                                "BothRural_WE_SC" = "#FF8000",
                                "BothRural_SC_WE" = "#FF8000", 
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanSC_RuralBI" = "#99FF66",
                                "UrbanSC_RuralWE" = "#FFCC66"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

Significance_Grouped_urbrur <-  Significance_Grouped %>% 
  filter(!(UrbanRural %in% Significance_Grouped_bothurb$UrbanRural)) %>% 
  filter(!(UrbanRural %in% Significance_Grouped_bothrur$UrbanRural))

ggplot(Significance_Grouped_urbrur, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("UrbanRural_BI"= "blue", "UrbanRural_WE"= "red",
                                "UrbanRural_SC"= "yellow",
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "purple3",
                                "UrbanSC_RuralBI" = "#99FF66",
                                "UrbanBI_RuralSC" = "limegreen",
                                "UrbanSC_RuralWE" = "#FFCC66",
                                "UrbanWE_RuralSC" = "goldenrod3"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()



# Make region-specific comparisons legible
Significance_Grouped_BI_urb <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "BI")) %>% filter(str_detect(UrbanRural, "Urban")) %>% 
  filter(!(str_detect(UrbanRural, "RuralBI")))

ggplot(Significance_Grouped_BI_urb, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrbanWE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothUrban_BI_WE" = "#8000FF",
                                "BothUrban_SC_BI" = "forestgreen",
                                "BothUrban_WE_SC" = "#FF8000",
                                "UrbanBI_RuralWE" = "#CC66FF",
                               "UrbanBI_RuralSC" = "#99FF66", "UrbanRural_BI"= "grey2"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

Significance_Grouped_BI_rur <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "BI")) %>% filter(str_detect(UrbanRural, "Rural")) %>% 
  filter(!(str_detect(UrbanRural, "UrbanBI")))

ggplot(Significance_Grouped_BI_rur, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrbanWE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothRural_BI" = "dodgerblue",
                                "BothRural_WE" = "pink3",
                                "BothRural_SC"= "goldenrod2",
                                "BothRural_BI_WE" = "#8000FF",
                                "BothRural_SC_BI" = "forestgreen",
                                "BothRural_WE_SC" = "#FF8000",
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanSC_RuralBI" = "#99FF66",
                                "UrbanBI_RuralSC" = "#99FF66", "UrbanRural_BI"= "grey2"))+
  labs(
       x = "Bin Number",
       y = "Percentage Significant",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

Significance_Grouped_WE_urb <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "WE")) %>% filter(str_detect(UrbanRural, "Urban")) %>% 
  filter(!(str_detect(UrbanRural, "RuralWE")))

ggplot(Significance_Grouped_WE_urb, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothUrban_BI_WE" = "#8000FF",
                                "BothUrban_SC_BI" = "forestgreen",
                                "BothUrban_SC_WE" = "#FF8000",
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanWE_RuralSC" = "#FFCC66",
                                "UrbanBI_RuralSC" = "#99FF66", "UrbanRural_WE"= "grey2"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()


Significance_Grouped_WE_rur <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "WE")) %>% filter(str_detect(UrbanRural, "Rural")) %>% 
  filter(!(str_detect(UrbanRural, "UrbanWE")))

ggplot(Significance_Grouped_WE_rur, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothRural_BI" = "dodgerblue",
                                "BothRural_WE" = "pink3",
                                "BothRural_SC"= "goldenrod2",
                                "BothRural_BI_WE" = "#8000FF",
                                "BothRural_SC_BI" = "forestgreen",
                                "BothRural_WE_SC" = "#FF8000",
                                "BothRural_SC_WE" = "#FF8000", 
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanSC_RuralBI" = "#99FF66",
                                "UrbanSC_RuralWE" = "#FFCC66", "UrbanRural_WE"= "grey2"))+
  labs(
    x = "Bin Number",
    y = "Percentage Significant",
    color = "Pairwise Urban/Rural Comparison") +
  theme_bw()

Significance_Grouped_SC_urb <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "SC")) %>% filter(str_detect(UrbanRural, "Urban")) %>% 
  filter(!(str_detect(UrbanRural, "RuralSC")))

ggplot(Significance_Grouped_SC_urb, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothUrban_BI_WE" = "#8000FF",
                                "BothUrban_SC_BI" = "forestgreen",
                                "BothUrban_SC_WE" = "#FF8000",
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanSC_RuralWE" = "#FFCC66",
                                "UrbanSC_RuralBI" = "#99FF66", "UrbanRural_SC"= "grey2"))+
  labs(y = "Percentage Significant",
       x = "Bin Number",
       #y = "Ratio",
       color = "Pairwise Urban/Rural Comparison") +
  theme_bw()


Significance_Grouped_SC_rur <-  Significance_Grouped %>% 
  filter(str_detect(UrbanRural, "SC")) %>% filter(str_detect(UrbanRural, "Rural")) %>% 
  filter(!(str_detect(UrbanRural, "UrbanSC")))

ggplot(Significance_Grouped_SC_rur, aes(x = TimeBin, y = Percentage, color = UrbanRural, group = UrbanRural)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BothUrban_BI" = "blue2", "BothUrban_WE" = "red",
                                "BothUrban_SC" = "#FFD700",
                                "BothRural_BI" = "dodgerblue",
                                "BothRural_WE" = "pink3",
                                "BothRural_SC"= "goldenrod2",
                                "BothRural_BI_WE" = "#8000FF",
                                "BothRural_SC_BI" = "forestgreen",
                                "BothRural_WE_SC" = "#FF8000",
                                "BothRural_SC_WE" = "#FF8000", 
                                "UrbanBI_RuralWE" = "#CC66FF",
                                "UrbanWE_RuralBI" = "#CC66FF",
                                "UrbanBI_RuralSC" = "#99FF66",
                                "UrbanWE_RuralSC" = "#FFCC66", "UrbanRural_SC"= "grey2"))+
  labs(
    x = "Bin Number",
    y = "Percentage Significant",
    color = "Pairwise Urban/Rural Comparison") +
  theme_bw()



