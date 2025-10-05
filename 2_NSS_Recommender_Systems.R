# Recommender Systems run using matrix factorization
# see Žliobaitė (2022), doi: 10.1111/2041-210X.13916, adapted from https://github.com/zliobaite/fossilrec
library(tidyverse)
library(pROC)

# Read in data: site (x) vs taxon (y) matrix. Can have additional columns with site information (here, age, lat/long, etc.)
data_all <- read.csv("Species_PA_noCEE_19Sep25.csv", header=T, row.names = 1)
# Read in matrices with fraction of NISP for each taxa by site, only for sites where sieving was employed (for evaluation)
occ_nisp<- read.csv("Sieved_Species_NISPfraction_noCEE_19Sep25.csv", header=T, row.names = 1)  
occ_sieve <- as.matrix(occ_nisp)
t <- length(which(!(is.na(occ_sieve))))
zs <- length(which(occ_sieve==0))  #how many 0s?
zs/t # fraction of sieve site possible occurrences that are negative
ps <- length(which(occ_sieve>0.2)) #how many species make up more than 20% of occurrences at a site?
ps/t  

# If using specific time bins, define here
sitefilt <- read.csv("Site_Metadata_SpeciesLevel_noCEE_190925.csv")
bins <- unique(sitefilt$time_bins2)
bins
thisbin <- sitefilt[sitefilt$time_bins2 %in% bins[6],]    #8 bins
data_all_bin <- data_all %>% filter(rownames(data_all) %in% thisbin$DB_Assemblage_ID)
occ_nisp_bin <- occ_nisp %>% filter(rownames(occ_nisp) %in% thisbin$DB_Assemblage_ID)

# setwd to folder to save results
write.csv(data_all_bin, "Species_PA_1500_1700.csv")
write.csv(occ_nisp_bin, "Sieved_Species_NISPfraction_1500_1700.csv")

# Set up recommender system run
# if using filtered bins, otherwise use data_all and occ_nisp as it is
data_all <- data_all_bin
occ_nisp <-  occ_nisp_bin

n_users <- dim(data_all)[1]  # number of rows or sites
n_columns <- dim(data_all)[2]  #number of columns 
n_items <- n_columns  # number of taxon columns (if added locality data columns at end, subtract that number)

do_auc_cv <- TRUE #includes cross validation
rec_mat <- TRUE   # saves recommender score values

#parameters-- these values are my initial suggestions based on some trials
alpha <- 5  #higher = more confidence in presences vs. absences
n_features <- 10   #how many topics in inner matrix
lambda <- 18  #regularization, higher links co-occurring taxa more strongly in predictions
cik <- 10  #iternations of matrix factorization (leave at 10 here)

n_remove_test <- 10   #number of occurrences to remove for cross-validation tests
cik_cv <- 10  #how many iterations of cross-validation (lower for speed)

set.seed(1981)
#to save outputs
setwd("")

results_all <- c()
 # In this case, I loop across parameters for alpha, n_features, and lambda 
for (alpha in c(35)){  #can choose just one value to start
  print(paste0("alpha =", alpha))
  for (n_features in c(3,4,5,6,8,10, 12)){
    for (lambda  in c(15, 18, 20, 25, 30, 40)){   #all of these parameters can be varied further
      print(lambda)
      data_now <- as.matrix(data_all)  #data_now matrix -- check this is only taxon columns
      rownames(data_now) <- paste(rownames(data_all) )
      input_marginalT <- colSums(data_now) #sum of each species' occurrences
      input_marginalS <- rowSums(data_now)  #occurrence count for sites
    
      
      factorize <- function(data_occ,n_features,alpha,lambda,cik,rec_mat){  #main RS function
        n_sites <- dim(data_occ)[1]
        n_genera <- dim(data_occ)[2]
        C <- 1 + alpha*data_occ  # alpha value multiplies occ's, gives more weight to presences
        P <- as.matrix((data_occ > 0) + 0)
        X <- matrix(rnorm(n_sites*n_features), ncol = n_features)
        Y <- matrix(rnorm(n_genera*n_features), ncol = n_features)
        Ix <- diag(n_sites)
        Iy <- diag(n_genera)
        Ilambda <- lambda*diag(n_features)
        for (sk in 1:cik){ #iterates
          #Y-transpose-Y and X-transpose-X
          xTx <- t(X)%*%X
          yTy <- t(Y)%*%Y
          # Loop through all sites
          for (uu in 1:n_sites){
            Cu <- diag(C[uu,])
            Cu_I <- Cu - Iy
            for_inversion <- yTy + t(Y)%*%Cu_I%*%Y + Ilambda  #introduces lambda regularization
            X[uu,] <- solve(for_inversion)%*%t(Y)%*%Cu%*%P[uu,] #computes inverse of yTy matrix, multiplies
          }
          # Loop through all genera
          for (ii in 1:n_genera){
            Ci <- diag(C[,ii])
            Ci_I <- Ci - Ix
            for_inversion <- xTx + t(X)%*%Ci_I%*%X + Ilambda
            Y[ii,] <- solve(for_inversion)%*%t(X)%*%Ci%*%P[,ii]
          }
        }
        rec <- X%*%t(Y)
        #print(mean(rec))
        #print(sd(rec))
        #print(max(rec))
        #print(min(rec))
        if (rec_mat){  #saves inner matrices
          rownames(X) <- rownames(data_all)
          rownames(Y) <- colnames(data_all) 
          write.csv(X, paste0("Xmat",alpha,"_", n_features, "_",lambda, "looprun.csv"))
          write.csv(Y, paste0("Ymat",alpha, "_",n_features, "_",lambda, "looprun.csv"))
        }
        return(rec)
      }
      
      #full fit on train data
      REC <- round(factorize(data_now,n_features,alpha,lambda,cik,TRUE),digits = 2)
      colnames(REC) <- colnames(data_all)
      rownames(REC) <- rownames(data_all) 
      REC_all <- REC
      
      write.csv(REC_all, paste0("REC_",alpha, "_", n_features, "_",lambda, "looprun.csv"))
      
      
      # ROC on training data, defined true positives and negatives
      ind_sieve <- which(!is.na(occ_nisp))  #sieved sites with NISP fraction
      occ_sieve <- as.matrix(occ_nisp)
                    
      #cross-validation
      res_cv_raw <- c()
      res_cv_auc <- c()
      
      for (sk in 1:cik_cv){
        ind_occ <- which(occ_sieve>0.2) #where sieve bones from taxa are over 20% of assemblage, "true positive"
        ind_neg <- which(occ_sieve==0)  # where sieve sites have no bones from taxon, "true negative"
        ind_true_remove <- sample(ind_occ, n_remove_test) # select 10 true positives
        ind_false_remove <- sample(ind_neg, n_remove_test)  # select 10 true negatives
        
        data_now_test <- data_now
        data_now_test[ind_true_remove] <- 0  #for given true pos, make observation absence
        REC_cv <- round(factorize(data_now_test,n_features,alpha,lambda,cik,FALSE),digits = 2) #run RS
        res_cv_raw <- rbind(res_cv_raw,c(mean(REC_cv[ind_true_remove]),mean(REC_cv[ind_false_remove])))
        
        if (do_auc_cv){
          data_test <- REC_cv[c(ind_false_remove,ind_true_remove)] #for 10 cv 0s and 10 true 0s
          labels_test <- c(rep(0,length(ind_false_remove)),rep(1,length(ind_true_remove)))
          data_frame_test <- as.data.frame(cbind(data_test,labels_test))
          suppressMessages({
          roc_test <- roc(labels_test,data_test)  # calculate AUC
          auc_test <-  auc(roc_test)
          })
        }else{
          auc_test <- 0
        }
        res_cv_auc <- c(res_cv_auc,auc_test)  # add auc score for this iteration to rest
      }
      cv_mean_predictions <- round(apply(res_cv_raw,2,mean),digits = 3)
      cv_sd_predictions <- round(apply(res_cv_raw,2,sd),digits = 3)
      cv_auc <- round(mean(res_cv_auc),digits = 3)
      cv_auc_sd <- round(sd(res_cv_auc),digits = 3)
      
      #add marginals
      REC_marginalT <- colSums(REC)
      marginalT_increase <- ((REC_marginalT-input_marginalT)/input_marginalT)*100  #vs input
      taxa_marginal <- mean(marginalT_increase, na.rm=T)  # percent increase in occurrences of each taxon at all sites
      
      REC_marginalS <- rowSums(REC)
      marginalS_increase <- ((REC_marginalS-input_marginalS)/input_marginalS)*100  # vs input
      site_marginal <- mean(marginalS_increase, na.rm=T) # percent increase of occurrence counts at each site
      
      mean_data_now <- round(mean(data_now),digits = 3)  #mean input value
      mean_predictions_now <- round(mean(REC),digits = 3)  #mean REC score for whole output
      mean_predictions_pos <- round(mean(REC[ind_occ]),digits = 3) # mean REC for trues positives
      mean_predictions_neg <- round(mean(REC[ind_neg]),digits = 3)  # mean REC score for true negatives
      animals_predict <- apply(round(REC),2,sum)  #sum REC values per taxon
      animals_true <- apply(round(data_now),2,sum)  #input occ sum per taxon
      mean_error_animals <- round(mean(abs(animals_predict-animals_true)),digits = 3)
      localities_predict <- apply(round(REC),1,sum)   #sum REC values per site
      localities_true <- apply(round(data_now),1,sum)  # input occ sum per site
      mean_error_localities <- round(mean(abs(localities_predict-localities_true)),digits = 3)
      correlation_everything <- round(cor(as.vector(REC),as.vector(data_now)),digits = 3) # REC score correlated to observed
      err <- REC-data_now
      mean_absolute_error <- round(mean(abs(err)),digits = 3) #absolute error of REC vs. input
      
      
      #Summary of results-- compare parameter combinations for performance in cross-validation
      results_all <- rbind(results_all,cbind(alpha,n_features,lambda,cik,n_remove_test,cik_cv,
                                             cv_mean_predictions[1],cv_sd_predictions[1],cv_mean_predictions[2],
                                             cv_sd_predictions[2],cv_auc,cv_auc_sd, 
                                             mean_predictions,
                                             mean_data, mean_predictions_pos, mean_predictions_neg,
                                             mean_error_animals,mean_error_localities,
                                             correlation_everything,mean_absolute_error, taxa_marginal,
                                             site_marginal))
      colnames(results_all)[7:10] <- c("cv_occ_REC_mean", "cv_neg_REC_mean",
                                       "cv_occ_REC_sd", "cv_neg_REC_sd")
      write.csv(results_all, "output_results_loop.csv")      
      
      
    }
  }
}
