library(randomForest)
library(caret)


ic50_data <- read.csv("/media/rosner-lab/ubuntu_hdd1/Ali_MLinC_HWs/Ali_ML_HW4/MLiC_HW4_data/GDSC/GDSC_dose_response_fitted_cosmic_mimick.csv")
descriptors <- read.delim("/media/rosner-lab/ubuntu_hdd1/Ali_MLinC_HWs/Ali_ML_HW4/MLiC_HW4_data/GDSC/GDSC-drug-descriptors-and-fingerprints/descriptors.txt")
ECFPs <- read.delim("/media/rosner-lab/ubuntu_hdd1/Ali_MLinC_HWs/Ali_ML_HW4/MLiC_HW4_data/GDSC/GDSC-drug-descriptors-and-fingerprints/ECFP_1024_m0-2_b2_c_0.txt", header=FALSE)

# USE ONLY THE DESCRIPTORS FIRST

# Prepare the data frames for the merge() function
descriptors <- descriptors[,2:ncol(descriptors)]
colnames(descriptors)[1] <- "NSC"

ic50_data <- ic50_data[, c(1:2,4)]

Dat1 <- merge(ic50_data, descriptors, by = "NSC")

# As an example, work on MDA-MB-231 cell line (COSMIC.905960)
MB231_temp1 <- Dat1[Dat1$CELLNAME == "COSMIC.905960",]

# Convert descriptors to numeric type
data.char <- apply(MB231_temp1[, 4:ncol(MB231_temp1)], 2, as.character)
data.num <- apply(data.char, 2, as.numeric)
MB231_temp1 <- cbind(MB231_temp1[,1:3], data.num)

# NA cleanup
table(is.na(MB231_temp1$GROWTH)) # there are no rows with a missing label (or outcome data)
MB231_des <- MB231_temp1[, -which(colMeans(is.na(MB231_temp1)) > 0.5)] # remove columns with more than 50% NA values
dim(MB231_des)

set.seed(123)
MB231_des_imputed <- rfImpute(x = MB231_des[,4:ncol(MB231_des)], y = MB231_des[,3], iter=5, ntree=300)
MB231_des_imputed <- cbind(MB231_des[,1:2], MB231_des_imputed)
colnames(MB231_des_imputed)[3] <- "LN_IC50"
save(MB231_des_imputed, file = "MB231_des_imputed.rda")

# We also need to remove the constant variables before PCA analysis as they have zero variance
MB231_des_imputed2 <- MB231_des_imputed[, sapply(MB231_des_imputed, function(x) length(unique(x))) > 1] # Remove the constant columns


# Perform PCA to reduce dimensionality using prcomp() function
MB231_des_pca <- prcomp(MB231_des_imputed2[,3:ncol(MB231_des_imputed2)], scale. = T, center = T)
MB231_des.train <- cbind(MB231_des_imputed2[,1:2], MB231_des_pca$x[,1:50]) # PC50 is when we reach 90% cumulative variance

# Perform Random Forest using ranger() from "caret" with 5-fold cross-validation
tc <- trainControl(method = "cv", number = 5)

set.seed(123)
rf_MB231_des <- train(x = MB231_des.train[, 3:ncol(MB231_des.train)], 
                      y = MB231_des.train[, 2], 
                      trControl = tc, 
                      method = "ranger", 
                      importance = "impurity")  # R2: -0.0256 (womp womp!)
rf_MB231_des$finalModel 
                     
set.seed(123)
glm_MB231_des <- train(x = MB231_des.train[, 3:ncol(MB231_des.train)], 
                      y = MB231_des.train[, 2], 
                      trControl = tc, 
                      method = "glm") # R2: 0.0105 (womp womp!)
glm_MB231_des$results


# USE ONLY THE FINGERPRINT DATA

# Prepare the data frames for the merge() function
colnames(ECFPs)[1] <- "NSC"
Dat2 <- merge(ic50_data, ECFPs, by = "NSC")
# As an example, work on MDA-MB-231 cell line (COSMIC.905960)
MB231_temp2 <- Dat2[Dat2$CELLNAME == "COSMIC.905960",]

table(apply(MB231_temp2, 2, function(x) any(is.na(x)))) # There are no NA values in this data frame.
MB231_fp <- MB231_temp2[, sapply(MB231_temp2, function(x) length(unique(x))) > 1] # Remove the constant columns (there were 19 of them)
MB231_fp_pca <- prcomp(MB231_fp[,3:ncol(MB231_fp)], scale. = T, center = T) # reaches 90% by PC130
MB231_fp.train <- cbind(MB231_fp[,1:2], MB231_fp_pca$x[,1:130])
# Perform Random Forest using ranger() from "caret" with 5-fold cross-validation (also glm)
tc <- trainControl(method = "cv", number = 5)

set.seed(123)
rf_MB231_fp <- train(x = MB231_fp.train[, 3:ncol(MB231_fp.train)], 
                      y = MB231_fp.train[, 2], 
                      trControl = tc, 
                      method = "ranger", 
                      importance = "impurity")  # R2: -0.0156 (womp womp!)
rf_MB231_fp$finalModel

# USE BOTH DESCRIPTOR AND FINGERPRINT DATA

# Merge the imputed desriptor data with the fingerprint data for MB231 cells
MB231_fp1 <- MB231_fp[, c(1,3:ncol(MB231_fp))] 
MB231_comb <- merge(MB231_des_imputed2, MB231_fp1, by = "NSC")

MB231_comb_pca <- prcomp(MB231_comb[,3:ncol(MB231_comb)], scale. = T, center = T) # reaches 90% by PC100
MB231_comb.train <- cbind(MB231_comb[,1:2], MB231_comb_pca$x[,1:100])
# Perform Random Forest using ranger() from "caret" with 5-fold cross-validation
tc <- trainControl(method = "cv", number = 5)

set.seed(123)
rf_MB231_comb <- train(x = MB231_comb.train[, 3:ncol(MB231_comb.train)], 
                     y = MB231_comb.train[, 2], 
                     trControl = tc, 
                     method = "ranger", 
                     importance = "impurity")  # R2: -0.0356 (womp womp!)
rf_MB231_comb$finalModel



# INSTEAD OF IC50 VALUES, TRY ADDING DIFFERENT DOSES AND GROWTH VALUES FOR TRAINING THE MODEL
GDSC_dose_response <- read.csv("/media/rosner-lab/ubuntu_hdd1/Ali_MLinC_HWs/Ali_ML_HW4/MLiC_HW4_data/GDSC/GDSC_dose_response.csv")
# Extract the dose-response curve for MB231 cells 
DR_MB231 <- subset(GDSC_dose_response, GDSC_dose_response$CELLNAME == "MDA-MB-231", select = c("NSC", "CELLNAME", "GROWTH", "LOG_CONCENTRATION")) #217 unique drug names
# Merge the dose response data with discriptor and fingerprint data
Dat3 <- merge(DR_MB231, descriptors, by = "NSC")
length(unique(Dat3$NSC)) # 215 drugs out of 217 had descriptor data
table(unique(Dat3$NSC) %in% MB231_des_imputed2$NSC) # These 215 drugs included all of 214 drugs that we imputed the NA values for earlier.
# For the sake of saving time, I will just use the imputed data frame
Dat3 <- merge(DR_MB231, MB231_des_imputed2[,c(1,3:ncol(MB231_des_imputed2))], by = "NSC")
Dat3 <- merge(Dat3, ECFPs, by = "NSC")
Dat3 <- Dat3[, sapply(Dat3, function(x) length(unique(x))) > 1] # Remove the constant columns
# PCA 
DR_MB231_pca <- prcomp(Dat3[,3:ncol(Dat3)], scale. = T, center = T)
plot(DR_MB231_pca, type = "l")
## compute standard deviation of each principal component
std_dev <- DR_MB231_pca$sdev
## compute variance
pr_var <- std_dev^2
## proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]
## scree plot
plot(prop_varex[1:15], xlab = "Principal Component",
       ylab = "Proportion of Variance Explained",
       type = "b")
## cumulative scree plot
plot(cumsum(prop_varex[1:150]), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       type = "b") # 150PC captures more than 95%

DR_MB231.train <- cbind(Dat3[,1:2], DR_MB231_pca$x[,1:150]) 
# Random Forest
set.seed(123)
rf_DR_MB231 <- train(x = DR_MB231.train[, 3:ncol(DR_MB231.train)], 
                       y = DR_MB231.train[, 2], 
                       trControl = tc, 
                       method = "ranger", 
                       importance = "impurity")  # R2: 0.701 (WOOHOOO!~)
rf_DR_MB231$finalModel

# INVESTIGATE DESCRIPTORS AND FINGERPRINTS INDIVIDUALLY (ALONG WITH THE DOSAGE DATA)
Dat4 <- Dat3[, 1:2227] # Descriptors
Dat5 <- Dat3[, c(1:3, 2228:ncol(Dat3))] # Fingerprints
## Dat 4 - Descriptors
DR_MB231_pca2 <- prcomp(Dat4[,3:ncol(Dat4)], scale. = T, center = T)
std_dev <- DR_MB231_pca2$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(prop_varex[1:15], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex[1:150]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # roughly 80PC captures more than 95%
DR_MB231.train2 <- cbind(Dat4[,1:2], DR_MB231_pca2$x[,1:80]) 

set.seed(123)
rf2_DR_MB231 <- train(x = DR_MB231.train2[, 3:ncol(DR_MB231.train2)], 
                     y = DR_MB231.train2[, 2], 
                     trControl = tc, 
                     method = "ranger", 
                     importance = "impurity")  # R2: 0.702 (Didn't change at all)
rf2_DR_MB231$finalModel

## Dat 5 - Fingerprints
DR_MB231_pca3 <- prcomp(Dat5[,3:ncol(Dat5)], scale. = T, center = T)
std_dev <- DR_MB231_pca3$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(prop_varex[1:15], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex[1:175]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # roughly 160PC captures more than 95%
DR_MB231.train3 <- cbind(Dat5[,1:2], DR_MB231_pca3$x[,1:160]) 

set.seed(123)
rf3_DR_MB231 <- train(x = DR_MB231.train3[, 3:ncol(DR_MB231.train3)], 
                      y = DR_MB231.train3[, 2], 
                      trControl = tc, 
                      method = "ranger", 
                      importance = "impurity")  # R2: 0.704 (Didn't change at all)
rf3_DR_MB231$finalModel



# PERFORM THE SAME ANALYSIS ON MORE CELLL LINES
### Since adding the descriptors vs. fingerprints vs. both did not make a change on R2 values at all, I will use fingerprints only (along with dose info)
### This way, I do not have to impute missing values in the descriptors for each cell line.

# Extract the dose-response curve for Hs578T cells 
DR_HS587T <- subset(GDSC_dose_response, GDSC_dose_response$CELLNAME == "Hs-578-T", select = c("NSC", "CELLNAME", "GROWTH", "LOG_CONCENTRATION")) 
Dat6 <- merge(DR_HS587T, ECFPs, by = "NSC")
Dat6 <- Dat6[, sapply(Dat6, function(x) length(unique(x))) > 1] # Remove the constant columns
# PCA 
DR_HS587T_pca <- prcomp(Dat6[,3:ncol(Dat6)], scale. = T, center = T)
plot(DR_HS587T_pca, type = "l")
std_dev <- DR_HS587T_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(cumsum(prop_varex[1:175]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # 160PC captures more than 95%

DR_HS587T.train <- cbind(Dat6[,1:2], DR_HS587T_pca$x[,1:160]) 
# Random Forest
set.seed(123)
rf_DR_HS578T <- train(x = DR_HS587T.train[, 3:ncol(DR_HS587T.train)], 
                     y = DR_HS587T.train[, 2], 
                     trControl = tc, 
                     method = "ranger", 
                     importance = "impurity")  
rf_DR_HS578T$finalModel # R2: 0.749 

# Extract the dose-response curve for MCF7 cells 
DR_MCF7 <- subset(GDSC_dose_response, GDSC_dose_response$CELLNAME == "MCF7", select = c("NSC", "CELLNAME", "GROWTH", "LOG_CONCENTRATION")) 
Dat7 <- merge(DR_MCF7, ECFPs, by = "NSC")
Dat7 <- Dat7[, sapply(Dat7, function(x) length(unique(x))) > 1] # Remove the constant columns
# PCA 
DR_MCF7_pca <- prcomp(Dat7[,3:ncol(Dat7)], scale. = T, center = T)
plot(DR_MCF7_pca, type = "l")
std_dev <- DR_MCF7_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(cumsum(prop_varex[1:175]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # 160PC captures more than 95%

DR_MCF7.train <- cbind(Dat7[,1:2], DR_MCF7_pca$x[,1:160]) 
# Random Forest
set.seed(123)
rf_DR_MCF7<- train(x = DR_MCF7.train[, 3:ncol(DR_MCF7.train)], 
                      y = DR_MCF7.train[, 2], 
                      trControl = tc, 
                      method = "ranger", 
                      importance = "impurity")  
rf_DR_MCF7$finalModel # R2: 0.724



# Extract the dose-response curve for HeLa cells 
DR_HeLa <- subset(GDSC_dose_response, GDSC_dose_response$CELLNAME == "HeLa", select = c("NSC", "CELLNAME", "GROWTH", "LOG_CONCENTRATION")) 
Dat8 <- merge(DR_HeLa, ECFPs, by = "NSC")
Dat8 <- Dat8[, sapply(Dat8, function(x) length(unique(x))) > 1] # Remove the constant columns
# PCA 
DR_HeLa_pca <- prcomp(Dat8[,3:ncol(Dat8)], scale. = T, center = T)
plot(DR_HeLa_pca, type = "l")
std_dev <- DR_HeLa_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(cumsum(prop_varex[1:175]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # 160PC captures more than 95%

DR_HeLa.train <- cbind(Dat8[,1:2], DR_HeLa_pca$x[,1:160]) 
# Random Forest
set.seed(123)
rf_DR_HeLa<- train(x = DR_HeLa.train[, 3:ncol(DR_HeLa.train)], 
                   y = DR_HeLa.train[, 2], 
                   trControl = tc, 
                   method = "ranger", 
                   importance = "impurity")  
rf_DR_HeLa$finalModel # R2: 0.771


# Extract the dose-response curve for MIA-PaCa-2 cells 
DR_panc <- subset(GDSC_dose_response, GDSC_dose_response$CELLNAME == "MIA-PaCa-2", select = c("NSC", "CELLNAME", "GROWTH", "LOG_CONCENTRATION")) 
Dat9 <- merge(DR_panc, ECFPs, by = "NSC")
Dat9 <- Dat9[, sapply(Dat9, function(x) length(unique(x))) > 1] # Remove the constant columns
# PCA 
DR_panc_pca <- prcomp(Dat9[,3:ncol(Dat9)], scale. = T, center = T)
plot(DR_panc_pca, type = "l")
std_dev <- DR_panc_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(cumsum(prop_varex[1:175]), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # 160PC captures more than 95%

DR_panc.train <- cbind(Dat9[,1:2], DR_panc_pca$x[,1:160]) 
# Random Forest
set.seed(123)
rf_DR_panc<- train(x = DR_panc.train[, 3:ncol(DR_panc.train)], 
                   y = DR_panc.train[, 2], 
                   trControl = tc, 
                   method = "ranger", 
                   importance = "impurity")  
rf_DR_panc$finalModel # R2: 0.846


# rank the models based on R2 values of fingerprint training
models <- c("rf3_DR_MB231", "rf_DR_HeLa", "rf_DR_MCF7", "rf_DR_HS578T", "rf_DR_panc")
r2s <- c(0.704, 0.771, 0.724, 0.749, 0.846)
compare.table <- as.data.frame(cbind(models, r2s))
colnames(compare.table) <- c("rf_models", "Rsquared")
