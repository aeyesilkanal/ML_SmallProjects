library(caret)
library(MLmetrics)
library(pROC)
library(dplyr)
library(regexr)
library(stringr)
library(ggplot2)

# Load the data frames
cl_exp_clean <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/Cell_line_COSMIC_ID_gene_expression_transposed_clean.tsv", header = FALSE)
cl_exp_t <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/Cell_line_RMA_proc_basalExp_transposed.tsv")
cl_details <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/Cell_Lines_Details.tsv")
cl_details_COSMIC <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/Cell_Lines_Details_COSMIC.tsv")
cl_RACS <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/RACS_in_cell_lines.tsv")
drug_SMILES_alias <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/GDSC-Drugs-with-SMILES-Aliases.tsv")
drugs_secreened <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/GDSC_Screened_Compounds.tsv")
gene_CN <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/Gene_level_CN.tsv")
dose_response <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/v17_fitted_dose_response.tsv")
WES_variants <- read.delim("~/Desktop/Ali_ML_HW2/Data_frames/WES_variants.tsv")

# Filter drug-cell line pairs to be used for training based on their AUC (AUC > 0.9)
dose_response_AUCorder <- dose_response[order(-dose_response$AUC),]
plot(dose_response_AUCorder$AUC, main = "AUC distribution across drug-cell line pairs", xlab = "index of drug - cell line pairs", ylab = "AUC")
doseRes_filt <- subset(dose_response_AUCorder, AUC >= 0.9, select = -Dataset_version)

# Group by drug id and isolate the list of drugs to be considered after filtering
doseRes_filt <- doseRes_filt[order(doseRes_filt$DRUG_ID),]
nrow(doseRes_filt[doseRes_filt$DRUG_ID == 1,]) # e.g. number of cell lines on which DRUG 1 has been tested - 353
drugs_filt <- unique(doseRes_filt$DRUG_ID)
length(drugs_filt) # total drug number was not reduced 

# Add column names to the clean data frame
col_names <- colnames(cl_exp_t)
colnames(cl_exp_clean) <- c("COSMIC_ID",col_names[2:length(col_names)])

# Create feature matrix and outcome vector for each drug 
for(i in 1:length(drugs_filt)) {
  drug_id <- drugs_filt[i]
  doseRes_sub <- subset(doseRes_filt, doseRes_filt$DRUG_ID == drug_id)
  cosmic_ids <- doseRes_sub$COSMIC_ID
  feat_mat <- subset(cl_exp_clean, cl_exp_clean$COSMIC_ID %in% cosmic_ids)
  feat_mat <- feat_mat[order(feat_mat$COSMIC_ID),]
  feat_name <- paste0("feat_mat_", drug_id)
  assign(feat_name, feat_mat)
  # Turns out not all cell lines that a particular drug has been tested on have expression data
  doseRes_sub2 <- subset(doseRes_sub, doseRes_sub$COSMIC_ID %in% feat_mat$COSMIC_ID, select = c(COSMIC_ID, LN_IC50))
  doseRes_sub2 <- doseRes_sub2[order(doseRes_sub2$COSMIC_ID),]
  outcome_name <- paste0("outcome_", drug_id)
  assign(outcome_name, doseRes_sub2)
} 

# Index the feature matrices for each drug
feat_mat_index <- ls(pattern = 'feat_mat_')
outcome_index <- ls(pattern = 'outcome_')[1:265]

# Filter features down to 30 by selecting features that are most highly correlated with the outcome
## First create matrices that have the COSMIC_ID in the first column, outcome in the second column, followed by all the features.
for (i in 1:length(outcome_index)) {
  feat <- get(feat_mat_index[i])
  outcome <- get(outcome_index[i])
  mat <- cbind(feat[,1],outcome[,2],feat[,2:ncol(feat)])
  colnames(mat)[c(1,2)] <- c("COSMIC_ID", "LN_IC50")
  drug_id <- str_split(outcome_index[i], "_")[[1]][2]
  mat_name <- paste0("fullMat_", drug_id)
  assign(mat_name, mat)
}
fullMat_index <- ls(pattern = 'fullMat_')

for (i in 1:length(fullMat_index)) {
  matx <- get(fullMat_index[i])
  write.table(matx, fullMat_index[i], sep = "\t", row.names = FALSE, quote = FALSE)
}

## For each feature matrix, extract the top 30 features that correlate with the outcome (positively or negatively)
for (i in 1:length(fullMat_index)) {
  mat <- get(fullMat_index[i])
  corr_mat <- cor(mat[,2:ncol(mat)], method = "spearman")
  outcome_corr <- corr_mat[1,]
  names <- names(sort(abs(outcome_corr), decreasing = TRUE)[2:31])
  drug_id <- str_split(fullMat_index[i], "_")[[1]][2]
  var_name <- paste0("topVar_", drug_id)
  assign(var_name, names)
}
topVar_index <- ls(pattern = 'topVar_')

## Subset the feature matrices based on the top 30 features 
for (i in 1:length(topVar_index)){
  mat <- get(fullMat_index[i])
  topVars <- get(topVar_index[i])
  subMat <- cbind(mat[,1:2], mat[,topVars])
  drug_id <- str_split(topVar_index[i], "_")[[1]][2]
  subMat_name <- paste0("subMat_", drug_id)
  assign(subMat_name, subMat)
}
subMat_index <- ls(pattern = 'subMat_[0-9]')

# Set up 5 fold cross-validation to be used in all training methods
# set.seed(100)
tc <- trainControl(method = "cv", number = 5)

# Run randomForest using "ranger" on all the drugs 
for (i in 1:length(subMat_index)) {
  matx <- get(subMat_index[i])
  ranger_temp <- train(x = matx[,3:ncol(matx)], y = matx[,2], trControl = tc, method = "ranger", importance = "impurity")
  drug_id <- str_split(subMat_index[i], "_")[[1]][2]
  model_name <- paste0("ranger_", drug_id)
  assign(model_name, ranger_temp)
}
ranger_index <- ls(pattern = 'ranger_[0-9]')

# for (i in 1:length(subMat_index)) {
#   matx <- get(subMat_index[i])
#   ranger_temp <- train(x = matx[,3:ncol(matx)], y = matx[,2], trControl = tc, method = "knn")
#   drug_id <- str_split(subMat_index[i], "_")[[1]][2]
#   model_name <- paste0("ranger_", drug_id)
#   assign(model_name, ranger_temp)
# }

# Aproaching the data set as a three-way categorical problem
## Z-tranform the LN_IC50 column and assign a class to each COSMIC_ID based on a -1/+1 Z-score threshold. 
## (-2/+2 was too stringent, some drugs did not have "sensitive" or "resistant)

cosmicClass <- function(x) {
  if(!is.data.frame(x)) stop("x must be a data frame")

  ic50_z <- scale(x[,2], center = TRUE, scale = TRUE)
  x1 <- cbind(ic50_z, x)
  colnames(x1)[1] <- "LN_IC50_z"

  cosmic_class <- c()
  for (i in 1:nrow(x1)){
    if (x1[i,1] <= -1) {
      cosmic_class <- c(cosmic_class, "sensitive")
    } else if (x1[i,1] >= 1) {
      cosmic_class <- c(cosmic_class, "resistant")
    } else {
      cosmic_class <- c(cosmic_class, "intermediate")
    }
  }
  cosmic_class <- as.factor(cosmic_class)
  x1 <- cbind(cosmic_class, x1)
  colnames(x1)[1] <- "COSMIC_class"
  x <- x1
  return(x)
}


for (i in 1:length(subMat_index)){
  matx <- get(subMat_index[i])
  matx1 <- cosmicClass(matx)
  drug_id <- str_split(subMat_index[i], "_")[[1]][2]
  subMat_name <- paste0("subMat_", drug_id)
  assign(subMat_name, matx1)
}

# for (i in 1:length(subMat_index)){
#   matx <- get(subMat_index[i])
#   matx1 <- matx[,-c(1:2)]
#   drug_id <- str_split(subMat_index[i], "_")[[1]][2]
#   subMat_name <- paste0("subMat_", drug_id)
#   assign(subMat_name, matx1)
# }

for (i in 1:length(subMat_index)) {
  matx <- get(subMat_index[i])
  ranger_temp <- train(x = matx[,5:ncol(matx)], y = matx[,1], trControl = tc, method = "ranger", importance = "impurity")
  drug_id <- str_split(subMat_index[i], "_")[[1]][2]
  model_name <- paste0("ranger2_", drug_id)
  assign(model_name, ranger_temp)
}
ranger2_index <- ls(pattern = 'ranger2_[0-9]')

# Isolate the R2 values and rank the regression models
r2_scores <- c()
for (i in 1:length(ranger_index)) {
  results <- get(ranger_index[i])
  Score <- results$finalModel$r.squared
  r2_scores <- c(r2_scores, Score)
}
regPerform <- cbind(ranger_index, r2_scores)
regPerform <- as.data.frame(regPerform)
regPerform$r2_scores <- sapply(regPerform$r2_scores, as.character)
regPerform$r2_scores <- sapply(regPerform$r2_scores, as.numeric)
regPerform_ranked <-  regPerform[order(regPerform$r2_scores, decreasing = TRUE),]
head(regPerform_ranked, n = 10)
tail(regPerform_ranked, n = 10)
write.table(regPerform_ranked, "rf_regres_ranked.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Isolate the F1 values and rank the classification models
F1_scores <- c()
for (i in 1:length(ranger2_index)) {
  results <- get(ranger2_index[i])
  subMat <- get(subMat_index[i])
  Score <- F1_Score(y_pred = results$finalModel$predictions, y_true = subMat$COSMIC_class)
  F1_scores <- c(F1_scores, Score)
}
classPerform <- cbind(ranger2_index, F1_scores)
classPerform <- as.data.frame(classPerform)
classPerform$F1_scores <- sapply(classPerform$F1_scores, as.character)
classPerform$F1_scores <- sapply(classPerform$F1_scores, as.numeric)
classPerform_ranked <-  classPerform[order(classPerform$F1_scores, decreasing = TRUE),]
head(classPerform_ranked, n = 10)
tail(classPerform_ranked, n = 10)
write.table(classPerform_ranked, "rf_class_ranked.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# PART 2:  Try to improve the predictions by looking incorporating other information about the cell lines and the drugs. 
### idea: add mutational profiles of the cell lines to predict drug response. filter the genes based on mutation frequency (# of rows per gene across all the variants identified)

mut_counts <- table(WES_variants$Gene)
mut_counts_sorted <- sort(mut_counts, decreasing = TRUE)
plot(mut_counts_sorted, main = "Ranking of mutated genes by frequency across cell lines", ylab = "Mutation count")
topMuts <- mut_counts_sorted[mut_counts_sorted >= 100] # 403 genes mutated over 100 times
topMuts_names <- names(topMuts)

WES_cosmic <- unique(WES_variants$COSMIC_ID)
clines <- c()
mut_profile <- c()
for (i in 1:length(WES_cosmic)) {
  cline <- WES_cosmic[i]
  genes <- subset(WES_variants, WES_variants$COSMIC_ID == cline, select = "Gene")[,1]
  cl_muts <- topMuts_names %in% genes
  mut_profile <- rbind(mut_profile, cl_muts)
  clines <- c(clines, cline)
}
mutMat <- cbind(clines, mut_profile)
cols <- paste0(topMuts_names,"_mut")
colnames(mutMat) <- c("COSMIC_ID", cols)
rownames(mutMat) <- NULL
write.table(mutMat, "mutMat.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# add the IC50 values to the mutation matrix for those cell lines that have it in order to make the new feature matrix.
for (i in 1:length(subMat_index)) {
  featMatx <- get(subMat_index[i])
  cosmic_overlap <- intersect(featMatx$COSMIC_ID, mutMat[,1])
  newFeatMatx <- subset(featMatx, featMatx$COSMIC_ID %in% cosmic_overlap)
  newMutMatx <- subset(mutMat, mutMat[,1] %in% cosmic_overlap)
  newFeatMatx <- newFeatMatx[order(newFeatMatx$COSMIC_ID),]
  newMutMatx <- newMutMatx[order(newMutMatx[,1]),]
  newSubMat <- cbind(newFeatMatx,newMutMatx[,2:ncol(newMutMatx)])
  drug_id <- str_split(subMat_index[i], "_")[[1]][2]
  subMat_name <- paste0("subMat2_", drug_id)
  assign(subMat_name, newSubMat)
}
subMat2_index <- ls(pattern = 'subMat2_[0-9]')

# Since the mutation features were added to the feature matrix columnwise in a ranked order,
# I can simply shrink the feature matrix by removing columns from the end of the matrix.
# This would select for most highly mutated/altered genes

# Pick the top 50 mutated genes (up to 84th column of the subMat2 format)
for (i in 1:length(subMat2_index)) {
  mat <- get(subMat2_index[i])
  mat <- mat[,1:84]
  drug_id <- str_split(subMat2_index[i], "_")[[1]][2]
  newMat_name <- paste0("geneNmut_", drug_id)
  assign(newMat_name, mat)
}
geneNmut_index <- ls(pattern = 'geneNmut_[0-9]')

# Train the drug predictors on the new matrix that has expression and highly mutated gene information
set.seed(100)
## regression
for (i in 1:length(geneNmut_index)) {
  matx <- get(geneNmut_index[i])
  ranger_temp <- train(x = matx[,5:ncol(matx)], y = matx[,4], trControl = tc, method = "ranger", importance = "impurity")
  drug_id <- str_split(geneNmut_index[i], "_")[[1]][2]
  model_name <- paste0("ranger3_", drug_id)
  assign(model_name, ranger_temp)
}
ranger3_index <- ls(pattern = 'ranger3_[0-9]')
## classification
for (i in 1:length(geneNmut_index)) {
  matx <- get(geneNmut_index[i])
  ranger_temp <- train(x = matx[,5:ncol(matx)], y = matx[,1], trControl = tc, method = "ranger", importance = "impurity")
  drug_id <- str_split(geneNmut_index[i], "_")[[1]][2]
  model_name <- paste0("ranger4_", drug_id)
  assign(model_name, ranger_temp)
}
ranger4_index <- ls(pattern = 'ranger4_[0-9]')

## Isolate the R2 values and rank the new regression models
r2_scores <- c()
for (i in 1:length(ranger3_index)) {
  results <- get(ranger3_index[i])
  Score <- results$finalModel$r.squared
  r2_scores <- c(r2_scores, Score)
}
regPerform <- cbind(ranger3_index, r2_scores)
regPerform <- as.data.frame(regPerform)
regPerform$r2_scores <- sapply(regPerform$r2_scores, as.character)
regPerform$r2_scores <- sapply(regPerform$r2_scores, as.numeric)
regPerform_ranked <-  regPerform[order(regPerform$r2_scores, decreasing = TRUE),]
head(regPerform_ranked, n = 10)
tail(regPerform_ranked, n = 10)
write.table(regPerform_ranked, "rf_regres2_ranked.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## Isolate the F1 values and rank the new classification models
F1_scores <- c()
for (i in 1:length(ranger4_index)) {
  results <- get(ranger4_index[i])
  subMat <- get(subMat2_index[i])
  Score <- F1_Score(y_pred = results$finalModel$predictions, y_true = subMat$COSMIC_class)
  F1_scores <- c(F1_scores, Score)
}
classPerform <- cbind(ranger4_index, F1_scores)
classPerform <- as.data.frame(classPerform)
classPerform$F1_scores <- sapply(classPerform$F1_scores, as.character)
classPerform$F1_scores <- sapply(classPerform$F1_scores, as.numeric)
classPerform_ranked <-  classPerform[order(classPerform$F1_scores, decreasing = TRUE),]
head(classPerform_ranked, n = 10)
tail(classPerform_ranked, n = 10)
write.table(classPerform_ranked, "rf_class2_ranked.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## Merge the first and second model for each drug to compare side by side.
regres1 <- read.delim("rf_regres_ranked.txt")
regres2 <- read.delim("rf_regres2_ranked.txt")
drugs_regres1 <- str_split(regres1$ranger_index, "_", simplify = TRUE)[,2]
reg1 <- cbind(regres1, drugs_regres1)
drugs_regres2 <- str_split(regres2$ranger3_index, "_", simplify = TRUE)[,2]
reg2 <- cbind(regres2, drugs_regres2)
colnames(reg1)[3] <- "DRUG_ID"
colnames(reg2)[3] <- "DRUG_ID"
regs <- merge(reg1, reg2, by = "DRUG_ID")
colnames(regs)[3] <- "r2_score_1"
colnames(regs)[5] <- "r2_score_2"
regs <- regs[order(regs$r2_score_1, decreasing = TRUE),]
regs <- subset(regs, select = c(DRUG_ID, r2_score_1, r2_score_2))
write.table(regs, "rf_reg1reg2_compare.txt", sep = "\t", quote = FALSE, row.names = FALSE)

class1 <- read.delim("rf_class_ranked.txt")
class2 <- read.delim("rf_class2_ranked.txt")
drugs_class1 <- str_split(class1$ranger2_index, "_", simplify = TRUE)[,2]
cla1 <- cbind(class1, drugs_class1)
drugs_class2 <- str_split(class2$ranger4_index, "_", simplify = TRUE)[,2]
cla2 <- cbind(class2, drugs_class2)
colnames(cla1)[3] <- "DRUG_ID"
colnames(cla2)[3] <- "DRUG_ID"
clas <- merge(cla1, cla2, by = "DRUG_ID")
colnames(clas)[3] <- "F1_score_1"
colnames(clas)[5] <- "F1_score_2"
clas <- clas[order(clas$F1_score_1, decreasing = TRUE),]
clas <- subset(clas, select = c(DRUG_ID, F1_score_1, F1_score_2))
write.table(clas, "rf_cla1cla2_compare.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# TRY OTHER MODELS
# glmboost, regression
for (i in 1:length(geneNmut_index)) {
  matx <- get(geneNmut_index[i])
  glmboost_temp <- train(x = matx[,5:ncol(matx)], y = matx[,4], trControl = tc, method = "glmboost")
  drug_id <- str_split(geneNmut_index[i], "_")[[1]][2]
  model_name <- paste0("glmboost_", drug_id)
  assign(model_name, glmboost_temp)
}
glmboost_index <- ls(pattern = 'glmboost_[0-9]')
# knn, regression
for (i in 1:length(geneNmut_index)) {
  matx <- get(geneNmut_index[i])
  knn_temp <- train(x = matx[,5:ncol(matx)], y = matx[,4], trControl = tc, method = "knn")
  drug_id <- str_split(geneNmut_index[i], "_")[[1]][2]
  model_name <- paste0("knn_", drug_id)
  assign(model_name, knn_temp)
}
knn_index <- ls(pattern = 'knn_[0-9]')
# compare these models to random forest from part 1
resamps_1 <- resamples(list(rf = ranger_1, glmboost = glmboost_1, knn = knn_1))
summary(resamps_1)

resamps_182 <- resamples(list(rf = ranger_182, glmboost = glmboost_182, knn = knn_182))
summary(resamps_182)

# PART 3: Pick a cell line and apply the predictors to it. 
## MDA-MB-231 (cosmic: 905960)
mb231_exp <- subset(cl_exp_clean, cl_exp_clean$COSMIC_ID == "905960")
### regression predictions
set.seed(100)
ic50_pred <- c()
for (i in 1:length(ranger_index)){
  model <- get(ranger_index[i])
  pred <- predict(model, mb231_exp)
  ic50_pred <- c(ic50_pred, pred)
}
predict.table <- as.data.frame(cbind(ranger_index, ic50_pred))
predict.table$ic50_pred <- as.numeric(as.character(predict.table$ic50_pred))

predict.table.ranked <- as.data.frame(predict.table[order(ic50_pred),])
predict.table.ranked$ic50_pred <- as.numeric(as.character(predict.table.ranked$ic50_pred))
rownames(predict.table.ranked) <- as.character(predict.table.ranked$ranger_index)
predict.table.ranked$z_scores <- as.numeric(as.character(scale(predict.table.ranked$ic50_pred, center = TRUE, scale = TRUE)))

p <- ggplot(predict.table.ranked, aes(x = rank(predict.table.ranked$z_scores), y = z_scores)) 
p + geom_point(stat = "identity") + ylim(-4,4) + labs(x = "drug rank", title = "Ranking of predicted LN_IC50 z scores on cell line MDA-MB-231 (COSMIC_ID:905960)")


write.table(predict.table, "predict_ranger_mb231.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(predict.table.ranked, "predict_ranger_ranked_mb231.txt", sep = '\t', quote = FALSE, row.names = FALSE)

### Rank the experimental IC_50s for MDA-MB-231
mb231_ic50s <- subset(doseRes_filt, doseRes_filt$COSMIC_ID == "905960", select = c("COSMIC_ID", "DRUG_ID", "LN_IC50"))
mb231_ic50s$z_scores <- as.numeric(as.character(scale(mb231_ic50s$LN_IC50, center = TRUE, scale = TRUE)))
mb231_ic50s_ranked <- mb231_ic50s[order(mb231_ic50s$z_scores),]
write.table(mb231_ic50s_ranked, "LN_IC50s_ranked_mb231.txt", sep = '\t', quote = FALSE, row.names = FALSE)
#### Compare the drug rankings between prediction and experimental (but first, filter the prediction ranking for those drugs that have been tested for the 231 cell line)
predict.table.ranked$DRUG_ID <- as.factor(str_split(predict.table.ranked$ranger_index, "_", simplify = TRUE)[,2])
drugs <- as.character(intersect(predict.table.ranked$DRUG_ID, mb231_ic50s_ranked$DRUG_ID))
subtable <- subset(predict.table.ranked, predict.table.ranked$DRUG_ID %in% drugs, select = c("DRUG_ID", "ic50_pred"))
newTable <- merge.data.frame(subtable, mb231_ic50s_ranked, by = "DRUG_ID")
mb231_compare <- newTable[order(newTable$ic50_pred),]
cor.test(mb231_compare$ic50_pred, mb231_compare$LN_IC50, method = "spearman")

## Automate this for all the cell lines tested with drugs. 
set.seed(100)

for (i in 1:nrow(cl_exp_clean)) {
  cl <- cl_exp_clean$COSMIC_ID[i]
  expression <- subset(cl_exp_clean, cl_exp_clean$COSMIC_ID == cl)
  
  ic50_pred <- c()
  for (j in 1:length(ranger_index)){
    model <- get(ranger_index[j])
    pred <- predict(model, expression)
    ic50_pred <- c(ic50_pred, pred)
  }
  
  predict.table <- as.data.frame(cbind(ranger_index, ic50_pred))
  predict.table$ic50_pred <- as.numeric(as.character(predict.table$ic50_pred))
  predict.table.ranked <- as.data.frame(predict.table[order(ic50_pred),])
  predict.table.ranked$ic50_pred <- as.numeric(as.character(predict.table.ranked$ic50_pred))
  predict.table.ranked$z_scores <- as.numeric(as.character(scale(predict.table.ranked$ic50_pred, center = TRUE, scale = TRUE)))
  
  if (cl %in% doseRes_filt$COSMIC_ID == TRUE) {
    cl_ic50s <- subset(doseRes_filt, doseRes_filt$COSMIC_ID == cl, select = c("COSMIC_ID", "DRUG_ID", "LN_IC50"))
    predict.table.ranked$DRUG_ID <- as.factor(str_split(predict.table.ranked$ranger_index, "_", simplify = TRUE)[,2])
    drugs <- as.character(intersect(predict.table.ranked$DRUG_ID, cl_ic50s$DRUG_ID))
    subtable <- subset(predict.table.ranked, predict.table.ranked$DRUG_ID %in% drugs, select = c("DRUG_ID", "ic50_pred"))
    newTable <- merge.data.frame(subtable, cl_ic50s, by = "DRUG_ID")
    cl_compare <- newTable[order(newTable$ic50_pred),]
    filename <- paste0("compare_", cl, ".txt")
    write.table(cl_compare, file = filename, sep = '\t', quote = FALSE, row.names = FALSE)
  }
}

### Extract the spearman rho score and the p-val.
setwd("~/Desktop/Ali_ML_HW2/compare")
compare_list <- list.files()

cellLine <- c()
rho <- c()
pval <- c()
for (i in 1:length(compare_list)) {
  filename <- compare_list[i]
  ic50_pred <- read.delim(filename)[,2]
  LNic50 <- read.delim(filename)[,4]
  results <- cor.test(ic50_pred, LNic50, method = "spearman")
  RHO <- results$estimate
  PVAL <- results$p.value
  rho <- c(rho, RHO)
  pval <- c(pval, PVAL)
  
  filename <- gsub(".txt", "", filename)
  cl <- gsub("compare_", "", filename)
  cellLine <- c(cellLine, cl)
}

rows <- c()
for (i in 1:length(compare_list)) {
  filename <- compare_list[i]
  nrows <- nrow(read.delim(filename))
  rows <- c(rows, nrows)
}
rows  
  
spearman_summary <- as.data.frame(cbind(cellLine, rho, pval, rows))
rownames(spearman_summary) <- NULL
colnames(spearman_summary) <- c("Cell_line", "rho", "pval", "num__of_drugs")
spearman_summary$rho <- as.numeric(as.character(spearman_summary$rho))
spearman_summary$pval <- as.numeric(as.character(spearman_summary$pval))

write.table(spearman_summary, "spearman_summary.txt", sep = '\t', quote = FALSE, row.names = FALSE) 
 
mb231_comparison <- get("compare_905960.txt")

# INDEPENDENT QUERIES
## table showing number of cell lines per drug 
geneNmut_index <- ls(pattern = 'geneNmut_[0-9]') # 265 drugs
cl_per_drug <- c()
drug_name <- c()
for (i in 1:length(geneNmut_index)) {
  cl_num <- nrow(get(geneNmut_index[i]))
  cl_per_drug <- c(cl_per_drug, cl_num)
  drug_id <- str_split(geneNmut_index[i], "_")[[1]][2]
  drug_name <- c(drug_name, drug_id)
}
table <- cbind(drug_name, cl_per_drug)
table <- as.data.frame(table)
colnames(table) <- c("DRUG_ID", "Num_CL_per_Drug")
write.table(table, "Num_CL_per_Drug.txt", sep = '\t', quote = FALSE, row.names = FALSE) 

## random forest models (regression and classification) with respect to the number of cell lines they were trained on. 
b <- read.delim("rf_regres_ranked.txt")
drug_ids_ranked <- str_split(b$ranger_index, "_", simplify = TRUE)[,2]
b <- cbind(b, drug_ids_ranked)
colnames(b)[3] <- "DRUG_ID"
c <- merge(b, table, by = "DRUG_ID")
c <- c[order(c$r2_scores, decreasing = TRUE),]
write.table(c, "rf_regres_ranked2.txt", sep = '\t', quote = FALSE, row.names = FALSE)
c$Num_CL_per_Drug <- as.numeric(as.character(c$Num_CL_per_Drug))
plot(x = rank(-c$r2_scores), y = c$Num_CL_per_Drug, main = "Correlation between r2 score of the model \n and the number of cell lines (per drug)", xlab = "ranking (highest to lowest r2 score)", ylab = "# of cell lines per drug")

d <- read.delim("rf_class_ranked.txt")
drug_ids_ranked2 <- str_split(d$ranger2_index, "_", simplify = TRUE)[,2]
d <- cbind(d, drug_ids_ranked2)
colnames(d)[3] <- "DRUG_ID"
e <- merge(d, table, by = "DRUG_ID")
e <- e[order(e$F1_scores, decreasing = TRUE),]
write.table(e, "rf_class_ranked2.txt", sep = '\t', quote = FALSE, row.names = FALSE)
e$Num_CL_per_Drug <- as.numeric(as.character(e$Num_CL_per_Drug))
plot(x = rank(-e$F1_scores), y = e$Num_CL_per_Drug, main = "Correlation between F1 score of the model \n and the number of cell lines (per drug)", xlab = "ranking (highest to lowest F1 score)", ylab = "# of cell lines per drug")


## remove data frames from the environment
rm(list = feat_mat_index)
rm(list = outcome_index)
rm(list = fullMat_index)
rm(list = subMat_index)

