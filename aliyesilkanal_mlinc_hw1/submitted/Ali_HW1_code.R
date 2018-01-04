#####################################################################################################################
########### Machine Learning in Cancer - HW1 - Ali Ekrem Yesilkanal #################################################
#####################################################################################################################

# install and load the necessary packages
library(jsonlite)
library(regexr)
library(stringr)
library(randomForest)
library(e1071)
library(caret)
library(doParallel)
library(foreach)
library(ranger)
library(gplots)

# parse the JSON
rawdata <- fromJSON("out.json")
data_init <- rawdata$data$hits
ssm <- rawdata$data$hits$ssm
cases <- rawdata$data$hits$case

# flatten the data table
mutation_subtype <- ssm$mutation_subtype #all "single base substitution"
chromosome <- ssm$chromosome
mut_position <- ssm$end_position
change <- sapply(ssm$genomic_dna_change, function(x) FUN = str_sub(x, -3))
case_id <- cases$case_id
primary_site <- cases$primary_site
mut_id <- data_init$id

## unlist the gene names
genes_temp <- unlist(ssm$consequence, recursive = FALSE, use.names = FALSE)
gene_symbol <- c()
for (i in 1:length(genes_temp)) {
  name <- unlist(genes_temp[[i]][1,])
  gene_symbol <- c(gene_symbol,name)
}
gene_symbol <- as.vector(gene_symbol)

# Construct the clean data frame
clean_data <- cbind(mut_id,case_id,primary_site, gene_symbol, chromosome, mut_position, change)
colnames(clean_data) <- c("mut_id", "case_id", "primary_site", "gene_symbol", "chromosome", "mut_position", "change")
rownames(clean_data) <- NULL
clean_data <- as.data.frame(clean_data)
write.table(clean_data, "clean_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Explore the data set
## Group by case_id (patient)
clean_case_sorted <- clean_data[order(clean_data$case_id),]
write.table(clean_case_sorted, "clean_case_sorted.txt", sep = "\t", quote = FALSE, row.names = FALSE)
## Example case of a gene that is mutated multiple times within the same case_id
ex_case <- clean_case_sorted[clean_case_sorted$case_id == "0011a67b-1ba9-4a32-a6b8-7850759a38cf",]
ex_case[ex_case$gene_symbol == "ZNF714",]


############################################################################################################################################
# Build features

## Feature 1:  Most frequent DNA change in a case  #########################################################################################
most_freq_change <- c()
for (i in 1:nrow(case_site_index)) {
  change_temp <- names(sort(table(clean_data[clean_data$case_id == case_site_index[i,1],"change"]), decreasing = TRUE)[1])
  most_freq_change <- c(most_freq_change, change_temp)
}
case_site_index_feat1 <- cbind(case_site_index, most_freq_change)

varNames <- names(table(case_site_index_feat1$most_freq_change)) #"A>C" "A>G" "A>T" "C>A" "C>G" "C>T" "G>A" "G>C" "G>T" "T>A" "T>C" "T>G"

changeHits <- c()
for (i in 1:nrow(case_site_index_feat1)){
  change_temp <- varNames %in% case_site_index_feat1$most_freq_change[i]
  changeHits <- rbind(changeHits, change_temp)
}
rownames(changeHits) <- NULL
varNames1 <- gsub(">","to", varNames)
colnames(changeHits) <- varNames1
case_site_index_feat1 <- cbind(case_site_index_feat1, changeHits)
write.table(case_site_index_feat1, "case_site_index_feat1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## For the random forest, we need one column with the classes ("primary_site") and the feature matrix 
case_site_index_feat1 <- case_site_index_feat1[, c(2, 4:ncol(case_site_index_feat1))]

## Split up the data set at 80% training and 20% test
set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat1), replace = T, prob = c(0.8,0.2))
case.train <- case_site_index_feat1[sample.ind==1,]
case.test <- case_site_index_feat1[sample.ind==2,]

## Primary site distribution in the Training vs. Test
table(case.train$primary_site)/nrow(case.train)*100
table(case.test$primary_site)/nrow(case.test)*100

## Run RandomForest with 'ranger'
primarySite.ranger.rf <- ranger(formula = NULL,
                                dependent.variable.name = "primary_site",
                                case.train,
                                num.trees = 500,
                                importance= 'impurity')

importance(primarySite.ranger.rf)

primarySite.ranger.rf$predictions
primarySite.ranger.rf$prediction.error # 0.84 error suggests 16% accuracy
primarySite.ranger.rf$variable.importance
barplot(sort(primarySite.ranger.rf$variable.importance, decreasing = TRUE), ylab = "Mean Decerease Gini", xlab = "1st feature set variables", main = "Variable importance")

predict.rf <- predict(primarySite.ranger.rf, data = case.test)
predict.rf$predictions
confusion.matrix.test <- table(case.test$primary_site, predict.rf$predictions)
heatmap.2(x = confusion.matrix.test, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix.test,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "1st Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
dev.off()

# Calculate the accuracy of the prediction in the test set
accuracy <- function(x) {
  total.matched <- sum(diag(x))
  total <- sum(x)
  percent.match <- (total.matched/total)*100
  return(percent.match)
}

accuracy(confusion.matrix.test) # 17.36% accuracy


## Feature 2: Genes mutated ####################################################################################################
genes_all <- unique(gene_symbol) #18177 genes mutated total

mutHits <- c()
for (i in 1:nrow(case_site_index)) {
  case_id_temp <- case_site_index[i,1]
  mutated_genes <- clean_data[case_id == case_id_temp,"gene_symbol"]
  mut_hits <- genes_all %in% mutated_genes
  mutHits <- rbind(mutHits, mut_hits)
}
colnames(mutHits) <- genes_all
rownames(mutHits) <- NULL
case_site_index_feat2 <- cbind(case_site_index, mutHits)
write.table(case_site_index_feat2, "case_site_index_feat2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Drop "case_id" and convert to a data.frame
case_site_index_feat2 <- case_site_index_feat2[,2:ncol(case_site_index_feat2)]

## Split up the data set at 80% training and 20% test
set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat2), replace = T, prob = c(0.8,0.2))
case.train2 <- case_site_index_feat2[sample.ind==1,]
case.test2 <- case_site_index_feat2[sample.ind==2,]
## Primary site distribution in the Training vs. Test
table(case.train2$primary_site)/nrow(case.train)*100
table(case.test2$primary_site)/nrow(case.test)*100

primarySite.ranger.rf2 <- ranger(formula = NULL, 
                               dependent.variable.name = "primary_site",
                               case.train2,
                               num.trees = 500,
                               importance= 'impurity')

primarySite.ranger.rf2$predictions
primarySite.ranger.rf2$prediction.error # 0.74 error suggests 26% accuracy
primarySite.ranger.rf2$variable.importance
barplot(sort(primarySite.ranger.rf2$variable.importance, decreasing = TRUE), ylab = "Mean Decerease Gini", xlab = "2nd feature set variables", main = "Variable importance")

# top 20 important genes mutated
barplot(sort(primarySite.ranger.rf2$variable.importance, decreasing = TRUE)[1:20], ylab = "Mean Decerease Gini", xlab = "2nd feature set variables", main = "Variable importance - Top 20", las = 2, cex.names = 0.8)

# Test the forest on the test set
predict.rf2 <- predict(primarySite.ranger.rf2, data = case.test2)
predict.rf2$predictions
confusion.matrix.test2 <- table(case.test2$primary_site, predict.rf2$predictions)
heatmap.2(x = confusion.matrix.test2, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix.test2,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "2nd Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
dev.off()

accuracy(confusion.matrix.test2) # 25.23%

# Select genes that have higher importance value
var_imp <- primarySite.ranger.rf2$variable.importance
length(var_imp[var_imp > 1]) # 815 genes
sort(var_imp[var_imp > 1]) # APC, ATRX, VHL, PBRM1, PTEN, CDH1, CDKN2A, TP53 are among the top 
genes_important <- names(var_imp[var_imp > 1])
case_site_index_feat2.1 <- case_site_index_feat2[,colnames(case_site_index_feat2) %in% genes_important]

# subset the feature matrix with the 815 important features
case.train2.impGenes <- case.train2[,colnames(case.train2) %in% genes_important]
dim(case.train2.impGenes)
case.train2.impGenes <- cbind(case.train2$primary_site, case.train2.impGenes)
colnames(case.train2.impGenes)[1] <- "primary_site"

case.test2.impGenes <- case.test2[,colnames(case.test2) %in% genes_important]
dim(case.test2.impGenes)
case.test2.impGenes <- cbind(case.test2$primary_site, case.test2.impGenes)
colnames(case.test2.impGenes)[1] <- "primary_site"

# Rerun ranger with the important features
primarySite.ranger.rf2.imp <- ranger(formula = NULL, 
                                 dependent.variable.name = "primary_site",
                                 case.train2.impGenes,
                                 num.trees = 500,
                                 importance= 'impurity')

primarySite.ranger.rf2.imp$prediction.error # 0.731 error in the training (stays the same)

predict.rf2.imp <- predict(primarySite.ranger.rf2.imp, data = case.test2.impGenes)
predict.rf2.imp$predictions
confusion.matrix.test2.imp <- table(case.test2.impGenes$primary_site, predict.rf2.imp$predictions)
heatmap.2(x = confusion.matrix.test2.imp, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix.test2,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "2nd Feature Set (shrunk) Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
dev.off()
accuracy(confusion.matrix.test2.imp) # 25.5% did not improve much (stays the same)

## Feature 3: Total number of mutations per case  ##############################################################################
num.mut <- numeric()
for (i in 1:nrow(case_site_index)) {
  case <- case_site_index$case_id[i]
  num <- nrow(clean_data[clean_data$case_id == case,])
  num.mut <- c(num.mut, num)
}

case_site_index_feat3 <- cbind(case_site_index$primary_site, case_site_index_feat2.1, num.mut)
colnames(case_site_index_feat3)[1] <- "primary_site"

set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat3), replace = T, prob = c(0.8,0.2))
case.train3 <- case_site_index_feat3[sample.ind==1,]
case.test3 <- case_site_index_feat3[sample.ind==2,]
## Primary site distribution in the Training vs. Test
table(case.train3$primary_site)/nrow(case.train)*100
table(case.test3$primary_site)/nrow(case.test)*100

primarySite.ranger.rf3 <- ranger(formula = NULL, 
                                dependent.variable.name = "primary_site",
                                case.train3,
                                num.trees = 500,
                                importance= 'impurity')

primarySite.ranger.rf3$predictions
primarySite.ranger.rf3$prediction.error # 0.715 error suggests 26% accuracy
primarySite.ranger.rf3$variable.importance
barplot(sort(primarySite.ranger.rf3$variable.importance, decreasing = TRUE)[1:20], ylab = "Mean Decerease Gini", xlab = "3rd feature set variables - Top 20", main = "Variable importance",las = 2, cex.names = 0.8)

predict.rf3 <- predict(primarySite.ranger.rf3, data = case.test3)
predict.rf3$predictions
confusion.matrix.test3 <- table(case.test3$primary_site, predict.rf3$predictions)
heatmap.2(x = confusion.matrix.test3, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix.test3,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "3rd Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
accuracy(confusion.matrix.test3) # 28.32%


## Feature 4: Combine feat3 with feat1 ########################################################################################
case_site_index_feat4 <- cbind(case_site_index_feat3, case_site_index_feat1[,4:15])

set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat4), replace = T, prob = c(0.8,0.2))
case.train4 <- case_site_index_feat4[sample.ind==1,]
case.test4 <- case_site_index_feat4[sample.ind==2,]
## Primary site distribution in the Training vs. Test
table(case.train4$primary_site)/nrow(case.train)*100
table(case.test4$primary_site)/nrow(case.test)*100

primarySite.ranger.rf4 <- ranger(formula = NULL, 
                                 dependent.variable.name = "primary_site",
                                 case.train4,
                                 num.trees = 500,
                                 importance= 'impurity')

primarySite.ranger.rf4$prediction.error # 0.693 error
barplot(sort(primarySite.ranger.rf4$variable.importance, decreasing = TRUE)[1:20], ylab = "Mean Decerease Gini", xlab = "4th feature set variables - Top 20", main = "Variable importance",las = 2, cex.names = 0.8)

predict.rf4 <- predict(primarySite.ranger.rf4, data = case.test4)
predict.rf4$predictions
confusion.matrix.test4 <- table(case.test4$primary_site, predict.rf4$predictions)
heatmap.2(x = confusion.matrix.test4, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'row', 
          col = bluered(100),
          cellnote = confusion.matrix.test4,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "4th Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
accuracy(confusion.matrix.test4) # 30.41%

## Feature 5: Chromosomes ####################################################################################################
chrom_all <- unique(clean_data$chromosome)

chromosomeHits <- c()
for (i in 1:nrow(case_site_index)) {
  case_id_temp <- case_site_index[i,1]
  chroms <- clean_data[case_id == case_id_temp,"chromosome"]
  chrom_hits <- chrom_all %in% chroms
  chromosomeHits <- rbind(chromosomeHits, chrom_hits)
}
colnames(chromosomeHits) <- chrom_all
rownames(chromosomeHits) <- NULL
chromosomeHits <- as.data.frame(chromosomeHits)
case_site_index_feat5 <- cbind(case_site_index$primary_site, chromosomeHits)
colnames(case_site_index_feat5)[1] <- "primary_site"

set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat5), replace = T, prob = c(0.8,0.2))
case.train5 <- case_site_index_feat5[sample.ind==1,]
case.test5 <- case_site_index_feat5[sample.ind==2,]


primarySite.ranger.rf5 <- ranger(formula = NULL, 
                                 dependent.variable.name = "primary_site",
                                 case.train5,
                                 num.trees = 500,
                                 importance= 'impurity')

primarySite.ranger.rf5$prediction.error # 0.821 error
barplot(sort(primarySite.ranger.rf5$variable.importance, decreasing = TRUE), ylab = "Mean Decerease Gini", xlab = "5th feature set variables", main = "Variable importance",las = 2, cex.names = 0.8)


predict.rf5 <- predict(primarySite.ranger.rf5, data = case.test5)
predict.rf5$predictions
confusion.matrix.test5 <- table(case.test5$primary_site, predict.rf5$predictions)
heatmap.2(x = confusion.matrix.test5, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix.test5,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "5th Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
accuracy(confusion.matrix.test5) # 17.796% very poor by itself. What if we combine it with the rest of the features?

## Feature 6: Combine chromosome hits with feat 4 ############################################################################
case_site_index_feat6 <- cbind(case_site_index_feat4, chromosomeHits)

set.seed(111)
sample.ind <- sample(2, nrow(case_site_index_feat6), replace = T, prob = c(0.8,0.2))
case.train6 <- case_site_index_feat6[sample.ind==1,]
case.test6 <- case_site_index_feat6[sample.ind==2,]


primarySite.ranger.rf6 <- ranger(formula = NULL, 
                                 dependent.variable.name = "primary_site",
                                 case.train6,
                                 num.trees = 500,
                                 importance= 'impurity')

primarySite.ranger.rf6$prediction.error # 0.712 error
barplot(sort(primarySite.ranger.rf6$variable.importance, decreasing = TRUE)[1:40], ylab = "Mean Decerease Gini", xlab = "6th feature set variables - Top 40", main = "Variable importance",las = 2, cex.names = 0.8)

predict.rf6 <- predict(primarySite.ranger.rf6, data = case.test6)
predict.rf6$predictions
confusion.matrix.test6 <- table(case.test6$primary_site, predict.rf6$predictions)
heatmap.2(x = confusion.matrix.test6, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'row', 
          col = bluered(100),
          cellnote = confusion.matrix.test6,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10),
          main = "6th Feature Set Predictions - Confusion Matrix",
          ylab = "True class",
          xlab = "Predicted class")
accuracy(confusion.matrix.test6) # 29.6% adding the chromosome mutation information did not add much


## Bin low frequency cancer types as "other" ##################################################################################33
sort(table(case_site_index$primary_site), decreasing = TRUE)
other <- c("Lymph Nodes", "Bile Duct", "Eye", "Thymus", "Pleura", "Bone Marrow", "Testis")
case_site_index2 <- case_site_index
levels(case_site_index2$primary_site)[c(2,4,10,15,18,23,24)] <- "other"
table(case_site_index2$primary_site)

case_site_index2_feat6 <- cbind(case_site_index2$primary_site, case_site_index_feat6[,2:853])
colnames(case_site_index2_feat6)[1] <- "primary_site"

set.seed(111)
sample.ind <- sample(2, nrow(case_site_index2_feat6), replace = T, prob = c(0.8,0.2))
case2.train6 <- case_site_index2_feat6[sample.ind==1,]
case2.test6 <- case_site_index2_feat6[sample.ind==2,]
table(case2.train6$primary_site)/nrow(case.train)*100
table(case2.test6$primary_site)/nrow(case.test)*100


primarySite2.ranger.rf6 <- ranger(formula = NULL, 
                                 dependent.variable.name = "primary_site",
                                 case2.train6,
                                 num.trees = 500,
                                 importance= 'impurity')

primarySite2.ranger.rf6$prediction.error # 0.712 error
barplot(primarySite2.ranger.rf6$variable.importance)

predict2.rf6 <- predict(primarySite2.ranger.rf6, data = case2.test6)
predict2.rf6$predictions
confusion.matrix2.test6 <- table(case2.test6$primary_site, predict2.rf6$predictions)
heatmap.2(x = confusion.matrix2.test6, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = 'none', 
          symm = TRUE, 
          trace = 'none', 
          scale = 'column', 
          col = bluered(100),
          cellnote = confusion.matrix2.test6,
          notecol = "black",
          keysize = 1.0,
          density.info = 'none',
          margins = c(10,10))
accuracy(confusion.matrix2.test6) # 29.42% Accuracy
