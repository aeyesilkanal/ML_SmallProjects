library(caret)
library(MLmetrics)
library(pROC)
library(dplyr)
library(regexr)
library(stringr)
library(ggplot2)
library(readxl)
library(zoo)

# Load the necessary data sets
census_clean <- read_csv("C:/Users/Ali/Desktop/MLiC - HW3/census_clean.csv")
melanomav08 <- read_delim("C:/Users/Ali/Desktop/MLiC - HW3/melanomav08.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Avg_Daily_Sunlight <- read_excel("C:/Users/Ali/Desktop/MLiC - HW3/NLDAS Daily Sunlight 1979-2011_clean.xlsx")

colnames(Avg_Daily_Sunlight)[2] <- "FIPS"

# Merge the data frames
census_clean_df <- as.data.frame(census_clean)
melanomav08_df <- as.data.frame(melanomav08)
Avg_Daily_Sunlight_df <- as.data.frame(Avg_Daily_Sunlight)

census_clean_df$FIPS <- as.factor(census_clean_df$FIPS)
melanomav08_df$FIPS <- as.factor(melanomav08_df$FIPS)
Avg_Daily_Sunlight_df$FIPS <- as.factor(Avg_Daily_Sunlight_df$FIPS)

data <- merge(melanomav08_df, census_clean_df, by = "FIPS")
data <- merge(data, Avg_Daily_Sunlight_df, by = "FIPS")

# Select the feature matrix
feat_mat <- subset(data, select = c("FIPS", "Age-Adjusted Incidence Rate", "Value (Dollars)", "Avg Daily Sunlight (KJ/m²)"))

# Omit the rows with missing data
feat_mat.naOmit <- na.omit(feat_mat)

# Train the model using random forest, linear model, or generalized linear model with 5-fold cross validation
tc <- trainControl(method = "cv", number = 5)

set.seed(100)
ranger1.1 <- train(x = feat_mat.naOmit[,3:ncol(feat_mat.naOmit)], y = feat_mat.naOmit[,2], trControl = tc, method = "ranger", importance = "impurity") # random forest
lm1.1 <- train(x = feat_mat.naOmit[,3:ncol(feat_mat.naOmit)], y = feat_mat.naOmit[,2], trControl = tc, method = "lm") # linear model
glm1.1 <- train(x = feat_mat.naOmit[,3:ncol(feat_mat.naOmit)], y = feat_mat.naOmit[,2], trControl = tc, method = "glm") # generalized linear model

summary(resamples(list(ranger = ranger1.1, lm = lm1.1, glm = glm1.1)))

# Increase the number of features by adding:
### Average UV exposure per county
### Min, Max, and average incidence rate in the neighboring counties
### Min, Max, and average UV exposure in the neighboring counties
uv_county <- read_excel("uv-county.xlsx")
county_adjacency <- read.delim2("/media/rosner-lab/ubuntu_hdd1/Ali_MLinC_HWs/Ali_ML_HW3/MLiC - HW3/county_adjacency.txt", header=FALSE)

county_adj <- na.locf(county_adjacency)
colnames(county_adj) <- c("County_center", "FIPS_center", "County_neighbor", "FIPS")
county_adj$FIPS <- gsub(" ","", county_adj$FIPS)
county_adj <- subset(county_adj, select = c("FIPS_center", "FIPS"))

uv_county$COUNTY_FIPS <- substr(uv_county$COUNTY_FIPS,regexpr("[^0]",uv_county$COUNTY_FIPS),nchar(uv_county$COUNTY_FIPS))
colnames(uv_county)[3] <- "FIPS"
uv_county$FIPS <- as.factor(uv_county$FIPS)
uv_county <- subset(uv_county, select = c("FIPS", "UV_ Wh/m²"))

## subset the adjacency data for counties that are in the intial feature matrix (NAs omitted)
county_adj <- subset(county_adj, county_adj$FIPS %in% feat_mat.naOmit$FIPS)

Neighbor_incidence <- c()
for (i in 1:nrow(county_adj)) {
  fips <- county_adj$FIPS[i]
  rate <- feat_mat.naOmit[feat_mat.naOmit$FIPS == fips, "Age-Adjusted Incidence Rate"]
  Neighbor_incidence <- c(Neighbor_incidence,rate)
}
new_feats <- cbind(county_adj, Neighbor_incidence)
neighbor_ave_incedence <- aggregate(new_feats$Neighbor_incidence, list(new_feats$FIPS_center), mean)
neighbor_min_incedence <- aggregate(new_feats$Neighbor_incidence, list(new_feats$FIPS_center), min)
neighbor_max_incedence <- aggregate(new_feats$Neighbor_incidence, list(new_feats$FIPS_center), max)

colnames(neighbor_ave_incedence) <- c("FIPS", "neighbor_ave_incedence")
colnames(neighbor_min_incedence) <- c("FIPS", "neighbor_min_incedence")
colnames(neighbor_max_incedence) <- c("FIPS", "neighbor_max_incedence")

new_feat_mat <- merge(feat_mat.naOmit, neighbor_ave_incedence, by = "FIPS")
new_feat_mat <- merge(new_feat_mat, neighbor_min_incedence, by = "FIPS")
new_feat_mat <- merge(new_feat_mat, neighbor_max_incedence, by = "FIPS")

## add the UV data from each county and the neighboring counties
new_feat_mat <- merge(new_feat_mat, uv_county, by = "FIPS")

Neighbor_uv <- c()
for (i in 1:nrow(county_adj)) {
  fips <- county_adj$FIPS[i]
  uv <- uv_county[uv_county$FIPS == fips, "UV_ Wh/m²"]
  Neighbor_uv <- c(Neighbor_uv,uv)
}
Neighbor_uv <- unlist(Neighbor_uv, use.names = FALSE)

uv_feats <- cbind(county_adj, Neighbor_uv)
neighbor_ave_uv <- aggregate(uv_feats$Neighbor_uv, list(uv_feats$FIPS_center), mean)
neighbor_min_uv <- aggregate(uv_feats$Neighbor_uv, list(uv_feats$FIPS_center), min)
neighbor_max_uv <- aggregate(uv_feats$Neighbor_uv, list(uv_feats$FIPS_center), max)

colnames(neighbor_ave_uv) <- c("FIPS", "neighbor_ave_uv")
colnames(neighbor_min_uv) <- c("FIPS", "neighbor_min_uv")
colnames(neighbor_max_uv) <- c("FIPS", "neighbor_max_uv")

new_feat_mat <- merge(new_feat_mat, neighbor_ave_uv, by = "FIPS")
new_feat_mat <- merge(new_feat_mat, neighbor_min_uv, by = "FIPS")
new_feat_mat <- merge(new_feat_mat, neighbor_max_uv, by = "FIPS")

## Train the model using random forest with 5-fold cross validation
set.seed(100)
ranger2.1 <- train(x = new_feat_mat[,3:ncol(new_feat_mat)], y = new_feat_mat[,2], trControl = tc, method = "ranger", importance = "impurity") # random forest
lm2.1 <- train(x = new_feat_mat[,3:ncol(new_feat_mat)], y = new_feat_mat[,2], trControl = tc, method = "lm") # linear model
glm2.1 <- train(x = new_feat_mat[,3:ncol(new_feat_mat)], y = new_feat_mat[,2], trControl = tc, method = "glm") # generalized linear model

summary(resamples(list(ranger = ranger2.1, lm = lm2.1, glm = glm2.1)))

## Plot the variable importance based on the random forest model
vars <- ranger2.1$finalModel$variable.importance
par(mar=c(12,16,5,16))
a <- barplot(sort(vars, decreasing = TRUE), las=1, names.arg="", ylab = "Variable importance", cex.axis = 0.8)
text(a[,1], -10, srt = 45, adj= 1.1, xpd = TRUE, labels = names(sort(vars, decreasing = TRUE)) , cex=0.8)

## Export the feature matrix
write.table(new_feat_mat, "new_feature_mat.txt", sep = "\t", quote = FALSE, row.names = FALSE)
