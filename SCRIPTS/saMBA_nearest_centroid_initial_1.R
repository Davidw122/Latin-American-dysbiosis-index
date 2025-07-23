install.packages("class", "randomForest", "caret")
install.packages("caret", dependencies = c("Depends", "Suggests"))
library("tidyverse")
library("vegan")
library("ggplot2")
library("MASS")
library(class)
library(caret)
library(randomForest)
library(dplyr)

# random forest - nearest centroid

metadata <- metadata_filtered

table(metadata$instrument_model)
table(metadata$broad_scale_environmental_context)

# finds "HiSeq" and gets rid of it
keep_samples <- metadata$run_accession[!grepl("HiSeq", metadata$instrument_model)]

keep_samples

metadata_filtered <- metadata %>% filter(run_accession %in% keep_samples)

table(metadata_filtered$instrument_model)
table(metadata_filtered$broad_scale_environmental_context)

prevalence_threshold <- 0.1

prevalence <- rowSums(count_table_filtered > 0) / ncol(count_table_filtered)
keep_genera <- prevalence >= prevalence_threshold

count_table_final <- count_table_filtered[keep_genera, ]
count_table_rel <- sweep(x = count_table_final, 
                         MARGIN = 2, 
                         STATS = colSums(count_table_final),
                         FUN = "/")
colSums(count_table_rel)

pseudocount <- 1e-6
count_table_log <- log10(count_table_rel + pseudocount)

data_matrix <- t(count_table_log)

# Make sure rownames match
common_samples <- intersect(rownames(data_matrix), metadata_filtered$run_accession)

# Filter both to keep only common samples
data_matrix <- data_matrix[common_samples, ]
metadata_filtered <- metadata_filtered[metadata_filtered$run_accession %in% common_samples, ]

# order metadata
metadata_ordered <- metadata_filtered %>% 
  mutate(order = match(rownames(data_matrix), run_accession)) %>%
  arrange(order)


#data_matrix <- readRDS("data_matrix.RDS")
#metadata_ordered <- readRDS("metadata_ordered.RDS")

  
  industrialization_labels <- ifelse(
    metadata_ordered$broad_scale_environmental_context %in% c("dense settlement biome", "urban biome", "urban_location"),
    "More_Industrialized",
    "Less_Industrialized"
  )

industrialization_labels <- factor(industrialization_labels, levels = c("Less_Industrialized", "More_Industrialized"))

  

#number of samples
nrow(data_matrix)
# number of features
ncol(data_matrix)

table(industrialization_labels)


set.seed(123)
train_indices <- createDataPartition(industrialization_labels, p=0.7, list=FALSE)
train_data <- data_matrix[train_indices, ]
test_data <- data_matrix[-train_indices, ]
train_labels <- industrialization_labels[train_indices]
test_labels <- industrialization_labels[-train_indices]

length(train_labels)
length(test_labels)
table(train_labels)
table(test_labels)

# train random forest model
rf_model <- randomForest(x = train_data, y=train_labels, ntree = 800, importance = TRUE)

rf_predictions <- predict(rf_model, test_data)

confusion_matrix <- confusionMatrix(rf_predictions, test_labels)
print(confusion_matrix)

importance_scores <- importance(rf_model)
importance_scores

importance_df <- data.frame(
  Genus = rownames(importance_scores),
  MeanDecreaseAccuracy = importance_scores[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance_scores[, "MeanDecreaseGini"]
)


top_genera <- importance_df[order(-importance_df$MeanDecreaseAccuracy), ][1:20, ]

top_genera

plot1 <- ggplot(top_genera, 
                aes(x = reorder(Genus, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
                geom_bar(stat = "identity", fill = "steelblue") +
                coord_flip() +
                labs(x = "Genus", y = "Mean Decrease in Accuracy") +
                theme_minimal()
                
             

plot1



centroid_more <- colMeans(data_matrix[industrialization_labels == "More_Industrialized", , drop = FALSE])
centroid_less <- colMeans(data_matrix[industrialization_labels == "Less_Industrialized", , drop = FALSE])

euclid_dist <- function(sample_vec, centroid_vec) {
  sqrt(sum((sample_vec - centroid_vec) ^ 2))
  
}

DistMore <- apply(data_matrix, 1, euclid_dist, centroid_vec = centroid_more)
DistLess <- apply(data_matrix, 1, euclid_dist, centroid_vec = centroid_less)

industrialization_index <- 1 - (DistMore / (DistMore + DistLess))

metadata_with_index <- metadata_ordered
metadata_with_index$industrialization_index <- industrialization_index

metadata_with_index$pedict_class <- 
  ifelse(DistMore < DistLess, "More_Industrialized", "Less_Industrialized")

plot2 <- ggplot(metadata_with_index, 
                aes(x = broad_scale_environmental_context,
                    y = industrialization_index)) +
  geom_boxplot(aes(fill = broad_scale_environmental_context)) +
  geom_jitter(width = .2)

plot2

anova_result <- aov(industrialization_index ~ broad_scale_environmental_context, data = metadata_with_index)
summary(anova_result)
TukeyHSD(anova_result)
confusion_matrix
adonis2(bray_dist ~ industrialization_labels, data = metadata_ordered)


