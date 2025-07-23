install.packages("tidyverse", "vegan", "ggplot2", "plotly")
install.packages("BiocManager")
BiocManager::install("sva") 

library("tidyverse")
library("vegan")
library("ggplot2")
library("plotly")
library("sva")

# note about accessing charts
# [left of comma is rows, right of comma is columns]
# [,] this accesses all rows & columns
# [,column 1] accesses all rows but only column 1, and same rule works for rows or certain columns


# Load data

# count_table, metadata
count_table <- read.delim("clean_genus_count_table_v0.1.0.tsv",
                          sep = "\t",
                          header = TRUE,
                          row.names = 1,
                          check.names = FALSE
                          )
metadata <- read.delim("saMBA_metadata_expanded_filtering_cleaned_with_country_project.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE
                       )
nrow(count_table) # number of genera identified in count table
ncol(count_table) # number of samples in count table

nrow(metadata) # number of genera identified in count table
ncol(metadata) # 413 samples in count table

metadata <- metadata[,-1] # gets rid of 1st column
metadata_samples <- metadata$run_accession # metadata_samples includes names of all the samples

# Gets rid of potentially ambiguous values
metadata <- metadata %>%
  filter(project_name != "Characterize the gut microbiome of indigenous and 
         urban individuals from the Brazilian Amazon, a region of enormous cultural, 
         ethnic, and biological diversity.")# %>%
  #filter(!grepl("urban biome", broad_scale_environmental_context, ignore.case = TRUE))

# colnames(count_table) converts the names of the columns in count_table into a list-like structure
common_samples <- intersect(colnames(count_table), metadata_samples) # 407 in common

length(common_samples)

# filter count table so it only includes columns in common w/ metadata samples
count_table_filtered <- count_table[, colnames(count_table) %in% common_samples] 

# filter count table so it only includes columns in common w/ count_table samples
metadata_filtered <- metadata[metadata$run_accession %in% common_samples, ]

# find and plot the sampling depth of each sample
## how much they sequenced

taxa_sums <- rowSums(count_table_filtered)
count_table_clean <- count_table_filtered[taxa_sums > 0, ] 
# filtered taxa that were no longer present

sample_depths <- colSums(count_table_clean) # depth of 10026 to 2198595

plot_data <- data.frame(sample = names(sample_depths), depth = sample_depths)

#ggplot(plot_data,
  #     aes(x = depth)) +
  #     geom_histogram(bins = 25, fill = "green")

# saves count_table file with raw counts
count_table_raw <- count_table_clean

# Important - converts count_table to relative abundances
count_table_clean <- sweep(count_table_clean, 2, colSums(count_table_clean), FUN = "/")  

count_table_t <- t(count_table_clean)

bray_dist <- vegdist(count_table_t, method = "bray")

pcoa_result <- cmdscale(bray_dist, eig = TRUE, k  = min(nrow(count_table_t)-1), 10)

pcoa_result

pcoa_df <- data.frame(
  sample_accession = rownames(pcoa_result$points),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2]
)

colnames(pcoa_df)[1] <- "run_accession"

pcoa_meta <- merge(pcoa_df, metadata_filtered, by = "run_accession")


ggplot(pcoa_meta,
       aes(x = PC1, y = PC2, color = broad_scale_environmental_context, shape = instrument_model)) +
  geom_point(size = 3) + 
  labs(
    x = "PC1",
    y = "PC2"
  ) +
theme_minimal()

batch <- metadata_filtered$instrument_model
mod <- model.matrix(~ broad_scale_environmental_context, data = metadata_filtered)
combat_data <- ComBat(dat = as.matrix(count_table_clean), batch = batch, mod = mod)

bray_dist_combat <- vegdist(t(combat_data), method = "bray")
pcoa_result_combat <- cmdscale(bray_dist_combat, eig = TRUE, k = 2)

pcoa_df_combat <- data.frame(
  run_accession = rownames(pcoa_result_combat$points),
  PC1 = pcoa_result_combat$points[, 1],
  PC2 = pcoa_result_combat$points[, 2]
)

pcoa_meta_combat <- merge(pcoa_df_combat, metadata_filtered, by = "run_accession")

ggplot(pcoa_meta_combat,
       aes(x = PC1, y = PC2, color = broad_scale_environmental_context, shape = instrument_model)) +
  geom_point(size = 3) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal()

write.table(count_table_clean, "filtered_genus_count_table.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(metadata_filtered, "filtered_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(pcoa_meta, "pcoa_metadata_merged.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


