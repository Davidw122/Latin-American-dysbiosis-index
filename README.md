# Latin American Dysbiosis Index: A Machine Learning & Cultural Insight Project
Welcome to the repository for the Latin American Dysbiosis Index (LADI) project — a human-centered research initiative that combines computational biology, machine learning, and cultural research to better understand the impact of industrialization and westernization on gut microbiota across Latin America.

### Project Summary

This project investigates how shifts in lifestyle and diet influence gut microbiome composition by leveraging:
- Publicly available microbiome sequencing data (saMBA)
- Random forest and nearest centroid classifiers (NCC)
- Dimensionality reduction and diversity metrics (e.g., Bray-Curtis)
- A novel Industrialization Index for visualization and diagnostic insights

The goal is to lay the foundation for a culturally grounded Latin American dysbiosis diagnostic tool, informed by machine learning, microbiome data, and public health perspectives.

### Motivation

Inspired by personal experiences and cultural heritage, this project explores how urbanization and dietary shifts in Latin American communities may be altering microbial health. By investigating these changes quantitatively, I aim to empower culturally specific health tools that respect traditional food practices and biodiversity.

### Methodology

1. Data Preprocessing
- Filtered raw genus-level count tables from public microbiome datasets
- Removed low-prevalence taxa
- Applied log-transformation and pseudocount smoothing

2. Industrialization Labeling
- Samples were categorized into "More Industrialized" and "Less Industrialized" based on the metadata's broad_scale_environmental_context.

4. Random Forest Classification
- Trained a 1000-tree model to distinguish industrialization classes
- Achieved ~79.5% accuracy on test set
- Identified top microbial genera contributing to classification

4. Nearest Centroid Classification (NCC)
- Computed class centroids in microbiome space
- Calculated Euclidean distances for each sample to centroids
- Potential to create a continuous Industrialization Index ranging from 0 (rural-like) to 1 (urban-like)

5. Statistical Analysis
- PERMANOVA & ANOVA to assess community differences across biomes
- Tukey's HSD for post-hoc comparisons

### Key Findings
- Significant microbiome differences were detected across environmental contexts (p < 0.001).
- The NCC-based index offers an interpretable and intuitive visualization of "microbial urbanization".
- Top genera contributing to industrialization-related differences include Prevotella, Bacteroides, and Ruminococcus species.

### Next Steps
Refine index to account for specific regional dietary practices
Interview Latin American individuals about traditional food and health practices
Build a demo web interface for the Dysbiosis Index

### About Me
I’m a high school researcher passionate about blending STEM with cultural preservation. I built this project using R, public datasets, and a design-thinking approach to computational biology.


