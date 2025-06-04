# packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("rpx", "MSnbase", "Biobase", "pheatmap", "FactoMineR", "factoextra", "ggplot2", "dplyr"))

library(rpx)
library(MSnbase)
library(Biobase)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)

# data dowloading
library(rpx)
px1 <- PXDataset("PXD000001")  
mztab <- pxget(px1, "F063721.dat-mztab.txt")



# mztab file reading
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")
sampleNames(qnt) <- reporterNames(TMT6)
qnt <- filterNA(qnt)

# student test (non corrected)
group1 <- c(1, 2, 3)  #control 
group2 <- c(4, 5, 6)  # treated

t_test_results <- apply(exprs(qnt), 1, function(y) {
  t.test(y[group1], y[group2], var.equal = TRUE)$p.value
})

# pvalues adjumstment
t_test_results_adj <- p.adjust(t_test_results, method = "BH")

#Data frame & results creations
results_df <- data.frame(
  Protein = fData(qnt)$accession,
  p_value = t_test_results,
  adj_p_value = t_test_results_adj
)

# Filtring significant proteins
significant_proteins <- results_df %>%
  filter(p_value < 0.05)
head(significant_proteins)

#  log2 Fold Chang
log_fc <- function(data){
  vect_result <- rep(NA, length(rownames(exprs(qnt))))
  for (i in 1:length(rownames(exprs(qnt)))) {
    temp <- log2(mean(exprs(qnt)[i, group2])/mean(exprs(qnt)[i, group1]))
    vect_result[i] <- temp
  }
  return(vect_result)
}
logFC <- log_fc(exprs(qnt))
results_df$logFC <- logFC
results_df$logp_value <- -log10(results_df$p_value)

# Volcano plot
results_df$significance <- "Non significatif"
results_df$significance[results_df$logp_value > -log10(0.05) & abs(results_df$logFC) > 0.5] <- "Significatif"

volcano_plot <- ggplot(results_df, aes(x = logFC, y = logp_value, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot des Protéines", x = "Log2 Fold Change", y = "-Log10(p-value)") +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
print(volcano_plot)
ggsave("volcano_plot.png", plot = volcano_plot, width = 6, height = 5, dpi = 300)


#Heatmap
heatmap_data <- exprs(qnt)
pheatmap(heatmap_data,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "Heatmap des peptides/protéines")

# ACP analysis
pca_res <- PCA(t(exprs(qnt)), graph = FALSE)
fviz_pca_ind(pca_res, col.ind = "cos2", gradient.cols = c("blue", "red"))
ggsave("ACP_TMT.png", width = 6, height = 5, dpi = 300)
