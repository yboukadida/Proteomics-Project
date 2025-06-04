library(corrplot)

# keep sample names
sample_names <- rownames(pData(qnt))

# create condition column
conditions <- factor(c(rep("CTRL", 3), rep("TRT", 3)))

# replace pData(qnt) with a df with conditions and keep sample names
pData(qnt) <- data.frame(condition = conditions, row.names = sample_names)

# Verification
print(pData(qnt))






# Extraction des échantillons par condition
ctrl_samples <- colnames(exprs(qnt))[pData(qnt)$condition == "CTRL"]
trt_samples <- colnames(exprs(qnt))[pData(qnt)$condition == "TRT"]

# Calcul des matrices de corrélation
cor_ctrl <- cor(exprs(qnt)[, ctrl_samples], use = "pairwise.complete.obs")
cor_trt <- cor(exprs(qnt)[, trt_samples], use = "pairwise.complete.obs")

# Affichage des corrélations par heatmap
library(pheatmap)

pheatmap(cor_ctrl,
         main = "Corrélation réplicats CTRL",
         fontsize = 10)

pheatmap(cor_trt,
         main = "Corrélation réplicats TRT",
         fontsize = 10)

