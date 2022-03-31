library(corrplot)
library(RColorBrewer)
data <- read.table("CorrelationsData.tsv", header = T, sep = "\t", row.names = 1)
data_corr <- cor(data, method = c("pearson"))
res <- cor.mtest(data_corr, conf.level = .95)
png(file = "Correlations.png", height=1000, width=900, type="cairo")

corrplot(data_corr, type = "lower", order = "hclust", tl.col = "black", tl.cex =1.5, cl.cex=1.5, tl.srt = 25, col = brewer.pal(n = 8, name = "RdBu"), p.mat = res$p, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 4, pch.col = "black")

dev.off()

