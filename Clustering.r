library(ComplexHeatmap)
data = as.matrix(read.table("GeneExpData.tsv", sep="\t", header = T, row.names = 1))
dist = dist(data)
hclust(dist)
hclust = hclust(dist)
Heatmap(data, show_row_names = FALSE, width = ncol(data)*unit(20, "mm"))
