#data <- read.csv('/Desktop/codex1.csv', row.names = 1) #Moneera directory
data <- read.csv('data/codex1.csv.gz', row.names=1) #Kalen directory

# - If you remove cells, why did you remove them? 

gexp <- data[,4:ncol(data)]
totalgexp <- rowSums(gexp)
good.cells <- names(totalgexp)[totalgexp > 0]
gexp <- gexp[good.cells,]
pos <- data[,1:2]


# variance of each protein? To know which one to be normalized 
variance <- apply(gexp, 2, var)

# Sort proteins by variance in descending order
sorted_variances <- sort(variance, decreasing = TRUE)
sorted_variances

# Choose the top proteins with highest variance
high_var_proteins <- names(sorted_variances[1:ncol(gexp)])

# Apply normalization to the selected proteins
norm_gexp <- gexp
for (protein in high_var_proteins) {
  norm_gexp[, protein] <- norm_gexp[, protein] / sum(norm_gexp[, protein]) * median(rowSums(norm_gexp))
}

norm_gexp[, -which(colnames(norm_gexp) %in% high_var_proteins)] <- 0
norm_gexp <- log10(norm_gexp + 1)

# data was normalized based on protein variance

# Reduced dimensional representation of cell clusters.
library(Rtsne)
library(ggplot2)


## PCA analysis
pcs <- prcomp(norm_gexp)

## tSNE analysis
set.seed(0)
emb <- Rtsne(norm_gexp)


# - Kmeans clustering, justify the choice of K

library(segmented)
ks <- 1:28
set.seed(0)
out <- do.call(rbind, lapply(ks, function(k) {
  com <- kmeans(norm_gexp, centers = k)
  c(within = com$tot.withinss, between = com$betweenss)
}))
plot(ks,out[,1],type='b',main='WSS')
plot(ks,out[,2],type='b',main='BSS')

fit <- segmented(lm(out[,1] ~ ks))

breakpoint <- fit$psi[1]

plot(ks, out[,1], type='b', main='WSS')
plot(ks,out[,2],type='b',main='BSS')
abline(fit, col='red')

abline(v=breakpoint, col='blue')

# in this data, the choice of k is 5, because that's the value where the graph seems to plateau 
# this code works by doing a loop through ks (1 -28) to return wss and bss 
# as matrix, then fits a segmented regression line to the plot, then plots the WSS and  segmented regression 
# line, and adds a vertical line at the estimated elbow point

#  Tissue representation of cell clusters

library(gridExtra)
num_center = 5
set.seed(0)
com <- kmeans(norm_gexp, centers=num_center)
# plot tSNE analysis results
df <- data.frame(pos, emb$Y, pcs$x, celltype = as.factor(com$cluster))
head(df)
p1 <- ggplot(df, aes(x = X1, y = X2, col=celltype)) +
  geom_point(size = 0.01) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")) +
  theme_classic()
p2 <- ggplot(df, aes(x = PC1, y = PC2, col=celltype)) +
  geom_point(size = 0.01) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")) +
  theme_classic()
grid.arrange(p1,p2)

p3 <- ggplot(df, aes(x=x, y=y, col=celltype)) +
  geom_point(size=0.01)

p3



# Summary of differentially expressed genes for cell-types of interest.
library(ggrepel)
## pick a cluster
cluster_number = 4
cluster.of.interest <- names(which(com$cluster == cluster_number))
cluster.other <- names(which(com$cluster != cluster_number))
genes <- colnames(norm_gexp)

## loop through each cluster to find differentially expressed genes for each cluster vs. remaining clusters
volcano_plts <- lapply(seq_len(num_center), function(i) {
  ## pick a cluster
  cluster.of.interest <- names(which(com$cluster == i))
  cluster.other <- names(which(com$cluster != i))
  

## loop through my genes and test each one
out <- sapply(genes, function(g) {
  a <- norm_gexp[cluster.of.interest, g]
  b <- norm_gexp[cluster.other, g]
  pvs <- wilcox.test(a,b,alternative='two.sided')$p.val
  #pvs <- as.data.frame(pvs)
  log2fc <- log2(mean(a)/mean(b))
  #log2fc <- as.data.frame(log2fc)
  c(pvs=pvs,log2fc=log2fc)
})

## volcano plot
# run Bonferroni correction determine a p value cutoff (conservative method with low false-positive rate but high false-negative rate as a trade-off)
cutoff.pvs <- 0.05/ncol(norm_gexp) # false-positive rate would be 1-(1-0.05/ncol(norm_gexp))^ncol(norm_gexp)
cutoff.log2fc <- 1
df_temp <- data.frame(pvs=out[1,], log2fc=out[2,])
## prepare a data frame for plotting
df2 <- df_temp
# add a column of NAs to df2 titled diffexpressed
df2$diffexpressed <- 'NO'
# if log2fc > cutoff.log2fc and pvalue < cutoff.pvs, set as "Up" 
df2$diffexpressed[df2$log2fc > cutoff.log2fc & df2$pvs < cutoff.pvs] <- "Up"
# if log2fc < -cutoff.log2fc and pvalue < cutoff.pvs, set as "Down"
df2$diffexpressed[df2$log2fc < -cutoff.log2fc & df2$pvs < cutoff.pvs] <- "Down"
# add a column of NAs to df2 titled genelabel that will contain names of differentially expressed genes
df2$genelabel <- NA
df2$genelabel[df2$diffexpressed != 'NO'] <- rownames(df2)[df2$diffexpressed != 'NO']

## get the most significant genes (upregulated) from the chosen cluster
significant_genes <- df2$gene[df2$diffexpressed == 'Up']
#sort by p-value
significant_genes_sorted <- significant_genes[order(out[1,][significant_genes])]
#choose first one from list
most_significant_gene <- significant_genes_sorted[1]

volcano_plt <- ggplot(df2, aes(x=log2fc,y=-log10(pvs), col=diffexpressed, label=genelabel)) +
  geom_point() +
  scale_color_manual(values = c('blue', 'black', 'red'),
                     name = 'Expression',
                     breaks = c('Down', 'NO', 'Up'),
                     labels = c('Down', 'N.S.', 'Up')) +
  ggtitle(paste('Cluster ',i,' vs. Other Clusters')) +
  xlab('log2 fold change') +
  ylab('-log10(p value)') +
  geom_label_repel() +
  geom_vline(xintercept = c(-cutoff.log2fc, cutoff.log2fc), col='black',linetype='dashed') +
  geom_hline(yintercept = -log10(cutoff.pvs), col='black',linetype='dashed') +
  theme_classic()

list(volcano_plt, significant_genes, most_significant_gene)

})

grid.arrange(volcano_plts[[1]][[1]],volcano_plts[[2]][[1]],volcano_plts[[3]][[1]], volcano_plts[[4]][[1]], volcano_plts[[5]][[1]])

# For the cell type of interest, cluster 4 was chosen, using a log2fc cutoff value of 1,
# the upregulated genes are CD35, SMActin, HLDA.DR, CD20, ECAD and CD21, while the downregulted genes are:
# CD8, Lyve1, CD3e. Using a function to sort in order of ascending p value for upregulated genes, the most significant gene is CD20

