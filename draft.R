library(Rtsne)
library(ggplot2)
library(segmented)
library(gridExtra)
library(ggrepel)


data <- read.csv('/Users/moneeraalsuwailm/Desktop/codex1.csv', row.names = 1)
data <- read.csv('data/codex1.csv.gz', row.names=1) #Kalen directory


gexp <- data[,4:ncol(data)]
totalgexp <- rowSums(gexp)
good.cells <- names(totalgexp)[totalgexp > 0]
gexp <- gexp[good.cells,]
pos <- data[,1:2]

# cells are removed if they have no expression (zero counts) across
# all selected columns (genes), this code is calculating the sum of each row
# and returns result vector, then it select the names of the rows where
# the sum is greater than 0 (0 expression) and assigns them to good.cells
# then subsets the df to only include the rows of good.cells.


# variance of each protein? To know which one to be normalized 
variance <- apply(gexp, 2, var)
hist(variance)

df_var <- data.frame(variance)
p <- ggplot(data = df_var, mapping = aes(x=log10(variance))) +
  geom_histogram(breaks = seq(0,8), 
                 binwidth = 1, 
                 boundary=0, closed="right", 
                 color="white", fill="black") +
  scale_x_continuous(n.breaks=8) +
  scale_y_continuous(n.breaks=10) +
  labs(title="Variance before Normalization", 
       x="log10(variance)", y = "frequency") +
  theme_bw()

p

# Sort proteins by variance in descending order
sorted_variances <- sort(variance, decreasing = TRUE)
sorted_variances

# Choose the top proteins with highest variance
threshold <- 1.0e+06
high_var_proteins <- names(sorted_variances[sorted_variances > threshold])

# Subset the expression matrix to only include high variance proteins
gexp_high_var <- gexp[, high_var_proteins]

# Apply normalization to the selected proteins
norm_high_gexp <- gexp_high_var
for (protein in high_var_proteins) {
  norm_high_gexp[, protein] <- norm_high_gexp[, protein] / sum(norm_high_gexp[, protein]) * median(rowSums(norm_high_gexp))
}

norm_high_gexp[, -which(colnames(norm_high_gexp) %in% high_var_proteins)] <- 0
norm_high_gexp <- log10(norm_high_gexp + 1)

# variance after normalization
norm_high_variance <- apply(norm_high_gexp, 2, var)
hist(norm_high_variance)

df_norm_high_var <- data.frame(norm_high_variance)
p_norm_high <- ggplot(data = df_norm_high_var , mapping = aes(x=log10(norm_high_variance))) +
  geom_histogram(breaks = seq(-5,-1), 
                 #binwidth = 1, 
                 boundary=0, closed="right", 
                 color="white", fill="black") +
  #scale_x_continuous(n.breaks=8) +
  #scale_y_continuous(n.breaks=10) +
  labs(title="Variance after Normalization of Only Highest Variable Proteins", 
       x="log10(variance)", y = "frequency") +
  theme_bw()


p_norm_high

#variaince of protein in the y axess insted of expression 
# Apply normalization to all proteins
norm_gexp <- gexp
for (protein in colnames(norm_gexp)) {
  norm_gexp[, protein] <- norm_gexp[, protein] / sum(norm_gexp[, protein]) * median(rowSums(norm_gexp))
}

norm_gexp <- log10(norm_gexp + 1)

# variance after normalization
norm_variance <- apply(norm_gexp, 2, var)
hist(norm_variance)

df_norm_var <- data.frame(norm_variance)
p_norm <- ggplot(data = df_norm_var , mapping = aes(x=log10(norm_variance))) +
  geom_histogram(breaks = seq(-5,-1), 
                 binwidth = 1, 
                 boundary=0, closed="right", 
                 color="white", fill="black") +
  #scale_x_continuous(n.breaks=8) +
  #scale_y_continuous(n.breaks=10) +
  labs(title="Variance after Normalization of All Proteins", 
       x="log10(variance)", y = "frequency") +
  theme_bw()

p_norm

grid.arrange(p, p_norm_high, p_norm, ncol=3)

## PCA analysis
pca <- prcomp(gexp)

#get number of PCAs
cumulative_var <- cumsum(pca$sdev^2/sum(pca$sdev^2))

plot(cumulative_var, xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained", type = "b")

desired_var <- 0.85 # for example, let's say we want to explain 85% of the variance
num_pcs <- which(cumulative_var > desired_var)[1]

num_pcs

df <- gexp
df$PC1 <- pca$x[, 1]
df$PC2 <- pca$x[, 2]
df$PC3 <- pca$x[, 3]

#Plot first two principal components for CD15
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = CD15)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD15", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD15") +
  xlab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for Vimentin
p2 <- ggplot(df, aes(x = PC2, y = PC3, color = Vimentin)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "Vimentin", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for Vimentin") +
  xlab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC3 (", round(pca$sdev[3]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for CD8
p3 <- ggplot(df, aes(x = PC3, y = PC2, color = CD8)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD8", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD8") +
  xlab(paste0("PC3 (", round(pca$sdev[3]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Combine the two plots into one image
grid.arrange(p1, p2, p3, ncol = 3)



## PCA analysis
pca <- prcomp(norm_gexp)

#get number of PCAs
cumulative_var <- cumsum(pca$sdev^2/sum(pca$sdev^2))

plot(cumulative_var, xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained", type = "b")

desired_var <- 0.85 # for example, let's say we want to explain 85% of the variance
num_pcs <- which(cumulative_var > desired_var)[1]

num_pcs

df <- norm_gexp
df$PC1 <- pca$x[, 1]
df$PC2 <- pca$x[, 2]

#Plot first two principal components for CD15
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = CD15)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD15", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD15") +
  xlab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for Vimentin
p2 <- ggplot(df, aes(x = PC1, y = PC2, color = Vimentin)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "Vimentin", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for Vimentin") +
  xlab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for CD8
p3 <- ggplot(df, aes(x = PC1, y = PC2, color = CD8)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD8", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD8") +
  xlab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Combine the two plots into one image
grid.arrange(p1, p2, p3, ncol = 3)

## PCA analysis
pca <- prcomp(norm_high_gexp)

#get number of PCAs
cumulative_var <- cumsum(pca$sdev^2/sum(pca$sdev^2))

plot(cumulative_var, xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance Explained", type = "b")

desired_var <- 0.85 # for example, let's say we want to explain 85% of the variance
num_pcs <- which(cumulative_var > desired_var)[1]

num_pcs

df <- norm_high_gexp
df$PC1 <- pca$x[, 1]
df$PC2 <- pca$x[, 2]
df$PC3 <- pca$x[, 3]

#Plot first two principal components for CD15
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = CD15)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD15", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD15") +
  xlab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for Vimentin
p2 <- ggplot(df, aes(x = PC3, y = PC2, color = Vimentin)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "Vimentin", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for Vimentin") +
  xlab(paste0("PC3 (", round(pca$sdev[3]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Plot first two principal components for CD8
p3 <- ggplot(df, aes(x = PC2, y = PC1, color = CD8)) +
  geom_point(size = 3) +
  scale_color_gradient(name = "CD8", low = "green", high = "red") +
  ggtitle("PCA of Spatial Transcriptomics Data for CD8") +
  xlab(paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "% variance)")) +
  ylab(paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "% variance)"))

#Combine the two plots into one image
grid.arrange(p1, p2, p3, ncol = 3)


# Combine the normalized data for all proteins and high variance proteins
combined_data <- rbind(data.frame(protein = rep("All Proteins", ncol(norm_gexp)), 
                                  value = as.vector(t(norm_gexp))), 
                       data.frame(protein = rep("High Variance Proteins", ncol(norm_high_gexp)), 
                                  value = as.vector(t(norm_high_gexp))))

# Create a box plot of the normalized expression values for all proteins and high variance proteins
ggplot(combined_data, aes(x = protein, y = value, fill = protein)) +
  geom_boxplot() +
  scale_fill_manual(values = c("All Proteins" = "grey", "High Variance Proteins" = "red")) +
  labs(title = "Comparison of Normalized Expression of All Proteins and High Variance Proteins",
       y = "Normalized Expression (log10 scale)",
       x = "") +
  theme_bw()
#histogram insted of box 

ggplot(combined_data, aes(x = value, fill = protein)) +
  geom_histogram(binwidth = 0.2) +
  scale_fill_manual(values = c("All Proteins" = "grey", "High Variance Proteins" = "red")) +
  labs(title = "Comparison of Normalized Expression of All Proteins and High Variance Proteins",
       y = "Count",
       x = "Normalized Expression (log10 scale)") +
  theme_bw()
len

hist(log10(data$CD45RO+1))
hist(log10(data$Podoplanin+1))


# Reduced dimensional representation of cell clusters.


get_kmeans_clusters <- function(gexp, n_proteins){
  
  # Kmeans clustering, justify the choice of K
  ks <- 1:n_proteins
  out <- do.call(rbind, lapply(ks, function(k) {
    com <- kmeans(gexp, centers = k)
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
  
  
}

plot_tnse_analysis_results <- function(pos, emb, pcs, com){
  
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
  
  print(grid.arrange(p1,p2))
  
  p3 <- ggplot(df, aes(x=x, y=y, col=celltype)) +
    geom_point(size=0.01)
  
  print(p3)
  
}

get_significant_genes <- function(cluster_number, com, gexp){
  
  # Summary of differentially expressed genes for cell-types of interest.
  cluster.of.interest <- names(which(com$cluster == cluster_number))
  cluster.other <- names(which(com$cluster != cluster_number))
  genes <- colnames(gexp)
  
  ## loop through my genes and test each one
  out <- sapply(genes, function(g) {
    a <- gexp[cluster.of.interest, g]
    b <- gexp[cluster.other, g]
    pvs <- wilcox.test(a,b,alternative='two.sided')$p.val
    #pvs <- as.data.frame(pvs)
    log2fc <- log2(mean(a)/mean(b))
    #log2fc <- as.data.frame(log2fc)
    c(pvs=pvs,log2fc=log2fc)
  })
  
  ## volcano plot
  # run Bonferroni correction determine a p value cutoff (conservative method with low false-positive rate but high false-negative rate as a trade-off)
  cutoff.pvs <- 0.05/ncol(gexp) # false-positive rate would be 1-(1-0.05/ncol(norm_gexp))^ncol(norm_gexp)
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
  
  volcano_plt <- ggplot(df2, aes(x=log2fc,y=-log10(pvs), col=diffexpressed, label=genelabel)) +
    geom_point() +
    scale_color_manual(values = c('blue', 'black', 'red'),
                       name = 'Expression',
                       breaks = c('Down', 'NO', 'Up'),
                       labels = c('Down', 'N.S.', 'Up')) +
    ggtitle(paste('Cluster ',cluster_number,' vs. Other Clusters')) +
    xlab('log2 fold change') +
    ylab('-log10(p value)') +
    geom_label_repel() +
    geom_vline(xintercept = c(-cutoff.log2fc, cutoff.log2fc), col='black',linetype='dashed') +
    geom_hline(yintercept = -log10(cutoff.pvs), col='black',linetype='dashed') +
    theme_classic()
  
  print(volcano_plt)
  
  return(significant_genes_sorted)
  
}


######################################################
# Using all proteins normalized
######################################################

## PCA analysis
pcs <- prcomp(norm_gexp)

## tSNE analysis
set.seed(1)
emb <- Rtsne(norm_gexp)

#get number of clusters
set.seed(1)
get_kmeans_clusters(norm_gexp, ncol(norm_gexp))

#number of clusters
set.seed(1)
num_center = 9
com <- kmeans(norm_gexp, centers=num_center)

#plot tsne results
plot_tnse_analysis_results(pos, emb, pcs, com)


#pick a cluster
cluster_number = 3

#get significant genes
significant_genes <- get_significant_genes(cluster_number, com, norm_gexp)
most_significant_gene <- significant_genes[1]

result1 <- paste("Cluster chosen = ", paste(cluster_number),"\n\nCase 1: All proteins normalized\n\nSignificant genes = ", paste(significant_genes, collapse = ", "), "\nMost significant gene = ", paste(most_significant_gene))
cat(result1)

## loop through each cluster to find differentially expressed genes for each cluster vs. remaining clusters
volcano_plts <- lapply(seq_len(num_center), function(i) {
  ## pick a cluster
  cluster.of.interest <- names(which(com$cluster == i))
  cluster.other <- names(which(com$cluster != i))
  
  
  ## loop through my genes and test each one
  out <- sapply(colnames(norm_gexp), function(g) {
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

p <- grid.arrange(volcano_plts[[1]][[1]],
                  volcano_plts[[2]][[1]],
                  volcano_plts[[3]][[1]], 
                  volcano_plts[[4]][[1]], 
                  volcano_plts[[5]][[1]],
                  volcano_plts[[6]][[1]], 
                  volcano_plts[[7]][[1]], 
                  volcano_plts[[8]][[1]],
                  volcano_plts[[9]][[1]]
)

#### plot for a figure ####
## set custom theme for plots
plot_theme <- theme_classic() +
  theme(
    text = element_text(size = 25),
    legend.key.size = unit(0.75, 'cm')
  )

# set the green color palette
green_palette <- c("#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b")

# df <- data.frame(pos, emb$Y, pcs$x, celltype = as.factor(com$cluster))
# p3 <- ggplot(df, aes(x=x, y=y, col=celltype==8)) +
#   geom_point(size=0.01)
# p3
# p3 <- ggplot(df, aes(x=x, y=y, col=celltype==5)) +
#   geom_point(size=0.01)
# p3

cluster.of.interest <- 5
gene.of.interest <- c('CD21', 'CD20', 'CD35', 'HLADR')
df <- data.frame(pos, pcs$x[,1:2], emb$Y, celltype = as.factor(com$cluster), gene=norm_gexp[,gene.of.interest[1]])
p8 <- ggplot(df, aes(x = x, y = y, col = celltype == cluster.of.interest)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(green_palette[5], 'black'),
                     name = paste('Cluster', cluster.of.interest),
                     labels = c('FALSE', 'TRUE')) +
  ggtitle(paste0('Spatial (Cluster ', cluster.of.interest, ')')) +
  plot_theme

p8


######################################################
# Using only proteins with high variance normalized
######################################################

## PCA analysis
pcs_high <- prcomp(norm_high_gexp)

## tSNE analysis
set.seed(1)
emb_high <- Rtsne(norm_high_gexp)

#get number of clusters
set.seed(1)
get_kmeans_clusters(norm_high_gexp, ncol(norm_high_gexp))


#number of clusters
set.seed(1)
num_center_high = 5
com_high <- kmeans(norm_high_gexp, centers=num_center_high)


#plot tsne results
plot_tnse_analysis_results(pos, emb_high, pcs_high, com_high)

#pick a cluster
cluster_number = 1

#get significant genes
significant_genes <- get_significant_genes(cluster_number, com_high, norm_high_gexp)
most_significant_gene <- significant_genes[1]

result2 <- paste("Cluster chosen = ", paste(cluster_number),"\n\nCase 2: High variance proteins normalized\n\nSignificant genes = ", paste(significant_genes, collapse = ", "), "\nMost significant gene = ", paste(most_significant_gene))
cat(result2)

## loop through each cluster to find differentially expressed genes for each cluster vs. remaining clusters
volcano_plts <- lapply(seq_len(num_center_high), function(i) {
  ## pick a cluster
  cluster.of.interest <- names(which(com_high$cluster == i))
  cluster.other <- names(which(com_high$cluster != i))
  
  
  ## loop through my genes and test each one
  out <- sapply(colnames(norm_high_gexp), function(g) {
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
  cutoff.pvs <- 0.05/ncol(norm_high_gexp) # false-positive rate would be 1-(1-0.05/ncol(norm_gexp))^ncol(norm_gexp)
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

p <- grid.arrange(volcano_plts[[1]][[1]],
                  volcano_plts[[2]][[1]],
                  volcano_plts[[3]][[1]], 
                  volcano_plts[[4]][[1]], 
                  volcano_plts[[5]][[1]]
                  )


