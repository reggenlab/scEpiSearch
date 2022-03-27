library(rliger)
data.rna <- read.csv("/storage/vibhor/neeteshp/shreyas_tools_testing/MCA_reference.csv")
data.atac <- read.csv("/storage/vibhor/neeteshp/shreyas_tools_testing/liger_data_query_/liver-endothelial_liger.txt",sep="\t",header=F)
rownames(data.rna) = data.rna[1:nrow(data.rna),1]
data.rna = data.rna[1:nrow(data.rna),2:ncol(data.rna)]
row.names(data.rna) = toupper(row.names(data.rna))

data.atac = data.atac[!duplicated(data.atac$V1),] 
data.atac = data.atac[!(is.na(data.atac$V1) | data.atac$V1==""), ]
rownames(data.atac) = data.atac[1:nrow(data.atac),1]
data.atac = data.atac[1:nrow(data.atac),2:ncol(data.atac)]
row.names(data.atac) = toupper(row.names(data.atac))

endo = data.atac[,1:318]
macro = data.atac[,319:413]
bj = data.atac[,414:437]
bmmc.data <- list(atac_1 = endo,atac_2=macro,atac_3=bj, rna = pbmc.rna)

bmmc.data <- list(atac = data.atac, rna = pbmc.rna)
int.bmmc <- createLiger(bmmc.data)
int.bmmc <- normalize(int.bmmc)
int.bmmc <- selectGenes(int.bmmc, datasets.use = 2)
int.bmmc <- scaleNotCenter(int.bmmc, remove.missing = FALSE)
int.bmmc <- optimizeALS(int.bmmc, k = 20)
int.bmmc <- quantile_norm(int.bmmc,knn_k=5)
int.bmmc <- runUMAP(int.bmmc, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
p = plotByDatasetAndCluster(int.bmmc, axis.labels = c('UMAP 1', 'UMAP 2'),return.plots = TRUE)
#Get tsne coordinates from Liger :  
tsne_df <- data.frame(int.bmmc@tsne.coords)
colnames(tsne_df) <- c("Dim1", "Dim2")
tsne_df[['Dataset']] <- unlist(lapply(1:length(int.bmmc@H), function(x) {
  rep(names(int.bmmc@H)[x], nrow(int.bmmc@H[[x]]))
}))
