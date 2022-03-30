#LIGER :
######################################### Load reference dataset ############################################
data.rna <- read.csv("MCA_reference.csv")
########################################Load query atac seq dataset ######################################
data.atac <- read.csv("h1esc_liger.tsv",sep="\t")
rownames(data.rna) = data.rna[1:nrow(data.rna),1]
data.rna = data.rna[1:nrow(data.rna),2:ncol(data.rna)]
row.names(data.rna) = toupper(row.names(data.rna))

data.atac = data.atac[!duplicated(data.atac$V1),] 
data.atac = data.atac[!(is.na(data.atac$V1) | data.atac$V1==""), ]
rownames(data.atac) = data.atac[1:nrow(data.atac),1]
data.atac = data.atac[1:nrow(data.atac),2:ncol(data.atac)]

bmmc.data <- list(atac = data.atac, rna = data.rna)
int.bmmc <- createLiger(bmmc.data)
int.bmmc <- normalize(int.bmmc)
int.bmmc <- selectGenes(int.bmmc, datasets.use = 2)
int.bmmc <- scaleNotCenter(int.bmmc, remove.missing = FALSE)
int.bmmc <- optimizeALS(int.bmmc, k = 20)
int.bmmc <- quantile_norm(int.bmmc,knn_k=5)
int.bmmc <- runUMAP(int.bmmc, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
p = plotByDatasetAndCluster(int.bmmc, axis.labels = c('UMAP 1', 'UMAP 2'),return.plots = TRUE))
#Get tsne coordinates from Liger :  
tsne_df <- data.frame(int.bmmc@tsne.coords)
  colnames(tsne_df) <- c("Dim1", "Dim2")
  tsne_df[['Dataset']] <- unlist(lapply(1:length(int.bmmc@H), function(x) {
    rep(names(int.bmmc@H)[x], nrow(int.bmmc@H[[x]]))
  }))







