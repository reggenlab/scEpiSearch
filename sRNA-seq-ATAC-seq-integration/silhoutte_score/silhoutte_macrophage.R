library(corrplot)
# Clustering
library(cluster) 
library(factoextra)
library(dplyr)

data_seurat <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/seurat-embed/mouse/embeds-MCA-cell-Mouse-macrophages.txt', header = TRUE,sep="\t")
head(data_seurat,3)
drops <- c("X")
data_seurat = data_seurat[ , !(names(data_seurat) %in% drops)]
colnames(data_seurat) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_seurat, grepl('Macro|macrophage_mouse', X))
data2 = dplyr::filter(data_seurat, !grepl('Macro|macrophage_mouse', X))
data_seurat_l = rbind(data1,data2)
labels_seurat = data_seurat_l$X
data1$X = "bcell"
data_seurat = rbind(data1,data2)
colnames(data_seurat) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_seurat,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_seurat$X = d1$X
head(data_seurat,5)

data_conos <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/macrophage_mouse_conos_tsne.txt.csv',header=TRUE,sep=",")
head(data_conos,3)
data_conos <- data_conos[-c(1,2)]
colnames(data_conos) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_conos, grepl('Macro|macrophage_mouse', X))
data2 = dplyr::filter(data_conos, !grepl('Macro|macrophage_mouse', X))
data_conos_l = rbind(data1,data2)
labels_conos = data_conos_l$X
data1$X = "bcell"
data_conos = rbind(data1,data2)
#data_conos = data1
colnames(data_conos) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_conos,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_conos$X = d1$X
head(data_conos,5)

data_liger <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/liger_macrophae_mouse_tsne.txt', header = TRUE,sep="\t")
head(data_liger,3)
#data_liger <- data_liger[-c(1)]
colnames(data_liger) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_liger, grepl('Macro|macrophage_mouse', X))
data2 = dplyr::filter(data_liger, !grepl('Macro|macrophage_mouse', X))
data_liger_l = rbind(data1,data2)
labels_liger = data_liger_l$X
data1$X = "bcell"
data_liger = rbind(data1,data2)
#data_liger = data1
colnames(data_liger) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_liger,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_liger$X = d1$X
head(data_liger,5)

data_scepi <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/macrophage_liver_scepi_tsne.txt', header = TRUE,sep=" ")
data1 = dplyr::filter(data_scepi, grepl('Macro|macrophage_mouse', X))
data2 = dplyr::filter(data_scepi, !grepl('Macro|macrophage_mouse', X))
data_scepi_l = rbind(data1,data2)
labels_scepi = data_scepi_l$X
data1$X = "bcell"
data_scepi = rbind(data1,data2)
#data_scepi = data1
colnames(data_scepi) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scepi,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scepi$X = d1$X
head(data_scepi,5)

df_seurat <- data_seurat[-c(3)]
df_conos <- data_conos[-c(3)]
df_liger <- data_liger[-c(3)]
df_scepi <- data_scepi[-c(3)]

silhouette_score_seurat <- function(k){
  #km <- kmeans(df_seurat, centers = k)
  ss <- silhouette(data_seurat$X, dist(df_seurat))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_conos <- function(k){
  #km <- kmeans(df_conos, centers = k)
  ss <- silhouette(data_conos$X, dist(df_conos))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_liger <- function(k){
  #km <- kmeans(df_liger, centers = k)
  ss <- silhouette(data_liger$X, dist(df_liger))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_scepi <- function(k){
  #km <- kmeans(df_scepi, centers = k)
  ss <- silhouette(data_scepi$X, dist(df_scepi))
  #mean(ss[, 3])
  ss[,3]
}
k=5
avg_sil_seurat <- sapply(k, silhouette_score_seurat)
avg_sil_conos <- sapply(k, silhouette_score_conos)
avg_sil_liger <- sapply(k, silhouette_score_liger)
avg_sil_scepi <- sapply(k, silhouette_score_scepi)
#plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

avg_sil_scepi = as.data.frame(avg_sil_scepi)
avg_sil_scepi$X = labels_scepi
sil_scepi = dplyr::filter(avg_sil_scepi, grepl('Macrophage_Mouse', X))

avg_sil_seurat = as.data.frame(avg_sil_seurat)
avg_sil_seurat$X = labels_seurat
sil_seurat = dplyr::filter(avg_sil_seurat, grepl('Macrophage_Mouse', X))

avg_sil_conos = as.data.frame(avg_sil_conos)
avg_sil_conos$X = labels_conos
sil_conos = dplyr::filter(avg_sil_conos, grepl('macrophage_mouse', X))

avg_sil_liger = as.data.frame(avg_sil_liger)
avg_sil_liger$X = labels_liger
sil_liger = dplyr::filter(avg_sil_liger, grepl('Macrophage_mouse', X))

boxplot(c(sil_seurat$V1),c(sil_scepi$V1),c(sil_liger$V1),c(sil_conos$V1),names = c("Seurat", "ScEpiSearch", "Liger", "Conos"),
        col = c("orange","red",'green','blue'),
        main = "Macrophage-Mouse",
        border = "brown",
        notch = FALSE)
