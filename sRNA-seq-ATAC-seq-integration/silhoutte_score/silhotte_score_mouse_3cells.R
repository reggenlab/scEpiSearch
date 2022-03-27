library(corrplot)
# Clustering
library(cluster) 
library(factoextra)
library(dplyr)

data_seurat <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/seurat-embed/mouse/embeds-MCA-cell-Mouse-liver-3cells.txt', header = TRUE,sep="\t")
drops <- c("X")
data_seurat = data_seurat[ , !(names(data_seurat) %in% drops)]
colnames(data_seurat) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_seurat, grepl('Macro|Macrophage_Mouse', X))
data2 = dplyr::filter(data_seurat, grepl('B cell|Bcell_Mouse', X))
data3 = dplyr::filter(data_seurat, grepl('Endo|Endothelial_Mouse', X))
data4 = dplyr::filter(data_seurat, !grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))
data_seurat_l = rbind(data1,data2,data3,data4)
labels_seurat = data_seurat_l$X
data1$X = "macrophage"
data2$X = "bcell"
data3$X = "endothelial"
data_seurat = rbind(data1,data2,data3,data4)
colnames(data_seurat) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_seurat,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_seurat$X = d1$X
head(data_seurat,5)

data_conos <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/3cells_mouse_conos_tsne.txt.csv',sep=",")
data_conos <- data_conos[-c(1,5,6,7,8,9,10,11)]
colnames(data_conos) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_conos, grepl('Macro|Macrophage_mouse', X))
data2 = dplyr::filter(data_conos, grepl('B cell|Bcell_mouse', X))
data3 = dplyr::filter(data_conos, grepl('Endo|Endothelial_mouse', X))
data4 = dplyr::filter(data_conos, !grepl('Macro|Macrophage_mouse|B cell|Bcell_mouse|Endo|Endothelial_mouse', X))
data_conos_l = rbind(data1,data2,data3,data4)
labels_conos = data_conos_l$X
data1$X = "macrophage"
data2$X = "bcell"
data3$X = "endothelial"
data_conos = rbind(data1,data2,data3,data4)
#data_conos = data1
colnames(data_conos) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_conos,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_conos$X = d1$X
head(data_conos,5)

data_liger <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/liger_3cell_mouse_tsne.txt', header = TRUE,sep="\t")
#data_liger <- data_liger[-c(1)]
colnames(data_liger) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_liger, grepl('Macro|Macrophage_mouse', X))
data2 = dplyr::filter(data_liger, grepl('B cell|Bcell_mouse', X))
data3 = dplyr::filter(data_liger, grepl('Endo|endothelial_mouse', X))
data4 = dplyr::filter(data_liger, !grepl('Macro|Macrophage_mouse|B cell|Bcell_mouse|Endo|Endothelial_mouse', X))
data_liger_l = rbind(data1,data2,data3,data4)
labels_liger = data_liger_l$X
data1$X = "macrophage"
data2$X = "bcell"
data3$X = "endothelial"
data_liger = rbind(data1,data2,data3,data4)
#data_liger = data1
colnames(data_liger) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_liger,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_liger$X = d1$X
head(data_liger,5)

data_scepi <- read.csv('/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/3cells_liver_scepi_tsne.txt', header = TRUE,sep=" ")
data1 = dplyr::filter(data_scepi, grepl('Macro|Macrophage_Mouse', X))
data2 = dplyr::filter(data_scepi, grepl('B cell|Bcell_Mouse', X))
data3 = dplyr::filter(data_scepi, grepl('Endo|Endothelial_Mouse', X))
data4 = dplyr::filter(data_scepi, !grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))
data_scepi_l = rbind(data1,data2,data3,data4)
labels_scepi = data_scepi_l$X
data1$X = "macrophage"
data2$X = "bcell"
data3$X = "endothelial"
data_scepi = rbind(data1,data2,data3,data4)
#data_scepi = data1
colnames(data_scepi) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scepi,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scepi$X = d1$X
head(data_scepi,5)

df_seurat <- data_seurat[-c(3)]
df_conos <- data_conos[-c(3)]
df_liger <- data_liger[-c(3)]
df_scepi <- data_scepi[-c(3)]
k=5
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
avg_sil_seurat <- sapply(k, silhouette_score_seurat)
avg_sil_conos <- sapply(k, silhouette_score_conos)
avg_sil_liger <- sapply(k, silhouette_score_liger)
avg_sil_scepi <- sapply(k, silhouette_score_scepi)
#plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

avg_sil_scepi = as.data.frame(avg_sil_scepi)
avg_sil_scepi$X = labels_scepi
sil_scepi = dplyr::filter(avg_sil_scepi, grepl('Macrophage_Mouse|Bcell_Mouse|Endothelial_Mouse', X))

avg_sil_seurat = as.data.frame(avg_sil_seurat)
avg_sil_seurat$X = labels_seurat
sil_seurat = dplyr::filter(avg_sil_seurat, grepl('Macrophage_Mouse|Bcell_Mouse|Endothelial_Mouse', X))

avg_sil_conos = as.data.frame(avg_sil_conos)
avg_sil_conos$X = labels_conos
sil_conos = dplyr::filter(avg_sil_conos, grepl('Macrophage_mouse|Bcell_mouse|Endothelial_mouse', X))

avg_sil_liger = as.data.frame(avg_sil_liger)
avg_sil_liger$X = labels_liger
sil_liger = dplyr::filter(avg_sil_liger, grepl('Macrophage_mouse|Bcell_mouse|Endothelial_mouse', X))

boxplot(c(sil_seurat$V1),c(sil_scepi$V1),c(sil_liger$V1),c(sil_conos$V1),names = c("Seurat", "ScEpiSearch", "Liger", "Conos"),
        col = c("orange","red",'green','blue'),
        main = "3cells-Mouse",
        border = "brown",
        notch = TRUE)
