library(corrplot)
# Clustering
library(cluster) 
library(factoextra)
library(dplyr)

data_mint <- read.csv('tsne_silhoutte/case_4_mint.csv', header = FALSE,sep=",")
drops <- c("V1","V3","V5")
data_mint = data_mint[ , !(names(data_mint) %in% drops)]
colnames(data_mint) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_mint, grepl('Human-GM|Mouse-Bcell', X))
data2 = dplyr::filter(data_mint, grepl('Human-Tcell|Mouse-Tcell', X))
data1$X = "bcell"
data2$X = "tcell"
data_mint = rbind(data1,data2)
colnames(data_mint) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_mint,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_mint$X = d1$X
dim(data_mint)

data_scanorama <- read.csv('tsne_silhoutte/scanorama_case4.csv',sep=",",header=FALSE)
data_scanorama <- data_scanorama[-c(1)]
colnames(data_scanorama) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_scanorama, grepl('Human-GM|Mouse-Bcell', X))
data2 = dplyr::filter(data_scanorama, grepl('Human-Tcell|Mouse-Tcell', X))
data1$X = "bcell"
data2$X = "tcell"
data_scanorama = rbind(data1,data2)
colnames(data_scanorama) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scanorama,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scanorama$X = d1$X
dim(data_scanorama)

data_scale <- read.csv('tsne_silhoutte/scale_case4.csv',sep=",",header=FALSE)
data_scale <- data_scale[-c(1)]
colnames(data_scale) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_scale, grepl('Human-GM|Mouse-Bcell', X))
data2 = dplyr::filter(data_scale, grepl('Human-Tcell|Mouse-Tcell', X))
data1$X = "bcell"
data2$X = "tcell"
data_scale = rbind(data1,data2)
colnames(data_scale) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scale,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scale$X = d1$X
dim(data_scale)

data_scvi <- read.csv('tsne_silhoutte/scvi_case4.csv',sep=",",header=FALSE)
data_scvi <- data_scvi[-c(1)]
colnames(data_scvi) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_scvi, grepl('Human-GM|Mouse-Bcell', X))
data2 = dplyr::filter(data_scvi, grepl('Human-Tcell|Mouse-Tcell', X))
data1$X = "bcell"
data2$X = "tcell"
data_scvi = rbind(data1,data2)
colnames(data_scvi) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scvi,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scvi$X = d1$X
dim(data_scvi)

data_scepi <- read.csv('tsne_silhoutte/scepi_pos_case4.csv',sep=",",header=FALSE)
data_scepi <- data_scepi[-c(1)]
colnames(data_scepi) = c("TSNE_1","TSNE_2","X")
data1 = dplyr::filter(data_scepi, grepl('Human-GM|Mouse-Bcell', X))
data2 = dplyr::filter(data_scepi, grepl('Human-Tcell|Mouse-Tcell', X))
data1$X = "bcell"
data2$X = "tcell"
data_scepi = rbind(data1,data2)
colnames(data_scepi) = c("TSNE_1","TSNE_2","X")
d1 = as.data.frame(apply(data_scepi,2,function(x) {x<-as.numeric(factor(x,levels = unique(x)))}))
data_scepi$X = d1$X
dim(data_scepi)

df_mint <- data_mint[-c(3)]
df_scale <- data_scale[-c(3)]
df_scanorama <- data_scanorama[-c(3)]
df_scepi <- data_scepi[-c(3)]
df_scvi <- data_scvi[-c(3)]

silhouette_score_mint <- function(k){
  #km <- kmeans(df_seurat, centers = k)
  ss <- silhouette(data_mint$X, dist(df_mint))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_scale <- function(k){
  #km <- kmeans(df_conos, centers = k)
  ss <- silhouette(data_scale$X, dist(df_scale))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_scanorama <- function(k){
  #km <- kmeans(df_liger, centers = k)
  ss <- silhouette(data_scanorama$X, dist(df_scanorama))
  #mean(ss[, 3])
  ss[,3]
}
silhouette_score_scvi <- function(k){
  #km <- kmeans(df_liger, centers = k)
  ss <- silhouette(data_scvi$X, dist(df_scvi))
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
avg_sil_scale <- sapply(k, silhouette_score_scale)
avg_sil_scanorama <- sapply(k, silhouette_score_scanorama)
avg_sil_scvi <- sapply(k, silhouette_score_scvi)
avg_sil_mint <- sapply(k, silhouette_score_mint)
avg_sil_scepi <- sapply(k, silhouette_score_scepi)
#plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

boxplot(c(avg_sil_scale),c(avg_sil_scepi),c(avg_sil_scanorama),c(avg_sil_scvi),c(avg_sil_mint), names = c("SCALE", "ScEpiSearch", "SCANORAMA", "SCVI","MINT"),
        col = c("orange","red",'green','blue','yellow'),
        main = "Case-2",
        border = "brown",
        notch = TRUE)

