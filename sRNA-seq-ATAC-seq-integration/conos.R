library(conos)
library(pagoda2)
library(parallel)
library(ggplot2)
library(Matrix)
pbmc.rna <- read.csv("/storage/vibhor/neeteshp/shreyas_tools_testing/MCA_reference.csv",row.names=1)
row.names(pbmc.rna) = toupper(row.names(pbmc.rna))
gene.activities <- read.table("/storage/vibhor/neeteshp/shreyas_tools_testing/seurat-embed/activity_matrix_new/3cells_activity_matrix_mouse.txt",row.names=1,sep=" ")
row.names(gene.activities) = toupper(row.names(gene.activities))
endo = gene.activities[,1:318]
macro = gene.activities[,319:413]
bj = gene.activities[,414:437]
gene.activities <- as(gene.activities,"dgCMatrix") 
endo <- as(endo,"dgCMatrix") 
macro <- as(macro,"dgCMatrix") 
bj <- as(bj,"dgCMatrix") 
rna_dgc =  as(pbmc.rna, "dgCMatrix")
data <- list(atac = gene.activities, rna = rna_dgc)
data <- list(atac_1 = endo,atac_2=macro,atac_3=bj, rna = rna_dgc)

p2l <- mclapply(data,basicP2proc,n.odgenes=3e3,nPcs=30,make.geneknn=F,n.cores=30,mc.cores=1)
metadata.rna = read.csv("/storage/vibhor/neeteshp/shreyas_tools_testing/tsne_coords/MCA_reference_labels.csv",header=F)
rownames(metadata.rna) = colnames(rna_dgc)

l.con <- Conos$new(p2l,n.cores=30)
l.con$buildGraph(k=15,k.self=5,k.self.weigh=0.01,ncomps=30,n.odgenes=5e3,space='PCA') 

l.con$findCommunities(resolution=1.5)
l.con$embedGraph(alpha=1/2);

p1 <- l.con$plotGraph(font.size=c(3,5),title='conos clusters',alpha=0.2) #+ annotate("text",  x=-Inf, y = Inf, label = "clusters", vjust=1, hjust=0)
p2 <- l.con$plotGraph(groups=as.factor(metadata.rna),mark.groups=T,alpha=0.2,plot.na=F,title='annotation: RNA',font.size=c(3,5))+xlim(range(l.con$embedding[,1]))+ylim(range(l.con$embedding[,2]));
p2c <- l.con$plotGraph(groups=atac.annotation,mark.groups=T,alpha=0.2,plot.na=F,title='annotation: ATAC',font.size=c(3,5))+xlim(range(l.con$embedding[,1]))+ylim(range(l.con$embedding[,2]));
p3 <- l.con$plotGraph(color.by='sample',mark.groups=T,alpha=0.1,show.legend=T,title='platform',raster=T)+theme(legend.position=c(1,1),legend.justification = c(1,1))+guides(color=guide_legend(ncol=2,override.aes = list(size=3,alpha=0.8)))
tsne = p3$data
tsne$Group <- as.character(tsne$Group)
metadata.rna$V2 = rownames(metadata.rna)
for(i in 1:nrow(tsne)){
  if (tsne$CellName[i] %in% metadata.rna$V2)
  {
    ind = metadata.rna$V2 == tsne$CellName[i]
    tsne[i,4] = metadata.rna$V1[ind][1]
  }
}
r1 = data.frame(l.con$samples$atac$embeddings$PCA$largeVis)
r1$CellName = rownames(r1)
r1$Group = 'Myoblast_human'
colnames(r1) = c("x","y","CellName","Group")
r1$x = tsne$x[tsne$Group=='atac']
r1$y = tsne$y[tsne$Group=='atac']
#
r1$x = r1$x+rnorm(3*10, sd=0.1)
r1$y = r1$y+rnorm(3*10, sd=0.1)
tsne = rbind(tsne,r1)
tsne$Group[tsne$Group == 'atac'] = 'macrophage_mouse'
tsne$Group[tsne$Group == 'atac_1'] = 'Endothelial_mouse'
tsne$Group[tsne$Group == 'atac_2'] = 'Macrophage_mouse'
tsne$Group[tsne$Group == 'atac_3'] = 'Bcell_mouse'
library(dplyr)
colnames(tsne) = c("a","TSNE_1","TSNE_2","X")
data1 = dplyr::filter(tsne, grepl('Macro|Macrophage_mouse|B cell|Bcell_mouse|Endo|Endothelial_mouse', X))
data2 = dplyr::filter(tsne, !grepl('Macro|Macrophage_mouse|B cell|Bcell_mouse|Endo|Endothelial_mouse', X))

data1 = dplyr::filter(tsne, grepl('Macro|macrophage_mouse', X))
data2 = dplyr::filter(tsne, !grepl('Macro|macrophage_mouse', X))
data1 = dplyr::filter(tsne, grepl('Endo|endothelial_mouse', X))
data2 = dplyr::filter(tsne, !grepl('Endo|endothelial_mouse', X))
data1 = dplyr::filter(tsne, grepl('T cell|B cell|mono|dendritic|PBMC_human', X))
data2 = dplyr::filter(tsne, !grepl('T cell|B cell|mono|dendritic|PBMC_human', X))
data1 = dplyr::filter(data, grepl('Embryonic-Stem-Cell|H1ESC|GM_human|B cell|BJ_human|fibro', X))
data2 = dplyr::filter(data, !grepl('Embryonic-Stem-Cell|H1ESC|GM_human|B cell|BJ_human|fibro', X))
data1 = dplyr::filter(tsne, grepl('Embryonic-Stem-Cell|H1ESC', X))
data2 = dplyr::filter(tsne, !grepl('Embryonic-Stem-Cell|H1ESC', X))
data1 = dplyr::filter(tsne, grepl('GM_human|B cell', X))
data2 = dplyr::filter(tsne, !grepl('GM_human|B cell', X))
data1 = dplyr::filter(tsne, grepl('Myoblast_human|myocyte', X))
data2 = dplyr::filter(tsne, !grepl('Myoblast_human|myocyte', X))
data1 = dplyr::filter(tsne, grepl('neuron|neuron_human', tolower(X)))
data2 = dplyr::filter(tsne, !grepl('neuron|neuron_human', tolower(X)))
data2$X = "other"
f_data = rbind(data1,data2)
colnames(f_data) = c("a","TSNE_1","TSNE_2","X")
f_data= f_data[order(-as.numeric(factor(f_data$X))),]
library(Polychrome)
library(ggpubr)
seed <- c("#ff0000", "#00ff00")
mycolors <- createPalette(length(unique(f_data$X)), seed, prefix="mine")
names(mycolors) <- levels(as.factor(f_data$X))
mycolors['B cell_Vpreb3 high(Peripheral_Blood)'] = 'black'
mycolors['other'] = '#DEDEDE33'
mycolors['H1ESC'] = 'purple'
colScale <- scale_colour_manual(name = "X",values = mycolors)
g=ggscatter(f_data, x = "TSNE_1", y = "TSNE_2",
            color = "X",size=1)+colScale+font("legend.text",color = "black",size = 10)+theme(legend.position='right')+ guides(colour = guide_legend(override.aes = list(size=10)))
g


drops <- c("Group")
tsne = tsne[ , !(names(tsne) %in% drops)]

