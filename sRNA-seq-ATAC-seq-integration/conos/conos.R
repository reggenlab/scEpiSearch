library(conos)
library(pagoda2)
library(parallel)
library(ggplot2)
library(Matrix)
########################Load reference dataset file###########################################
pbmc.rna <- read.csv("scRNA-scATAC-integration/MCA_reference.csv",row.names=1)
row.names(pbmc.rna) = toupper(row.names(pbmc.rna))

######################Load query atac-seq (gene activity scores) ###################################
#gene.activities <- read.table("scRNA-scATAC-integration/seurat-embed/activity_matrix_new/macropahe_activity_matrix_mouse.txt",row.names=1,sep=" ")
#gene.activities <- read.table("scRNA-scATAC-integration/seurat-embed/activity_matrix_new/endothelial_activity_matrix_mouse.txt",row.names=1,sep=" ")
#gene.activities <- read.table("scRNA-scATAC-integration/seurat-embed/activity_matrix_new/endothelial_activity_matrix_mouse.txt",row.names=1,sep=" ")
row.names(gene.activities) = toupper(row.names(gene.activities))

gene.activities <- read.table("scRNA-scATAC-integration/seurat-embed/activity_matrix_new/PBMC_activity_matrix.txt",row.names=1,sep=" ")
#gene.activities <- read.table("scRNA-scATAC-integration/gene_activity_GM.txt",row.names=1,sep=",")
#gene.activities <- read.table("scRNA-scATAC-integration/gene_activity_h1esc.txt",row.names=1,sep=",")
#gene.activities <- read.table("scRNA-scATAC-integration/gene_activity_neuron.txt",row.names=1,sep=",")
#gene.activities <- read.table("scRNA-scATAC-integration/gene_activity_MYOBLAST.txt",row.names=1,sep=",")

##################subset query cells based on celltype if required (For Liver 3 cells)#################################
#gene.activities <- read.table("scRNA-scATAC-integration/seurat-embed/activity_matrix_new/3cells_activity_matrix_mouse.txt",row.names=1,sep=" ")
#row.names(gene.activities) = toupper(row.names(gene.activities))
#endo = gene.activities[,1:318]
#macro = gene.activities[,319:413]
#bcell = gene.activities[,414:437]
#endo <- as(endo,"dgCMatrix") 
#macro <- as(macro,"dgCMatrix") 
#bcell <- as(bcell,"dgCMatrix") 
#data <- list(atac_1 = endo,atac_2=macro,atac_3=bcell, rna = rna_dgc)

################convert queries to dgcMatrix#############################
gene.activities <- as(gene.activities,"dgCMatrix") 
rna_dgc =  as(pbmc.rna, "dgCMatrix")
data <- list(atac = gene.activities, rna = rna_dgc)

#####################combine reference and query matrices#########################

p2l <- mclapply(data,basicP2proc,n.odgenes=3e3,nPcs=30,make.geneknn=F,n.cores=30,mc.cores=1)
metadata.rna = read.csv("scRNA-scATAC-integration/MCA_reference_labels.csv",header=F)
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
#assign names to atac groups according to their cell types
tsne$Group[tsne$Group == 'atac'] = 'PBMC_human'
#tsne$Group[tsne$Group == 'atac'] = 'macrophage_mouse'
#tsne$Group[tsne$Group == 'atac_1'] = 'Endothelial_mouse'
#tsne$Group[tsne$Group == 'atac_2'] = 'Macrophage_mouse'
#tsne$Group[tsne$Group == 'atac_3'] = 'Bcell_mouse'

