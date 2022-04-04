library('Seurat')
library('Matrix')

###  using seurat 3.2
###################################################

##  This file contains codes for co-embedding single-cell ATAC-seq data with reference mouse cell expression profiles 
#   You can download needed files from 
#   
#####################################################################################################
#   First part is for co-embedding Human cell scATAC-seq profile with mouse single-cell expression profile
#  
#
#
########## loading mouse cell atlas expression profile 

rna = read.csv('scRNA-scATAC-integration/MCA_reference_raw_count.tsv', header=T, sep=" ") ;
data = rna ;
rownames(data) = toupper(rownames(data)) ;

label=read.csv('scRNA-scATAC-integration/MCA_reference_labels.csv', header=F) ;

#########

atacf = c('scRNA-scATAC-integration/seurat data/new-h1esc_seurat_hg19.csv' , 'scRNA-scATAC-integration/seurat data/new-neuron_seurat_hg19.csv', 'scRNA-scATAC-integration/seurat data/new-GM_seurat_hg19.csv' , 'scRNA-scATAC-integration/seurat data/pbmc-good-cells.csv') ;

for (i in 1:4) 
{

##  for every iteration we define object for RNA to reduce artfact from previous iteration calculation

mca.expr <- CreateSeuratObject(counts = as.matrix(data), project = "myRNA", min.cells = 3, min.features = 200)
mca.expr <- NormalizeData(mca.expr, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(mca.expr)
mca.expr$celltype = as.matrix(label) ;
mca.expr <- ScaleData(mca.expr, features = all.genes)
mca.expr@assays$RNA@data@x[is.na(mca.expr@assays$RNA@data@x)] <- 0
mca.expr <- FindVariableFeatures(mca.expr, selection.method = "vst", nfeatures = 2000)
mca.expr <- RunPCA(mca.expr)
mca.expr <- RunTSNE(mca.expr, dims = 1:30)

#mca.expr <- FindNeighbors(mca.expr, dims = 1:30)
#mca.expr <- FindClusters(mca.expr, resolution = 0.5)


#####  Loading scATAC-seq profile and making gene-activity matrix ###########

atac = read.csv(atacf[i], header=F);
atacd = atac[, 2:ncol(atac)] ;
rownames(atacd) = atac[,1] ;
activity.matrix <- CreateGeneActivityMatrix(peak.matrix =atacd, annotation.file ="Homo_sapiens.GRCh37.75.gtf", seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)


human.atac <- CreateSeuratObject(counts = atacd, assay = "ATAC", project = "my_ATAC")
human.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
human.atac$tech <- "atac"
DefaultAssay(human.atac) <- "ACTIVITY"
#human.atac <- FindVariableFeatures(human.atac)
#human.atac <- NormalizeData(human.atac)
#human.atac <- ScaleData(human.atac)

DefaultAssay(human.atac) <- "ATAC"
#  here we have use peaks with sum of count atleast 10 acorss the cells as our data-size is small 
thr = 0.05*ncol(atacd) ;
VariableFeatures(human.atac) <- names(which(Matrix::rowSums(human.atac) > 2))
human.atac <- RunLSI(human.atac, n = 30, scale.max = NULL)


## integrating expression and scATAC-seq profile
##  our query data-size is small to we use k.filter=80 ##### 


transfer.anchors <- FindTransferAnchors(reference = mca.expr, query = human.atac, features =VariableFeatures(object = mca.expr),  reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca", k.anchor=5, k.filter=80 )

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = mca.expr$celltype,  weight.reduction = human.atac[["lsi"]])

human.atac <- AddMetaData(human.atac, metadata = celltype.predictions)


genes.use <- VariableFeatures(mca.expr )
refdata <- GetAssayData( mca.expr, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = "cca")

# this line adds the imputed data matrix to the human.atac object
human.atac[["RNA"]] <- imputation
coembed <- merge(x = mca.expr, y = human.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
#p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
#CombinePlots(list(p1, p2))


embeds = Embeddings(coembed[["tsne"]])

fname = paste("embedsHuman-MCA-cell-", atacf[i],  sep="") ;
write.table(embeds, file=fname) ;

}


######################################################################################

#    
#
#
################################################################
#  The portion below contains the code for co-embedding mouse single-cell ATAC-seq with reference mouse expression profile
#
#
###################################################


atacf = c('scRNA-scATAC-integration/seurat data/mouse-liver-endothelial-readcount.csv', 'scRNA-scATAC-integration/seurat data/mouse-liver-macrophage-readcount.csv', 'scRNA-scATAC-integration/seurat data/mouse-liver-endo-1:318-macrp-319:413-B-414:437.csv','scRNA-scATAC-integration/seurat data/mouse-liver-Bdata-readcount.csv' ) ;

peakf = 'mouse-liver-peaks.txt' ;

for (i in 1:4) 
{

## we create RNA-seq object in order to avoid any bias from previous iteration 

mca.expr <- CreateSeuratObject(counts = as.matrix(data), project = "myRNA", min.cells = 3, min.features = 200)
mca.expr <- NormalizeData(mca.expr, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(mca.expr)
mca.expr$celltype = as.matrix(label) ;
mca.expr <- ScaleData(mca.expr, features = all.genes)
mca.expr@assays$RNA@data@x[is.na(mca.expr@assays$RNA@data@x)] <- 0
mca.expr <- FindVariableFeatures(mca.expr, selection.method = "vst", nfeatures = 2000)
mca.expr <- RunPCA(mca.expr)
mca.expr <- RunTSNE(mca.expr, dims = 1:30)


#######  Loading scATAC-seq profile and making gene-activity matrix ###########


atac = read.table(atacf[i]);
peak = read.table(peakf) ;
rownames(atac) = as.matrix(peak) ; 

activity.matrix <- CreateGeneActivityMatrix(peak.matrix =atac, annotation.file ="gencode.vM1.annotation.gtf", seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)

mouse.atac <- CreateSeuratObject(counts = atac, assay = "ATAC", project = "my_ATAC")
mouse.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
mouse.atac$tech <- "atac"
DefaultAssay(mouse.atac) <- "ACTIVITY"
#mouse.atac <- FindVariableFeatures(mouse.atac)
mouse.atac <- NormalizeData(mouse.atac)
mouse.atac <- ScaleData(mouse.atac)

DefaultAssay(mouse.atac) <- "ATAC"
thr = 0.05 * ncol(atac) ;
VariableFeatures(mouse.atac) <- names(which(Matrix::rowSums(mouse.atac) > thr))
mouse.atac <- RunLSI(mouse.atac, n = 30, scale.max = NULL)


############ integrating expression and scATAC-seq profile

transfer.anchors <- FindTransferAnchors(reference = mca.expr, query = mouse.atac, features =VariableFeatures(object = mca.expr),  reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca",  k.anchor = 5, k.filter = 80 , npcs=30)

# transfer.anchors <- FindTransferAnchors(reference = mca.expr, query = mouse.atac, features =VariableFeatures(object = mca.expr),  reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca"  )

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = mca.expr$celltype,  weight.reduction = mouse.atac[["lsi"]])

mouse.atac <- AddMetaData(mouse.atac, metadata = celltype.predictions)


genes.use <- VariableFeatures(mca.expr )
refdata <- GetAssayData( mca.expr, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = mouse.atac[["lsi"]])

# this line adds the imputed data matrix to the mouse.atac object
mouse.atac[["RNA"]] <- imputation
coembed <- merge(x = mca.expr, y = mouse.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
#p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
#CombinePlots(list(p1, p2))


embeds = Embeddings(coembed[["tsne"]])

fname = paste("embedsMouse-MCA-cell-", atacf[i],  sep="\t") ;
write.table(embeds, file=fname) ;

}
