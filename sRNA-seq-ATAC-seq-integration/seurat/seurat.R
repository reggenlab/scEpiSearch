atac = read.csv('new-h1esc_seurat_hg19.csv')
atacd = atac[, 2:ncol(atac)] ;
rownames(atacd) = atac[,1] ;
M = Matrix(as.matrix(atacd), byrow=TRUE ,nrow=nrow(atacd),sparse=TRUE);
activity.matrix <- CreateGeneActivityMatrix(peak.matrix =M,
annotation.file ="Homo_sapiens.GRCh37.75.gtf", seq.levels = c(1:22,
"X", "Y"), upstream = 2000, verbose = TRUE)

pbmc.atac <- CreateSeuratObject(counts = M, assay = "ATAC", project = "my_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
pbmc.atac$tech <- "atac"

#for RNA-seq
rna = read.csv('/home/cellsearch/cellatlassearch_shreya/epigenome_search/compare_methods/MCA_reference.csv')
;
data = rna[,2:ncol(rna)] ;
rownames(data) = toupper(rna[,1]) ;
label=read.csv('/home/cellsearch/cellatlassearch_shreya/epigenome_search/compare_methods/MCA_reference_labels.csv')
;
label = c(as.matrix(label) , "tri" ) ;


pbmc <- CreateSeuratObject(counts = as.matrix(data), project =
"myFit", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA@data@x[is.na(pbmc@assays$RNA@data@x)] <- 0
pbmc$celltype = as.matrix(label) ;


#reu = VariableFeatures(object = pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#pdf('trial.pdf') ;
#pbmc <- RunTSNE(pbmc, dims = 1:20)
#dev.off() ;


transfer.anchors <- FindTransferAnchors(reference = pbmc, query =
pbmc.atac, features =VariableFeatures(object = pbmc),  reference.assay
= "RNA", query.assay = "ACTIVITY", reduction = "cca", k.anchor=3,
k.filter=20, k.score=20)

celltype.predictions <- TransferData(anchorset = transfer.anchors,
refdata = pbmc$celltype,  weight.reduction = "cca")


pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)


genes.use <- VariableFeatures(pbmc )
refdata <- GetAssayData(pbmc, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata =
refdata, weight.reduction = "cca")

# this line adds the imputed data matrix to the pbmc.atac object
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc, y = pbmc.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype,
coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))