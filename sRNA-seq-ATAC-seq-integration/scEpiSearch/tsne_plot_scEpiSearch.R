#####################################visualization of result#####################################
data = read.table('scepisearch_integration/3cells_liver_scepi_tsne.txt',sep=" ")
#data = read.csv("scepisearch_integration/pbmc_scepi_tsne.txt",sep=" ")
#data = read.csv("scepisearch_integration/macrophage_liver_scepi_tsne.txt",sep=" ")
#data = read.csv("scepisearch_integration/endothelial_liver_scepi_tsne.txt",sep=" ")
#data = read.csv("scepisearch_integration/h1esc_scepi_tsne.txt",sep=" ")
#data = read.csv("scepisearch_integration/neuron_scepi_tsne.txt",sep=" ")
#data = read.csv("scepisearch_integration/GM_scepi_tsne.txt",sep=" ")

#data = data[,c(2,3,4)]

#####################################grep celltypes to be colored separately #############################

library(dplyr)
colnames(data) = c("TSNE_1","TSNE_2","X")
#Un-comment the 2 lines for which dataset was loaded.
#data1 = dplyr::filter(data, grepl('T cell|B cell|mono|dendritic|PBMC_human', X))
#data2 = dplyr::filter(data, !grepl('T cell|B cell|mono|dendritic|PBMC_human', X))
#data1 = dplyr::filter(data, grepl('Embryonic-Stem-Cell|H1ESC', X))
#data2 = dplyr::filter(data, !grepl('Embryonic-Stem-Cell|H1ESC', X))
#data1 = dplyr::filter(data, grepl('GM_human|B cell', X))
#data2 = dplyr::filter(data, !grepl('GM_human|B cell', X))
#data1 = dplyr::filter(data, grepl('Myoblast_Human|myocyte', X))
#data2 = dplyr::filter(data, !grepl('Myoblast_Human|myocyte', X))
#data1 = dplyr::filter(data, grepl('neuron|neuron_human', X))
#data2 = dplyr::filter(data, !grepl('neuron|neuron_human', X))
#data1 = dplyr::filter(data, grepl('Endo|Endothelial_Mouse', X))
#data2 = dplyr::filter(data, !grepl('Endo|Endothelial_Mouse', X))
#data1 = dplyr::filter(data, grepl('Macro|Macrophage_Mouse', X))
#data2 = dplyr::filter(data, !grepl('Macro|Macrophage_Mouse', X))
data1 = dplyr::filter(data, grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))
data2 = dplyr::filter(data, !grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))

data2$X = "other"
f_data = rbind(data1,data2)
#colnames(f_data) = c("a","TSNE_1","TSNE_2","X")
colnames(f_data) = c("TSNE_1","TSNE_2","X")
f_data= f_data[order(-as.numeric(factor(f_data$X))),]

#########################################################################ggplot based
library(Polychrome)
library(ggpubr)
seed <- c("#ff0000", "#00ff00")
mycolors <- createPalette(length(unique(f_data$X)), seed, prefix="mine")
names(mycolors) <- levels(as.factor(f_data$X))
mycolors['other'] = '#DEDEDE33'
#mycolors['Marginal zone B cell(Spleen)'] = '#66CC00'
#mycolors['H1ESC'] = 'purple'
colScale <- scale_colour_manual(name = "X",values = mycolors)
g=ggscatter(f_data, x = "TSNE_1", y = "TSNE_2",
            color = "X",size=1)+colScale+font("legend.text",color = "black",size = 10)+theme(legend.position='right')+ guides(colour = guide_legend(override.aes = list(size=10)))
g


