#####################################visualization of result#####################################
data = read.table(,'macrophage_liver_scepi_tsne.txt')
#####################################grep celltypes to be colored separately #############################
data1 = dplyr::filter(data, grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))
data2 = dplyr::filter(data, !grepl('Macro|Macrophage_Mouse|B cell|Bcell_Mouse|Endo|Endothelial_Mouse', X))
data2$X = "other"
f_data = rbind(data1,data2)
colnames(f_data) = c("TSNE_1","TSNE_2","X")
f_data= f_data[order(-as.numeric(factor(f_data$X))),]

library(Polychrome)
library(ggpubr)
seed <- c("#ff0000", "#00ff00")
mycolors <- createPalette(length(unique(f_data$X)), seed, prefix="mine")
names(mycolors) <- levels(as.factor(f_data$X))
mycolors['other'] = 'grey73'
#mycolors['H1ESC'] = 'purple'
colScale <- scale_colour_manual(name = "X",values = mycolors)
g=ggscatter(f_data, x = "TSNE_1", y = "TSNE_2",
            color = "X",size=1)+colScale+font("legend.text",color = "black",size = 20)+theme(legend.position='right')+ guides(colour = guide_legend(override.aes = list(size=10)))
g

