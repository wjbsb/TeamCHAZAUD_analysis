


library(Seurat)
library(qs)

data = qread("/media/mna_bioinfo/DD_linux/references/Mmu/singlecell/mckellard/downsample_1k.qs")

mtdt = merge(data@meta.data,data@reductions$umap_harmony@cell.embeddings,by=0) %>% column_to_rownames(var="Row.names")
# 
FeaturePlot(data,features=c("Fshr","Cga"),reduction="umap_harmony",cols=c("blue","orange","red"),order = T) +
  DimPlot(data,reduction="umap_harmony",group.by = "harmony_factorIDs",label=T)
# 
# unique(data$harmony_factorIDs)
# 
# sdata <- subset(x = data, )
# VlnPlot(data,features="Fshr",group.by = "harmony_factorIDs",split.by = "sex",log = T)



library(tidyverse)
Fshr_data = FetchData(data,vars = c("Fshr","Cga")) %>% merge(mtdt,by=0)

table(Fshr_data$Cga>0)
Fshr_data_f = filter(Fshr_data, Fshr>0)
MIN=min(na.omit(Fshr_data_f$Fshr))
MID=quantile(na.omit(Fshr_data_f$Fshr),.5)
MAX=max(na.omit(Fshr_data_f$Fshr))

p1 = ggplot() +
  geom_point(data=dplyr::select(Fshr_data,-Fshr),aes(x=UMAPHarmony_1,y=UMAPHarmony_2),color="grey80") +
  geom_point(data=arrange(Fshr_data_f,Fshr), aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=Fshr)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(0.1,1,2.9),
                        limits = c(0.1,2.9)) +
  theme_bw()+
  theme(legend.position="bottom")

Cga_data_f = filter(Fshr_data, Cga>0)
MIN=min(na.omit(Cga_data_f$Fshr))
MID=quantile(na.omit(Cga_data_f$Fshr),.5)
MAX=max(na.omit(Cga_data_f$Fshr))

p2 = ggplot() +
  geom_point(data=dplyr::select(Fshr_data,-Cga),aes(x=UMAPHarmony_1,y=UMAPHarmony_2),color="grey80") +
  geom_point(data=arrange(Cga_data_f,Cga), aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=Cga)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(0.08,1.01,2.9),
                        limits =  c(0.08,2.9)) +
  theme_bw() +
  theme(legend.position="bottom")

library(patchwork)
p12 = p1 + p2



library(ggplot2)
p4 = ggplot(filter(Fshr_data,Fshr>0 & injury.agent=="cardiotoxin"), aes(x=harmony_factorIDs,y=Fshr,fill=sex)) + 
  geom_boxplot() +
  geom_jitter(alpha=.2)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1)) + 
  facet_grid(source.label~age.months) +
  coord_flip() +
  theme(legend.position = "bottom") +
  plot_spacer()
# ggplot(filter(Fshr_data,Cga>0 & injury.agent=="cardiotoxin"), aes(x=harmony_factorIDs,y=Cga,fill=sex)) + 
#   geom_boxplot() +
#   geom_jitter(alpha=.2)+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90,hjust=1)) + 
#   facet_grid(source.label~age.months) +
#   coord_flip() +
#   theme(legend.position = "bottom")

ggsave(plot=DimPlot(data,reduction="umap_harmony",group.by = "harmony_factorIDs",label=T), path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="global_view_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot=p12, path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="view_expression_Fshr_Cga_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot=p4, path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="boxplot_Fshr_Cga_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)



