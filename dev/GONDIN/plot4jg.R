


library(Seurat)
library(qs)

data = qread("/media/mna_bioinfo/DD_linux/references/Mmu/singlecell/mckellard/downsample_1k.qs")

mtdt = merge(data@meta.data,data@reductions$umap_harmony@cell.embeddings,by=0) %>% column_to_rownames(var="Row.names")


FeaturePlot(data,features=c("Esr1","Esr2"),reduction="umap_harmony",cols=c("blue","orange","red"),order = T) +
  DimPlot(data,reduction="umap_harmony",group.by = "harmony_factorIDs",label=T)

unique(data$harmony_factorIDs)

sdata <- subset(x = data, )


VlnPlot(data,features="Esr1",group.by = "harmony_factorIDs",split.by = "sex",log = T)



library(tidyverse)
# esr1_data = FetchData(data,vars = c("Esr1","Esr2")) %>% merge(mtdt,by=0) %>%
esr1_data = FetchData(data,vars = c("Esr1","Esr2")) %>% merge(mtdt,by=0) %>%
  mutate(age.months=factor(age.months,levels=order))

esr1_data_f = filter(esr1_data, Esr1>0)
MIN=min(na.omit(esr1_data_f$Esr1))
MID=quantile(na.omit(esr1_data_f$Esr1),.5)
MAX=max(na.omit(esr1_data_f$Esr1))

p1 = ggplot() +
  geom_point(data=dplyr::select(esr1_data,-Esr1),aes(x=UMAPHarmony_1,y=UMAPHarmony_2),color="grey80") +
  geom_point(data=arrange(esr1_data_f,Esr1), aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=Esr1)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(0.1,1,2.9),
                        limits = c(0.1,2.9)) +
  theme_bw()+
  theme(legend.position="bottom")

esr2_data_f = filter(esr1_data, Esr2>0)
MIN=min(na.omit(esr2_data_f$Esr1))
MID=quantile(na.omit(esr2_data_f$Esr1),.5)
MAX=max(na.omit(esr2_data_f$Esr1))

p2 = ggplot() +
  geom_point(data=dplyr::select(esr1_data,-Esr2),aes(x=UMAPHarmony_1,y=UMAPHarmony_2),color="grey80") +
  geom_point(data=arrange(esr2_data_f,Esr2), aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=Esr2)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(0.08,1.01,2.9),
                        limits =  c(0.08,2.9)) +
  theme_bw() +
  theme(legend.position="bottom")

library(patchwork)
p12 = p1 + p2



library(ggplot2)
p4 = ggplot(filter(esr1_data,Esr1>0 & injury.agent=="cardiotoxin"), aes(x=harmony_factorIDs,y=Esr1,fill=sex)) + 
  geom_boxplot() +
  geom_jitter(alpha=.2)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1)) + 
  facet_grid(source.label~age.months) +
  coord_flip() +
  theme(legend.position = "bottom") +
ggplot(filter(esr1_data,Esr2>0 & injury.agent=="cardiotoxin"), aes(x=harmony_factorIDs,y=Esr2,fill=sex)) + 
  geom_boxplot() +
  geom_jitter(alpha=.2)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1)) + 
  facet_grid(source.label~age.months) +
  coord_flip() +
  theme(legend.position = "bottom")

ggsave(plot=DimPlot(data,reduction="umap_harmony",group.by = "harmony_factorIDs",label=T), path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="global_view_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot=p12, path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="view_expression_esr1_esr2_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)
ggsave(plot=p4, path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="boxplot_esr1_esr2_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)



