


library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)
data = qread("/media/mna_bioinfo/DD_linux/references/Mmu/singlecell/mckellard/downsample_1k.qs")

mtdt = merge(data@meta.data,data@reductions$umap_harmony@cell.embeddings,by=0) %>% column_to_rownames(var="Row.names")


library(tidyverse)
Fshr_data = FetchData(data,vars = c("Bdnf")) %>% merge(mtdt,by=0)

Fshr_data_f = filter(Fshr_data, Bdnf>0)
MIN=min(na.omit(Fshr_data_f$Bdnf))
MID=quantile(na.omit(Fshr_data_f$Bdnf),.5)
MAX=max(na.omit(Fshr_data_f$Bdnf))

p1 = ggplot() +
  geom_point(data=dplyr::select(Fshr_data,-Bdnf),aes(x=UMAPHarmony_1,y=UMAPHarmony_2),color="grey80") +
  geom_point(data=arrange(Fshr_data_f,Bdnf), aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=Bdnf)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1"),
                        breaks = c(0.1,1,2.9),
                        limits = c(0.1,2.9)) +
  theme_bw()+
  theme(legend.position="bottom")

p2 = ggplot() +
  geom_point(data=mtdt, aes(x=UMAPHarmony_1,y=UMAPHarmony_2,color=harmony_factorIDs)) +
  theme_bw() +
  theme(legend.position="bottom")

library(ggplot2)
p4 = ggplot(filter(Fshr_data,Bdnf>0 & injury.agent=="cardiotoxin"), aes(x=harmony_factorIDs,y=Bdnf,fill=sex)) + 
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
ggsave(plot= (p1 + p4), path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/GONDIN",
       filename="view_expression_Bdnf_fromMckellarCosgrove.pdf",width = 20, height = 10, dpi=1000,units="in",
       device = cairo_pdf)


