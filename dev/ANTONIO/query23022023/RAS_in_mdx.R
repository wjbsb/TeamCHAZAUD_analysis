
path = "~/Bureau/development/github/single_cell_RNA_seq"
setwd(dir = path)

source("~/Bureau/development/github/single_cell_RNA_seq/functions_v2.R")
library(latex2exp)
rm(path,pathData)

# configuration_scrnaseq@paths@output_archives = "/home/legrand-lab/Bureau/workflow/scrnaseq/save_annotation_v3_pipeline/Rds_Rdata"
# configuration_scrnaseq@paths@output_plots = "/home/legrand-lab/Bureau/workflow/scrnaseq/save_annotation_v3_pipeline/Rds_Rdata"


integrated = qread(file = paste0("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/FABIEN_LEGRAND/scMDX_serieA/new_thresholds/Rds_Rdata/integrated_ccbr_annotation.qs"))

# FeaturePlot(integrated,reduction="fitsne_by_harmony",features=c("Rab10","Dmd"),blend = T,cols=c("grey80","red","blue"),pt.size = 1,ncol = 2,split.by = "sample")
# integrated$group = paste0(integrated$group,"_",integrated$Hypergroup)
# RidgePlot(integrated,features=c("Rab10","Dmd"),group.by = "group",sort = "increasing")

###DEMANDE 23022023
genesRAS = c("Cdc42","Dock10","Grb2","Mapk1","Mapk14","Mapk3","Mylk","Mylk2","Rab10","Rab11b","Rab12","Rab14","Rab18","Rab1a","Rab1b","Rab21","Rab22b",
"Rab2a","Rab35","Rab3a","Rab5a","Rab5b","Rab5c","Rab6a","Rab7a","Rac1","R-Ras")
###DEMANDE 28062023
genes = c("Snca","Ctsb","Ctsd","Dcn","Lamp2","Hspg2","Fga","F10","Bgn","CD74","Gpc1","S100a8","S100a9","Plg","Alb")
###DEMANDE 23102023
genes = c("Anxa1","Anxa2","Anxa3","Anxa4","Anxa5")
#wrkshop
genes = c("Nuak1","Serbf1")
plot_dea = function(dataset,genes,filename,path) {
  fdata = FetchData(dataset,assay="RNA",vars=genes)
  dataset = ScaleData(dataset,features = genes,assay = "RNA")
  
  genesRAS_present=colnames(fdata) %>% str_replace_all("rna_","")
  genesRAS_present2=paste0("rna_",genesRAS_present)
  # p = lapply(genesRAS_present, function(x) plot_imbalance_score_gene(integrated,"fitsne_by_harmony",x))
  # names(p)=genesRAS_present
  
  
  # scaling_data = dataset@assays$RNA@scale.data %>% data.frame()
  
  # for (i in 1:length(genesRAS_present)) {
  #   # is = plot_imbalance_score_v3(data_seurat = integrated,
  #   #                              reduceDimension_selected = toupper("fitsne_by_harmony"),
  #   #                              cl=genesRAS_present[i],
  #   #                              genes = T)
  #   tmp = ggarrange(
  #     FeaturePlot(integrated,features = genesRAS_present[i],reduction="fitsne_by_harmony",order=T,pt.size = 1.5,min.cutoff = "q5",cols = c("grey70","red")),
  #     # plot_density(integrated,genesRAS_present2[i],reduction="fitsne_by_harmony",slot="counts"),
  #     # is[[1]],
  #     # is[[2]],
  #     ncol=2,nrow=2)
  #   dir.create("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/ANTONIO/query28062023")
  #   ggsave(filename = paste0(genesRAS_present[i],".pdf"),device = cairo_pdf,path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/ANTONIO/query28062023",width = 15,height = 20,dpi=1000,plot = tmp)
  # }
  
  # colnames(fdata) = colnames(fdata) %>% str_replace_all("rna_","")
  # mtdt = extract_metadata_reduction_seurat(integrated) %>% merge(fdata,by=0) %>% pivot_longer(cols=colnames(fdata))
  
  Idents(dataset) = "predicted.id"
  dataset$grp = paste0(dataset$sample,"_",dataset$predicted.id)
  
  Idents(dataset) = "grp"
  dea = FindAllMarkers(dataset,features=genesRAS_present,assay = "RNA",test.use = "MAST",latent.vars = c("sample","lineage","condition"),logfc.threshold = -Inf,min.pct = 0)
  dea_m <- dea %>% 
    separate(col="cluster",into=c("sample","predicted.id","part2"),sep = "_",remove=F) %>%
    unite(predicted.id, predicted.id:part2,sep="_",remove=F,na.rm=T) %>%
    dplyr::select(-part2)
  MIN=min(na.omit(dea$avg_log2FC)) %>% round()
  MID=quantile(na.omit(dea$avg_log2FC),.5) %>% round()
  MAX=max(na.omit(dea$avg_log2FC)) %>% round()
  
  library(ggh4x)
  
  subtype_order = dplyr::select(dataset@meta.data,hypergroup,predicted.id) %>% 
    unique() %>%
    arrange(hypergroup) %>%       
    pull(predicted.id)
  
  
  dea_n = dea_m %>%
    mutate(
      predicted.id = factor(predicted.id, levels=subtype_order)) %>%
    full_join(dplyr::select(dataset@meta.data,hypergroup,predicted.id) %>% 
                unique(),
              by="predicted.id") %>%
    na.omit() %>%
    mutate(predicted.id = factor(predicted.id, levels=subtype_order)) %>%
    group_by(hypergroup,predicted.id,sample)
  #   typehypergroup = case_when(
  #          cluster == "Red Blood Cells" ~ "Red Blood Cells",
  #          cluster == "Schwann Cells" ~ "Schwann Cells",
  #          cluster == "Regulatory T Cells" ~ "Regulatory T Cells",
  #          cluster %in% c("Tenocytes_1","Tenocytes_2") ~ "Tenocytes",
  #          cluster %in% c("Pro-inflammatory MΦ","Resolution-Phase MΦ","Dendritic MΦ") ~ "Immun Cells",
  #          cluster %in% c("Quiescent MuSCs","Myocytes","Profilerating MuSCs","Myonuclei") ~ "Skeletal Muscle Cells",
  #          cluster %in% c("Sprouting ECs","Capillary ECs_2","Veinous ECs","Lympathic ECs","Capillary ECs_1","Arterial ECs") ~ "Endothelial Cells",
  #          cluster %in% c("Fibroblasts","Activated FAPs","Adventitial FAPs","Parenchymal FAPs") ~ "Fibroblasts",
  #          cluster %in% c("Pericytes","Smooth Muscle Cells") ~ "Mural Cells"
  #        )) %>%
  # mutate(typehypergroup=factor(typehypergroup,
  #                              levels=c("Skeletal Muscle Cells","Regulatory T Cells","Endothelial Cells","Schwann Cells","Mural Cells","Immun Cells","Fibroblasts","Tenocytes","Red Blood Cells")))
  # dotplot1 = ggplot(dea_n, aes(x=gene,y=predicted.id,color=avg_log2FC,size=pct.1*100)) + 
  #   geom_point() +
  #   scale_size_continuous(labels = scales::percent_format(scale = 1)) +
  #   scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1"),
  #                         na.value="grey80",
  #                         breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
  #                         values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
  #   labs(x="",y="",color="Average Log2 Fold Change",size="% expression in cluster") +
  #   facet_wrap(~sample) +
  #   # facet_grid2(hypergroup~cluster,scales = "free",space = "free",shrink = T,switch = "both",
  #   #             strip = strip_themed()) +
  #   theme_legrand2 +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
  #         axis.text.y = element_text(family="Helvetica",size=11,face="bold.italic"),
  #         strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
  #         strip.background.x = element_rect(),
  #         strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
  #         legend.position = "bottom")
  # dotplot1
  
  
  # 
  # scaling_data_m = extract_metadata_reduction_seurat(integrated) %>%
  #   merge(t(scaling_data) %>% data.frame(),by=0) %>%
  #   pivot_longer(cols=rownames(scaling_data)) %>%
  #   # mutate(
  #   #   Subtype = factor(Subtype, levels=subtype_order),
  #   #   typehypergroup = case_when(
  #   #     Subtype == "Red Blood Cells" ~ "Red Blood Cells",
  #   #     Subtype == "Schwann Cells" ~ "Schwann Cells",
  #   #     Subtype == "Regulatory T Cells" ~ "Regulatory T Cells",
  #   #     Subtype %in% c("Tenocytes_1","Tenocytes_2") ~ "Tenocytes",
  #   #     Subtype %in% c("Pro-inflammatory MΦ","Resolution-Phase MΦ","Dendritic MΦ") ~ "Immun Cells",
  #   #     Subtype %in% c("Quiescent MuSCs","Myocytes","Profilerating MuSCs","Myonuclei") ~ "Skeletal Muscle Cells",
  #   #     Subtype %in% c("Sprouting ECs","Capillary ECs_2","Veinous ECs","Lympathic ECs","Capillary ECs_1","Arterial ECs") ~ "Endothelial Cells",
  #   #     Subtype %in% c("Fibroblasts","Activated FAPs","Adventitial FAPs","Parenchymal FAPs") ~ "Fibroblasts",
  #   #     Subtype %in% c("Pericytes","Smooth Muscle Cells") ~ "Mural Cells"
  #   #   )) %>%
  #   # mutate(typehypergroup=factor(typehypergroup,
  #   #                              levels=c("Skeletal Muscle Cells","Regulatory T Cells","Endothelial Cells","Schwann Cells","Mural Cells","Immun Cells","Fibroblasts","Tenocytes","Red Blood Cells"))) %>%
  #   group_by(predicted.id,name) %>%
  #   summarise(stats=median(value))
  # 
  # head(dea_m)
  # head(scaling_data_m)
  # 
  # dea_scaling = inner_join(dea_m,scaling_data_m, by=c("cluster"="predicted.id","gene"="name"))
  # 
  # MIN=min(na.omit(dea_scaling$avg_log2FC)) %>% round()
  # MID=quantile(na.omit(dea_scaling$avg_log2FC),.5) %>% round()
  # MAX=max(na.omit(dea_scaling$avg_log2FC)) %>% round()
  # 
  # dotplot1 = ggplot(dea_scaling, aes(x=gene,y=cluster,color=avg_log2FC,size=pct.1*100)) + 
  #   geom_point() +
  #   scale_size_continuous(labels = scales::percent_format(scale = 1)) +
  #   scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1"),
  #                         na.value="grey80",
  #                         breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
  #                         values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
  #   labs(x="",y="",color="Average Log2 Fold Change",size="% expression in cluster") +
  #   # facet_grid2(typehypergroup~.,scales = "free",space = "free",shrink = T,switch = "both",
  #   #             strip = strip_themed()) +
  #   theme_legrand2 +
  #   theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
  #         axis.text.y = element_text(family="Helvetica",size=11,face="bold.italic"),
  #         strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
  #         strip.background.x = element_rect(),
  #         strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
  #         legend.position = "bottom")
  # dotplot1
  # 
  # MIN=min(na.omit(dea_scaling$stats)) %>% round()
  # MID=quantile(na.omit(dea_scaling$stats),.5) %>% round()
  # MAX=max(na.omit(dea_scaling$stats)) %>% round()
  # 
  # dotplot2 = ggplot(dea_scaling, aes(x=gene,y=cluster,color=stats,size=pct.1*100)) + 
  #   geom_point() +
  #   scale_size_continuous(labels = scales::percent_format(scale = 1)) +
  #   scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1"),
  #                         na.value="grey80",
  #                         breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
  #                         values=scales::rescale(c(MIN,MID,MAX)),limits=c(-1,1)) +
  #   labs(x="",y="",color="Scaling Log-Normalized Count",size="% expression in cluster") +
  #   # facet_grid2(typehypergroup~.,scales = "free",space = "free",shrink = T,switch = "both",
  #   #             strip = strip_themed()) +
  #   theme_legrand2 +
  #   theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
  #         axis.text.y = element_text(family="Helvetica",size=11,face="bold.italic"),
  #         strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
  #         strip.background.x = element_rect(),
  #         strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
  #         legend.position = "bottom")
  # 
  # ggsave(filename = "DEA_Scaling_RAS_genes_in_mouse.pdf",plot=dotplot1/dotplot2,path = "/home/legrand-lab/Bureau/demande_antionio/",
  #        width = 20,height = 15,dpi=1000,device = cairo_pdf)
  # 
  # 
  # scaling_data_m2 = extract_metadata_reduction_seurat(integrated) %>% 
  #   merge(t(scaling_data) %>% data.frame(),by=0) %>% 
  #   pivot_longer(cols=rownames(scaling_data)) %>%
  #   group_by(predicted.id,sample,name) %>%
  #   summarise(stats=median(value)) %>%
  #   mutate(
  #     Subtype = factor(predicted.id, levels=subtype_order))
  #   #   typehypergroup = case_when(
  #   #     Subtype == "Red Blood Cells" ~ "Red Blood Cells",
  #   #     Subtype == "Schwann Cells" ~ "Schwann Cells",
  #   #     Subtype == "Regulatory T Cells" ~ "Regulatory T Cells",
  #   #     Subtype %in% c("Tenocytes_1","Tenocytes_2") ~ "Tenocytes",
  #   #     Subtype %in% c("Pro-inflammatory MΦ","Resolution-Phase MΦ","Dendritic MΦ") ~ "Immun Cells",
  #   #     Subtype %in% c("Quiescent MuSCs","Myocytes","Profilerating MuSCs","Myonuclei") ~ "Skeletal Muscle Cells",
  #   #     Subtype %in% c("Sprouting ECs","Capillary ECs_2","Veinous ECs","Lympathic ECs","Capillary ECs_1","Arterial ECs") ~ "Endothelial Cells",
  #   #     Subtype %in% c("Fibroblasts","Activated FAPs","Adventitial FAPs","Parenchymal FAPs") ~ "Fibroblasts",
  #   #     Subtype %in% c("Pericytes","Smooth Muscle Cells") ~ "Mural Cells"
  #   #   )) %>%
  #   # mutate(typehypergroup=factor(typehypergroup,
  #   #                              levels=c("Skeletal Muscle Cells","Regulatory T Cells","Endothelial Cells","Schwann Cells","Mural Cells","Immun Cells","Fibroblasts","Tenocytes","Red Blood Cells")))
  # 
  # dea_scaling_2 = inner_join(dea_m,scaling_data_m2, by=c("cluster"="predicted.id","gene"="name"))
  # MIN=min(na.omit(dea_scaling_2$stats)) %>% round()
  # MID=quantile(na.omit(dea_scaling_2$stats),.5) %>% round()
  # MAX=max(na.omit(dea_scaling_2$stats)) %>% round()
  
  # dotplot3 = ggplot(dea_scaling_2, aes(x=gene,y=cluster,color=stats,size=pct.1*100)) + 
  dotplot3 = ggplot(dea_n, aes(x=gene,y=predicted.id,color=avg_log2FC,size=pct.1*100)) + 
    geom_point() +
    scale_size_continuous(labels = scales::percent_format(scale = 1)) +
    scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1"),
                          na.value="grey80",
                          breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
    labs(x="",y="",color="Scaling Log-Normalized Count",size="% expression in cluster") +
    # facet_wrap(sample,nrow = 2,ncol = 2) +
    facet_grid(hypergroup~sample,scales = "free",space = "free",shrink = T,switch = "both") +
    # facet_grid2(vars(hypergroup),vars(sample),scales = "free",space = "free",shrink = T,switch = "both",
    #             strip = strip_themed()) +
    theme_legrand2 +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
          axis.text.y = element_text(family="Helvetica",size=11,face="bold.italic"),
          strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
          strip.background.x = element_rect(),
          strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
          legend.position = "bottom")
  
  # ggsave(filename = "Scaling_RAS_genes_in_mouse_by_sample.pdf",plot=dotplot3,path = "/home/legrand-lab/Bureau/demande_antionio/",
  #        width = 20,height = 15,dpi=1000,device = cairo_pdf)
  
  dataset = AddMetaData(dataset,fdata)
  fdata_knn = SmoothKNN(dataset,signature.names = setNames(as.list(genes),genes),reduction = "fitsne_by_harmony") %>%
    extract_metadata_reduction_seurat()
  p = list()
  for(g in paste0(genes,"_kNN")) {
    
    tmp = filter(fdata_knn, fdata_knn[[g]]>0) %>% arrange(.data[[g]])
    
    MIN=min(na.omit(tmp[[g]])) %>% round()
    MID=quantile(na.omit(tmp[[g]],.5)) %>% round()
    MAX=max(na.omit(tmp[[g]])) %>% round()
    
    p[[g]] = ggplot() + 
      geom_point(data=dplyr::select(fdata_knn,-g), aes(x=FItSNE_1,y=FItSNE_2),color="grey80") + 
      geom_point(data= tmp, aes(x=FItSNE_1,y=FItSNE_2,color=.data[[g]])) +
      scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1"),
                          na.value="grey80",
                          breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                          values=scales::rescale(c(MIN,MID,MAX)),limits=c(MIN,MAX)) +
      theme_legrand2 +
      labs(color=g)
  }
  
  if (!dir.exists(path)) {
    dir.create(path)
  }
  
  ggsave(filename,plot = dotplot3,device = cairo_pdf,path = path,
         width = 25,height = 15,dpi=1000)
  ggsave(paste0(filename,"_view.pdf"),plot = ggarrange(plotlist=p,ncol=2,nrow=3) + DimPlot(dataset,group.by="predicted.id",label=T,reduction="fitsne_by_harmony") ,device = cairo_pdf,path = path,
         width = 25,height = 15,dpi=1000)
  ^return(dotplot3)
}
plot_dea(dataset =integrated,
         genes=genes, 
         filename="query23102023_allgenes.pdf",
         path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/ANTONIO/query23102023/")
plot_dea(dataset =integrated,
         genes=genes, 
         # filename="query28062023_allgenes.pdf",
         filename= "workshopt.pdf",
         path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/ANTONIO/query28062023/")
