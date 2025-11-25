library(readr)

search = read_csv("/media/legrand-lab/Stockage/InstitutNeuroMyoGene_DATA/collaboration/MNA2/Minor_Request/projet_Rmounier_1/demande_RMOUNIER.csv",num_threads = 32,col_names = T) %>%
  column_to_rownames(var="...1") %>%
  dplyr::select(-`...10`)


data = readxl::read_xlsx(path="/home/legrand-lab/Téléchargements/GO_term_summary_20230307_092327.xlsx") %>%
  rename("Annotated_Term"="Annotated Term") %>% 
  mutate(Annotated_Term = factor(Annotated_Term,
                                 levels=unique(Annotated_Term),
                                 labels=c("ESPM","CS","ESAPM","ECESPM","VSC")))

data_f = filter(data, Symbol %in% rownames(integrated)) 
p1 = data_f %>% 
  ggplot() +
  geom_bar(aes(x=Annotated_Term,..count../sum(..count..),fill=Evidence)) +
  labs(x="Annotated Term",y="Relative proportion") +
  scale_y_continuous(limits=c(0,1),n.breaks=6,labels = c("0"="0%","0.2"="20%","0.4"="40%","0.6"="60%","0.8"="80%","1"="100%")) +
  theme_legrand2 +
  theme(axis.text.x = element_text(angle=45,hjust = 1))
p1


make_references = function(dataset, evidence) {
  tmp_1 = filter(dataset, Evidence == evidence) %>% dplyr::select(`MGI Gene/Marker ID`,Annotated_Term,Symbol,Evidence) %>% unique() 
  tmp_2 = pivot_wider(data = tmp_1, names_from = "Annotated_Term",values_from="Symbol")
  tmp_3 = tmp_1 %>% group_by(Symbol,`MGI Gene/Marker ID`) %>% summarise(pop=n())
  
  return(inner_join(tmp_2,tmp_3))
}
IDA = make_references(data_f,"IDA")

# ida = read_delim("/home/legrand-lab/Bureau/references/Mmu/biomarkers/MGImarkerQuery_20210723_101150.txt")
integrated = qs::qread(file="/home/legrand-lab/Bureau/workflow/MNA2/scMdx/MacrophagesInMuscles/integratedData_withoutvisium_v6.1.qs",nthreads = 32)

Idents(integrated) = "group"
dea_conditionnelle = FindAllMarkers(integrated,features=IDA$Symbol %>% unique(),logfc.threshold = -Inf,test.use="MAST",min.pct=0)

qsave(dea_conditionnelle, file="/home/legrand-lab/Bureau/demande_antionio/IDAonscMacrophages_v3.qs",nthreads = 32)
unique(dea_conditionnelle$cluster)
dea_slides_top2 = dea_conditionnelle %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) >=1) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n=2) %>%
  separate(col=cluster,sep="_",into=c("Original_Louvain","a","b","c","d","e","f")) %>%
  mutate(sample=ifelse(a=="MT",paste(a,b,c,d,e,f,sep="_"),a),
         diff = pct.1 - pct.2) %>%
  dplyr::select(-a,-b,-c,-d,-e,-f) %>%
  filter(abs(diff)>=0.25)

Idents(integrated)="sample"
DotPlot(integrated,features=dea_slides_top2$gene %>% unique(),cols = c("blue","yellow","red")) + theme(axis.text.x = element_text(angle=45,hjust = 1))

MIN=min(na.omit(dea_slides_top2$avg_log2FC))
MID=quantile(na.omit(dea_slides_top2$avg_log2FC),.5)
MAX=max(na.omit(dea_slides_top2$avg_log2FC))
library(ggh4x)
stripx <- strip_themed(background_x = elem_list_rect(fill = rainbow(2)))

p1 = make_plots(integrated,groupby = "Original_Louvain")

library(latex2exp)
dea_slides_top2 = dea_slides_top2 %>%
  mutate(newsample =  factor(sample, 
                             levels=c("B6","BM","BW","D2","DM","DW","MT_CD45_CD11b_fib_MDX_S2]","MT_CD45_CD11b_old_MDX_S3]","MT_CD45_CD11b_young_MDX_S1]","MT_CD45_CD11b_young_MDX_S2]"),
                             labels=TeX(c("$C57BL/6J-Dmd^{mdx} 12w$ ",
                                          "$C57BL/6J-Dmd^{mdx} 8w$ ",
                                          "$C57BL/6J-Dmd^{wt} 8w$ ",
                                          "$DBA/2-Dmd^{mdx}$ 12w",
                                          "$DBA/2-Dmd^{mdx}$ 8w",
                                          "$DBA/2-Dmd^{wt}$ 8w",
                                          "$C57BL/6J-Dmd^{mdx} Fib 8w$ ",
                                          "$C57BL/6J-Dmd^{mdx} old$ ",
                                          "$C57BL/6J-Dmd^{mdx} young 2$ ",
                                          "$C57BL/6J-Dmd^{mdx} young 1$ "))
                             ))
  
l = lapply(p1$palettes$sample, function(x) element_rect(fill=x))
dotplot = ggplot(dea_slides_top2, aes(x=gene,y=1,fill=avg_log2FC,size=pct.1*100)) + 
  geom_point(shape=21) +
  scale_size_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_gradientn(colours=c("blue","white","red"),
                        na.value="grey80",
                        # breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                        values=scales::rescale(c(MIN,0,MAX)),limits=c(MIN,MAX)) +
  labs(x="",y="",color="Average Log2 Fold Change",size="% expression in cluster") +
  theme_bw() +
  facet_grid2(newsample~Original_Louvain,scales = "free",space = "free",shrink = T,switch = "both",
              labeller = label_parsed,
              strip = strip_themed(
                # background_x = list(element_rect(fill = "magenta"),
                #                     element_rect(fill = "orange")),
                background_y = l)) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,family = "Helvetica",size=11,face = "bold.italic"),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.x.bottom = element_text(family = "Helvetica",size=11,face="bold.italic"),
        strip.background = element_blank(),
        # strip.background.x = element_rect(),
        strip.text.y.left = element_text(angle=0,family="Helvetica",size=11,face="bold.italic"),
        legend.position = "bottom")

ggsave(filename="DEA_scMacrophages.pdf",device=cairo_pdf,width=30,height = 7,dpi=1000,path="/home/legrand-lab/Bureau/demande_antionio/",plot = dotplot)

integrated@meta.data$new_clusters = paste0(integrated@meta.data$Original_Louvain,"_",integrated@meta.data$sample)
Idents(integrated)="new_clusters"
dea_conditionnelle = FindAllMarkers(integrated,features=IDA$Symbol %>% unique(),logfc.threshold = -Inf,test.use="MAST",min.pct=0,latent.vars = "condition")
qsave(dea_conditionnelle, file="/home/legrand-lab/Bureau/demande_antionio/IDAonscMacrophages_v4.qs",nthreads = 32)

dea_slides_top2 = dea_conditionnelle %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n=2)

DotPlot(integrated,features=dea_slides_top2$gene %>% unique(),cols = c("blue","yellow","red")) + theme(axis.text.x = element_text(angle=45,hjust = 1))



integrated=ScaleData(integrated,features = dea_slides_top2$gene %>% unique())
mtdt = merge(extract_metadata_reduction_seurat(integrated),
             FetchData(integrated,vars = dea_slides_top2$gene %>% unique(),slot = "scale.data"),
             by=0) %>% 
  pivot_longer(cols=dea_slides_top2$gene %>% unique())

a = mtdt$name %>% unique()

ggplot(filter(mtdt, name %in% a[1:10]),aes(x=name,y=value,fill=Original_Louvain)) + geom_boxplot() + facet_wrap(~group)
