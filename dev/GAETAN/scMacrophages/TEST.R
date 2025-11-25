

scMacrophages =  qs::qread(paste0(configuration_scrnaseq@paths@output_archives,"/integratedData.qs"),nthreads = 32) %>%
  SplitObject(split.by="sample")
scMDX_macrophages = qs::qread("/media/legrand-lab/Stockage/InstitutNeuroMyoGene_DATA/EqLegrand/scRNA-seq_doc_travail/versioning/version_5/Rds_Rdata/Macrophages.qs",nthreads=32)  %>% 
  SplitObject(split.by="sample")

sc_list = c(scMacrophages,scMDX_macrophages)
sc_list_filtered_features = list()
for (cso in names(sc_list)) {
  print(cso)
  DefaultAssay(sc_list[[cso]]) = "RNA"
  nonzero <- GetAssayData(object = sc_list[[cso]], slot = "counts") > 1
  keep_genes <- Matrix::rowSums(nonzero) >= configuration_scrnaseq@parameters@keep_genes_threshold_min
  print(table(keep_genes))
  sc_list_filtered_features[[cso]] <- 
    CreateSeuratObject(counts = GetAssayData(object = sc_list[[cso]], slot = "counts")[keep_genes,])
  
  tmp = data.frame(
    new_nCount_RNA = colSums(x = sc_list_filtered_features[[cso]], slot = "counts"), # nCount_RNA
    new_nFeature_RNA = colSums(x = GetAssayData(object = sc_list_filtered_features[[cso]], slot = "counts") > 0),
    sample = cso)
  sc_list_filtered_features[[cso]]=AddMetaData(sc_list_filtered_features[[cso]],metadata = tmp)
  sc_list_filtered_features[[cso]]@meta.data = dplyr::select(sc_list_filtered_features[[cso]]@meta.data,-orig.ident)
}  

sc_list_filtered_features = lapply(sc_list_filtered_features, function(x) {
  x = NormalizeData(x)
  x = FindVariableFeatures(x, 
                           selection.method=configuration_scrnaseq@parameters@selection.method_value,
                           nfeatures = configuration_scrnaseq@parameters@nfeatures_value)
})

genes = lapply(sc_list_filtered_features, function(x) rownames(x))
upset(fromList(genes),nintersects = NA,nsets = 8)

features = SelectIntegrationFeatures(object.list = sc_list_filtered_features)
anchors = FindIntegrationAnchors(object.list = sc_list_filtered_features, anchor.features = features)
nbcells = lapply(sc_list_filtered_features, function(x) dim(x)) 
nbcells = do.call("rbind",nbcells) %>% data.frame()
colnames(nbcells) = c("features","cells")
k.weight = ifelse(min(nbcells$cells)<100,50,100)
integrated = IntegrateData(anchorset = anchors,k.weight = k.weight)  

################ 
integrated = ScaleData(integrated)
integrated = RunPCA(object = integrated, seed.use=42,verbose = T,weight.by.var = T,npcs = 50)
EP_pca = elbowplot(integrated,"pca")
integrated <- RunHarmony(integrated, 
                         group.by.vars = "sample", 
                         assay.use = 'integrated',
                         plot_convergence = TRUE, 
                         dims.use = 1:EP_pca$lastpoint_pctvarchangeto0.1)
EP_harmony = elbowplot(integrated,"harmony")
integrated = reduction_plot(integrated,
                            dims = EP_harmony$lastpoint_pctvarchangeto0.1,
                            ncells = configuration_scrnaseq@parameters@ncells,
                            reduction = "harmony")
# configuration_scrnaseq@parameters@resolution = found_resolution(integrated,EP_harmony)
configuration_scrnaseq@parameters@resolution = 0.51
integrated = neighbors_clustering(data = integrated, 
                                  dims = EP_harmony$lastpoint_pctvarchangeto0.1, 
                                  reduction = "harmony", 
                                  resolution = 0.51)

scMDX = qs::qread("/media/legrand-lab/Stockage/InstitutNeuroMyoGene_DATA/EqLegrand/scRNA-seq_doc_travail/versioning/version_5/Rds_Rdata/muscle.integrated_harmonycorrection_annotated.qs",nthreads=32)
tmp =  data.frame(Subtype=scMDX$Subtype)
integrated = AddMetaData(integrated, metadata=tmp)

deg = FindAllMarkers(object = integrated,
                     assay = "RNA",
                     # features =  c("Folr2","Lyve1","Lyz1","Retnla","Thbs1","Plac8","Bcl2a1b","Tnf","Spp1","Gpnmb","Il1b","S100a9","S100a8","Top2a","Pclaf","Igfbp7","Sparc"),
                     logfc.threshold = -Inf,
                     test.use = "MAST",
                     min.pct = 0,only.pos = F,
                     latent.vars = "sample",densify = T)
f = c("Folr2","Lyve1","Lyz1","Retnla","Thbs1","Plac8","Bcl2a1b","Tnf","Spp1","Gpnmb","Il1b","S100a9","S100a8","Top2a","Pclaf","Igfbp7","Sparc")
deg2 = deg %>% mutate(gene=factor(gene,levels=f))

MIN=min(na.omit(deg2$avg_log2FC))
MID=quantile(na.omit(deg2$avg_log2FC),.5)
MAX=max(na.omit(deg2$avg_log2FC))
# col_scales_scores_max = deg2$col_scales_scores[deg2$avg_log2FC == MAX] %>% unique()


data = dplyr::select(deg2, avg_log2FC, gene, cluster) %>% 
  pivot_wider(id_cols="cluster",names_from="gene",values_from="avg_log2FC") %>%
  column_to_rownames(var="cluster")
features =  c("Folr2","Lyve1","Lyz1","Retnla","Thbs1","Plac8","Bcl2a1b","Tnf","Spp1","Gpnmb","Il1b","S100a9","S100a8","Top2a","Pclaf","Igfbp7","Sparc")
data = dplyr::select(data,features)

integrated = ScaleData(integrated,assay="RNA")
DefaultAssay(integrated)="RNA"
fetchdata = FetchData(integrated,vars =c("Ptprc","Itgam","Lyve1","H2-Ab1","H2-Eb1"),slot = "scale.data") %>%
  rownames_to_column(var="barcode") %>%
  pivot_longer(cols=c("Ptprc","Itgam","Lyve1","H2-Ab1","H2-Eb1"), names_to="gene",values_to="scale_data") %>%
  mutate(scale_data=ifelse(scale_data<0,0,scale_data)) %>%
  pivot_wider(id_cols = "barcode",names_from = "gene",values_from="scale_data") %>%
  column_to_rownames(var="barcode")

library(ComplexHeatmap)
tmp = FetchData(integrated,vars =c("Ptprc","Itgam","Lyve1","H2-Ab1","H2-Eb1"),slot = "scale.data") %>%
  rownames_to_column(var="barcode") %>%
  pivot_longer(cols=c("Ptprc","Itgam","Lyve1","H2-Ab1","H2-Eb1"), names_to="gene",values_to="scale_data") %>%
  data.frame()
tmp2 = dplyr::select(integrated@meta.data,sample,seurat_clusters) %>%
  rownames_to_column(var="barcode")

tmp = inner_join(tmp,tmp2,by="barcode")

ggplot(filter(tmp,scale_data>0), aes(x=seurat_clusters,y=scale_data,color=gene)) +
  geom_boxplot() +
  facet_grid(~sample)

Heatmap(tmp,
        show_row_dend = T,show_column_dend = T,
        show_row_names = F,show_column_names = T,
        row_km = 3,
        column_km = 3,use_raster = F)



mtdt = extract_metadata_reduction_seurat(integrated) %>% 
  merge(fetchdata,by=0) %>% 
  column_to_rownames(var="Row.names")

# Define the dot product function
dot <- function(v1, v2) {
  return(sum(v1 * v2))
}
calculate_project_ortho=function(x,y,type) {
  point=c(x,y)
  line=c(1,1)
  
  projection <- (dot(point, line) / dot(line, line)) * line
  distancePoint2Diag <- sqrt(sum((point - projection)^2))
  distancePoint2Orig <- sqrt(sum((projection - line)^2))
  
  if (type == "D") {
    return(distancePoint2Diag)
  } else if (type == "O") {
    return(distancePoint2Orig)
  }
}

xy = dplyr::select(mtdt,sample,FItSNE_1,FItSNE_2,Ptprc,Itgam) %>%
  mutate(Ptprc = ifelse(Ptprc<0,0,Ptprc),
         Itgam = ifelse(Itgam<0,0,Itgam)) %>%
  rowwise() %>%
  mutate(distD=calculate_project_ortho(Ptprc,Itgam,"D"),
         sens=ifelse(Ptprc<Itgam,-1,1),
         rounddist=round(distD*10),
         distO=calculate_project_ortho(Ptprc,Itgam,"O")) %>%
  mutate(distsens=rounddist*sens) %>%
  data.frame()
inpcol = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple")
cInp = strsplit(inpcol, "; ")[[2]]
if(cInp[1] == "Red (Gene1)"){
  c10 = c(255,0,0)
} else if(cInp[1] == "Orange (Gene1)"){
  c10 = c(255,140,0)
} else {
  c10 = c(0,255,0)
}
c01 = c(0,0,255)
c00=c(217,217,217); c11=c10+c01; nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2
gg = data.table::data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1)))
gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid
gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22) {
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22
  oup = oup / (xy*xy)
  return(oup)
}
gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
gg = gg[, c("v1", "v2", "cMix")]

xy$v1 = round(nTot*xy$Ptprc/max(xy$Ptprc))
xy$v2 = round(nTot*xy$Itgam/max(xy$Itgam))
xy$v0 = xy$v1 + xy$v2
xy = gg[xy, on = c("v1", "v2")]
xy=xy[order(v0)]
xy_f = dplyr::select(xy,FItSNE_1,FItSNE_2)
xy$state = xy$distsens*xy$distO
p1=ggplot(xy,aes(x=FItSNE_1,y=FItSNE_2)) +
  geom_point(data=xy_f, color="snow2") +
  geom_point(color=xy$cMix) +
  scale_color_gradientn(cInp, colours=cList[[1]]) +
  facet_grid(distsens<0*distO~sample) +
  theme_bw()

BiocManager::install("changepoint")
library(changepoint)
res = c(cp_ptprc=cpt.var(xy$Ptprc)@param.est$mean,
        cp_itgam=cpt.var(xy$Itgam)@param.est$mean,
        cp_dist0_high=cpt.meanvar(xy$distO)@param.est$mean[1],
        cp_dist0_low=cpt.meanvar(xy$distO)@param.est$mean[2])

xy = xy %>%
  mutate(stateGene=case_when(
    Ptprc <= res$cp_ptprc & Itgam <= res$cp_itgam & dist0 <= cp_dist0_low~ "Low-Low",
    Ptprc > res$cp_ptprc & Itgam > res$cp_itgam & dist0 > cp_dist0_low~ "Low-Low",
    Ptprc <= res$cp_ptprc & Itgam <= res$cp_itgam & dist0 <= cp_dist0_low~ "Low-Low",
    Ptprc <= res$cp_ptprc & Itgam <= res$cp_itgam & dist0 <= cp_dist0_low~ "Low-Low",
  ))


p2=ggplot(gg,aes(v1,v2)) +
  geom_tile(fill=gg$cMix)

p3=ggplot(xy,aes(Ptprc,Itgam)) +
  # geom_tile(fill=gg$cMix) +
  geom_jitter() +
  geom_hline(yintercept = c(cpt.var(xy$Ptprc)@param.est$mean,
                            cpt.meanvar(xy$distO)@param.est$mean[1]),linetype="dashed",color=c("red","green")) +
  geom_vline(xintercept = c(cpt.var(xy$Ptprc)@param.est$mean,
                            cpt.meanvar(xy$distO)@param.est$mean[2]),linetype="dashed",color=c("red","green")) +
  theme_bw()

p1|p2

MIN=min(na.omit(xy$distO))
MID=quantile(na.omit(xy$distO),.5)
MAX=max(na.omit(xy$distO))
xy=xy[order(xy$distO),]
ggplot(arrange(xy,abs(distO)),aes(x=FItSNE_1,y=FItSNE_2,color=distO,alpha=distsens)) +
  geom_point() +
  scale_color_gradientn(colours=c("blue","orange","red"),
                        na.value="grey80",
                        breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                        values=scales::rescale(c(MIN,MID,MAX)),limits=round(c(MIN,MAX),1)) +
  theme_bw() +
  facet_grid(distsens<0~sample)

ggplot(xy,aes(distsens)) + geom_density()



library(ggplot2)

# Define the point and the line
point <- c(3, 4)
line <- c(1, 1)



# Calculate the projection of the point onto the line
projection <- (dot(point, line) / dot(line, line)) * line
# Calculate the distance between the point and its projection
distance <- sqrt(sum((point - projection)^2))

# Create a data frame for plotting
data <- data.frame(x = c(point[1], projection[1]), y = c(point[2], projection[2]), label = c("Point", "Projection"))

# Plot the point, the line, and the projection
ggplot(data, aes(x, y, label = label)) +
  geom_point() +
  geom_segment(aes(x = 0, y = 0, xend = 4, yend = 4), color = "red") +
  ggtitle(paste("Distance: ", round(distance, 2)))




mtdt$score_mc = ((mtdt$Ptprc-mtdt$Itgam)/(mtdt$Ptprc+mtdt$Itgam))+mtdt$Ptprc*mtdt$Itgam

MIN=min(na.omit(mtdt$score_mc))
MID=quantile(na.omit(mtdt$score_mc),.5)
MAX=max(na.omit(mtdt$score_mc))
mtdt = mtdt %>%
  mutate(stateMC = case_when(
    cut_score_mc == "(-3.97,-1.13]" ~ "Igam+",
    cut_score_mc == "(-1.13,1.7]" ~ "Both",
    cut_score_mc == "(1.7,4.53]" ~ "Ptprc+"
  ))
p1=ggplot(mtdt,aes(x=FItSNE_1,y=FItSNE_2)) +
  geom_point(data = dplyr::select(mtdt, -sample,-score_mc), color = "grey70", alpha=alphaP,size=sizeP) +
  geom_point(aes(color = score_mc),alpha=1, show.legend = T) +
  scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1", type = "discrete"),
                        na.value="grey80",
                        breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                        values=scales::rescale(c(MIN,MID,MAX)),limits=round(c(MIN,MAX),1)) +
  theme_bw() +
  # scale_color_manual(values=plots$palettes$sample)
  facet_wrap(~sample)
p2=ggplot(mtdt,aes(x=FItSNE_1,y=FItSNE_2)) +
  geom_point(data = dplyr::select(mtdt, -sample,-score_mc), color = "grey70", alpha=alphaP,size=sizeP) +
  geom_point(aes(color = stateMC),alpha=1, show.legend = T) +
  scale_color_manual(values = c())
  # scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1", type = "discrete"),
  #                       na.value="grey80",
  #                       breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
  #                       values=scales::rescale(c(MIN,MID,MAX)),
  #                       limits=round(c(MIN,MAX),1)) +
  theme_bw() +
  # scale_color_manual(values=plots$palettes$sample)
  facet_grid(stateMC~sample)

  
  p1=ggplot(mtdt,aes(x=Ptprc)) + geom_density()
  p2=ggplot(mtdt,aes(x=Itgam)) + geom_density()
  p1/p2

library(ComplexHeatmap)
Heatmap(matrix = data,
        cluster_rows = T,
        cluster_columns = F)

dend = as.dendrogram(hclust(dist(data)))
plot(dend)
dend <- as.dendrogram(hclust(dist(with(deg2, tapply(gene, cluster)))))
deg2$new_cluster <- factor(deg2$cluster, levels = labels(dend))
# Stacked Percent
library(dendextend)
p1 <- ggplot(dend, horiz = T) 
p2 <- ggplot(deg, aes(size=pct.1*100, y=cluster, x=gene)) + 
  geom_point()

library(cowplot)
plot_grid(p1, p2, align = "h")


ggplot(deg2, aes(x=gene,y=cluster,size=pct.1*100,color=avg_log2FC)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_gradientn(colours=wesanderson::wes_palette("Zissou1", type = "discrete"),
                        na.value="grey80",
                        breaks=round(as.vector(sort(c(MIN,MID,MAX))),1),
                        values=scales::rescale(c(MIN,MID,MAX)),limits=round(c(MIN,MAX),1))
head(deg)
