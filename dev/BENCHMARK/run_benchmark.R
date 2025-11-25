
############# FUNCTIONS
library(SPOTlight)
# library(BayesPrism)
library(patchwork)
library(scuttle)
library(SpatialFeatureExperiment)
library(org.Mm.eg.db)
library(EBImage)
library(terra)
library(sf)
source("/home/mna_bioinfo/Bureau/development/github/single_cell_RNA_seq/functions_v2.R")
source("/home/mna_bioinfo/Bureau/development/github/single_cell_RNA_seq/in_development/deconvolution/bayesprism_adapted_functions.R")
source("~/Bureau/development/github/single_cell_RNA_seq/in_development/advanced_analysis/final_annotation_process_functions.R")

# sfe <- SFEData::McKellarMuscleData("full")


seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  if (class(seu@images[[img_id]])[1] == "VisiumV2") {
    spatialCoords <- as.matrix(GetTissueCoordinates(seu)[,1:2])
  } else {
    spatialCoords <- as.matrix(
      seu@images[[img_id]]@coordinates[, c("imagerow", "imagecol")])
  }
  
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}
elbowplot_visium <- function(visium, reduction) {
  
  pct = apply(Embeddings(visium,reduction),2,sd)
  
  
  
  # pct <- data_cso[[paste(reduction)]]@stdev / sum(data_cso[[paste(reduction)]]@stdev) *100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  if (is.na(co1)) {
    co1 <- which(cumu > 50 & pct < 5)[1]
  }
  
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  # if (reduction=="pca")  {reduction <- "PCA"} else if (reduction=="harmony") {reduction <- "Harmony"}
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  return(list(
    plot = ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
      geom_text() + 
      labs(x="Cumulative Sums of % Standard Deviation",
           y="% Standard Deviation",
           title=reduction,
           subtitle=paste0("Lim=",co2),
           color="Usable Dimension class") +
      theme(legend.position = "bottom",
            text = element_text(face="bold")),
    lastpoint_pctvarchangeto0.1=co2)
  )
}
raster2polygon <- function(seg, keep = 0.2) {
  r <- rast(as.array(seg)) |> 
    terra::trans() |> flip() |> trans()
  r[r < 1] <- NA
  contours <- st_as_sf(as.polygons(r, dissolve = TRUE))
  simplified <- rmapshaper::ms_simplify(contours, keep = keep)
  list(full = contours,
       simplified = simplified)
}
############# PARAMETERS INPUT
#single cell
outdir="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/result/FABIEN_LEGRAND/serieA/downstream/single_cell_0.5.27052024/"
#spatial
output_path="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/BENCHMARK/"
parameters_output=outdir
nthreads = 20

gse_repertories = setNames(list.dirs("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/fastq/VISIUM/benchmark/",full.names = T,recursive = F)[1:3],
                basename(list.dirs("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/fastq/VISIUM/benchmark/",full.names = T,recursive = F)[1:3])
)

sfe_tissue_final = list()
for (gse in 1:length(gse_repertories)) {
  print(gse)
  outs = setNames(list.dirs(gse_repertories[gse],full.names = T,recursive = F),
                  basename(list.dirs(gse_repertories[gse],full.names = T,recursive = F)))
  
  for (sample in names(outs)) {
    print(sample)
    # if (!file.exists(paste0(outs[[sample]],"/",sample,"_scalefactors_json.json"))) {
    #   print("decompression gz...")
    #   lf = list.files(path=outs[[sample]], pattern=".gz")
    #   lapply(lf, function(x) R.utils::gunzip(paste0(outs[[sample]],"/",x)))
    # }
    print("move files in subdirectories")
    # if (!dir.exists(paste0(outs[[sample]],"/outs"))) {
      # if (gse != 1) {
      #   dir.create(paste0(outs[[sample]],"/outs"))
      #   require(fs)
      #   if (file.exists(paste0(outs[[sample]],"/",sample,"_filtered_feature_bc_matrix.h5"))) {
      #     file_move(paste0(outs[[sample]],"/",sample,"_filtered_feature_bc_matrix.h5"),paste0(outs[[sample]],"/outs/filtered_feature_bc_matrix.h5"))
      #   } else {
      #     print("no h5file")
      #   }
      #   dir.create(paste0(outs[[sample]],"/outs/spatial"))
      #   if (file.exists(paste0(outs[[sample]],"/",sample,"_scalefactors_json.json"))) {
      #     file_move(paste0(outs[[sample]],"/",sample,"_scalefactors_json.json"),paste0(outs[[sample]],"/outs/spatial/scalefactors_json.json"))
      #     file_move(paste0(outs[[sample]],"/",sample,"_tissue_lowres_image.png"),paste0(outs[[sample]],"/outs/spatial/tissue_lowres_image.png"))
      #     file_move(paste0(outs[[sample]],"/",sample,"_tissue_hires_image.png"),paste0(outs[[sample]],"/outs/spatial/tissue_hires_image.png"))
      #     file_move(paste0(outs[[sample]],"/",sample,"_tissue_positions_list.csv"),paste0(outs[[sample]],"/outs/spatial/tissue_positions_list.csv"))
      #   } else {
      #     print("no classical data")
      #     file_move(paste0(outs[[sample]],"/",sample,"_Stitched_Image_0.tif"),paste0(outs[[sample]],"/outs/spatial/Stitched_Image_0.tiff"))
      #   }
      # }
      # if (gse != 3) {
      # file_move(paste0(outs[[sample]],"/",sample, "_metrics_summary.csv"),paste0(outs[[sample]],"/outs/metrics_summary.csv"))
      # file_move(paste0(outs[[sample]],"/",sample, "_molecule_info.h5"),paste0(outs[[sample]],"/outs/molecule_info.h5"))
      # file_move(paste0(outs[[sample]],"/",sample,"_aligned_fiducials.jpg"), paste0(outs[[sample]],"/outs/spatial/aligned_fiducials.jpg"))
      # file_move(paste0(outs[[sample]],"/",sample,"_detected_tissue_image.jpg"),paste0(outs[[sample]],"/outs/spatial/detected_tissue_image.jpg"))
      # 
      # dir.create(paste0(outs[[sample]],"/outs/filtered_feature_bc_matrix"))
      # file_move(paste0(outs[[sample]],"/",sample, "_barcodes.tsv"),paste0(outs[[sample]],"/outs/filtered_feature_bc_matrix/barcodes.tsv"))
      # file_move(paste0(outs[[sample]],"/",sample, "_features.tsv"),paste0(outs[[sample]],"/outs/filtered_feature_bc_matrix/features.tsv"))
      # file_move(paste0(outs[[sample]],"/",sample, "_matrix.mtx"),paste0(outs[[sample]],"/outs/filtered_feature_bc_matrix/matrix.mtx"))
      # 
    
      if (gse == 3) {
        print("3")
        sfe <- read10xVisiumSFE(dirs=paste0(outs[[sample]],"/outs"),
                                sample_id = sample,
                                type="HDF5",
                                images = "hires",
                                unit = "full_res_image_pixel",zero.policy = T,
                                data="filtered")
      } 
      if (gse %in% 1) {
        print("1")
        sfe <- read10xVisiumSFE(dirs=paste0(outs[[sample]],"/outs"),
                                sample_id = sample,
                                type="HDF5",
                                data="filtered",
                                load=F)
      }
      if (gse == 2) {
        print("2")
        sfe <- read10xVisiumSFE(dirs=paste0(outs[[sample]],"/outs"),
                                sample_id = sample,
                                type="sparse",
                                data="filtered",
                                load=F)
        
      }
    scale_factors = jsonlite::fromJSON(txt = paste0(outs[[sample]],"/outs/spatial/scalefactors_json.json"))
    print("annotation")
    bitr_res = clusterProfiler::bitr(geneID = rownames(sfe),
                                     fromType="ENSEMBL",
                                     toType=c("ENTREZID","SYMBOL","GENETYPE"),
                                     OrgDb = org.Mm.eg.db,
                                     drop=T)
    sfe = sfe[rownames(sfe) %in% bitr_res$ENSEMBL[bitr_res$GENETYPE == "protein-coding"], ]
    is_mt <- str_detect(rowData(sfe)$symbol, "^mt-")
    sfe <- addPerCellQCMetrics(sfe, subsets = list(mito = is_mt))
    print("image treatment")
    img <- readImage(paste0(outs[[sample]],"/outs/spatial/tissue_hires_image.png"))
    img2 <- img
    colorMode(img2) <- Grayscale
    display(img2)
    img_gray <- normalize(img2)
    thresh <- otsu(img_gray)
    nuc_th = combine( mapply(function(frame, th) frame > th, getFrames(img_gray), thresh, SIMPLIFY=FALSE) )
    display(nuc_th, all=TRUE)
    
    disc = makeBrush(21,shape = "disc",step = F)^2
    disc = disc / sum(disc)
    offset = 0.05
    nuc_bg = filter2( img_gray, disc )
    nuc_th = img_gray > nuc_bg + offset
    display(nuc_th, all=TRUE)
    
    
    img_bin <- img_gray < thresh
    display(img_bin)
    img_fh <- fillHull(img_bin)
    display(img_fh)
    
    mask_clean <- opening(img_fh, makeBrush(3, shape = 'disc'))
    display(mask_clean)
    
    mask_label <- bwlabel(mask_clean)
    display(mask_label)
    
    sizes <- table(mask_label)
    # main_obj <- which.max(sizes[-1]) + 1
    tissue_mask <- mask_label == 0
    display(tissue_mask)
    
    mask_flip <- EBImage::flip(x = tissue_mask)
    mask_fliprot <- EBImage::rotate(mask_flip,angle = 270)
    
    display(mask_fliprot,all=F,frame=1,method="raster")
    
    edge_mask <- dilate(mask_fliprot, makeBrush(3, shape='disc')) & !mask_fliprot
    display(edge_mask,all=F,frame=1,method="raster")
    
    tb <- raster2polygon(mask_fliprot)
    tb$simplified$geometry <- tb$simplified$geometry / scale_factors$tissue_hires_scalef
    tissueBoundary(sfe) <- tb$simplified
    sfe$int_tissue <- annotPred(sfe, colGeometryName = "spotPoly", 
                                annotGeometryName = "tissueBoundary",
                                pred = st_intersects)
    sfe$cov_tissue <- annotPred(sfe, colGeometryName = "spotPoly", 
                                annotGeometryName = "tissueBoundary",
                                pred = st_covered_by)
    sfe$diff_sr <- case_when(sfe$in_tissue == sfe$int_tissue ~ "same",
                             sfe$in_tissue & !sfe$int_tissue ~ "Space Ranger",
                             sfe$int_tissue & !sfe$in_tissue ~ "segmentation") |> 
      factor(levels = c("Space Ranger", "same", "segmentation"))
    sfe$diff_int_cov <- sfe$int_tissue != sfe$cov_tissue
    spot_ints <- annotOp(sfe, colGeometryName = "spotPoly", annotGeometryName = "tissueBoundary", op = st_intersection)
    sfe$pct_tissue <- st_area(spot_ints) / st_area(spotPoly(sfe)) * 100
    
    colData(sfe)$notgoodspot = ifelse(colData(sfe)$detected<=130,F,T)
    colData(sfe)$verif = factor(paste0(colData(sfe)$cov_tissue,"_",colData(sfe)$notgoodspot),
                                levels=c("TRUE_TRUE","FALSE_TRUE","TRUE_FALSE","FALSE_FALSE")
    )
    colData(sfe)$verif = ifelse(colData(sfe)$verif %in% c("TRUE_TRUE","FALSE_TRUE"),T,F)
    colData(sfe)$allpass = colData(sfe)$verif
    sfe_tissue <- sfe[,sfe$allpass]
    
    # plotSpatialFeature(sfe_tissue, "allpass",
    #                    annotGeometryName = "tissueBoundary",
    #                    annot_fixed = list(fill = NA, size = 0.5, color = "black"))
    
    colGraph(sfe_tissue, "visium") <- findVisiumGraph(sfe_tissue,sample_id = sample,zero.policy = T)
    sfe_tissue_final[[sample]] = sfe_tissue

  }
}
dir.create(paste0(output_path,"/downstream"))
qs::qsave(x = sfe_tissue_final,file = paste0(output_path,"/downstream/segmentation_done_v2.qs"))

sfe_tissue_final = qs::qread(file = paste0(output_path,"/downstream/segmentation_done_v2.qs"))
# source("~/Bureau/development/github/visium_rnaseq/VisiumSD/3_spatialPCA.R")
# SpatialPCA=qs::qread(file = paste0(output_path,"/downstream/spatial_pca.qs"))

# outs = c(setNames(list.dirs(gse_repertories[[1]],full.names = T,recursive = F),
#                   basename(list.dirs(gse_repertories[[1]],full.names = T,recursive = F))),
#          setNames(list.dirs(gse_repertories[[2]],full.names = T,recursive = F),
#                   basename(list.dirs(gse_repertories[[2]],full.names = T,recursive = F)))
# )
outs = c(setNames(list.dirs(gse_repertories,full.names = T,recursive = F),
                  basename(list.dirs(gse_repertories,full.names = T,recursive = F)))
)

sfe_tissue_seurat = list()
for (sample in names(sfe_tissue_final)) {
  print(sample)
  cso = Load10X_Spatial(data.dir = paste0(outs[[sample]],"/outs"), #"/media/mna_bioinfo/MNA2_Stockage/VISIUM/B6lineage_singlecell/D-GA-B6mdx/outs",
                        assay = "spatial",
                        filter.matrix = F,
                        to.upper = F)
  
  library(magrittr)
  cso_f = subset(cso, features = rowData(sfe_tissue_final[[sample]])$symbol)
  cso_f = subset(cso_f, cells = colnames(sfe_tissue_final[[sample]]),Update.slots = T,Update.object = T)
  cso_f = SCTransform(cso_f, assay = "spatial", return.only.var.genes = TRUE, variable.features.rv.th = 1.3,method="glmGamPoi")
  cso_f = RunPCA(cso_f)
  # cso_f@reductions$SPCA = CreateDimReducObject(embeddings = SpatialPCA[[sample]][colnames(cso_f),],key = "SPCA_", assay = "SCT_vst")
  pca_elbow = elbowplot_visium(cso_f,"pca")
  # spca_elbow = elbowplot_visium(cso_f,"SPCA")
  cso_f <- FindNeighbors(cso_f, reduction = "pca", dims = 1:pca_elbow$lastpoint_pctvarchangeto0.1,k.param = 30,graph.name = "knn_pca")
  # cso_f <- FindNeighbors(cso_f, reduction = "SPCA", dims = 1:spca_elbow$lastpoint_pctvarchangeto0.1,k.param = 30,graph.name = "knn_spca")
  
  for (x in seq(0.1,2,0.1)) {
    a = FindClusters(cso_f, verbose = FALSE,graph.name = "knn_pca",cluster.name = paste0("pca_leiden_",x),method = 4,resolution = x)
    cso_f[[paste0("pca_leiden_",x)]] = a[[paste0("pca_leiden_",x)]]
    # b = FindClusters(cso_f, verbose = FALSE,graph.name = "knn_spca",cluster.name = paste0("spca_leiden_",x),method = 4,resolution = x)
    # cso_f[[paste0("spca_leiden_",x)]] = b[[paste0("spca_leiden_",x)]]
  }
  # cso_f <- RunUMAP(cso_f, reduction = "SPCA", metric="euclidean", dims = 1:spca_elbow$lastpoint_pctvarchangeto0.1,reduction.key="UMAPSPCA_",reduction.name="UMAP_SPCA")
  cso_f <- RunUMAP(cso_f, reduction = "pca",  metric="euclidean", dims = 1:pca_elbow$lastpoint_pctvarchangeto0.1,reduction.key="UMAPPCA_",reduction.name="UMAP_PCA")
  
  # DimPlot(cso_f,group.by = c("pca_leiden_1","spca_leiden_1")) /
  #   SpatialDimPlot(cso_f,group.by = c("pca_leiden_1","spca_leiden_1"))
  
  
  sfe_tissue_seurat[[sample]] = cso_f
}

qs::qsave(x=sfe_tissue_seurat, paste0(output_path,"/downstream/spatial_data.qs"))

sfe_tissue_seurat = qs::qread(paste0(output_path,"/downstream/spatial_data.qs"))

library(scran)
library(BiocParallel)
nthreads = 20
############# FUNCTIONS
library(SPOTlight)
# library(BayesPrism)
library(InstaPrism)
library(patchwork)
# # source("/home/mna_bioinfo/Bureau/development/github/single_cell_RNA_seq/functions_v2.R")
# source("/home/mna_bioinfo/Bureau/development/github/single_cell_RNA_seq/in_development/deconvolution/bayesprism_adapted_functions.R")
# source("~/Bureau/development/github/single_cell_RNA_seq/in_development/advanced_analysis/final_annotation_process_functions.R")
getLocalMinMax <- function(x){
  
  inflect <- function(x, threshold = 1){
    up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
    down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
    a    <- cbind(x,up,down)
    list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  }  
  
  n <- 2
  d <- density(x)
  bottoms <- lapply(1:n, function(x) inflect(d$y, threshold = x)$minima)
  tops <- lapply(1:n, function(x) inflect(d$y, threshold = x)$maxima)
  
  cf.1 <- grDevices::colorRampPalette(c("pink","red"))
  cf.2 <- grDevices::colorRampPalette(c("cyan","blue"))
  
  MIN <- unlist(bottoms) %>% unique()
  if(isEmpty(MIN)) {MIN <- 1}
  MAX <- unlist(tops) %>% unique()
  
  DATA <- data.frame(index = MIN, x=d$x[MIN], y=d$y[MIN], val="minimum") %>%
    rbind(data.frame(index = MAX, x=d$x[MAX], y=d$y[MAX], val="maximum")) %>% arrange(index)
  
  DATA <- DATA %>% arrange(x)
  
  return(DATA)
}
seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  
  # seurat_to_spe(seu = sfe_tissue_seurat[[sample]],sample_id = "S1",img_id = "S1")
  # seu = sfe_tissue_seurat[[sample]]
  # sample_id = "S1"
  # img_id = "S1"
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  spatialCoords = GetTissueCoordinates(seu)[,1:2] %>% as.matrix()
  colnames(spatialCoords) = c("imagerow","imagecol")
  
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  spe
}

refPhi_obj_hypergroup=list()
refPhi_obj_subcluster=list()
for (comparison in c("lineage_b6","lineage_d2")) {
  refPhi_obj_hypergroup[[comparison]] = qs::qread(file = paste0("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/deconvolution/hypergroup/requestJanuary2025/modeleprism_",comparison,".qs"))
}
subset = list()
spatial_obj = list()
deconvolute_matrix_m = list()
deconvolute_matrix = list()
p=list()
for (comparison in c("lineage_b6","lineage_d2")) {
  print(paste0("* ",comparison))
  spatial_obj[[comparison]] = list()
  if (comparison == "lineage_b6") {
    data_space = list(GSM5979972_C57BL10=GetAssayData(sfe_tissue_seurat$GSM5979972_C57BL10 ,assay="spatial") %>% t(),  
                      GSM5979974_MDX=GetAssayData(sfe_tissue_seurat$GSM5979974_MDX ,assay="spatial") %>% t(),
                      Vis5A_2dpi=GetAssayData(sfe_tissue_seurat$Vis5A_2dpi ,assay="spatial") %>% t(),
                      Vis7B_5dpi=GetAssayData(sfe_tissue_seurat$Vis7B_5dpi ,assay="spatial") %>% t(),
                      Vis9A_7dpi=GetAssayData(sfe_tissue_seurat$Vis9A_7dpi ,assay="spatial") %>% t() )
  } else {
    data_space = list(GSM5979973_DBA2J=GetAssayData(sfe_tissue_seurat$GSM5979973_DBA2J ,assay="spatial") %>% t(),  
                      GSM5979975_D2MDX=GetAssayData(sfe_tissue_seurat$GSM5979975_D2MDX ,assay="spatial") %>% t(),
                      GSM7055901_WT1=GetAssayData(sfe_tissue_seurat$GSM7055901_WT1 ,assay="spatial") %>% t(),
                      GSM7055902_WT2=GetAssayData(sfe_tissue_seurat$GSM7055902_WT2 ,assay="spatial") %>% t(),
                      GSM7055903_MDX1=GetAssayData(sfe_tissue_seurat$GSM7055903_MDX1 ,assay="spatial") %>% t(),
                      GSM7055904_MDX2=GetAssayData(sfe_tissue_seurat$GSM7055904_MDX2 ,assay="spatial") %>% t(),
                      GSM7055905_WT3=GetAssayData(sfe_tissue_seurat$GSM7055905_WT3 ,assay="spatial") %>% t(),
                      GSM7055906_MDX3=GetAssayData(sfe_tissue_seurat$GSM7055906_MDX3 ,assay="spatial") %>% t(),
                      GSM7055907_MDX4=GetAssayData(sfe_tissue_seurat$GSM7055907_MDX4 ,assay="spatial") %>% t(),
                      GSM7055908_MDX5=GetAssayData(sfe_tissue_seurat$GSM7055908_MDX5 ,assay="spatial") %>% t(),
                      GSM7230943_Day1_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230943_Day1_PostInjury ,assay="spatial") %>% t(),
                      GSM7230944_Day3_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230944_Day3_PostInjury ,assay="spatial") %>% t(),
                      GSM7230945_Day5_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230945_Day5_PostInjury ,assay="spatial") %>% t())
  }
  for (sample in names(data_space)) {
    set.seed(42)
    
    if (!file.exists(paste0(output_path,"/downstream/result_instaprism_",sample,".qs"))) {
      print(" go instaprism")
      spatial_obj[[comparison]][[sample]] = InstaPrism(bulk_Expr = data_space[[sample]] %>% t() %>% as.matrix(),
                                                       refPhi_cs = refPhi_obj_hypergroup[[comparison]][[1]])
      thetaH = spatial_obj[[comparison]][[sample]]@Post.ini.ct@theta %>% t() %>% data.frame()
      
      qs::qsave(spatial_obj[[comparison]][[sample]],file = paste0(output_path,"/downstream/result_instaprism_",sample,".qs"))
    } else {
      print(" read qs")
      spatial_obj[[comparison]][[sample]] = qs::qread(paste0(output_path,"/downstream/result_instaprism_",sample,".qs"))
      thetaH = spatial_obj[[comparison]][[sample]]@Post.ini.ct@theta %>% t() %>% data.frame()
    }
    
    print("fix threshold")
    thetaH_m = thetaH
    thetaH_m[thetaH_m < 0.1] <- 0
    
    # Z = get_Z_array(thetaH)
    
    tmp = list()
    for (col in colnames(thetaH_m)) {
      if (col == "MuSCs") {
        tmp[[col]] = scater::isOutlier(metric= thetaH_m[,col],type="lower") %>%
          as.character() %>%
          factor(levels=c("TRUE","FALSE"),labels=c("nokeep","keep")) %>%
          as.character()
      } else {
        # print(col)
        tmp[[col]] = scater::isOutlier(metric= thetaH_m[,col],type="higher") %>%
          as.character() %>%
          factor(levels=c("TRUE","FALSE"),labels=c("keep","nokeep")) %>%
          as.character()
      }
      print("lower")
      print(scater::isOutlier(metric= thetaH_m[,col],type="lower"))
      print("higher")
      print(scater::isOutlier(metric= thetaH_m[,col],type="higher"))
    }
    tmp = do.call("cbind",tmp) %>% data.frame()
    colnames(tmp)=paste0("keep_",colnames(tmp))
    rownames(tmp)=rownames(thetaH_m)
    thetaH_m = merge(thetaH_m,tmp,by=0) %>% column_to_rownames(var="Row.names")
    
    sfe_tissue_seurat[[sample]] = AddMetaData(sfe_tissue_seurat[[sample]], thetaH_m)
    
    subset[[sample]] = list()
    p[[sample]] = list()
    for (hypergroup in colnames(thetaH)[1:7]) {
      #for subset
      if (length(rownames(sfe_tissue_seurat[[sample]]@meta.data[sfe_tissue_seurat[[sample]][[paste0("keep_",hypergroup)]]=="keep",])) != 0) {
        
      
        subset[[sample]][[hypergroup]] = subset(sfe_tissue_seurat[[sample]],
                                                 cells = rownames(sfe_tissue_seurat[[sample]]@meta.data[sfe_tissue_seurat[[sample]][[paste0("keep_",hypergroup)]]=="keep",]))
        subset[[sample]][[hypergroup]]@images = sfe_tissue_seurat[[sample]]@images
        
        if (length(unique(c(0.1,round(max(subset[[sample]][[hypergroup]][[hypergroup]]),1))))==1) {
          custom_breaks = c(0.1,0.2)
        } else  {
          custom_breaks = c(0.1,round(max(subset[[sample]][[hypergroup]][[hypergroup]]),1))
        }
        p[[sample]][[hypergroup]] = SpatialFeaturePlot( subset[[sample]][[hypergroup]], 
                                                        features=hypergroup,pt.size.factor = 2) +
          if (custom_breaks[[2]] != 0) {
            scale_fill_gradient2(low = "white", high = "darkred",
                               # color_hypergroup[[hypergroup]],
                               midpoint = .1,
                               breaks=seq(custom_breaks[[1]],custom_breaks[[2]],0.1),
                               # labels=c(0.1,1),
                               na.value="white",
                               limits=custom_breaks)}
        p[[sample]][[hypergroup]] = p[[sample]][[hypergroup]] +
          SpatialDimPlot(sfe_tissue_seurat[[sample]],group.by=paste0("keep_",hypergroup),pt.size.factor = 4) +
          ggplot(sfe_tissue_seurat[[sample]]@meta.data,aes(x=.data[[paste0("keep_",hypergroup)]],y=.data[[hypergroup]])) + 
          geom_boxplot(outlier.size = -1) + 
          geom_jitter(shape=1) +
          geom_hline(yintercept = 0.1,linetype="dashed",color="red") +
          theme_bw() +
          ggplot(sfe_tissue_seurat[[sample]]@meta.data, aes(x=.data[[hypergroup]],fill=.data[[paste0("keep_",hypergroup)]])) + 
          geom_histogram(bins=100,alpha=0.7) +
          geom_density(alpha=0.5) +
          geom_vline(xintercept = 0.1,linetype="dashed",color="red") +
          labs(x=paste0("Deconvolute Score ",hypergroup)) +
          theme_bw()
      } else {
        p[[sample]][[hypergroup]] = NULL
      }
    }
    print("add thetaH_m")
    sfe_tissue_seurat[[sample]] = AddMetaData(sfe_tissue_seurat[[sample]],thetaH_m)
    # colGraph(sfe_tissue_final[[sample]], "visium_B") <- findVisiumGraph(sfe_tissue_final[[sample]],style = "B")
    colGraph(sfe_tissue_final[[sample]], "visium_B") <- findVisiumGraph(sfe_tissue_final[[sample]],style = "B",zero.policy = T)
    
    deconvolute_matrix_m[[sample]] = thetaH_m
    deconvolute_matrix[[sample]] = thetaH
    
    print("calculate getis-ord gi*")
    library(Voyager)
    for (col in colnames(thetaH_m)[!grepl(x=colnames(thetaH_m),pattern = "keep_")]) {
      print(col)
      colData(sfe_tissue_final[[sample]])[[col]] = thetaH_m[,col]
      sfe_tissue_final[[sample]] = colDataUnivariate(sfe_tissue_final[[sample]],features=col,type="localG_perm",colGraphName="visium")
    }
    SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
    require(patchwork)
    SpatialFeaturePlot(sfe_tissue_seurat[[sample]],names(p[[sample]]),pt.size.factor = 3) &
      scale_fill_gradientn(limits=c(0, 1), colours=SpatialColors(n=100)) |
    plotLocalResult(sfe_tissue_final[[sample]],
                    "localG_perm",
                    colGeometryName = "spotPoly",
                    features = names(p[[sample]]), divergent=T, diverge_center = 0)
  }
}

qs::qsave(x=list(sfe_tissue_seurat=sfe_tissue_seurat,
                 sfe_tissue_final=sfe_tissue_final,
                 deconvolute_m=deconvolute_matrix_m,
                 deconvolute_raw=deconvolute_matrix), paste0(output_path,"/downstream/spatial_data_after_deconvolute_seurat.qs"))

# hypergroup_data = qs::qread(paste0(output_path,"/downstream/spatial_data_after_deconvolute_seurat.qs"),nthreads = 20)

# library(SpatialExperiment)
pieplot = list()
colhypergroup = c("FAPs","Endothelials","MuSCs","Immunes","Tenocytes","Murals","Neurals")
for (sample in names(sfe_tissue_seurat)) {
  print(sample)
  pieplot[[sample]] = SPOTlight::plotSpatialScatterpie(
    x = seurat_to_spe(seu = sfe_tissue_seurat[[sample]],sample_id = "slice1",img_id = "slice1"),
    y = sfe_tissue_seurat[[sample]]@meta.data[,colhypergroup],
    slice <- "slice1",
    cell_types = colhypergroup,
    img = T,
    scatterpie_alpha = 1,
    pie_scale = 0.4)
}

list_cso_ff_merged=qs::qread("/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/result/FABIEN_LEGRAND/serieA/downstream/single_cell_0.5.27052024/allcomparison_merged_v3.qs",nthreads = 10)  
reso=list(lineage_b6=c("All_Harmony_hypergroup_Immunes"="subharmony_0.2","All_Harmony_hypergroup_Endothelials"="subharmony_0.2","All_Harmony_hypergroup_FAPs"="subharmony_0.4","All_Harmony_hypergroup_MuSCs"="subharmony_0.5"),
          lineage_d2=c("All_Harmony_hypergroup_Immunes"="subharmony_0.5","All_Harmony_hypergroup_Endothelials"="subharmony_0.4","All_Harmony_hypergroup_FAPs"="subharmony_0.7","All_Harmony_hypergroup_MuSCs"="subharmony_0.2")
)
list_cso_fff_merged = qs::qread("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/VISIUM/archives/scupa/global_data.qs",nthreads = 10)[[2]]
if (!file.exists("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/prism/prismobject_subcluster7/prepare_subclustering.qs")) {
  for (comparison in names(list_cso_fff_merged)) {
    dt=list()
    for (hypergroup_tested in names(list_cso_fff_merged[[comparison]])) {
      
      if (comparison == "lineage_b6" & hypergroup_tested == "All_Harmony_hypergroup_FAPs") {
        print(hypergroup_tested)
        name_hypergroup = hypergroup_tested %>% stringr::str_replace_all(pattern = "All_Harmony_hypergroup_",replacement = "")
        
        # Idents(list_cso_fff_merged[[comparison]][[hypergroup_tested]]) = reso[[comparison]][[hypergroup_tested]]
        # DimPlot(list_cso_fff_merged[[comparison]][[hypergroup_tested]])
        
        dt[[hypergroup_tested]] = dplyr::select(list_cso_fff_merged[[comparison]][[hypergroup_tested]]@meta.data,
                                                reso[[comparison]][[hypergroup_tested]]) %>%
          mutate(subcluster=paste0("subclustering_",name_hypergroup,"_",.data[[reso[[comparison]][[hypergroup_tested]]]])) %>%
          dplyr::select(subcluster) %>% tibble::rownames_to_column(var="barcode") %>%
          filter(!subcluster %in% c("subclustering_FAPs_4","subclustering_FAPs_5"))
        remove = dplyr::select(list_cso_fff_merged[[comparison]][[hypergroup_tested]]@meta.data,
                               reso[[comparison]][[hypergroup_tested]]) %>%
          mutate(subcluster=paste0("subclustering_",name_hypergroup,"_",.data[[reso[[comparison]][[hypergroup_tested]]]])) %>%
          dplyr::select(subcluster) %>% tibble::rownames_to_column(var="barcode") %>%
          filter(subcluster %in% c("subclustering_FAPs_4","subclustering_FAPs_5"))
      } else {
        print(hypergroup_tested)

        Idents(list_cso_fff_merged[[comparison]][[hypergroup_tested]]) = reso[[comparison]][[hypergroup_tested]]
        DimPlot(list_cso_fff_merged[[comparison]][[hypergroup_tested]])
        
        name_hypergroup = hypergroup_tested %>% stringr::str_replace_all(pattern = "All_Harmony_hypergroup_",replacement = "")
        dt[[hypergroup_tested]] = dplyr::select(list_cso_fff_merged[[comparison]][[hypergroup_tested]]@meta.data,
                                                reso[[comparison]][[hypergroup_tested]]) %>%
          mutate(subcluster=paste0("subclustering_",name_hypergroup,"_",.data[[reso[[comparison]][[hypergroup_tested]]]])) %>%
          dplyr::select(subcluster) %>% tibble::rownames_to_column(var="barcode")
      }
    }
    
    dt[["nohypergroupdeconvolute"]] = dplyr::select(list_cso_ff_merged[[comparison]]@meta.data, All_Harmony_0.8_hypergroup) %>%
      filter(!All_Harmony_0.8_hypergroup %in% c("Endothelials","FAPs","Immunes","MuSCs")) %>%
      mutate(subcluster=All_Harmony_0.8_hypergroup) %>% dplyr::select(subcluster) %>% tibble::rownames_to_column(var="barcode")
    
    final = data.table::rbindlist(dt) %>%
      data.frame() %>%
      column_to_rownames(var="barcode")
    
    if (comparison == "lineage_b6") {
      list_cso_ff_merged[[comparison]] = subset(list_cso_ff_merged[[comparison]], cells = rownames(final))
      list_cso_ff_merged[[comparison]] = AddMetaData(list_cso_ff_merged[[comparison]],final)
    } else {
      list_cso_ff_merged[[comparison]] = AddMetaData(list_cso_ff_merged[[comparison]],final)
      DimPlot(list_cso_ff_merged[[comparison]],group.by = c("hypergroup","subcluster"))
    }
  }
  
  qs::qsave(x = list_cso_ff_merged, file="/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/prism/prismobject7/prepare_subclustering.qs")
} else {
  list_cso_ff_merged = qs::qread("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/prism/prismobject7/prepare_subclustering.qs")
}
clustername = lapply(list_cso_ff_merged, function(x) x@meta.data$subcluster %>% unique())

refPhi_obj_subcluster = list()
spatial_obj = list()
boxplot = list()
deconvolute_s_matrix_m = list()
deconvolute_s_matrix = list()
subsetcluster = list()
pcluster = list()
for (comparison in c("lineage_b6","lineage_d2")) {
  print(paste0("* ",comparison))
  spatial_obj[[comparison]] = list()
  if (comparison == "lineage_b6") {
    data_space = list(GSM5979972_C57BL10=GetAssayData(sfe_tissue_seurat$GSM5979972_C57BL10 ,assay="spatial") %>% t(),  
                      GSM5979974_MDX=GetAssayData(sfe_tissue_seurat$GSM5979974_MDX ,assay="spatial") %>% t(),
                      Vis5A_2dpi=GetAssayData(sfe_tissue_seurat$Vis5A_2dpi ,assay="spatial") %>% t(),
                      Vis7B_5dpi=GetAssayData(sfe_tissue_seurat$Vis7B_5dpi ,assay="spatial") %>% t(),
                      Vis9A_7dpi=GetAssayData(sfe_tissue_seurat$Vis9A_7dpi ,assay="spatial") %>% t() )
  } else {
    data_space = list(GSM5979973_DBA2J=GetAssayData(sfe_tissue_seurat$GSM5979973_DBA2J ,assay="spatial") %>% t(),  
                      GSM5979975_D2MDX=GetAssayData(sfe_tissue_seurat$GSM5979975_D2MDX ,assay="spatial") %>% t(),
                      GSM7055901_WT1=GetAssayData(sfe_tissue_seurat$GSM7055901_WT1 ,assay="spatial") %>% t(),
                      GSM7055902_WT2=GetAssayData(sfe_tissue_seurat$GSM7055902_WT2 ,assay="spatial") %>% t(),
                      GSM7055903_MDX1=GetAssayData(sfe_tissue_seurat$GSM7055903_MDX1 ,assay="spatial") %>% t(),
                      GSM7055904_MDX2=GetAssayData(sfe_tissue_seurat$GSM7055904_MDX2 ,assay="spatial") %>% t(),
                      GSM7055905_WT3=GetAssayData(sfe_tissue_seurat$GSM7055905_WT3 ,assay="spatial") %>% t(),
                      GSM7055906_MDX3=GetAssayData(sfe_tissue_seurat$GSM7055906_MDX3 ,assay="spatial") %>% t(),
                      GSM7055907_MDX4=GetAssayData(sfe_tissue_seurat$GSM7055907_MDX4 ,assay="spatial") %>% t(),
                      GSM7055908_MDX5=GetAssayData(sfe_tissue_seurat$GSM7055908_MDX5 ,assay="spatial") %>% t(),
                      GSM7230943_Day1_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230943_Day1_PostInjury ,assay="spatial") %>% t(),
                      GSM7230944_Day3_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230944_Day3_PostInjury ,assay="spatial") %>% t(),
                      GSM7230945_Day5_PostInjury=GetAssayData(sfe_tissue_seurat$GSM7230945_Day5_PostInjury ,assay="spatial") %>% t())
  }
  for (sample in names(data_space)) {
    print(toupper(sample))
    set.seed(42)
    if (!file.exists(paste0("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/deconvolution/subcluster/requestJanuary2025/modeleprism_",comparison,".qs"))) {
      refPhi_obj_subcluster <- refPrepare(sc_Expr = list_cso_ff_merged[[comparison]]@assays[["SCT"]]@counts,
                               cell.type.labels = list_cso_ff_merged[[comparison]]$All_Harmony_0.8_hypergroup,
                               cell.state.labels = list_cso_ff_merged[[comparison]]$subcluster)
      qs::qsave(refPhi_obj_subcluster,file = paste0("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/deconvolution/subcluster/requestJanuary2025/modeleprism_",comparison,".qs"))
    } else {
      refPhi_obj_subcluster[[comparison]] = qs::qread(file = paste0("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/deconvolution/subcluster/requestJanuary2025/modeleprism_",comparison,".qs"))
    }
    if (!file.exists(paste0(output_path,"/downstream/result_instaprism_subcluster",sample,".qs"))) {
      print(" go instaprism")
      
      spatial_obj[[comparison]][[sample]] = InstaPrism(bulk_Expr = data_space[[sample]] %>% t() %>% as.matrix(),
                                                       refPhi_cs = refPhi_obj_subcluster[[comparison]])
      thetaS = spatial_obj[[comparison]][[sample]]@Post.ini.cs@theta %>% t() %>% data.frame()
      
      qs::qsave(spatial_obj[[comparison]][[sample]],file = paste0(output_path,"/downstream/result_instaprism_subcluster",sample,".qs"))
    } else {
      print(" read qs")
      spatial_obj[[comparison]][[sample]] = qs::qread(paste0(output_path,"/downstream/result_instaprism_subcluster",sample,".qs"))
      thetaS = spatial_obj[[comparison]][[sample]]@Post.ini.cs@theta %>% t() %>% data.frame()
    }
    print("fix threshold")
    thetaS_m = thetaS
    thetaS_m[thetaS_m < 0.1] <- 0
    
    
    tmp = list()
    for (x in c("FAPs","Immunes","Endothelials","MuSCs")) {
      print(x)
      for (col in colnames(thetaS_m)[grepl(x=colnames(thetaS_m),pattern = x,ignore.case = "keep_")]) {
        print(paste0(" ",col))
        if (x == "MuSCs") {
          tmp[[col]] = scater::isOutlier(metric= thetaS_m[,col],type="lower") %>%
            as.character() %>%
            factor(levels=c("TRUE","FALSE"),labels=c("nokeep","keep")) %>%
            as.character()
        } else {
          print(paste0(" ",col))
          tmp[[col]] = scater::isOutlier(metric= thetaS_m[,col],type="higher") %>%
            as.character() %>%
            factor(levels=c("TRUE","FALSE"),labels=c("keep","nokeep")) %>%
            as.character()
        }
      }
    }
    
    tmp = do.call("cbind",tmp) %>% data.frame()
    colnames(tmp)=paste0("keep_",colnames(tmp))
    rownames(tmp)=rownames(thetaS_m)
    thetaS_m = merge(thetaS_m,tmp,by=0) %>% column_to_rownames(var="Row.names")
    
    sfe_tissue_seurat[[sample]] = AddMetaData(sfe_tissue_seurat[[sample]], thetaS_m)
    
    subsetcluster[[sample]] = list()
    pcluster[[sample]] = list()
    for (hypergroup in colnames(thetaH)[1:7]) {
      if (length(rownames(sfe_tissue_seurat[[sample]]@meta.data[sfe_tissue_seurat[[sample]][[paste0("keep_",hypergroup)]]=="keep",])) != 0) {
        colsubcluster = clustername[[comparison]][grepl(x=clustername[[comparison]],
                                                        pattern=hypergroup)]
        for (subcluster in colsubcluster) {
          
          if (length(unique(c(0.1,round(max(sfe_tissue_seurat[[sample]][[subcluster]]),1))))==1) {
            custom_breaks = c(0.1,0.2)
          } else  {
            custom_breaks = c(0.1,round(max(sfe_tissue_seurat[[sample]][[subcluster]]),1))
          }
          
          pcluster[[sample]][[subcluster]] = SpatialFeaturePlot(sfe_tissue_seurat[[sample]], 
                                                                features=subcluster,pt.size.factor = 2) +
            if (custom_breaks[[2]] != 0) {
              scale_fill_gradient2(low = "white", high = "darkred",
                                   # color_hypergroup[[hypergroup]],
                                   midpoint = .1,
                                   breaks=seq(custom_breaks[[1]],custom_breaks[[2]],0.1),
                                   # labels=c(0.1,1),
                                   na.value="white",
                                   limits=custom_breaks)}
          pcluster[[sample]][[subcluster]] = pcluster[[sample]][[subcluster]] +
            SpatialDimPlot(sfe_tissue_seurat[[sample]],group.by=paste0("keep_",subcluster),pt.size.factor = 4) +
            ggplot(sfe_tissue_seurat[[sample]]@meta.data,aes(x=.data[[paste0("keep_",subcluster)]],y=.data[[subcluster]])) + 
            geom_boxplot(outlier.size = -1) + 
            geom_jitter(shape=1) +
            geom_hline(yintercept = 0.1,linetype="dashed",color="red") +
            theme_bw() +
            ggplot(sfe_tissue_seurat[[sample]]@meta.data, aes(x=.data[[subcluster]],fill=.data[[paste0("keep_",subcluster)]])) + 
            geom_histogram(bins=100,alpha=0.7) +
            geom_density(alpha=0.5) +
            geom_vline(xintercept = 0.1,linetype="dashed",color="red") +
            labs(x=paste0("Deconvolute Score ",subcluster)) +
            theme_bw()
        }
      } else {
        pcluster[[sample]][[subcluster]] = NULL
      }
    }
    
    library(Voyager)
    for (col in colnames(thetaS_m)[!grepl(x=colnames(thetaS_m),pattern = "keep_")]) {
      print(col)
      colData(sfe_tissue_final[[sample]])[[col]] = thetaS_m[,col]
      sfe_tissue_final[[sample]] = colDataUnivariate(sfe_tissue_final[[sample]],features=col,type="localG_perm",colGraphName="visium")
    } 
    
    deconvolute_name = c(colhypergroup,clustername[[comparison]]) %>% unique()
    boxplot[[sample]] = merge(sfe_tissue_seurat[[sample]]@meta.data[,deconvolute_name],
                    GetTissueCoordinates(sfe_tissue_seurat[[sample]]),
                    by=0) %>% dplyr::select(-Row.names) %>%
      pivot_longer(cols=deconvolute_name,names_to = "deconvolute",values_to="DS") %>%
      # ggplot(aes(x=reorder(deconvolute,-DS,FUN = median),y=DS,fill=deconvolute)) +
      ggplot(aes(x=deconvolute,y=DS,fill=deconvolute)) +
      geom_boxplot(outlier.size = -1) + 
      geom_jitter(alpha=0.3) +
      theme_bw() +
      scale_y_log10() +
      theme(axis.text.x = element_text(angle=90,hjust = 1))
    
    deconvolute_s_matrix_m[[sample]] = thetaS_m
    deconvolute_s_matrix[[sample]] = thetaS
  }
}
qs::qsave(x=list(sfe_tissue_seurat=sfe_tissue_seurat,
                 sfe_tissue_final=sfe_tissue_final,
                 deconvolute_m=deconvolute_s_matrix_m,
                 deconvolute_raw=deconvolute_s_matrix,
                 boxplot=boxplot), paste0(output_path,"/downstream/spatial_data_after_deconvolute_subcluster_seurat.qs"))
data = qs::qread(paste0(output_path,"/downstream/spatial_data_after_deconvolute_subcluster_seurat.qs"),nthreads=20)





sfe <- SFEData::McKellarMuscleData("full")
sfe_tissue <- sfe[,colData(sfe)$in_tissue]
sfe_tissue <- sfe_tissue[rowSums(counts(sfe_tissue)) > 0,]
sfe_tissue <- logNormCounts(sfe_tissue)
colGraph(sfe_tissue, "visium") <- findVisiumGraph(sfe_tissue)
annotGraph(sfe_tissue, "myofiber_poly2nb") <- 
  findSpatialNeighbors(sfe_tissue, type = "myofiber_simplified", MARGIN = 3,
                       method = "poly2nb", zero.policy = TRUE)

annotGraph(sfe_tissue, "myofiber_poly2nb_B") <- 
  findSpatialNeighbors(sfe_tissue, type = "myofiber_simplified", MARGIN = 3,
                       method = "poly2nb", zero.policy = TRUE, style = "B")
sfe_tissue <- annotGeometryUnivariate(sfe_tissue, "localG_perm", "area", 
                                      annotGeometryName = "myofiber_simplified",
                                      annotGraphName = "myofiber_poly2nb_B",
                                      include_self = TRUE, zero.policy = TRUE)

annotGraph(sfe_tissue, "myofiber_poly2nb_C") <- 
  findSpatialNeighbors(sfe_tissue, type = "nuclei", MARGIN = 3,
                       method = "poly2nb", zero.policy = TRUE, style = "B")
sfe_tissue <- annotGeometryUnivariate(sfe_tissue, "localG_perm", "area", 
                                      annotGeometryName = "nuclei",
                                      annotGraphName = "myofiber_poly2nb_C",
                                      include_self = TRUE, zero.policy = TRUE)
plotLocalResult(sfe_tissue, "localG_perm", "area", 
                annotGeometryName = "nuclei",
                divergent = TRUE, diverge_center = 0)







Vis5A_2dpi_deconvolute = data$deconvolute_m$Vis5A_2dpi %>%
  rownames_to_column(var="barcode") %>%
  mutate(barcode=str_replace_all(string=barcode,pattern="-1",replacement="")) %>%
  filter(barcode %in% colnames(sfe_tissue)) %>%
  column_to_rownames(var="barcode")

absence = data.frame(barcode=colnames(sfe_tissue)[!colnames(sfe_tissue) %in% Vis5A_2dpi_deconvolute$barcode],
                     subclustering_FAPs_1=0,subclustering_Endothelials_2=0,subclustering_MuSCs_3=0,            
                     subclustering_MuSCs_4=0,subclustering_Immunes_2=0,subclustering_FAPs_2=0,            
                     subclustering_Immunes_3=0,subclustering_MuSCs_2=0,Tenocytes=0,                       
                     subclustering_Immunes_1=0,Murals=0,subclustering_Endothelials_3=0,   
                     subclustering_Endothelials_1=0,subclustering_MuSCs_1=0,subclustering_MuSCs_5=0,            
                     subclustering_FAPs_3=0,
                     Neurals=0,
                     subclustering_Immunes_4=0) %>%
  column_to_rownames(var="barcode")

cols = colnames(data$deconvolute_m$Vis5A_2dpi)[!grepl(x=colnames(data$deconvolute_m$Vis5A_2dpi),pattern = "keep_")]
Vis5A_2dpi_deconvolute2 = data.frame(rbind(Vis5A_2dpi_deconvolute[,cols],absence),
                                     check.names = F)
colData(sfe_tissue) = merge(colData(sfe_tissue), Vis5A_2dpi_deconvolute2,by=0)




sfe_tissue <- annotGeometryUnivariate(sfe_tissue, "localG_perm", "area", 
                                      annotGeometryName = "nuclei",
                                      annotGraphName = "myofiber_poly2nb_B",
                                      include_self = TRUE, zero.policy = TRUE)


for (col in colnames(Vis5A_2dpi_deconvolute2)) {
  print(col)
  sfe_tissue = colDataUnivariate(sfe_tissue,features=col,type="localG_perm",colGraphName="visium")
} 


plotLocalResult(sfe_tissue, "localG_perm", "nuclei", 
                annotGeometryName = "myofiber_simplified",
                divergent = TRUE, diverge_center = 0)

# plotLocalResult(data$sfe_tissue_final$GSM7055906_MDX3, "localG_perm", 
#                 features = colhypergroup,
#                 # features = clustername[[2]], 
#                 attribute="-log10p_adj Sim",
#                 colGeometryName = "spotPoly", divergent = TRUE,
#                 color=-log(0.05),
#                 diverge_center = 0, swap_rownames = "symbol")

pieplot = list()
pieplots = list()
cols=setNames(randomcoloR::distinctColorPalette(k = 7),
              c("FAPs","Endothelials","Immunes","Tenocytes","Murals","Neurals","MuSCs"))
cols_hypergroup=setNames(randomcoloR::distinctColorPalette(k = 6),
              c("FAPs","Endothelials","Immunes","Tenocytes","Murals","Neurals"))
# cols_d2_sub =setNames(randomcoloR::distinctColorPalette(10),
#                       c(clustername[[2]][c(1,4,5,7,8,9,11,14,15,16)]))
dir.create("/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/BENCHMARK/downstream/pieplot/")

for (lineage in c("lineage_b6","lineage_d2")) {
  if (lineage == "lineage_b6") {
    data_space = c("GSM5979972_C57BL10",   "GSM5979974_MDX")
    subcluster_col = clustername[[1]][c(1,3,5,6,7,8,10,14,15,16,17)]
  } else {
    data_space = c("GSM5979973_DBA2J","GSM5979975_D2MDX",
                    "GSM7055901_WT1","GSM7055902_WT2",
                    "GSM7055903_MDX1","GSM7055904_MDX2",
                    "GSM7055905_WT3","GSM7055906_MDX3",
                    "GSM7055907_MDX4","GSM7055908_MDX5",
                    "GSM7230943_Day1_PostInjury","GSM7230944_Day3_PostInjury",
                    "GSM7230945_Day5_PostInjury")
    subcluster_col = clustername[[2]][c(1,4,5,7,8,9,11,14,15,16)]
  }
  
  cols =setNames(randomcoloR::distinctColorPalette(length(subcluster_col)),
                 subcluster_col)
  
  for (sample in data_space) {
    print(sample)
    pieplot[[sample]] = SPOTlight::plotSpatialScatterpie(
      x = seurat_to_spe(seu = data$sfe_tissue_seurat[[sample]],sample_id = "slice1",img_id = "slice1"),
      # y = data$sfe_tissue_seurat[[sample]]@meta.data[,clustername[[2]][c(1,4,5,7,8,9,11,14,15,16)]],
      y = data$sfe_tissue_seurat[[sample]]@meta.data[,names(cols_hypergroup)],
      slice <- "slice1",
      cell_types = names(cols_hypergroup),
      # cell_types = clustername[[2]][c(1,4,5,7,8,9,11,14,15,16)],
      img = T,
      scatterpie_alpha = 1,
      pie_scale = 0.4) + scale_fill_manual(values=cols_hypergroup)
    
    ggsave(plot = pieplot[[sample]],
           filename = paste0("pieplot_hypergroup_",sample,".pdf"),
           device = cairo_pdf,width = 15,height = 10,dpi = 1000,
           path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/BENCHMARK/downstream/pieplot/"
    )
    pieplots[[sample]] = SPOTlight::plotSpatialScatterpie(
      x = seurat_to_spe(seu = data$sfe_tissue_seurat[[sample]],sample_id = "slice1",img_id = "slice1"),
      y = data$sfe_tissue_seurat[[sample]]@meta.data[,subcluster_col],
      slice <- "slice1",
      cell_types = subcluster_col,
      img = T,
      scatterpie_alpha = 1,
      pie_scale = 0.4) + scale_fill_manual(values=cols)
    
    ggsave(plot = pieplots[[sample]],
           filename = paste0("pieplot_subcluster_",sample,".pdf"),
           device = cairo_pdf,width = 15,height = 10,dpi = 1000,
           path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/BENCHMARK/downstream/pieplot/"
    )
  }
}

# ggarrange(plotlist=pieplot,labels=names(sfe_tissue_seurat)[c(7:8,10:12)])
# pieplot$GSM7055903_MDX1
# pieplot$GSM7055904_MDX2
# pieplot$GSM7055906_MDX3
# pieplot$GSM7055907_MDX4
# pieplot$GSM7055908_MDX5

ggsave(plot = p[[feature]][[sample]],
       filename = paste0(sample,"_",feature,".pdf"),
       device = cairo_pdf,width = 15,height = 10,dpi = 1000,
       path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/genes/"
)

p = list()
for (feature in c("Cacna1a","Mtmr10","Cxcl16")) {
  p[[feature]] = list()
  for (sample in names(data$sfe_tissue_seurat)) {
    data$sfe_tissue_final[[sample]] <- scuttle::logNormCounts(data$sfe_tissue_final[[sample]])
    p[[feature]][[sample]]=plotSpatialFeature(data$sfe_tissue_final[[sample]],
                                   colGeometryName = "spotPoly",
                                   annotGeometryName = "tissueBoundary", 
                                   features = feature,
                                   swap_rownames = "symbol",
                                   annot_fixed = list(fill = NA, color = "black"),
                                   shape=1) & 
      viridis::scale_fill_viridis() &
      guides(colour = guide_legend(override.aes = list(size=2), ncol = 2))
    
    ggsave(plot = p[[feature]][[sample]],
           filename = paste0(sample,"_",feature,".pdf"),
           device = cairo_pdf,width = 15,height = 10,dpi = 1000,
           path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED/genes/"
    )
  }
}
