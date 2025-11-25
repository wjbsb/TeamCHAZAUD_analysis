# sfe_tissue_final = qs::qread(file = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/VISIUM/FIXED/save/segmentation_done_v2.qs")
# sfe_tissue_seurat = qs::qread(file = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/VISIUM/FIXED/save/spatial_data.qs")
save = qs::qread(file = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/VISIUM/FIXED/save/save_deconvoluteONslide_v3.qs")
# scmdx = qs::qread(file = "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/result/FABIEN_LEGRAND/serieA/downstream/single_cell_0.5.27052024/allcomparison_merged_v3.qs")
results=save$dataDeconvolute
sfe_tissue_seurat=save$seurat_slides

colb6s=setNames(randomcoloR::distinctColorPalette(k = 10),
              c("FAPs","Endothelials","Immunes","Tenocytes","Murals","Neurals",
                "subclustering_MuSCs_1","subclustering_MuSCs_2","subclustering_MuSCs_3","subclustering_MuSCs_5"))
cold2s=setNames(randomcoloR::distinctColorPalette(k = 8),
              c("FAPs","Endothelials","Immunes","Tenocytes","Murals","Neurals",
                "subclustering_MuSCs_1","subclustering_MuSCs_2"))
p = (SPOTlight::plotSpatialScatterpie(
  x = seurat_to_spe(seu = sfe_tissue_seurat$B6mdx_rep1,sample_id = "S1",img_id = "S1"),
  y = sfe_tissue_seurat$B6mdx_rep1@meta.data[,names(colb6s)],
  slice <- "S1",
  cell_types = names(colb6s),
  img = T,
  scatterpie_alpha = 1,
  pie_scale = 0.4) + scale_fill_manual(values=colb6s) +
  SPOTlight::plotSpatialScatterpie(
    x = seurat_to_spe(seu = sfe_tissue_seurat$B6mdx_rep2,sample_id = "S1",img_id = "S1"),
    y = sfe_tissue_seurat$B6mdx_rep2@meta.data[,names(colb6s)],
    slice <- "S1",
    cell_types = names(colb6s),
    img = T,
    scatterpie_alpha = 1,
    pie_scale = 0.3) + scale_fill_manual(values=colb6s) ) /
(SPOTlight::plotSpatialScatterpie(
  x = seurat_to_spe(seu = sfe_tissue_seurat$D2mdx_rep1,sample_id = "S1",img_id = "S1"),
  y = sfe_tissue_seurat$D2mdx_rep1@meta.data[,names(cold2s)],
  slice <- "S1",
  cell_types = names(cold2s),
  img = T,
  scatterpie_alpha = 1,
  pie_scale = 0.5) + scale_fill_manual(values=cold2s) +
SPOTlight::plotSpatialScatterpie(
    x = seurat_to_spe(seu = sfe_tissue_seurat$D2mdx_rep2,sample_id = "S1",img_id = "S1"),
    y = sfe_tissue_seurat$D2mdx_rep2@meta.data[,names(cold2s)],
    slice <- "S1",
    cell_types = names(cold2s),
    img = T,
    scatterpie_alpha = 1,
    pie_scale = 0.5) + scale_fill_manual(values=cold2s)
)

ggsave("pieplot_global.pdf",
       plot = p,
       device = cairo_pdf,
       path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED",
       width = 15,height = 10,dpi = 1000
       )
#############################################

colnames(sfe_tissue_seurat$B6mdx_rep1@meta.data)[grepl(x=colnames(sfe_tissue_seurat$B6mdx_rep1@meta.data),
                                                       pattern="^subclustering")]
colb6s=setNames(randomcoloR::distinctColorPalette(k = 11),
                c("subclustering_FAPs_1","subclustering_FAPs_2","subclustering_FAPs_3",
                  "subclustering_Immunes_1","subclustering_Immunes_2","subclustering_Immunes_3","subclustering_Immunes_4",
                  "subclustering_MuSCs_1","subclustering_MuSCs_2","subclustering_MuSCs_3","subclustering_MuSCs_5"))

colnames(sfe_tissue_seurat$B6mdx_rep1@meta.data)[grepl(x=colnames(sfe_tissue_seurat$D2mdx_rep1@meta.data),
                                                       pattern="^subclustering")]
cold2s=setNames(randomcoloR::distinctColorPalette(k = 10),
                c("subclustering_FAPs_1","subclustering_FAPs_2","subclustering_FAPs_3",
                   "subclustering_Immunes_1","subclustering_Immunes_2","subclustering_Immunes_3","subclustering_Immunes_4","subclustering_Immunes_5",
                   "subclustering_MuSCs_1","subclustering_MuSCs_2"))

p2 = (SPOTlight::plotSpatialScatterpie(
  x = seurat_to_spe(seu = sfe_tissue_seurat$B6mdx_rep1,sample_id = "S1",img_id = "S1"),
  y = sfe_tissue_seurat$B6mdx_rep1@meta.data[,names(colb6s)],
  slice <- "S1",
  cell_types = names(colb6s),
  img = T,
  scatterpie_alpha = 1,
  pie_scale = 0.4) + scale_fill_manual(values=colb6s) +
    SPOTlight::plotSpatialScatterpie(
      x = seurat_to_spe(seu = sfe_tissue_seurat$B6mdx_rep2,sample_id = "S1",img_id = "S1"),
      y = sfe_tissue_seurat$B6mdx_rep2@meta.data[,names(colb6s)],
      slice <- "S1",
      cell_types = names(colb6s),
      img = T,
      scatterpie_alpha = 1,
      pie_scale = 0.3) + scale_fill_manual(values=colb6s) ) /
  (SPOTlight::plotSpatialScatterpie(
    x = seurat_to_spe(seu = sfe_tissue_seurat$D2mdx_rep1,sample_id = "S1",img_id = "S1"),
    y = sfe_tissue_seurat$D2mdx_rep1@meta.data[,names(cold2s)],
    slice <- "S1",
    cell_types = names(cold2s),
    img = T,
    scatterpie_alpha = 1,
    pie_scale = 0.5) + scale_fill_manual(values=cold2s) +
     SPOTlight::plotSpatialScatterpie(
       x = seurat_to_spe(seu = sfe_tissue_seurat$D2mdx_rep2,sample_id = "S1",img_id = "S1"),
       y = sfe_tissue_seurat$D2mdx_rep2@meta.data[,names(cold2s)],
       slice <- "S1",
       cell_types = names(cold2s),
       img = T,
       scatterpie_alpha = 1,
       pie_scale = 0.5) + scale_fill_manual(values=cold2s)
  )

ggsave("pieplot_global_subclusters.pdf",
       plot = p2,
       device = cairo_pdf,
       path = "/media/mna_bioinfo/MNA2_Stockage/CHAZAUD/result/VISIUM/downstream/FIXED",
       width = 15,height = 10,dpi = 1000
)
