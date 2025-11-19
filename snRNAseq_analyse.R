# tsi 2025/11/19

#### step 0 preparation ####
sample_name <- "sample"
main_dir <- "X:/test" 
output_dir <- file.path(main_dir, sample_name)  
if (!dir.exists(output_dir)) { 
  dir.create(output_dir, recursive = TRUE)       
}

#### step 1 upstream analysis ####
step_name = "step1_" 
PATH <- 'X:/GR/SC-seq/E20220094-1/0_Cellranger'
samples <- list.files(path = PATH)

sce_list <- lapply(X = samples, FUN = function(sample) {
  sce <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(PATH, sample)),
    project = sample,
    min.cells = 10,
    min.features = 10
  )
  sce[['percent.mt']] <- PercentageFeatureSet(sce, pattern = 'mt')
  return(sce)
})

pre_qc_cells <- sapply(sce_list, function(x) ncol(x))


for(i in 1:length(samples)) {
  cat(samples[i], ":", pre_qc_cells[i], "cells\n")
  sample_name <- ifelse(is.null(names(sce_list)[i]), samples[i], names(sce_list)[i])
  
}
#### step 2 QC ####
step_name = "step2_"
sce_list[[1]] <- subset(sce_list[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
sce_list[[2]] <- subset(sce_list[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
sce_list[[3]] <- subset(sce_list[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
post_qc_cells <- sapply(sce_list, function(x) ncol(x))

for(i in 1:length(samples)) {
  cat(samples[i], ":", post_qc_cells[i], "cells (", 
      round(post_qc_cells[i]/pre_qc_cells[i]*100, 1), "% retained)\n")
  
  sample_name <- ifelse(is.null(names(sce_list)[i]), samples[i], names(sce_list)[i])
}

sce_all <- merge(x = sce_list[[1]], 
                 y = sce_list[2:length(sce_list)], 
                 add.cell.ids = samples)
sce_all <- subset(sce_all, subset = percent.mt < 10)
saveRDS(sce_all, "results/final_filtered_seurat.rds")

#### step 3 find variable features and PCA ####

step_name = "step3_"
sce_all <- sce_all %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData()
sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all))

elbow_plot <- ElbowPlot(sce_all, ndims = 30)
print(elbow_plot)
#### step 4 Harmony & UMAP ####
step_name = "step4_"
var_explained <- sce_all@reductions$pca@stdev^2 / sum(sce_all@reductions$pca@stdev^2)
cum_var_explained <- cumsum(var_explained)
pc_80var <- which(cum_var_explained >= 0.8)[1]
pc_90var <- which(cum_var_explained >= 0.9)[1]

# t-SNE before harmony
selected_pcs <- 24  # from elbow plot
sce_all <- RunTSNE(sce_all, dims = 1:selected_pcs, reduction = "pca", reduction.name = "tsne_before")
p_tsne_before <- DimPlot(sce_all, reduction = "tsne_before", group.by = "orig.ident") + 
  ggtitle("t-SNE - Before Batch Correction") +
  theme(legend.position = "none")

# Harmony
sce_all_backup_1 = sce_all
sce_all <- RunHarmony(sce_all, 
                      group.by.vars = "orig.ident", 
                      dims.use = 1:selected_pcs)  # 使用与后续分析相同的PC数量

# t-SNE after harmony
sce_all <- RunTSNE(sce_all, reduction = "harmony", dims = 1:selected_pcs, reduction.name = "tsne_after")


#### step 5 clustering #### 
step_name = "step5_"
sce_all <- sce_all %>%
  RunTSNE(reduction = "harmony", dims = 1:selected_pcs, reduction.name = "harmony_tsne_v2") %>%
  FindNeighbors(reduction = "harmony", dims = 1:selected_pcs)

sce_all$unified_group <- "All Cells"

# clustering by res 0.1, 1, 0.1
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  cat(paste0("res = ", res, " clustering...\n"))
  sce_all <- FindClusters(sce_all, resolution = res)
}

# save Seurat after clustering
saveRDS(sce_all, file = file.path(output_dir, paste0(step_name, "_sce_clustered.rds")))

cluster_stats <- table(Idents(sce_all))
cluster_stats_filename <- file.path(output_dir, paste0(step_name, "_cluster_stats_selected_res.csv"))
write.csv(data.frame(Cluster = names(cluster_stats), Count = as.numeric(cluster_stats)), 
          file = cluster_stats_filename)


#### step 6 find marker and annotation #### 
step_name = "step6_"

sce_all <- JoinLayers(sce_all)
markers_all <- FindAllMarkers(
  object = sce_all,
  assay = DefaultAssay(sce_all),        
  only.pos = TRUE,                      
  min.pct = 0.1,                        
  logfc.threshold = 0.25,               
  test.use = "wilcox"                   
)

top10_per_cluster <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)


write.csv(markers_all, file.path(output_dir, paste0(step_name, "_markers_res0.2_all.csv")),
          row.names = FALSE)
write.csv(top10_per_cluster, file.path(output_dir, paste0(step_name, "_markers_res0.2_top10.csv")),
          row.names = FALSE)


cluster2type <- c(
  "0"="Oligodendrocyte",
  "1"="Astrocyte",
  "2"="Excitatory neuron",
  "3"="Oligodendrocyte",
  "4"="Excitatory neuron",
  "5"="Microglia",
  "6"="OPC",
  "7"="Inhibitory interneuron",
  "8"="Ependymal cell",
  "9"="Inhibitory interneuron",
  "10"="Endothelial cell",
  "11"="VLMC",
  "12"="Inhibitory interneuron",
  "13"="Pericyte"
)

Idents(sce_all) <- "RNA_snn_res.0.2"  
sce_all$celltype <- unname(cluster2type[as.character(Idents(sce_all))])

unique_celltypes <- unique(sce_all$celltype)


# ranking by biological significance
sce_all$celltype <- factor(
  sce_all$celltype,
  levels = c(
    "Excitatory neuron",
    "Inhibitory interneuron", 
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "Microglia",
    "Endothelial cell",
    "Pericyte",
    "VLMC",
    "Ependymal cell"
  )
)

celltype_colors <- c(
  "Excitatory neuron"      = "#F5B08D",
  "Inhibitory interneuron" = "#F6D7A8", 
  "Oligodendrocyte"        = "#A8D2EE",
  "OPC"                    = "#D0D0D0",
  "Astrocyte"              = "#CC79A7",
  "Microglia"              = "#009E73",
  "Endothelial cell"       = "#E5C2EB",
  "Pericyte"               = "#BEE7D0",
  "VLMC"                   = "#E2D1B5",
  "Ependymal cell"         = "#FFB6C1"
)

p_annot <- DimPlot(
  sce_all,
  reduction = "harmony_tsne_v2",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  cols = celltype_colors,
  label.size = 4
) + ggtitle(paste0("Cell Type Annotation (res = ", optimal_resolution, ")")) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_annot)

saveRDS(sce_all, file = file.path(output_dir, paste0(step_name, "_sce_total_annotated.rds")))

#### step 7 neuron #### 
step_name = "step7_"

# backup
sce_all$orig_celltype <- sce_all$celltype
sce_neuron <- subset(sce_all, celltype %in% c("Inhibitory interneuron","Excitatory neuron"))
DefaultAssay(sce_neuron) <- "RNA"
sce_neuron[["percent.mt"]] <- PercentageFeatureSet(sce_neuron, pattern = "^mt-|^MT-")
sce_neuron[["percent.rb"]] <- PercentageFeatureSet(sce_neuron, pattern = "^Rp[sl]|^RP[SL]")
sce_neuron[["percent.hb"]] <- PercentageFeatureSet(sce_neuron, pattern = "^Hb[^(p)]|^HB[^(P)]")
min_features <- 200    
max_features <- 8000   
max_mt_percent <- 25   
min_counts <- 500    

sce_neuron <- subset(sce_neuron, 
                     subset = nFeature_RNA >= min_features & 
                       nFeature_RNA <= max_features & 
                       percent.mt <= max_mt_percent &
                       nCount_RNA >= min_counts)

sce_neuron <- NormalizeData(sce_neuron, normalization.method = "LogNormalize", scale.factor = 10000)
sce_neuron <- FindVariableFeatures(sce_neuron, 
                                   selection.method = "vst", 
                                   nfeatures = 3000)

all_genes <- rownames(sce_neuron)
genes_to_exclude <- grep("^mt-|^MT-|^Rp[sl]|^RP[SL]|^Hb[^(p)]|^HB[^(P)]", all_genes, value = TRUE)
genes_to_scale <- setdiff(VariableFeatures(sce_neuron), genes_to_exclude)

sce_neuron <- ScaleData(sce_neuron, features = genes_to_scale)

# PCA
sce_neuron <- RunPCA(sce_neuron, features = genes_to_scale, npcs = 50)

elbow_plot <- ElbowPlot(sce_neuron, ndims = 50)
print(elbow_plot)

# save
saveRDS(sce_neuron, file = file.path(output_dir, paste0(step_name, "sce_neuron_filtered.rds")))

dims_use <- 1:30
sce_neuron <- FindNeighbors(sce_neuron, dims = dims_use)

resolutions <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0)

for(res in resolutions){
  sce_neuron <- FindClusters(
    sce_neuron,
    resolution = res,
    algorithm = 1, # 1 = Louvain
    cluster.name = paste0("neur_snn_res.", res)
  )
}

sce_neuron <- RunTSNE(sce_neuron, dims = dims_use)

# bar plot
resolution <- "neur_snn_res.1"  
cluster_info <- sce_neuron@meta.data[[resolution]]
sample_info <- sce_neuron@meta.data$orig.ident  
cluster_sample_table <- table(cluster_info, sample_info)
df_long <- as.data.frame(cluster_sample_table)
colnames(df_long) <- c("Cluster", "Sample", "Count")
df_long$Cluster <- factor(df_long$Cluster, 
                          levels = sort(as.numeric(levels(df_long$Cluster))))

write.csv(as.data.frame.matrix(cluster_sample_table), 
          file = file.path(output_dir, paste0(step_name, "_sample_cluster_table_", gsub("\\.", "_", resolution), ".csv")))


saveRDS(sce_neuron, file = file.path(output_dir, paste0(step_name, "_sce_neuron_v2.rds")))

Idents(sce_neuron) <- "neur_snn_res.1"
markers_neur <- FindAllMarkers(
  sce_neuron,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top5 <- markers_neur %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)
write.csv(markers_neur, file.path(output_dir, paste0(step_name, "_neuronal_markers_res1.0.csv")), row.names = FALSE)

#### step 8 proj neuron #### 
step_name = "step8_"

if(!exists("sce_neuron_backup")) {
  sce_neuron_backup <- sce_neuron
  cat("already backed up\n")
}

clusters_to_remove <- c(0, 2, 3)
all_clusters <- sort(unique(as.numeric(Idents(sce_neuron))))-1
clusters_to_keep <- setdiff(all_clusters, clusters_to_remove)
sce_neuron_only <- subset(sce_neuron, idents = clusters_to_keep)

old_ids <- c(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
new_ids <- 0:(length(old_ids)-1)
id_map <- data.frame(old_cluster = old_ids, new_cluster = new_ids)

meta <- sce_neuron_only@meta.data
orig <- as.numeric(as.character(meta$seurat_clusters))
meta$orig_cluster <- orig

meta$new_cluster <- id_map$new_cluster[ match(meta$orig_cluster, id_map$old_cluster) ]
sce_neuron_only@meta.data <- meta
Idents(sce_neuron_only) <- factor(meta$new_cluster, levels = new_ids)

# tsne
dims_use <- 1:30
perp <- max(5, min(50, floor((ncol(sce_neuron_only)-1)/3)))
set.seed(1234)
sce_neuron_only <- RunTSNE(sce_neuron_only, dims = dims_use,
                           perplexity = perp, check_duplicates = FALSE)

sce_neuron_only$new_cluster <- factor(sce_neuron_only$new_cluster,
                                      levels = sort(unique(sce_neuron_only$new_cluster)))
p_tsne <- DimPlot(sce_neuron_only, reduction = "tsne",
                  group.by = "new_cluster", label = TRUE, repel = TRUE)


neuron_map <- data.frame(
  new_cluster = 0:23,
  L1 = c("Ex", "Inh", "Inh", "Inh", "Inh", 
         "Ex", "Ex", "Inh", "Ex", "Ex", 
         "Ex", "Ex", "Inh", "Ex", "Ex", 
         "Ex", "Ex", "Ex", "Ex", "Ex", 
         "Ex", "Inh", "Ex", "Ex"),
  L2 = c("RelayMod", "GlyMixedInh", "GlyMixedInh", "GlyFastInh", "GlyFastInh",
         "RelayFast", "RelayMod", "GatingInh", "ModPep", "ModPep",
         "RelayFast", "RelayFast", "GatingInh", "RelayMod", "ModPep",
         "ModPep", "RelayFast", "StateMod", "ModPep", "StateMod",
         "ModPep", "GatingInh", "RelayFast", "RelayMod"),
  Subtype = c("Tshz2Ebf1Tenm2", "RelnNdnf", "Pax2Lhx1Penk", "FastGlyA", "FastGlyOtof",
              "FastRelayA", "SstGrm1", "Oprm1GlyGABA", "PenkGalOprm1", "SstGal",
              "FastRelayB", "FastTiming", "Oprm1Gpr176", "NtsCckPiezo2", "GalPACAPVGF",
              "PenkSstGrp", "FineTouchRelay", "NAmod", "NtsSstPtgds", "VGLUT23State",
              "PtgdsSlc4a4", "PenkEsr1Inh", "FastRelayC", "Calb1RelnGrp"),
  Short_CN = c("Tshz2/Ebf1黏附主传", "Reln/Ndnf抑制", "Pax2/Lhx1多肽抑制", "快速Gly抑制A", "快速Gly抑制Otof",
               "快速主传A", "Sst/mGluR1调制", "Oprm1门控GlyGABA", "Penk/Gal/Oprm1调制", "Sst/Gal调制",
               "快速主传B", "超快速定时", "Oprm1/Gpr176抑制", "Nts/Cck高频主传", "Gal/PACAP多肽",
               "Penk/Sst/Grp调制", "精细触觉主传", "去甲样状态", "Nts/Sst/Ptgds调制", "VGLUT2/3状态",
               "Ptgds整合", "Penk/Esr1抑制", "快速主传C", "Calb1/Reln/Grp主传"),
  major_class <- c(
    "Excit_Proj","Inhibitory","Inhibitory","Inhibitory","Inhibitory",
    "Excit_Proj","Excit_Proj","Inhibitory","Excit_Proj","Excit_Proj",
    "Excit_Proj","Excit_Proj","Inhibitory","Excit_Proj","Excit_Proj",
    "Excit_Proj","Excit_Proj","Excit_Proj","Excit_Proj","Excit_Proj",
    "Excit_Local","Excit_Proj","Inhibitory","Excit_Proj"),
  functional_layer <- c(
    "ModulatoryRelay","Inhibitory","Inhibitory","Inhibitory","Inhibitory",
    "FastRelay","ModulatoryRelay","Inhibitory","ModulatoryRelay","ModulatoryRelay",
    "FastRelay","FastRelay","Inhibitory","FastRelay","StateRelay",
    "ModulatoryRelay","FastRelay","StateRelay","ModulatoryRelay","StateRelay",
    "Local","Inhibitory","FastRelay","ModulatoryRelay"),
  KeyMarkers = c("Slc17a7,Satb2,Rora,Tshz2,Cntnap2",
                 "Gad1/2,Reln,Ndnf",
                 "Pax2,Lhx1,Slc6a5,Penk",
                 "Slc6a5,Syt2,Kcnc1",
                 "Slc6a5,Syt2,Otof,Kcnc1",
                 "Slc17a6,Syt2,Rorb,Rora",
                 "Slc17a6,Sst,Grm1/8",
                 "Slc6a5,Gad2,Oprm1",
                 "Slc17a6,Penk,Gal,Oprm1",
                 "Slc17a6,Sst,Gal",
                 "Slc17a6,Cux2,Grm1,Grid2",
                 "Slc17a6,Syt2,Unc13c",
                 "Gad2,Oprm1,Gpr176,Nalcn",
                 "Slc17a6,Nts,Cck,Piezo2",
                 "Slc17a6,Gal,Adcyap1,Vgf",
                 "Slc17a6,Penk,Sst,Grp",
                 "Slc17a6,Cux2,Etv1,Hs3st2",
                 "Slc17a6,Dbh,Ddc,Slc6a2,Clock",
                 "Nts,Sst,Ptgds,Ebf1-3",
                 "Slc17a6,Slc17a8,NPFF,Agrp,Pcsk1/2",
                 "Slc17a6,Ptgds,Slc4a4,Grin3a",
                 "Pax2,Lhx1,Penk,Esr1,Slc6a5",
                 "Slc17a6,Hs3st2,Rorb,Maf",
                 "Slc17a6,Calb1,Reln,Grp")
)


meta <- sce_neuron_only@meta.data
meta$new_cluster <- as.numeric(as.character(Idents(sce_neuron_only)))
idx <- match(meta$new_cluster, neuron_map$new_cluster)

meta$L1        <- neuron_map$L1[idx]
meta$L2        <- neuron_map$L2[idx]
meta$Subtype   <- neuron_map$Subtype[idx]
meta$Short_CN  <- neuron_map$Short_CN[idx]
meta$major_class<- neuron_map$major_class[idx]
meta$functional_layer<- neuron_map$functional_layer[idx]
meta$KeyMarkers<- neuron_map$KeyMarkers[idx]

meta$LabelCombined <- paste(meta$L1, meta$L2, meta$Subtype, sep = "|")
meta$Axis1_NT <- ifelse(meta$L1 == "Ex", "Glutamatergic",
                        ifelse(meta$L1 == "Inh", "Inhibitory", "Unknown"))
meta$Axis2_Function <- meta$L2

sce_neuron_only@meta.data <- meta

write.csv(neuron_map, file.path(output_dir, paste0(step_name, "_Neuron_cluster_annotation_map.csv")), row.names = FALSE)
saveRDS(sce_neuron_only, file.path(output_dir, paste0(step_name, "_sce_neuron_only_annotated.rds")))


gracile_excitatory_genes <- c(
  "Slc17a6",    
  "Camk2a",      
  "Bcl11b",     
  "Meis2"       
)

available_genes <- gracile_excitatory_genes[gracile_excitatory_genes %in% rownames(sce_neuron_only)]
missing_genes <- gracile_excitatory_genes[!gracile_excitatory_genes %in% rownames(sce_neuron_only)]

Idents(sce_neuron_only) <- sce_neuron_only@meta.data$new_cluster 

all_markers <- FindAllMarkers(sce_neuron_only, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              verbose = FALSE)

top_genes_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%  
  pull(gene) %>%
  unique()

final_genes <- unique(c(available_genes, top_genes_per_cluster))
final_genes <- final_genes[1:min(50, length(final_genes))]

if(!"scale.data" %in% names(sce_neuron_only@assays$RNA@layers) || 
   length(rownames(GetAssayData(sce_neuron_only, slot = "scale.data"))) == 0) {
  sce_neuron_only <- ScaleData(sce_neuron_only, features = rownames(sce_neuron_only))
}

expr_matrix <- GetAssayData(sce_neuron_only, slot = "scale.data")[final_genes, ]
expr_matrix <- as.matrix(expr_matrix)

sce_neuron_only@meta.data$new_cluster <- as.factor(sce_neuron_only@meta.data$new_cluster)
Idents(sce_neuron_only) <- sce_neuron_only@meta.data$new_cluster

cell_clusters <- Idents(sce_neuron_only)
cell_order <- order(cell_clusters)
expr_matrix_ordered <- expr_matrix[, cell_order]  # 使用标准化的矩阵
cluster_ordered <- cell_clusters[cell_order]


#### step 9 functional analysis #### 
step_name = "step9_"
set.seed(123)

exproj <- subset(sce_neuron_only, subset = major_class == "Excit_Proj")
exproj <- FindVariableFeatures(exproj, selection.method="vst", nfeatures=3000)

hvg_set <- VariableFeatures(exproj)

exproj <- ScaleData(exproj, features = hvg_set)
exproj <- RunPCA(exproj, features = hvg_set, npcs = 50, verbose=FALSE)

pca_mat <- Embeddings(exproj, "pca")[,1:30]
dm <- DiffusionMap(pca_mat, sigma = "local")
DCs <- eigenvectors(dm)
exproj@meta.data$DC1 <- DCs[,1]
exproj@meta.data$DC2 <- DCs[,2]

if("seurat_clusters" %in% colnames(plot_data)) {
  plot_data$group <- plot_data$seurat_clusters
} else {
  plot_data$group <- "All"
}

plot_data <- exproj@meta.data
plot_data$cell_id <- rownames(plot_data)

reductions <- names(exproj@reductions)
for(red in reductions) {
  coords <- Embeddings(exproj, reduction = red)[, 1:2]
  plot_data[paste0(red, "_1")] <- coords[, 1]
  plot_data[paste0(red, "_2")] <- coords[, 2]
}

if("seurat_clusters" %in% colnames(plot_data)) {
  plot_data$group <- plot_data$new_clusters
} else {
  plot_data$group <- "All"
}


if("pca_1" %in% colnames(plot_data)) {
  p1 <- ggplot(plot_data, aes(pca_1, pca_2, color = group)) + 
    geom_point(size = 0.5, alpha = 0.7) + theme_minimal() + labs(title = "PCA")
  print(p1)
}

if("tsne_1" %in% colnames(plot_data)) {
  p2 <- ggplot(plot_data, aes(tsne_1, tsne_2, color = group)) + 
    geom_point(size = 0.5, alpha = 0.7) + theme_minimal() + labs(title = "t-SNE")
  print(p2)
}

gene_sets <- list(
  OXPHOS = c("Ndufb1", "Ndufb2", "Cox4i1", "Cox5a", "Atp5f1a", "Atp5f1b", 
             "Uqcrb", "Uqcrc1", "Sdha", "Sdhb"),
  
  Synaptic = c("Syt2", "Rims1", "Rims2", "Snap25", "Vamp2", "Slc17a7", 
               "Slc17a6", "Cplx1", "Unc13a", "Kcnc1"),
  
  Plasticity = c("Arc", "Fos", "Fosb", "Egr1", "Egr2", "Egr4", "Junb", 
                 "Nr4a1", "Nr4a2", "Npas4", "Bdnf", "Homer1", "Dusp1", "Atf3"),
  
  AxonGuidance = c("Robo2", "Slit2", "Dcc", "Ntn1", "Epha4", "Epha5", 
                   "L1cam", "Ncam1", "Cdh2", "Cdh11", "Plxna2"),
  
  Cytoskeleton = c("Gap43", "Dcx", "Map1b", "Tubb3", "Dclk1", "Cald1", 
                   "Wasf1", "Baiap2", "Arpc3", "Gsn")
)

for(i in names(gene_sets)) {
  score_name <- paste0("TMP_", i, "1")
  exproj <- AddModuleScore(exproj, features = list(gene_sets[[i]]), 
                           name = score_name, seed = 1)
}

current_clusters <- levels(exproj$seurat_clusters)
new_names <- paste0("exproj", 0:(length(current_clusters)-1))
names(new_names) <- current_clusters

exproj$seurat_clusters <- plyr::mapvalues(exproj$seurat_clusters, 
                                          from = current_clusters, 
                                          to = new_names)

for(module in names(gene_sets)) {
  genes_found <- intersect(gene_sets[[module]], rownames(exproj))
  genes_missing <- setdiff(gene_sets[[module]], rownames(exproj))
  
  
  if(length(genes_missing) > 0) {
    cat(paste(genes_missing, collapse = ", "), "\n")
  }
}

actual_score_names <- c("TMP_OXPHOS11", "TMP_Synaptic11", "TMP_Plasticity11", 
                        "TMP_AxonGuidance11", "TMP_Cytoskeleton11")
score_data <- exproj@meta.data[, actual_score_names]
score_names <- actual_score_names


avg_scores <- exproj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise_at(actual_score_names, mean) %>%
  column_to_rownames("seurat_clusters")

colnames(avg_scores) <- module_names
saveRDS(exproj, file.path(output_dir, paste0(step_name, "_exproj.rds")))

#### step10 Endpoint differential gene analysis#### 
step_name = "step10_"

q_low  <- quantile(exproj$DC1, 0.1)
q_high <- quantile(exproj$DC1, 0.9)
exproj$DC1_state <- "Mid"
exproj$DC1_state[exproj$DC1 <= q_low]  <- "LowEnd"
exproj$DC1_state[exproj$DC1 >= q_high] <- "HighEnd"

meta_data <- exproj@meta.data
meta_data$DC1 <- exproj$DC1  
meta_data$DC1_state <- exproj$DC1_state

sample_dc1_summary <- exproj@meta.data %>%
  filter(orig.ident %in% c("C57-1", "C57-2", "C57GR")) %>%
  group_by(orig.ident) %>%
  summarise(
    mean_DC1 = mean(DC1, na.rm = TRUE),
    median_DC1 = median(DC1, na.rm = TRUE),
    sd_DC1 = sd(DC1, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  )

# scoring
ieg_genes <- c("Fos", "Jun", "Junb", "Jund", "Fosb", "Fosl1", "Fosl2",
               "Egr1", "Egr2", "Egr3", "Egr4", "Nr4a1", "Nr4a2", "Nr4a3",
               "Dusp1", "Dusp5", "Ier2", "Ier3", "Ier5", "Arc", "Npas4")
available_ieg <- intersect(ieg_genes, rownames(exproj))
exproj <- AddModuleScore(
  exproj,
  features = list(available_ieg),
  name = "IEG_score"
)

colnames(exproj@meta.data)[colnames(exproj@meta.data) == "IEG_score1"] <- "IEG_score"


# correlation between IEG & DC1 
SAMPLES <- c("C57-1", "C57-2", "C57GR")
plot_data <- exproj@meta.data %>%
  dplyr::filter(orig.ident %in% SAMPLES) %>%
  dplyr::select(orig.ident, DC1, DC2, IEG_score)

cor_overall <- cor(plot_data$DC1, plot_data$IEG_score, method = "spearman")

cor_by_sample <- plot_data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    correlation = cor(DC1, IEG_score, method = "spearman"),
    n_cells = n(),
    mean_IEG = mean(IEG_score),
    sd_IEG = sd(IEG_score),
    .groups = 'drop'
  )


# filter IEG
DefaultAssay(exproj) <- "RNA"
expr_matrix <- GetAssayData(exproj, slot = "data") 
ieg_expression <- expr_matrix[available_ieg, , drop = FALSE]
n_ieg_per_cell <- colSums(ieg_expression > 0)

exproj@meta.data$n_IEG_genes <- n_ieg_per_cell
exproj@meta.data$has_IEG <- n_ieg_per_cell > 0

ieg_stats_by_sample <- exproj@meta.data %>%
  filter(orig.ident %in% SAMPLES) %>%
  group_by(orig.ident) %>%
  summarise(
    n_total = n(),
    n_with_IEG = sum(has_IEG),
    n_without_IEG = sum(!has_IEG),
    pct_with_IEG = 100 * mean(has_IEG),
    mean_n_IEG_genes = mean(n_IEG_genes),
    .groups = 'drop'
  )

n_before <- ncol(exproj)
n_before_by_sample <- table(exproj@meta.data$orig.ident[exproj@meta.data$orig.ident %in% SAMPLES])

cells_no_ieg <- colnames(exproj)[!exproj@meta.data$has_IEG]
exproj_no_ieg <- subset(exproj, cells = cells_no_ieg)
exproj_filtered <- exproj_no_ieg

# filtered correlation between IEG & DC1
sample_data_filtered <- lapply(SAMPLES, function(s) {
  exproj_filtered@meta.data$DC1[exproj_filtered@meta.data$orig.ident == s]
})
names(sample_data_filtered) <- SAMPLES

cor_matrix_filtered <- outer(1:length(SAMPLES), 1:length(SAMPLES), Vectorize(function(i, j) {
  if(i == j) return(1.000)
  
  dc1_range <- range(c(sample_data_filtered[[i]], sample_data_filtered[[j]]), na.rm = TRUE)
  bins <- seq(dc1_range[1], dc1_range[2], length.out = 50)
  
  density1 <- hist(sample_data_filtered[[i]], breaks = bins, plot = FALSE)$density
  density2 <- hist(sample_data_filtered[[j]], breaks = bins, plot = FALSE)$density
  
  cor(density1, density2, use = "complete.obs")
}))

dimnames(cor_matrix_filtered) <- list(SAMPLES, SAMPLES)

# diff
cor_diff <- cor_matrix_filtered - cor_matrix
print(round(cor_diff, 3))


# DEG
Idents(exproj) <- exproj$DC1_state
markers_high_vs_low <- FindMarkers(exproj, 
                                   ident.1 = "HighEnd", 
                                   ident.2 = "LowEnd",
                                   logfc.threshold = 0.25, 
                                   min.pct = 0.1, 
                                   test.use = "wilcox")

deg_tab <- markers_high_vs_low %>% 
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05)

up_genes <- deg_tab %>% 
  filter(avg_log2FC > 0.5) %>%  
  pull(gene)

down_genes <- deg_tab %>% 
  filter(avg_log2FC < -0.5) %>%  
  pull(gene)

gene2id <- bitr(c(up_genes, down_genes), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Mm.eg.db)

up_ids <- gene2id$ENTREZID[gene2id$SYMBOL %in% up_genes]
down_ids <- gene2id$ENTREZID[gene2id$SYMBOL %in% down_genes]

# GO
ego_up_bp <- enrichGO(up_ids, 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

ego_down_bp <- enrichGO(down_ids, 
                        OrgDb = org.Mm.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE)

ego_up_cc <- enrichGO(up_ids, 
                      OrgDb = org.Mm.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

ego_down_cc <- enrichGO(down_ids, 
                        OrgDb = org.Mm.eg.db, 
                        ont = "CC", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE)
ego_up_mf <- enrichGO(up_ids,
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

ego_down_mf <- enrichGO(down_ids,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = TRUE)


# volcano
volcano_data <- markers_high_vs_low %>%
  rownames_to_column("gene") %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Upregulated",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    neg_log10_padj = -log10(p_val_adj + 1e-300)  
  )


# save

save_csv <- function(df, suffix, dir, step) {
  fn <- file.path(dir, paste0(step, "_", suffix, ".csv"))
  write.csv(df, fn, row.names = FALSE)
}

results_list <- list(
  "GO_BP_UP"   = ego_up_bp@result,
  "GO_BP_DOWN" = ego_down_bp@result,
  "GO_CC_UP"   = ego_up_cc@result,
  "GO_CC_DOWN" = ego_down_cc@result,
  "GO_MF_UP"   = ego_up_mf@result,
  "GO_MF_DOWN" = ego_down_mf@result
)

for (nm in names(results_list)) {
  save_csv(results_list[[nm]], nm, output_dir, step_name)
}

library(openxlsx)
wb <- createWorkbook()
for (sheet_name in names(results_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, results_list[[sheet_name]])
}
excel_filename <- file.path(output_dir, paste0(step_name, "_GO_enrichment_all_results.xlsx"))
saveWorkbook(wb, excel_filename, overwrite = TRUE)

write.csv(
  markers_high_vs_low,
  file = file.path(output_dir, paste0(step_name, "_DEG_HighEnd_vs_LowEnd.csv")),
  row.names = TRUE
)

