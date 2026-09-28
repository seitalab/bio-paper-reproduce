# A sustained type I IFN-neutrophil-IL-18 axis drives pathology during mucosal viral infection
# Lebratti et al. (2021): Figure 3A–D reanalysis in R - Seurat.

# ---- 1. Packages and paths ---------------------------------------------------
cat('\n[1/14] Loading packages and setting paths\n')
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
script_path <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
PROJECT_DIR <- dirname(dirname(normalizePath(script_path)))
DATA_DIR <- file.path(PROJECT_DIR, 'data', 'ifn-neutrophil', 'GSE161336', 'raw')
SOURCE_DIR <- file.path(PROJECT_DIR, 'data', 'ifn-neutrophil', 'paper-source')
FIGURE_DIR <- file.path(PROJECT_DIR, 'figures', 'ifn-neutrophil', 'r-seurat')
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
RANDOM_SEED <- 42
set.seed(RANDOM_SEED)
# Assay v3 stores a single full-gene counts/data matrix, avoiding split v5 layers.
options(Seurat.object.assay.version = 'v3')
cat('Seurat:', as.character(packageVersion('Seurat')), '\n')

# ---- 2. Load all seven 10x samples --------------------------------------------
cat('\n[2/14] Loading all seven samples\n')
samples <- data.frame(
  sample = paste0('GSM49050', 24:30),
  prefix = c('GSM4905024_mock_', 'GSM4905025_d1_HSV-1_',
             'GSM4905026_d1_HSV-2_', 'GSM4905027_d3_HSV-1_',
             'GSM4905028_d3_HSV-2_', 'GSM4905029_d5_HSV-1_',
             'GSM4905030_d5_HSV-2_'),
  condition = c('mock', 'HSV-1', 'HSV-2', 'HSV-1', 'HSV-2', 'HSV-1', 'HSV-2'),
  dpi = c(NA, 1, 1, 3, 3, 5, 5)
)
objects <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples[i, ]
  counts <- ReadMtx(
    mtx = file.path(DATA_DIR, paste0(s$prefix, 'matrix.mtx.gz')),
    cells = file.path(DATA_DIR, paste0(s$prefix, 'barcodes.tsv.gz')),
    features = file.path(DATA_DIR, paste0(s$prefix, 'features.tsv.gz')),
    feature.column = 2, unique.features = TRUE
  )
  # Seurat replaces underscores in feature names; strip mm10_ explicitly.
  rownames(counts) <- make.unique(sub('^mm10_', '', rownames(counts)))
  obj <- CreateSeuratObject(counts, min.cells = 0, min.features = 0,
                            project = s$sample)
  obj <- RenameCells(obj, new.names = paste0(s$sample, '_', colnames(obj)))
  obj$sample <- s$sample
  obj$condition <- s$condition
  obj$dpi <- s$dpi
  cat(s$sample, ':', ncol(obj), 'cells\n')
  obj
})
obj <- merge(objects[[1]], y = objects[-1], merge.data = FALSE)
obj$condition <- factor(obj$condition, levels = c('mock', 'HSV-1', 'HSV-2'))
rm(objects)
cat('Loaded:', ncol(obj), 'cells x', nrow(obj), 'genes\n')

# ---- 3. Mitochondrial QC -----------------------------------------------------
cat('\n[3/14] Removing cells with >5% mitochondrial counts\n')
mt_genes <- grep('^mt-', rownames(obj), ignore.case = TRUE, value = TRUE)
obj[['percent.mt']] <- PercentageFeatureSet(obj, features = mt_genes)
input_metadata <- obj[[]]
obj <- subset(obj, cells = colnames(obj)[obj$percent.mt <= 5])
cat(nrow(input_metadata), '->', ncol(obj), 'cells\n')

# ---- 4. Initial clustering for QC --------------------------------------------
cat('\n[4/14] Initial QC clustering\n')
# Used the same parameters as the notebook's initial clustering.
QC_N_HVGS <- 3000
QC_N_PCS <- 30
QC_N_NEIGHBORS <- 15
QC_RESOLUTION <- 0.5
qc <- NormalizeData(obj, scale.factor = 1e4, verbose = FALSE)
select_hvgs <- function(x, n) {
  x <- FindVariableFeatures(x, selection.method = 'mean.var.plot', verbose = FALSE)
  info <- x[['RNA']][[]]
  ranked <- order(info$mvp.dispersion.scaled, decreasing = TRUE, na.last = NA)
  VariableFeatures(x) <- head(rownames(info)[ranked], n)
  x
}
qc <- select_hvgs(qc, QC_N_HVGS)
qc <- ScaleData(qc, features = VariableFeatures(qc), scale.max = 10, verbose = FALSE)
qc <- RunPCA(qc, features = VariableFeatures(qc), npcs = 50,
             seed.use = RANDOM_SEED, verbose = FALSE)
qc <- FindNeighbors(qc, dims = seq_len(QC_N_PCS), k.param = QC_N_NEIGHBORS,
                    verbose = FALSE)
# Seurat clusters its SNN graph; Scanpy's neighbor graph is different.
qc <- FindClusters(qc, algorithm = 4, resolution = QC_RESOLUTION,
                   random.seed = RANDOM_SEED, n.iter = 2, verbose = FALSE)
obj$qc_cluster <- as.character(Idents(qc))
cat('QC clusters:', length(unique(obj$qc_cluster)), '\n')
rm(qc)

# ---- 5. Remove clusters with low median detected genes ------------------------
cat('\n[5/14] Cluster-level gene-count QC\n')
QC_MIN_MEDIAN_GENES <- 500 # align with notebook's median <500 rule
cluster_qc <- do.call(rbind, lapply(split(obj[[]], obj$qc_cluster), function(d) {
  data.frame(n_cells = nrow(d), median_genes = median(d$nFeature_RNA),
             mean_genes = mean(d$nFeature_RNA),
             fraction_below_500 = mean(d$nFeature_RNA < 500),
             median_total_counts = median(d$nCount_RNA),
             median_mt_pct = median(d$percent.mt))
}))
cluster_qc$remove_cluster <- cluster_qc$median_genes < QC_MIN_MEDIAN_GENES
print(cluster_qc)
low_quality <- rownames(cluster_qc)[cluster_qc$remove_cluster]
obj <- subset(obj, cells = colnames(obj)[!obj$qc_cluster %in% low_quality])
cat('Removed clusters:', paste(low_quality, collapse = ', '), '\n')
print(table(obj$sample))
cat('Retained:', ncol(obj), 'cells\n')

# ---- 6. Normalize and select variable genes again -----------------------------
cat('\n[6/14] Normalizing retained raw counts and selecting 3,000 HVGs\n')
obj <- NormalizeData(obj, scale.factor = 1e4, verbose = FALSE)
obj <- select_hvgs(obj, 3000)
cat('HVGs:', length(VariableFeatures(obj)), '\n')

# ---- 7. Regress total UMI, scale, and run PCA ----------------------------------
cat('\n[7/14] Negative-binomial UMI regression, scaling, and PCA\n')
# Paper p. 18: regress total UMI with a negative-binomial model.
# Fit raw counts; retain log-normalized data for marker plots and module scoring.
REGRESSION_MODEL <- 'negbinom'
obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = 'nCount_RNA',
                  model.use = REGRESSION_MODEL, use.umi = TRUE, scale.max = 10, verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50,
               seed.use = RANDOM_SEED, verbose = FALSE)
ggsave(file.path(FIGURE_DIR, 'pca_elbow.png'), ElbowPlot(obj, ndims = 50),
       width = 7, height = 4, dpi = 300)
cat('PCA input:', ncol(obj), 'cells x', length(VariableFeatures(obj)), 'HVGs\n')

# ---- 8. Final neighbor graph and clustering ----------------------------------
cat('\n[8/14] Final graph and Leiden clustering\n')
N_PCS <- 20
N_NEIGHBORS <- 15
LEIDEN_RESOLUTION <- 0.8
obj <- FindNeighbors(obj, dims = seq_len(N_PCS), k.param = N_NEIGHBORS, verbose = FALSE)
obj <- FindClusters(obj, algorithm = 4, resolution = LEIDEN_RESOLUTION,
                     random.seed = RANDOM_SEED, n.iter = 2, verbose = FALSE)
obj$cluster <- Idents(obj)
print(table(obj$cluster))

# ---- 9. Figure 3A: shared t-SNE and cluster IDs -------------------------------
cat('\n[9/14] Figure 3A: t-SNE\n')
obj <- RunTSNE(obj, dims = seq_len(N_PCS), seed.use = RANDOM_SEED, check_duplicates = FALSE)
xy <- Embeddings(obj, 'tsne')
plot_data <- data.frame(x = xy[, 1], y = xy[, 2], obj[[]], check.names = FALSE)
cluster_colors <- setNames(hcl.colors(nlevels(obj$cluster), 'Dynamic'), levels(obj$cluster))
centers <- function(d) aggregate(cbind(x, y) ~ cluster, d, median)
base_plot <- function(d) ggplot(d, aes(x, y)) +
  coord_cartesian(xlim = range(xy[, 1]), ylim = range(xy[, 2])) +
  labs(x = 't-SNE 1', y = 't-SNE 2') + theme_classic(base_size = 11)
label_clusters <- function(p, d) p + geom_label(
  data = centers(d), aes(x, y, label = cluster), inherit.aes = FALSE,
  size = 2.6, fill = 'white', linewidth = 0.1, label.padding = unit(0.08, 'lines'))
save_plot <- function(p, filename, width = 12, height = 5) {
  ggsave(file.path(FIGURE_DIR, filename), p, width = width, height = height, dpi = 300)
  cat('Saved:', filename, '\n')
}
p <- base_plot(plot_data) + geom_point(aes(color = cluster), size = 0.3) +
  scale_color_manual(values = cluster_colors, drop = FALSE) + labs(title = 'Figure 3A: Seurat clusters')
save_plot(label_clusters(p, plot_data), 'figure3A_tsne_clusters.png', 8, 6)

# ---- 10. Figure 3B: neutrophil markers ----------------------------------------
cat('\n[10/14] Figure 3B: S100a8 and Csf3r\n')
normalized <- GetAssayData(obj, assay = 'RNA', layer = 'data')
marker_types <- c(S100a8 = 'PMN', Csf3r = 'PMN', Csf1r = 'Monocytes / macrophages', Fcgr1 = 'Monocytes / macrophages', Cd74 = 'Antigen-presenting cells', Flt3 = 'Dendritic-cell marker', Epcam = 'Epithelial cells', Krt5 = 'Basal epithelial cells', Krt14 = 'Basal epithelial cells', Krt18 = 'Epithelial cells', Cd3d = 'T cells', Cd3e = 'T cells', Ncr1 = 'NK cells', Klrd1 = 'NK cells')
marker_plot <- function(gene, group) {
  d <- plot_data
  d$expression <- as.numeric(normalized[gene, ])
  d <- d[order(d$expression), ]
  p <- base_plot(d) + geom_point(aes(color = expression), size = 0.3) +
    scale_color_viridis_c(option = 'magma', limits = c(0, max(d$expression))) +
    labs(title = paste(marker_types[[gene]], gene, sep = ': '), color = 'Expression')
  label_clusters(p, d)
}
pmn_genes <- c('S100a8', 'Csf3r')
save_plot(wrap_plots(lapply(pmn_genes, marker_plot, group = 'PMN'), ncol = 2),
          'figure3B_neutrophil_markers.png')
marker_summary <- aggregate(t(as.matrix(normalized[pmn_genes, ])),
                             by = list(cluster = obj$cluster), FUN = mean)
marker_summary$mean_marker_expression <- rowMeans(marker_summary[, pmn_genes])
print(marker_summary[order(-marker_summary$mean_marker_expression), ])

# ---- 11. Extra cell-type marker panels and dot plot ---------------------------
cat('\n[11/14] Extra marker evidence\n')
# Only the PMN pair is paper-specified. The other markers are provisional evidence.
# Used the same markers from the notebook's marker panels.
group_markers <- list(PMN = pmn_genes,
  'Other myeloid' = c('Csf1r', 'Fcgr1', 'Cd74', 'Flt3'),
  Epithelial = c('Epcam', 'Krt5', 'Krt14', 'Krt18'),
  Lymphocytes = c('Cd3d', 'Cd3e', 'Ncr1', 'Klrd1'))
for (group in names(group_markers)) {
  cat(group, ':', paste(group_markers[[group]], collapse = ', '), '\n')
  p <- wrap_plots(lapply(group_markers[[group]], marker_plot, group = group), ncol = 2)
  save_plot(p, paste0('figure3A_', gsub(' ', '_', tolower(group)), '_markers.png'),
            height = 5 * ceiling(length(group_markers[[group]]) / 2))
}
save_plot(DotPlot(obj, features = group_markers, group.by = 'cluster') + RotatedAxis(),
          'figure3A_group_marker_dotplot.png', 12, 6)

# ---- 12. Figure 3C: mock and day-5 condition panels ----------------------------
cat('\n[12/14] Figure 3C: split the existing embedding\n')
# All seven samples contribute to the embedding; only mock/day 5 are displayed.
display_cells <- obj$condition == 'mock' | (!is.na(obj$dpi) & obj$dpi == 5)
display_data <- plot_data[display_cells, ]
conditions <- c('mock', 'HSV-1', 'HSV-2')
panels <- lapply(conditions, function(condition) {
  d <- display_data[display_data$condition == condition, ]
  cat(condition, ':', nrow(d), 'cells\n')
  p <- base_plot(d) + geom_point(aes(color = cluster), size = 0.3) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
    labs(title = paste0(condition, ' (n=', nrow(d), ')'))
  label_clusters(p, d)
})
save_plot(wrap_plots(panels, ncol = 3, guides = 'collect'), 'figure3C_tsne_by_condition.png', 16)

# ---- 13. Figure 3D: IFN-alpha response score ----------------------------------
cat('\n[13/14] Figure 3D: IFN-alpha response scoring\n')
gene_set <- fromJSON(file.path(SOURCE_DIR, 'HALLMARK_INTERFERON_ALPHA_RESPONSE.mouse.json'))
ifn_symbols <- gene_set$HALLMARK_INTERFERON_ALPHA_RESPONSE$geneSymbols
ifn_genes <- intersect(ifn_symbols, rownames(normalized))
cat('IFN-alpha genes matched:', length(ifn_genes), '/', length(ifn_symbols), '\n')
cat('Unmatched:', paste(setdiff(ifn_symbols, ifn_genes), collapse = ', '), '\n')
# AddModuleScore subtracts expression-matched control expression, like score_genes.
# Its binning/sampling differ: ctrl=50 is sampled PER GENE, not per occupied bin.
# nbin=25/ctrl=50 therefore do not reproduce Scanpy's exact control set or scores.
# The paper does not specify its scoring formula/settings or gene-set version.
obj <- AddModuleScore(obj, features = list(ifn_genes), pool = rownames(normalized),
                       nbin = 25, ctrl = 50, name = 'IFN_alpha_score',
                       seed = RANDOM_SEED, search = FALSE)
display_data$score <- obj$IFN_alpha_score1[display_cells]
score_limits <- range(display_data$score)
panels <- lapply(conditions, function(condition) {
  d <- display_data[display_data$condition == condition, ]
  d <- d[order(d$score), ]
  p <- base_plot(d) + geom_point(aes(color = score), size = 0.3) +
    scale_color_gradientn(colors = c('#3b4cc0', '#dddddd', '#b40426'), limits = score_limits) +
    labs(title = paste0(condition, ' (n=', nrow(d), ')'), color = 'IFN-alpha score')
  label_clusters(p, d)
})
save_plot(wrap_plots(panels, ncol = 3, guides = 'collect'), 'figure3D_ifn_alpha_response.png', 16)


cat('Figures:', FIGURE_DIR, '\n')
