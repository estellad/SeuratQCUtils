#' A generic histogram plotting function 
#' 
#' A generic histogram plotting function that will be used for gene-wise and 
#' cell-wise QC flags on their corresponding continuous variables. 
#'
#' @param plot_df A data frame contains one continuous variable and one binary 
#' flag for indicating which cells/spots/genes did not pass QC.
#' @param xvar One continuous variable for plotting, e.g. log of total library
#' size. 
#' @param yvar One binary flag with level TRUE/FALSE indicating if a cell/spot 
#' or a gene, e.g. libsize_drop. 
#'
#' @return A histogram, colored by QC flag status. 
#' @export
#'
#' @importFrom dplyr %>% 
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_manual labs xlab theme element_text
#' @importFrom hrbrthemes theme_ipsum
#'
#' @examples
#' \dontrun{
#' ## Chromium Example 
#' data("pbmc_small")
#' pbmc_small <- SeuratQCUtils::addQCMetricsPerCell_seu(pbmc_small)
#' plot_df <- data.frame(logtotal = log(as.numeric(pbmc_small$nCount_RNA)),
#'                       pbmc_small@meta.data)
#' p <- plotHist(plot_df, xvar = "logtotal", yvar = "libsize_drop")
#' p
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerCell_seu(vis)
#' plot_df <- data.frame(logtotal = log(as.numeric(vis$nCount_Spatial)),
#'                       vis@meta.data)
#' p <- plotHist(plot_df, xvar = "logtotal", yvar = "libsize_drop")
#' p
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerCell_seu(xe)
#' p <- SeuratQCUtils::plotHist(plot_df, xvar = "logtotal", yvar = "libsize_drop")
#' p
#' 
#' }
plotHist <- function(plot_df, xvar = "logtotal", yvar = "libsize_drop"){
  p <- plot_df %>%
    ggplot(aes(x = get(xvar), fill = get(yvar))) +
    geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 30) +
    scale_fill_manual(values = c("#404080", "#69b3a2")) +
    theme_ipsum() +
    labs(fill = "") + xlab(xvar) +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12))
  
  p
}

#' A per-cell library size QC histogram
#' 
#' A per-cell library size QC histogram that takes a Seurat object that has  
#' been preprocessed with `SeuratQCUtils::addQCMetrics()`, such that flag 
#' `libsize_drop` exists in meta.data. This function sources helper function to 
#' generate histogram with `SeuratQCUtils::plotHist()`.
#'
#' @param seu A Seurat object with per-cell QC flag `libsize_drop` stored in 
#' meta.data
#' @param yvar A binary flag name `libsize_drop` with level TRUE/FALSE  
#' indicating if a cell/spot will be removed.
#'
#' @return A histogram, colored by QC `libsize_drop` flag status. 
#' @export
#' 
#' @importFrom ggplot2 ggplot xlab ylab ggtitle
#' 
#' @examples
#' \dontrun{
#' ## Chromium Example
#' data_dir <- '/path/to/Chromium/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' chrom <- Seurat::Read10X_h5(paste0(data_dir, "filtered_feature_bc_matrix.h5"))
#' chrom <- CreateSeuratObject(counts = chrom)
#' chrom <- SeuratQCUtils::addQCMetricsPerCell_seu(chrom)
#' if(any(chrom$libsize_drop)){
#'   p <- SeuratQCUtils::plot_Hist_Low_Lib_Sizes(chrom, yvar = "libsize_drop")
#' }
#' p
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerCell_seu(vis)
#' if(any(vis$libsize_drop)){
#'   p <- SeuratQCUtils::plot_Hist_Low_Lib_Sizes(vis, yvar = "libsize_drop")
#' }
#' p
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerCell_seu(xe)
#' if(any(xe$libsize_drop)){
#'   p <- SeuratQCUtils::plot_Hist_Low_Lib_Sizes(xe, yvar = "libsize_drop")
#' }
#' p
#' 
#' }
plot_Hist_Low_Lib_Sizes <- function(seu, yvar = "libsize_drop"){
  CD <- seu@meta.data
  CD[[yvar]] <- factor(CD[[yvar]])
  CD[[yvar]] <- relevel(CD[[yvar]], "TRUE")
  
  plot_df <- data.frame(logtotal = log(as.numeric(seu@meta.data[, grepl(paste0("nCount_", names(seu@assays)[1]), 
                                                                        colnames(seu@meta.data))])),
                        CD)
  
  p <- plotHist(plot_df, xvar = "logtotal", yvar = yvar)
  
  p <- p + xlab("Log cell level total count") + ylab("Frequency") + 
    ggtitle('Cells with low library size')
  
  p
}

#' A per-cell mitochondria percentage QC histogram
#' 
#' A per-cell mitochondria percentage QC histogram that takes a Seurat object   
#' that has been preprocessed with `SeuratQCUtils::addQCMetrics()`, such that  
#' flag `mito_drop` and continuous variable `percent.mt` exist in meta.data.   
#' This function sources helper function `SeuratQCUtils::plotHist()` to 
#' generate histogram.
#'
#' @param seu A Seurat object with per-cell QC flag `mito_drop` and `percent.mt`
#' stored in meta.data. 
#' 
#' @param yvar A binary flag name `mito_drop` with level TRUE/FALSE  
#' indicating if a cell/spot will be removed.
#'
#' @return A histogram, colored by QC `mito_drop` flag status. 
#' @export
#' 
#' @importFrom ggplot2 xlab ylab ggtitle
#'
#' @examples
#' \dontrun{
#' ## Chromium Example
#' data_dir <- '/path/to/Chromium/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' chrom <- Seurat::Read10X_h5(paste0(data_dir, "filtered_feature_bc_matrix.h5"))
#' chrom <- CreateSeuratObject(counts = chrom)
#' chrom <- SeuratQCUtils::addQCMetricsPerCell_seu(chrom)
#' if(any(chrom$mito_drop)){
#'   p <- SeuratQCUtils::plot_Hist_High_Mito_Props(chrom, yvar = "mito_drop")
#' }
#' p
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerCell_seu(vis)
#' if(any(vis$mito_drop)){
#'   p <- SeuratQCUtils::plot_Hist_High_Mito_Props(vis, yvar = "mito_drop")
#' }
#' p
#' 
#' ## Xenium Example (Xenium has no mitochondria genes)
#' 
#' }
plot_Hist_High_Mito_Props <- function(seu, yvar = "mito_drop"){
  CD <- seu@meta.data
  CD[[yvar]] <- factor(CD[[yvar]])
  CD[[yvar]] <- relevel(CD[[yvar]], "TRUE")
  
  plot_df <- data.frame(log_mito_percent = log(seu$percent.mt),
                        CD)
  
  p <- plotHist(plot_df, xvar = "log_mito_percent", yvar = yvar)
  
  p <- p + xlab("Log cell level mitochondria percent") + ylab("Frequency") + 
    ggtitle('Cells with high mitochondria percentage')
  
  p
}

#' A per-gene abundance QC histogram
#' 
#' A per-gene abundance QC histogram that takes a Seurat object that has been 
#' preprocessed with `SeuratQCUtils::addQCMetrics()`, such that flag 
#' `lowgenecount_drop` and continuous variable gene `means` exist in row-wise
#' meta.data. This function sources helper function `SeuratQCUtils::plotHist()` 
#' to generate histogram.
#'
#' @param seu A Seurat object with per-gene QC flag `lowgenecount_drop` and 
#' `means` stored in row-wise meta.data. 
#'
#' @return A histogram, colored by QC `lowgenecount_drop` flag status. 
#' @export
#'
#' @importFrom ggplot2 xlab ylab ggtitle
#' 
#' @examples
#' \dontrun{
#' ## Chromium Example
#' data_dir <- '/path/to/Chromium/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' chrom <- Seurat::Read10X_h5(paste0(data_dir, "filtered_feature_bc_matrix.h5"))
#' chrom <- CreateSeuratObject(counts = chrom)
#' chrom <- SeuratQCUtils::addQCMetricsPerGene_seu(chrom)
#' if(any(unlist(chrom[["RNA"]][["lowgenecount_drop"]]))){
#'   p <- SeuratQCUtils::plot_Hist_Low_Abun_Genes(chrom)
#' }
#' p
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerGene_seu(vis)
#' if(any(unlist(vis[["Spatial"]][["lowgenecount_drop"]]))){
#'   p <- SeuratQCUtils::plot_Hist_Low_Abun_Genes(vis)
#' }
#' p
# 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerGene_seu(xe)
#' if(any(unlist(xe[["Xenium"]][["lowgenecount_drop"]]))){
#'   p <- SeuratQCUtils::plot_Hist_Low_Abun_Genes(xe)
#' }
#' p
#' }
plot_Hist_Low_Abun_Genes <- function(seu){
  plot_df <- data.frame(mean_genecount = log(unlist(seu[[names(seu@assays)[1]]][["means"]])),
                        lowgenecount_drop = factor(unlist(seu[[names(seu@assays)[1]]][["lowgenecount_drop"]])))
  
  plot_df$lowgenecount_drop <- relevel(plot_df$lowgenecount_drop, "TRUE")
  
  p <- plotHist(plot_df, xvar = "mean_genecount", yvar = "lowgenecount_drop")
  
  p <- p + xlab("Log mean count across all cells") + ylab("Frequency") + 
    ggtitle('Low abundance genes')
  
  p
}

