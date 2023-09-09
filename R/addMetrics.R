#' Inspired by OSCA book, the following function automatically adds gene-wise  
#' and cell-wise QC metrics on low library size, high mitochondria percentage, 
#' and low abundance genes to a Seurat object and store for later use. 
#' `addQCMetrics_seu()` is the combination of `addQCMetricsPerCell_seu()` and 
#' `addQCMetricsPerGene_seu()`.
#'
#' @param seu A Seurat object, can contain a raw assay of any name, such as 
#' "RNA" for chromium, "Spatial" for Visium, or "Xenium" for Xenium data.
#'
#' @return Returns a Seurat object with both per-cell and per-gene flags.
#' 1. Two per-cell flags are derived.  `libsize_drop` and `mito_drop` are stored 
#' in column meta data as per-cell low total count flag and high mitochondria  
#' gene percentage flag, respectively. Column-wise meta data of a Seurat object    
#' can be retrieved by `seu[["libsize_drop"]]` or `seu$libsize_drop`. 
#' 2. Two per-gene statistics and flag are derived. `means` and `lowgenecount_drop`
#' are stored in gene-wise meta data, which can be retrieved by, e.g. for a raw
#' assay with name "RNA", `seu[["RNA"]][["means"]]`. 
#' @export
#'
#' @examples 
#' \dontrun{
#' ## Chromium Example
#' data("pbmc_small")
#' pbmc_small <- addQCMetrics_seu(pbmc_small)
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetrics_seu(vis)
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetrics_seu(xe)
#' }
addQCMetrics_seu <- function(seu){
  ## Per cell
  seu$percent.mt <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  
  libsize_drop <- scuttle::isOutlier(
    metric = as.numeric(unlist(seu@meta.data[, grepl(paste0("nCount_", names(seu@assays)[1]), 
                                                     colnames(seu@meta.data))])), 
    type = "lower",
    log = TRUE) 
  
  mito_drop <- scuttle::isOutlier(
    metric = seu$percent.mt,
    type = "higher") 
  
  seu$libsize_drop <- libsize_drop
  seu$mito_drop <- mito_drop
  
  ## Per gene
  gene_means <- as.numeric(unlist(rowMeans(Seurat::GetAssayData(seu, "counts"), na.rm = TRUE)))
  lowgenecount_drop <- log(gene_means) < -5 | gene_means <= 0
  
  seu[[names(seu@assays)[1]]][["means"]] <- gene_means
  seu[[names(seu@assays)[1]]][["lowgenecount_drop"]] <- lowgenecount_drop
  
  return(seu)
}

#' Inspired by OSCA book, the following function automatically adds cell-wise 
#' QC metrics on low library size and high mitochondria percentage to a Seurat 
#' object and store for later use. 
#'
#' @param seu A Seurat object, can contain a assay of any name, such as 
#' "RNA" for chromium, "Spatial" for Visium, or "Xenium" for Xenium data.
#'
#' @return Returns a Seurat object with per-cell flags.
#' Two per-cell flags are derived.  `libsize_drop` and `mito_drop` are stored 
#' in column meta data as per-cell low total count flag and high mitochondria  
#' gene percentage flag, respectively. Column-wise meta data of a Seurat object    
#' can be retrieved by `seu[["libsize_drop"]]` or `seu$libsize_drop`. 
#' @export
#'
#' @examples
#' \dontrun{
#' ## Chromium Example
#' data("pbmc_small")
#' pbmc_small <- addQCMetricsPerCell_seu(pbmc_small)
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerCell_seu(vis)
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerCell_seu(xe)
#' }
addQCMetricsPerCell_seu <- function(seu){
  ## Per cell
  seu$percent.mt <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  
  libsize_drop <- scuttle::isOutlier(
    metric = as.numeric(unlist(seu@meta.data[, grepl(paste0("nCount_", names(seu@assays)[1]), 
                                                     colnames(seu@meta.data))])), 
    type = "lower",
    log = TRUE) 
  
  mito_drop <- scuttle::isOutlier(
    metric = seu$percent.mt,
    type = "higher") 
  
  seu$libsize_drop <- libsize_drop
  seu$mito_drop <- mito_drop
  
  return(seu)
}

#' Inspired by OSCA book, the following function automatically adds gene-wise  
#' QC metric on low abundance genes to a Seurat object and store for later use. 
#'
#' @param seu A Seurat object, can contain a assay of any name, such as 
#' "RNA" for chromium, "Spatial" for Visium, or "Xenium" for Xenium data.
#'
#' @return Returns a Seurat object with both per-cell and per-gene flags.
#' 2. Two per-gene statistics and flag are derived. `means` and `lowgenecount_drop`
#' are stored in gene-wise meta data, which can be retrieved by, e.g. for a raw
#' assay with name "RNA", `seu[["RNA"]][["means"]]`. 
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small <- addQCMetricsPerGene_seu(pbmc_small)
#' 
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerGene_seu(vis)
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerGene_seu(xe)
#' }
addQCMetricsPerGene_seu <- function(seu){
  ## Per gene
  gene_means <- as.numeric(unlist(rowMeans(Seurat::GetAssayData(seu, "counts"), na.rm = TRUE)))
  lowgenecount_drop <- log(gene_means) < -5 | gene_means <= 0
  
  seu[[names(seu@assays)[1]]][["means"]] <- gene_means
  seu[[names(seu@assays)[1]]][["lowgenecount_drop"]] <- lowgenecount_drop
  
  return(seu)
}
