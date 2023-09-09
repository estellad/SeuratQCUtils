#' Inspired by ggspavis in OSTA book, spatial QC plot for binary flags, such 
#' as `libsize_drop`, to see which spots or cells have not pass QC check in 
#' `addQCMetrics_seu()`. 
#'
#' @param seu_sp A Visium or Xenium Seurat object with spatial coordinates.
#' @param flag A binary flag contains two classes of TRUE or FALSE
#'
#' @return A spatial QC plot highlighting spots to be removed after QC in red.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Visium Example
#' data_dir <- '/path/to/Visium/outs/'
#' list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
#' vis <- Seurat::Load10X_Spatial(data.dir = data_dir)
#' vis <- SeuratQCUtils::addQCMetricsPerCell_seu(vis)
#' p <- SeuratQCUtils::plotQCSpatial_seu(seu_sp = vis, flag = "libsize_drop")
#' p
#' 
#' ## Xenium Example
#' data_dir <- '/path/to/Xenium/outs/'
#' list.files(data_dir) # Should show cell_feature_matrix.h5
#' xe <- Seurat::LoadXenium(data.dir = data_dir)
#' xe <- SeuratQCUtils::addQCMetricsPerCell_seu(xe)
#' p <- SeuratQCUtils::plotQCSpatial_seu(seu_sp = xe, flag = "libsize_drop")
#' p
#' }
plotQCSpatial_seu <- function(seu_sp, flag = "libsize_drop"){
  if(!("array_row" %in% colnames(seu_sp@meta.data))){
    seu_sp[[c("array_row", "array_col")]] <- Seurat::GetTissueCoordinates(seu_sp)[, 1:2]
    seu_sp[["cell_id"]] <- colnames(seu_sp)
  }
  
  df_cellmeta <- seu_sp@meta.data[, grepl(paste0("nCount_", names(seu_sp@assays)[1], 
                                                 "|nFeature_", names(seu_sp@assays)[1]), 
                                          colnames(seu_sp@meta.data))]
  df <- data.frame(seu_sp[[c("cell_id", "array_col", "array_row", flag)]])
  df <- cbind(df, df_cellmeta)
  ggplot(df, aes(x = array_row, y = array_col, color = get(flag))) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_color_manual(name = flag, values = c("gray85", "red")) + 
    ggtitle("QC dots") + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}


