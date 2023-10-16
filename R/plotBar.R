#' Create cell-type frequency bar plot in annotated Seurat object
#'
#' @param seu A Seurat object contains "Annotation" column in the meta data
#' @param decreasing logical that indicates if the bars should be arranged 
#' decreasing or not
#' @param ylabel default "Percentage", as we are plotting cell-type proportion
#' @param title default "Cell-type Frequency in Annotated Chromium"
#'
#' @return Returns a bar plot with cell-type percentages arranged in decreasing 
#' order for the input annotated Seurat object.
#' @export
#'
#' @importFrom ggplot2 ggplot aes after_stat geom_bar theme element_text ggtitle ylab
#' @importFrom hrbrthemes theme_ipsum
#' 
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' pbmc_small$Annotation <- pbmc_small$groups
#' p <- SeuratQCUtils::plotBar_seu(pbmc_small)
#' p
#' }
plotBar_seu <- function(seu, decreasing = TRUE, ylabel = "Percentage", 
                        title = "Cell-type Frequency in Annotated Chromium"){
  CD <- seu@meta.data
  cnt <- plyr::count(CD$Annotation)
  CD$Annotation <- factor(CD$Annotation, 
                          levels = cnt$x[order(cnt$freq, decreasing = decreasing)])
  
  p <- ggplot(data = CD, aes(x = Annotation)) + 
    geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    ylab(ylabel) + ggtitle(title)
  
  p
}
