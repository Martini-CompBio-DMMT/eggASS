#' @export
#'
attach_contrast_pvalues <- function(DMP, contrastFit) {
  padj <- apply(contrastFit$p.value, 2, p.adjust, method="BH")
  colnames(padj) <- paste0(colnames(padj), ".padjBH")
  row.names(padj) <- row.names(contrastFit$p.value)
  cbind(DMP, contrastFit$p.value[row.names(DMP),], padj[row.names(DMP),])
}
