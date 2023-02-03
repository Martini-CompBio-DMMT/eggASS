#' Prepare data for Manhattan plot
#'
#' @param dmpTable output table from DMP analysis (ChampLike)
#' @param adjPThr adjusted p-value threshold
#' @param pvalueThr p-value threshold
#' @param txtfile file to table. If NULL do not save.
#' @param dec decimal separator
#'
#' @export
#'
save_table_results <- function(dmpTable, adjPThr, pvalueThr, txtfile=NULL, dec=".") {
  top <- dmpTable[dmpTable$adj.P.Val <= adjPThr | dmpTable$P.Value<=pvalueThr, , drop=F]
  if (!is.null(txtfile)){
    write.table(top, file=txtfile, sep="\t", dec=dec)
  }
  top
}

#' Extract data with Pvalue or adjust PvalueThr
#'
#' @param dmpTable output table from DMP analysis (ChampLike)
#' @param adjPThr adjusted p-value threshold
#' @param pvalueThr p-value threshold
#'
#' @export
#'
extract_results <- function(dmpTable, pvalueThr=10^-5, adjPThr=0) {
  dmpTable[dmpTable$P.Value<=pvalueThr | dmpTable$adj.P.Val <= adjPThr, , drop=F]
}
