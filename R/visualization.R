#' Prepare data for Manhattan plot
#'
#' @param dmpTable output table from DMP analysis (ChampLike)
#' @param pngfile file to save PNG manhattan
#'
#' @export
#'
dmp_manhattan <- function(dmpTable, pngfile) {
  annotateTop=T
  mdata <- data.frame(cpg=row.names(dmpTable),
                      dmpTable[, c("CHR","MAPINFO", "P.Value", "feature", "cgi", "gene")],
                      check.names = F)
  mdata$CHR <- as.character(mdata$CHR)
  mdata$CHR  <- as.numeric(mdata$CHR)

  y.up <- max(ceiling(max(-log10(mdata$P.Value))),8)
  if (!sum(mdata$P.Value <= 10^-5))
    annotateTop=F
  png(pngfile, width = 1280, height = 980, res=150)
  qqman::manhattan(x=mdata, chr="CHR", bp = "MAPINFO", p="P.Value", snp="gene", annotatePval=10^-5,
                   annotateTop=annotateTop, ylim=c(0,y.up))
  dev.off()

}

#' Create qq plot
#'
#' @param dmpTable output table from DMP analysis (ChampLike)
#'
#' @export
#'
dmp_qq <- function(dmpTable) {
  require(ggplot2)
  pvector <- dmpTable$P.Value
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  xlab = paste("Expected -log P")
  ylab = paste("Observed -log P")
  qq <- data.frame(Expected=e, Observed=o)
  ggplot(qq, aes(x=Expected, y=Observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, col="red")
}

