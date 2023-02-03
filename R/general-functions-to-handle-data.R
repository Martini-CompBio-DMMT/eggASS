#' Function to compute PCA from control probes
#'
#' From RaMWAS vignette using minfi::getControlAddress()
#' @param controlProbeData beta values from control probes
#' @param npc number of PC
#'
#' @export
compute_pca <- function(controlProbeData, npc) {
  data = controlProbeData - rowMeans(controlProbeData)
  covmat = crossprod(data)
  eig = eigen(covmat)
  covariates.pca = eig$vectors[,seq_len(npc)]
  colnames(covariates.pca) = paste0('PC',seq_len(npc))
  row.names(covariates.pca) <- colnames(controlProbeData)
  covariates.pca
}

#' Function to compute PCA variance
#'
#' From RaMWAS vignette using minfi::getControlAddress()
#' @param controlProbeData beta values from control probes
#' @param npc number of PC
#'
#' @export
computePCAvariance <- function(ctrlProbe, npc) {
  data = ctrlProbe - rowMeans(ctrlProbe)
  covmat = crossprod(data)
  eig = eigen(covmat)
  head(eig$values, npc)/sum(eig$values) * 100
}


#' Function to M values
#'
#' @param beta normalized beta values
#' @param capEdges edges to compress distribution
#'
#' @export
convertToM <- function(beta, capEdges = c(0.001,0.999)) {
  beta<- replace(beta,which(beta <= capEdges[1]), capEdges[1])
  beta <- replace(beta,which(beta >= capEdges[2]),capEdges[2])
  beta <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values
  beta
}

#' Extract random sub-sumple
#'
#' @param beta normalized beta values
#' @return  sorted index of random probes
#'
#' @export
#'
extractRandomProbes <- function(beta, sample_size=100000, seeds=1234) {
  to_pick <- seq_len(nrow(beta))
  if (nrow(beta) > sample_size) {
    set.seed(seeds)
    to_pick <- sort(sample(nrow(beta), sample_size))
  }
  to_pick
}
