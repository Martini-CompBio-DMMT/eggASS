#' Smoothing function
#'
#' @param y beta matrix
#' @param x beta matrix
#' @param mycluster beta matrix
#' @param bpspan base pair span
#' @param cores number of cores
#'
methyl_smooth <- function(y, x, mycluster, bpspan=bpspan, cores=cores) {
  tmpPos <- split(x,mycluster)
  tmpLength <- unlist(lapply(tmpPos,function(x) length(x)))
  names(tmpLength) <- names(tmpPos)
  tmpSpan <- unlist(lapply(tmpPos,function(x) bpspan/median(diff(x))/length(x)))
  names(tmpSpan) <- names(tmpPos)
  tmpSpan[which(tmpSpan>1)] <- 1

  if(all(class(y)=="numeric")) {
    tmpCpG <- split(y,mycluster)
    tmpSmooth <- unlist(sapply(names(tmpCpG),function(m) limma::loessFit(tmpCpG[[m]],tmpPos[[m]],span = tmpSpan[m])$fitted))
  } else {
    smallfunction <- function(tmp_y,tmp_cluster,tmp_pos,tmp_span) {
      tmpCpG <- split(tmp_y,tmp_cluster)
      smallsmooth <- unlist(sapply(names(tmpCpG),function(m) limma::loessFit(tmpCpG[[m]],tmp_pos[[m]],span = tmp_span[m])$fitted))
      return(smallsmooth)
    }
    message("Parallel on Matrix.")
    parallelclusters <- doParallel::makeCluster(cores)
    doParallel::registerDoParallel(parallelclusters)
    foreach::getDoParWorkers()
    tmpSmooth <- foreach(i = 1:ncol(y), .combine = cbind) %dopar% smallfunction(y[,i],mycluster,tmpPos,tmpSpan)
  }
  return(tmpSmooth)
}


greaterOrEqual <- function(x,y) {
  precision <- sqrt(.Machine$double.eps)
  (x >= y) | (abs(x-y) <= precision)
}

set_cutoff <- function(permBeta, quantile=0.99) {
    quantile(abs(permBeta), quantile, na.rm = TRUE)
}


compute_fwer <- function(B) {
  # Calulcate p value and FWER.
  L <- V <- A <- as.list(rep(0, B))
  for (i in 1:B) {
    nulltab <- nulltabs[[i]]
    if (nrow(nulltab) > 0) {
      L[[i]] <- nulltab$L
      V[[i]] <- nulltab$value
      A[[i]] <- nulltab$area
    }
  }

  mytots <- apply(tab[,c("L","value")],1,function(i) {
    sapply(c(1:B),function(x) {sum(greaterOrEqual(L[[x]],i[1]) & greaterOrEqual(abs(V[[x]]),abs(i[2])))
    })
  })

  mytots2 <- t(sapply(tab$area, function(i) {sapply(c(1:B),function(x) {sum(greaterOrEqual(A[[x]],i))})}))

  rate1 <- apply(mytots,2,function(x) mean(x>0))
  pvalues1 <- apply(mytots,2,function(x) sum(x))/sum(sapply(nulltabs, nrow))
  rate2 <- rowMeans(mytots2 > 0)
  pvalues2 <- rowSums(mytots2)/sum(sapply(nulltabs, nrow))
}

