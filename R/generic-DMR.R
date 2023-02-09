#' Differentially methylated longitudinal test for regions
#'
#' General function to compute longitudinal test
#'
#' @param beta normalized beta values (methylation)
#' @param covariates the data.frame wit all covariates
#' @param formula the formula; last term is the considered phenotype
#' @param adjPvalDmr adjusted p-value cut off for DMR
#' @param adjust.method adjust pvalue method
#' @param probe_annotation for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @export
#'
computeDMR <- function(
    beta, covariates,
    formula="~ phenotype",
    probe_annotation = get("probe.features"),
    minProbes = 7,
    maxGap = 300,
    bpspan=1000,
    B=250,
    cutoff=NULL,
    pickCutoff=TRUE,
    cores = 3,
    adjPvalDmr=0.05,
    seeds=1234) {

  ori_beta <- beta
  check_samples_covariates(beta, covariates)
#   design <- model.matrix(as.formula(formula), covariates) #
#
#   interesting_var <- rev(all.vars(as.formula(formula)))[1] #
#   coef_idx <- grep(interesting_var, colnames(design)) #

  cpg.idx <- intersect(rownames(beta),rownames(probe.features)) #
  Anno <- probe.features[cpg.idx,]
  Anno <- Anno[order(Anno$CHR,Anno$MAPINFO),]
  cpg.idx <- rownames(Anno)
  cl <- bumphunter::clusterMaker(Anno$CHR,Anno$MAPINFO,maxGap=maxGap)
  names(cl) <- cpg.idx
  bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl) > minProbes)))]

  beta <- beta[bumphunter.idx,]
  chr <- Anno[bumphunter.idx,]$CHR
  pos <- Anno[bumphunter.idx,]$MAPINFO
  cluster <- cl[bumphunter.idx]

  beta <- beta[bumphunter.idx,]
  beta<- replace(beta,which(beta <= 0.001),0.001)
  beta <- replace(beta,which(beta >= 0.999),0.999)

  message("We used beta value / wanna switch to M?")
  Y <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values

  # Get all t statistic of the DMP according to the chosen (apply also for longitudinal)
  # rawBeta <- computeDMP(beta = ori_beta, covariates = covariates, formula = formula, adjPVal = 1)
  set.seed(seeds)
  rawBeta <- computeDMP(beta = ori_beta, covariates = covariates, formula = formula, adjPVal = 1)
  rawBeta <- rawBeta[rownames(beta),"t"]
  names(rawBeta) <- rownames(beta)

  # Smooth
  rawBeta_fitted <- methyl_smooth(y = rawBeta, x = pos, mycluster = cluster, bpspan=bpspan)

  # bootstrap or permuation
  set.seed(seeds)
  probe_permutation <- replicate(B,sample(rawBeta))

  # Smoothing on Permutation Beta, This process if very slow...
  permBeta <- methyl_smooth(probe_permutation,pos,cluster,bpspan=bpspan,cores=cores)

  # Set cutoff.
  message("[Paired DMR] Setting Cutoff.")
  if(pickCutoff == FALSE & is.null(cutoff)) {
    message("User setted cutoff value here.")
    cutoff <- cutoff
  }else{
    message("User did not set cutoff, 0.99 quantile of all permutation smooth value will be used as cutoff.")
    cutoff <- quantile(abs(permBeta), 0.99, na.rm = TRUE)
  }
  message(sprintf("[Paired DMR] cutoff: %s",round(cutoff, 3)))

  # Start to find bumps
  tab <- regionFinder(x = rawBeta_fitted, chr = chr, pos = pos, cluster = cluster,cutoff = cutoff, verbose = FALSE)
  message(sprintf("[Paired DMR] Found %s bumps.",nrow(tab)))

  chunksize <- ceiling(B/cores)
  subMat <- NULL
  require(doRNG)
  nulltabs <- foreach(subMat = iter(permBeta, by = "col", chunksize = chunksize),.combine = "c", .packages = "bumphunter") %dorng% {
    apply(subMat, 2, regionFinder, chr = chr, pos = pos,cluster = cluster, cutoff = cutoff,verbose = FALSE)}

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
  mytots <- apply(tab[,c("L","value")],1,function(i) sapply(c(1:B),function(x) sum(greaterOrEqual(L[[x]],i[1]) & greaterOrEqual(abs(V[[x]]),abs(i[2])))))
  mytots2 <- t(sapply(tab$area,function(i) sapply(c(1:B),function(x) sum(greaterOrEqual(A[[x]],i)))))
  rate1 <- apply(mytots,2,function(x) mean(x>0))
  pvalues1 <- apply(mytots,2,function(x) sum(x))/sum(sapply(nulltabs, nrow))
  rate2 <- rowMeans(mytots2 > 0)
  pvalues2 <- rowSums(mytots2)/sum(sapply(nulltabs, nrow))

  # Finally
  tab$p.value <- pvalues1
  tab$fwer <- rate1
  tab$p.valueArea <- pvalues2
  tab$fwerArea <- rate2
  tab <- tab[order(tab$fwer, -tab$area), ]

  PairedDMR <- tab[which(tab$p.valueArea <= adjPvalDmr),]
  message("You detected ",nrow(PairedDMR)," Paired DMRs with P value <= ",adjPvalDmr,".")

  if(nrow(PairedDMR) == 0) stop("No Paired DMR detected.")

  rownames(PairedDMR) <- paste("DMR",1:nrow(PairedDMR),sep="_")
  PairedDMR <- data.frame(PairedDMR[,1:3],width=PairedDMR[,3]-PairedDMR[,2],strand="*",PairedDMR[,4:14])
  colnames(PairedDMR)[1:3] <- c("seqnames","start","end")

  return(PairedDMR)
}

#' Compute Differentially methylated Regions
#'
#' General function to compute Differentially methylated Regions
#'
#' @param DMP full distribution of t from DMP test
#' @param minProbes min number of probes in cluster
#' @param maxGap max gap dimension
#' @param bpspan base pair span
#' @param B=250 number of Permutation for DMR
#' @param cutoff cut off number
#' @param pickCutoff set cutoff at give quantile. Default 0.99
#' @param cores multi core to use
#'
#' @export
#'
compute_DMR_from_DMP <- function(DMP,
    # probe_annotation = get("probe.features"),
    minProbes = 7,
    maxGap = 300,
    bpspan=1000,
    B=250,
    cutoff=NULL,
    pickCutoff=0.99,
    cores = 3,
    adjPvalDmr=0.05,
    seeds=1234) {

  cpg.idx <- intersect(rownames(DMP),rownames(probe.features)) #
  Anno <- probe.features[cpg.idx,]
  Anno <- Anno[order(Anno$CHR,Anno$MAPINFO),]
  cpg.idx <- rownames(Anno)
  cl <- bumphunter::clusterMaker(Anno$CHR,Anno$MAPINFO,maxGap=maxGap)
  names(cl) <- cpg.idx
  bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl) > minProbes)))]


  # beta <- beta[bumphunter.idx,]
  chr <- Anno[bumphunter.idx,]$CHR
  pos <- Anno[bumphunter.idx,]$MAPINFO
  cluster <- cl[bumphunter.idx]
  # Get all t statistic of the DMP according to the chosen (apply also for longitudinal)
  # rawBeta <- computeDMP(beta = ori_beta, covariates = covariates, formula = formula, adjPVal = 1)
  DMP <- DMP[bumphunter.idx,]
  rawBeta <- DMP[,"t"]
  names(rawBeta) <- rownames(DMP)

  # Smooth
  rawBeta_fitted <- methyl_smooth(y = rawBeta, x = pos, mycluster = cluster, bpspan=bpspan)

  # bootstrap or permuation
  set.seed(seeds)
  probe_permutation <- replicate(B,sample(rawBeta))
  # Smoothing on Permutation Beta, This process if very slow...
  permBeta <- methyl_smooth(probe_permutation,pos,cluster,bpspan=bpspan,cores=cores)

  # Set cutoff.
  if (is.null(cutoff)) {
    cutoff <- quantile(abs(permBeta), pickCutoff, na.rm = TRUE)
  }
  # Start to find bumps
  tab <- bumphunter::regionFinder(x = rawBeta_fitted, chr = chr, pos = pos,
                      cluster = cluster,cutoff = cutoff, verbose = FALSE)
  # message(sprintf("[Paired DMR] Found %s bumps.",nrow(tab)))

  chunksize <- ceiling(B/cores)
  subMat <- NULL
  require(doRNG)
  nulltabs <- foreach(subMat = iter(permBeta, by = "col", chunksize = chunksize),
                      .combine = "c", .packages = "bumphunter") %dorng% {
                        apply(subMat, 2, regionFinder, chr = chr, pos = pos,
                              cluster = cluster, cutoff = cutoff,verbose = FALSE)
                        }

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
  mytots <- apply(tab[,c("L","value")],1,function(i) sapply(c(1:B), function(x) sum(greaterOrEqual(L[[x]],i[1]) & greaterOrEqual(abs(V[[x]]),abs(i[2])))))
  mytots2 <- t(sapply(tab$area,function(i) sapply(c(1:B),function(x) sum(greaterOrEqual(A[[x]],i)))))
  rate1 <- apply(mytots,2,function(x) mean(x>0))
  pvalues1 <- apply(mytots,2,function(x) sum(x))/sum(sapply(nulltabs, nrow))
  rate2 <- rowMeans(mytots2 > 0)
  pvalues2 <- rowSums(mytots2)/sum(sapply(nulltabs, nrow))

  # Finally
  tab$p.value <- pvalues1
  tab$fwer <- rate1
  tab$p.valueArea <- pvalues2
  tab$fwerArea <- rate2
  tab <- tab[order(tab$fwer, -tab$area), ]

  PairedDMR <- tab[which(tab$p.valueArea <= adjPvalDmr),]
  # message("You detected ",nrow(PairedDMR)," Paired DMRs with P value <= ",adjPvalDmr,".")

  if(nrow(PairedDMR) == 0) stop("No Paired DMR detected.")

  rownames(PairedDMR) <- paste("DMR",1:nrow(PairedDMR),sep="_")
  PairedDMR <- data.frame(PairedDMR[,1:3],width=PairedDMR[,3]-PairedDMR[,2],strand="*",PairedDMR[,4:14])
  colnames(PairedDMR)[1:3] <- c("seqnames","start","end")

  return(PairedDMR)

}
