#' Differentially methylated longitudinal test
#'
#' General function to compute longitudinal test.
#' Probe annotation can be done with eggASS:::annotate_CpGs using 'probe_annotation'
#' object from ChAMP for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @param beta normalized beta values (methylation)
#' @param covariates the data.frame wit all covariates
#' @param formula the formula; last term is the considered phenotype
#' @param responseVar if NULL, use last var of the formula.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust pvalue method
#' @param useM use M values rather than beta
#' @param capEdges when using M values, beta are capped to limits
#'
#' @export
#'
computeDMP <- function(beta, covariates,
                       formula="~ phenotype",
                       responseVar=NULL,
                       adjPVal = 1,
                       adjust.method = "BH",
                       useM = T,
                       capEdges = c(0.001,0.999),
                       probe_annotation = get("probe.features")) {

  check_samples_covariates(beta, covariates)

  design <- model.matrix(as.formula(formula), covariates)

  if (is.null(responseVar)){
    interesting_var <- rev(all.vars(as.formula(formula)))[1]
    coef_idx <- grep(interesting_var, colnames(design))
  } else {
    if (length(responseVar)>1)
      stop("Only one variable can be specified.")
    formulaVars <- all.vars(as.formula(formula))
    checkv <- responseVar %in% formulaVars
    if (!all(checkv))
      stop(paste0(paste(responseVar, " var is missing from formula")))
    coef_idx <- grep(responseVar, colnames(design))
  }

  if (useM) {
    message("Using M-value")
    beta <- replace(beta,which(beta <= capEdges[1]), capEdges[1])
    beta <- replace(beta,which(beta >= capEdges[2]),capEdges[2])
    beta <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values
  }

  fit <- limma::eBayes(limma::lmFit(beta, design))
  # dmpFinder() # from minfi
  DMP <- limma::topTable(fit, number=nrow(beta),adjust.method=adjust.method,
                         p.value=adjPVal, coef=coef_idx)

  # annotations
  # https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html

  # message("You have found ",sum(DMP$adj.P.Val <= adjPVal),
  #         " significant paired DMPs with a ",adjust.method,
  #         " adjusted P-value below ", adjPVal,".")

  # require(ChAMP)
  # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  # data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  # probe.features is created as well

  # if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
  # com.idx <- intersect(rownames(DMP),rownames(probe_annotation))
  # DMP <- data.frame(DMP[com.idx,,drop=F], probe_annotation[com.idx,,drop=F])

  return(DMP)
}


#' Differentially methylated longitudinal test paired samples
#'
#' Probe annotation can be done with eggASS:::annotate_CpGs using 'probe_annotation'
#' object from ChAMP for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @param beta normalized beta values (methylation)
#' @param pair name of the variable with pairs
#' @param time name of time variable to be compared
#' @param covariates the data.frame wit all covariate
#' @param formula the formula; last term is the considered phenotype
#' @param coef the coefficient to be used for pvalue; NULL = last covariate.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust p-value method
#'
#' @export
#'
computeLongitudinalDMP <- function(beta, covariates,
                                   pair_var = "patient",
                                   time_var = "time",
                                   formula=NULL,
                                   coef = NULL,
                                   adjPVal = 1,
                                   adjust.method = "BH",
                                   probe_annotation = get("probe.features")) {
#
#   arraytype = arraytype[1]

  check_samples_covariates(beta, covariates)
  check_selected_covariates(all.vars(as.formula(formula)), covariates=covariates)

  missing <- setdiff(c(time_var, pair_var), colnames(covariates))
  if (length(missing))
    stop(paste(paste(missing, collapse=" and")), "variabels are missing")

  time <- covariates[[time_var]]
  pair <- covariates[[pair_var]]

  compare.group <- levels(time)

  if (is.null(compare.group))
    stop(paste0(time_var, "variable for time must be a factor"))
  if (length(compare.group)!=2)
    stop("Only time pairs allowed. Use general model to compare more that two times.")

  p <- time[which(time %in% compare.group)]
  compare.pair <- pair[which(time %in% compare.group)]

  if(!all(table(compare.pair)==2))
    stop("Valid Pairs for compare sampels are required. Odd numbers have been detected in your compared data. But in paired information, each patient's name should appear exactly twice.")
  if(!all(table(compare.pair,p)==1))
    stop("time and Pairs must corrsponding to each other. The match between your paired information and time is not correct.")

  beta <- beta[,which(time %in% compare.group)]

  beta_1 <- beta[,p == compare.group[1]][,order(compare.pair[p == compare.group[1]])]
  beta_2 <- beta[,p == compare.group[2]][,order(compare.pair[p == compare.group[2]])]
  beta.sub <- beta_2 - beta_1

  if (is.null(formula)){
    PairedDMP <- limma::topTable(limma::eBayes(limma::lmFit(beta.sub)),number=nrow(beta.sub),adjust.method=adjust.method,p.value=adjPVal)
  } else {
    covariates.sub <- covariates[colnames(beta.sub), ]
    design <- model.matrix(as.formula(formula), covariates.sub)
    varOfInterest <- rev(all.vars(as.formula(formula)))[1]
    coef_idx <- grep(varOfInterest, colnames(design))

    fit <- limma::eBayes(limma::lmFit(beta.sub, design))
    if (is.null(coef))
      coefs = coef_idx

    PairedDMP <- limma::topTable(fit, number=nrow(beta.sub),adjust.method=adjust.method, p.value=adjPVal, coef=coefs)
  }

  # message("You have found ",sum(PairedDMP$adj.P.Val <= adjPVal), " significant paired DMPs with a ",
  #         adjust.method," adjusted P-value below ", adjPVal,".")

  # com.idx <- intersect(rownames(PairedDMP),rownames(probe_annotation))
  # avg.substract <- rowMeans(beta.sub[com.idx,,drop=F])
  # PairedDMP <- data.frame(PairedDMP[com.idx,,drop=F],Ave_Sub=avg.substract,probe_annotation[com.idx,,drop=F])
  return(PairedDMP)
}

#' Differentially methylated longitudinal test with blocking factor
#'
#' General function to compute longitudinal tests with blocking factors
#' Be aware that we used a random sample of 100,000 probes to speed up analysis
#' Probe annotation can be done with eggASS:::annotate_CpGs using 'probe_annotation'
#' object from ChAMP for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @param beta normalized beta values (methylation)
#' @param covariates the data.frame wit all covariates
#' @param block blocking factor
#' @param formula the formula; last term is the considered phenotype
#' @param responseVar if NULL, use last var of the formula.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust pvalue method
#' @param useM use M rather than beta values
#' @param capEdges when using M values, beta are capped to limits
#'
#' @export
#'
compute_blocked_DMP <- function(beta, covariates,
                                block,
                                formula="~ phenotype",
                                responseVar=NULL,
                                adjPVal = 1,
                                adjust.method = "BH",
                                useM = FALSE,
                                capEdges = c(0.001,0.999),
                                seeds=1234,
                                probe_annotation = get("probe.features")) {

  check_samples_covariates(beta, covariates)

  if (!(block %in% colnames(covariates)))
    stop(paste0(block, " blocking factor is not present in covariates"))

  design <- model.matrix(as.formula(formula), covariates)

  if (is.null(responseVar)){
    interesting_var <- rev(all.vars(as.formula(formula)))[1]
    coef_idx <- grep(interesting_var, colnames(design))
    if (length(coef_idx) >1)
      warning(paste0("p-value computed using all coefficient containing ", interesting_var, " has been selected."))
  } else {
    if (length(responseVar)>1)
      stop("Only one variable can be specified.")
    formulaVars <- all.vars(as.formula(formula))
    checkv <- responseVar %in% formulaVars
    if (!all(checkv))
      stop(paste0(paste(responseVar, " var is missing from formula")))
    coef_idx <- grep(responseVar, colnames(design))
  }

  if (useM){
    message("Using M-value")
    beta <- convertToM(beta, capEdges=capEdges)
  }

  to_pick <- extractRandomProbes(beta, sample_size=100000, seeds=seeds)

  dupcor <- limma::duplicateCorrelation(beta[to_pick,], design, block=covariates[[block]])
  fitDupCor <- limma::lmFit(beta, design, block=covariates[[block]], correlation=dupcor$consensus)
  fitDupCor <- limma::eBayes(fitDupCor)

  total <- limma::topTable(fitDupCor, number=nrow(beta),adjust.method=adjust.method,
                          p.value=adjPVal, coef=coef_idx)
  DMP <- list(total=total)

  if (length(coef_idx)>1) {
    DMP_coefs  <- lapply(coef_idx, function(coefficient) {
      limma::topTable(fitDupCor, number=nrow(beta),adjust.method=adjust.method,
                             p.value=adjPVal, coef=coefficient)
    })
    names(DMP_coefs) <- colnames(design)[coef_idx]
    DMP <- c(DMP, DMP_coefs)
  }

  # message("You have found ",sum(DMP$adj.P.Val <= adjPVal),
  #         " significant paired DMPs with a ",adjust.method,
  #         " adjusted P-value below ", adjPVal,".")

  # require(ChAMP)
  # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  # data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  # probe.features is created as well

  # if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
  # com.idx <- intersect(rownames(DMP),rownames(probe_annotation))
  # DMP <- data.frame(DMP[com.idx,,drop=F], probe_annotation[com.idx,,drop=F])

  return(DMP)
}

#' Differentially methylated longitudinal test with blocking factor
#'
#' General function to compute longitudinal tests.
#' Probe annotation can be done with eggASS:::annotate_CpGs using 'probe_annotation'
#' object from ChAMP for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @param beta normalized beta values (methylation)
#' @param covariates the data.frame wit all covariates
#' @param block blocking factor
#' @param formula the formula; last term is the considered phenotype
#' @param responseVar if NULL, use last var of the formula.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust pvalue method
#' @param useM use M rather than beta values
#' @param capEdges when using M values, beta are capped to limits
#'
#' @export
#'
compute_blocked_dream_DMP <- function(beta, covariates,
                                formula_block = "~ time + (1|patient)",
                                formula = "~ time",
                                responseVar = NULL,
                                adjPVal = 1,
                                adjust.method = "BH",
                                useM = T,
                                capEdges = c(0.001,0.999),
                                probe_annotation = get("probe.features")) {

  require("variancePartition")
  check_samples_covariates(beta, covariates)
  check_selected_covariates(all.vars(as.formula(formula_block)), covariates=covariates)
  check_selected_covariates(all.vars(as.formula(formula)), covariates=covariates)

  if (useM){
    beta<- replace(beta,which(beta <= capEdges[1]), capEdges[1])
    beta <- replace(beta,which(beta >= capEdges[2]),capEdges[2])
    beta <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values
  }

  # design <- model.matrix(as.formula(formula), covariates)
  fitmm = variancePartition::dream(beta, formula_block, covariates)

  if (is.null(responseVar)){
    interesting_var <- rev(all.vars(as.formula(formula)))[1]
    coef_idx <- grep(interesting_var, colnames(fitmm$design))
  } else {
    if (length(responseVar)>1)
      stop("Only one variable can be specified.")

    formulaVars <- all.vars(as.formula(formula))
    checkv <- responseVar %in% formulaVars
    if (!all(checkv))
      stop(paste0(paste(responseVar, " var is missing from formula")))
    coef_idx <- grep(responseVar, colnames(fitmm$design))
  }

  fitmm = variancePartition::eBayes(fitmm)
  DMP = variancePartition::topTable( fitmm, coef=coef_idx, number=nrow(beta),
                                     adjust.method=adjust.method, p.value=adjPVal)

  # require(ChAMP)
  # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  # data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  # probe.features is created as well

  # if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
  # com.idx <- intersect(rownames(DMP),rownames(probe_annotation))
  # DMP <- data.frame(DMP[com.idx,,drop=F], probe_annotation[com.idx,,drop=F])

  return(DMP)
}


check_selected_covariates <- function(covs, covariates) {
  missing <- setdiff(covs, colnames(covariates))
  if (length(missing))
    stop(paste(paste(missing, collapse=" and ")), " variabel/s are missing")
}

check_samples_covariates <- function(beta, covariate) {
  if(!identical(colnames(beta), row.names(covariate)))
    stop("Miscmatch in sample names")
}

#'
#' @export
#'
annotate_CpGs <- function(DMP, probe_annotation = get("probe.features")) {
  if (!(nrow(DMP)))
    stop("No DMP provided.")

  com.idx <- intersect(rownames(DMP),rownames(probe_annotation))
  data.frame(DMP[com.idx,,drop=F], probe_annotation[com.idx,,drop=F])
}


#'
#' @export
#'
plotProbeMethylation <- function(probeId, y, covariates=covariates[, c("patient", "time", "DeltaMADRS.T8.T0")]) {
  require(ggplot2)

  library(ggplot2)
  data <- data.frame(methyl=y[probeId, ], covariates)
  data <- data[order(data$patient,data$time),]
  identical(as.character(data$patient[data$time=="T8"]), as.character(data$patient[data$time=="T0"]))
  diff <- data$methyl[data$time=="T8"]-data$methyl[data$time=="T0"]
  names(diff) <- as.character(data$patient[data$time=="T8"])
  diff <- diff[order(diff)]

  data$patient <- factor(as.character(data$patient), levels=names(diff))
  data <- data[order(data$patient,data$time),]
  data <- data[order(data$DeltaMADRS.T8.T0,data$patient,data$time),]
  data$patient <- factor(as.character(data$patient), levels=unique(data$patient))
  ggplot(data, aes(x=patient, y=methyl, color=time, size=DeltaMADRS.T8.T0)) +
    geom_point()


}


