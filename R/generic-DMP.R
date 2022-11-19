#' Differentially methylated longitudinal test
#'
#' General function to compute longitudinal test
#'
#' @param beta normalized beta values (methylation)
#' @param covariates the data.frame wit all covariates
#' @param formula the formula; last term is the considered phenotype
#' @param responseVar if NULL, use last var of the formula.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust pvalue method
#' @param useM use M values rather than beta
#' @param capEdges when using M values, beta are capped to limits
#' @param probe_annotation for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @export
#'
computeDMP <- function(beta, covariates,
                       formula="~ phenotype",
                       responseVar=NULL,
                       adjPVal = 1,
                       adjust.method = "BH",
                       useM = FALSE,
                       capEdges = c(0.001,0.999),
                       probe_annotation = get("probe.features")) {

  check_samples_covariates(beta, covariates)

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

  design <- model.matrix(as.formula(formula), covariates)
  fit <- limma::eBayes(limma::lmFit(beta, design))

  if (useM){
    beta<- replace(beta,which(beta <= capEdges[1]), capEdges[1])
    beta <- replace(beta,which(beta >= capEdges[2]),capEdges[2])
    beta <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values
  }
  # dmpFinder() # from minfi
  DMP <- limma::topTable(fit, number=nrow(beta),adjust.method=adjust.method,
                         p.value=adjPVal, coef=coef_idx)

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
#' @param beta normalized beta values (methylation)
#' @param pair name of the variable with pairs
#' @param time name of time variable to be compared
#' @param covariates the data.frame wit all covariate
#' @param formula the formula; last term is the considered phenotype
#' @param coef the coefficient to be used for pvalue; NULL = last covariate.
#' @param adjPVal adjusted p-value cut off
#' @param adjust.method adjust p-value method
#' @param probe_annotation for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @export
#'
computeLongitudinalDMP <- function(beta, covariates,
                                   pair_var = "patient",
                                   time_var = "time",
                                   formula=NULL,
                                   coef = NULL,
                                   adjPVal = 0.05,
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
  beta.sub <- beta_1 - beta_2

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

    PairedDMP <- limma::topTable(fit, number=nrow(beta),adjust.method=adjust.method, p.value=adjPVal, coef=coefs)
  }

  message("You have found ",sum(PairedDMP$adj.P.Val <= adjPVal), " significant paired DMPs with a ",
          adjust.method," adjusted P-value below ", adjPVal,".")

  com.idx <- intersect(rownames(PairedDMP),rownames(probe_annotation))
  avg.substract <- rowMeans(beta.sub[com.idx,,drop=F])
  PairedDMP <- data.frame(PairedDMP[com.idx,,drop=F],Ave_Sub=avg.substract,probe_annotation[com.idx,,drop=F])
  return(PairedDMP)
}

#' Differentially methylated longitudinal test with blocking factor
#'
#' General function to compute longitudinal tests
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
#' @param probe_annotation for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @export
#'
compute_blocked_DMP <- function(beta, covariates,
                                block,
                                formula="~ phenotype",
                                responseVar=NULL,
                                adjPVal = 0.05,
                                adjust.method = "BH",
                                useM = FALSE,
                                capEdges = c(0.001,0.999),
                                probe_annotation = get("probe.features")) {

  check_samples_covariates(beta, covariates)

  if (!(block %in% colnames(covariates)))
    stop(paste0(block, " blocking factor is not present in covariates"))

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

  if (useM){
    beta<- replace(beta,which(beta <= capEdges[1]), capEdges[1])
    beta <- replace(beta,which(beta >= capEdges[2]),capEdges[2])
    beta <- log((beta/(1-beta)),2) # M=log2(beta)/(1-beta) Convert to M values
  }

  dupcor <- limma::duplicateCorrelation(beta, design, block=covariates[[block]])
  fitDupCor <- limma::lmFit(beta, design, block=covariates[[block]], correlation=dupcor$consensus)
  fitDupCor <- limma::eBayes(fitDupCor)
  DMP <- limma::topTable(fitDupCor, number=nrow(beta),adjust.method=adjust.method,
                         p.value=adjPVal, coef=coef_idx)

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
#' General function to compute longitudinal tests
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
#' @param probe_annotation for EPIC array data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#'
#' @export
#'
compute_blocked_dream_DMP <- function(beta, covariates,
                                formula_block = "~ time + (1|patient)",
                                formula = "~ time",
                                responseVar = NULL,
                                adjPVal = 1,
                                adjust.method = "BH",
                                useM = FALSE,
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

annotate_CpGs <- function(DMP, probe_annotation = get("probe.features")) {
  if (!(nrow(DMP)))
    stop("No DMP provided.")

  com.idx <- intersect(rownames(DMP),rownames(probe_annotation))
  data.frame(DMP[com.idx,,drop=F], probe_annotation[com.idx,,drop=F])
}
