#' Differentially methylated longitudinal test
#'
#' @param beta normalized beta values (methylation)
#' @param pair the vector with paired samples.
#' @param pheno the phenotype of interest
#' @param covariates the covariates
#'
#' @export
#'
champ.PairedDMP <- function(beta, pair, pheno,
                            covariates = NULL,
                            adjPVal = 0.05,
                            adjust.method = "BH",
                            compare.group = NULL,
                            arraytype = c("EPIC", "450K"))
{

  arraytype = arraytype[1]

  if(length(unique(pheno))<=1) {
    stop("Invadil pheno: at least two class needed")
  } else {
    warning ("Remove this")
    message("<< Your pheno information contains following groups. >>")
    sapply(unique(pheno),function(x) message("<",x,">:",sum(pheno==x)," samples."))
    message("[The power of statistics analysis on groups contain very few samples may not strong.]\n")
  }

  if(is.null(compare.group)) {
    message("You did not assign compare groups. The first two groups: <",unique(pheno)[1],"> and <",unique(pheno)[2],">, will be compared automatically.")
    compare.group <- unique(pheno)[1:2]
  } else if(sum(compare.group %in% unique(pheno))==2) {
    message("As you assigned, champ.PairedDMP will compare ",compare.group[1]," and ",compare.group[2],".")
  }else{
    message("Seems you did not assign correst compare groups. The first two groups: <",unique(pheno)[1],"> and <",unique(pheno)[2],">, will be compared automatically.")
    compare.group <- unique(pheno)[1:2]
  }

  p <- pheno[which(pheno %in% compare.group)]
  compare.pair <- pair[which(pheno %in% compare.group)]

  if(!all(table(compare.pair)==2))
    stop("Valid Pairs for compare sampels are required. Odd numbers have been detected in your compared data. But in paired information, each patient's name should appear exactly twice.")
  if(!all(table(compare.pair,p)==1))
    stop("Pheno and Pairs must corrsponding to each other. The match between your paired information and pheno is not correct.")


  beta <- beta[,which(pheno %in% compare.group)]

  beta_1 <- beta[,p == compare.group[1]][,order(compare.pair[p == compare.group[1]])]
  beta_2 <- beta[,p == compare.group[2]][,order(compare.pair[p == compare.group[2]])]
  beta.sub <- beta_1 - beta_2
  PairedDMP <- limma::topTable(limma::eBayes(limma::lmFit(beta.sub)),number=nrow(beta.sub),adjust.method=adjust.method,p.value=adjPVal)

  message("You have found ",sum(PairedDMP$adj.P.Val <= adjPVal), " significant paired DMPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")

  if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
  com.idx <- intersect(rownames(PairedDMP),rownames(probe.features))
  avg.substract <- rowMeans(beta.sub[com.idx,])
  PairedDMP <- data.frame(PairedDMP[com.idx,],Ave_Sub=avg.substract,probe.features[com.idx,])

  return(PairedDMP)
}

