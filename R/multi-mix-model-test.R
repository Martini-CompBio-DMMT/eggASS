
test <- function() {

  load("../epigenetic-gwas-Baune/Psych-EMDR-cohort/smallMet_example.RData")
  beta <- beta_values
  row.names(covariates) <- covariates$Sample_Name

  devtools::load_all()

  require(lme4)
  require(lmerTest)

  check_samples_covariates(beta, covariates)

  formula = paste0("score ", formula)

  single_lmer_test_two_class(beta[1,], formula, covariates, "time")
  single_lmer_test_two_class_cmp(beta[1,], formula, covariates, "time")

  covariates$score <- beta[1,]
  results <- lmer(score~time*Smoker+(1|patient), data=covariates)
  # results <- lmer(Scores~Factor1*Factor2+(1|Subject), data=reg_data)
  r<- anova(results)
  as.data.frame(r)
  require(emmeans)
  emmeans(results, list(pairwise ~ time), adjust = "bonferroni")  #post hoc (two correction egs)
  emmeans(results, list(pairwise ~ time + Smoker), adjust = "bonferroni")  #post hoc (two correction egs)
}


single_lmer_test_two_class <- function(meth, formula, data, varOfInterest) {
  # formula="score~time+(1|patient)"
  # data=covariates
  data$score <- meth
  results <- lmer(as.formula(formula), data=data, REML = FALSE)
  r<- anova(results)
  as.data.frame(r)[varOfInterest, , drop=F]
}

single_lmer_test_two_class_cmp <- function(meth, formula, data, varOfInterest) {
  # formula="score~time+(1|patient)"
  # data=covariates
  data$score <- meth
  results <- lmer(as.formula(formula), data=data, REML = FALSE)
  emmeans::emmeans(results, list(as.formula(paste0("pairwise ~ ", varOfInterest))), adjust = "bonferroni")
}


