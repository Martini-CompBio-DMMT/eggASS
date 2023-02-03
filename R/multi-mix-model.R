# ### Adapted from
# ###  http://rcompanion.org/handbook/G_03.html
# ###  http://rcompanion.org/handbook/G_06.html
#
# library(lme4)
# library(lmerTest)
# library(lsmeans)
# library(multcompView)
# library(car)
# library(multcomp)
#
# Input = ("
#  Individual  Time   SBP    SL
#  A           T1     17.5    6
#  B           T1     18.4    7
#  C           T1     16.2    3
#  D           T1     14.5    4
#  E           T1     15.5    5
#  F           T1     18.9    6
#  G           T1     19.5    7
#  H           T1     21.1    8
#  I           T1     17.8    9
#  J           T1     16.8   10
#  K           T1     18.4    9
#  L           T1     17.3    8
#  M           T1     18.9    7
#  N           T1     16.4    6
#  O           T1     17.5    5
#  P           T1     15.0    4
#  A           T2     17.6    9
#  B           T2     18.5    10
#  C           T2     15.9    7
#  D           T2     14.9    5
#  E           T2     15.7    5
#  F           T2     18.9    10
#  G           T2     19.5    10
#  H           T2     21.5    10
#  I           T2     18.5   10
#  J           T2     17.1   10
#  K           T2     18.9   10
#  L           T2     17.5   10
#  M           T2     19.5   10
#  N           T2     16.5    8
#  O           T2     17.4    8
#  P           T2     15.6    5
# ")
#
# Data = read.table(textConnection(Input),header=TRUE)
#
# str(Data)
#
# ##############
#
# library(lme4)
#
# library(lmerTest)
#
# model = lmer(SBP ~ Time + SL  + (1|Individual),
#              data=Data,
#              REML=TRUE)
#
# anova(model)
#
# rand(model)
#
# ###############
#
# library(car)
# scatterplot(SBP ~ SL | Time, data = Data, smooth=F, reg.line=F)
#
# ###############
#
# library(multcompView)
# library(lsmeans)
#
# leastsquare = emmeans(model,
#                       pairwise ~ Time,
#                       adjust="tukey")
#
# CLD = emmeans:::cld.emm_list(leastsquare,
#                              alpha=0.05,
#                              Letters=letters)
#
# CLD
#
# ###############################
#
# library(ggplot2)
#
# qplot(x    = Time ,
#       y    = lsmean,
#       data = CLD) +
#
#   geom_errorbar(aes(
#     ymin  = lower.CL,
#     ymax  = upper.CL,
#     width = 0.15))
#
# ###############################
#
#
# T1  = Data$SBP[Data$Time=="T1"]
#
# T2  = Data$SBP[Data$Time=="T2"]
#
# Difference = T2 - T1
#
# X = Data$Individual[Data$Time=="T1"]
#
#
# barplot(Difference,
#         names.arg = X,
#         col="dark gray",
#         xlab="Individual",
#         ylab="Difference (T2 ñ T1)")
#
# #################################
#
# print(data.frame(Time = c("T1", "T2"), Mean = c(mean(T1), mean(T2))))
#
# ##################################################
#
#
# hist(residuals(model), col="darkgray")
#
# plot(predict(model), residuals(model))
