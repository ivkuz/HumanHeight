###
# Pesidual Predictor (RP, predF3) calculation
###

library(data.table)
library(lme4)

dfFull = fread(".tsv") # table with genotypes and phenotypes
local_path = ""

dfFull[, Centre:=as.factor(Centre)]
dfFull[, Gbatch:=as.factor(Gbatch)]
dfFull$year_visit=relevel(as.factor(dfFull$year_visit), ref = "2008")
dfFull$income=relevel(as.factor(dfFull$Income), ref = "3")



dim(dfFull)
dfFull = dfFull[!is.na(Income)]
dim(dfFull)
dfFull[, norm_score := (score - mean(score)) / sd(score), by = analysis_group]


fitFull = lmer(Height ~ Sex + norm_score + income + Year_of_birth + year_visit + as.factor(analysis_group) + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 +
                 (1 | Centre) +  (1 | Gbatch), na.action = na.exclude, data = dfFull)
summary(fitFull)



predH = predict(fitFull)
predF3 = predH - summary(fitFull)$coef[2,"Estimate"]*as.numeric(dfFull$Sex) - summary(fitFull)$coef[3,"Estimate"]*dfFull$norm_score
dfFull$predF3=scale(predF3)



if (current.user=="sslavskii"){
  write.table(dfFull, file = "/storage/projects/UKB-nonAdditive-41601/shared/genotypes_and_phenotypes_v4.tsv", 
              sep="\t", quote = F, row.names = F)
} else{
}



fit_log = lmer(log_height ~ Sex + norm_score + income + Year_of_birth + year_visit + as.factor(analysis_group) + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 +
                 (1 | Centre) +  (1 | Gbatch), na.action = na.exclude, data = dfFull)
summary(fit_log)

predLH = predict(fit_log)

predF3_log = predLH - summary(fit_log)$coef[2,"Estimate"]*dfFull$Sex-summary(fit_log)$coef[3,"Estimate"]*dfFull$norm_score
dfFull$predF3_log=scale(predF3_log)


dfFull$predH = scale(predH)
dfFull$predLH = scale(predLH)


sink("models.txt")

summary(lm(Height ~ Sex + norm_score, data = dfFull))
summary(lm(Height ~ Sex + norm_score + I(norm_score^2), data = dfFull))
summary(lm(Height ~ Sex + norm_score + I(norm_score^2) + I(norm_score^3), data = dfFull))

summary(lm(Height ~ predH, data = dfFull))
summary(lm(Height ~ predH + I(predH^2), data = dfFull))
summary(lm(Height ~ predH + I(predH^2) + I(predH^3), data = dfFull))

summary(lm(Height ~ Sex + norm_score, data = dfFull[analysis_group=="English"]))
summary(lm(Height ~ Sex + norm_score + I(norm_score^2), data = dfFull[analysis_group=="English"]))
summary(lm(Height ~ Sex + norm_score + I(norm_score^2) + I(norm_score^3), data = dfFull[analysis_group=="English"]))

summary(lm(Height ~ predH, data = dfFull[analysis_group=="English"]))
summary(lm(Height ~ predH + I(predH^2), data = dfFull[analysis_group=="English"]))
summary(lm(Height ~ predH + I(predH^2) + I(predH^3), data = dfFull[analysis_group=="English"]))





summary(lm(log_height ~ Sex + norm_score, data = dfFull))
summary(lm(log_height ~ Sex + norm_score + I(norm_score^2), data = dfFull))
summary(lm(log_height ~ Sex + norm_score + I(norm_score^2) + I(norm_score^3), data = dfFull))

summary(lm(log_height ~ predLH, data = dfFull))
summary(lm(log_height ~ predLH + I(predLH^2), data = dfFull))
summary(lm(log_height ~ predLH + I(predLH^2) + I(predLH^3), data = dfFull))

summary(lm(log_height ~ Sex + norm_score, data = dfFull[analysis_group=="English"]))
summary(lm(log_height ~ Sex + norm_score + I(norm_score^2), data = dfFull[analysis_group=="English"]))
summary(lm(log_height ~ Sex + norm_score + I(norm_score^2) + I(norm_score^3), data = dfFull[analysis_group=="English"]))

summary(lm(log_height ~ predLH, data = dfFull[analysis_group=="English"]))
summary(lm(log_height ~ predLH + I(predLH^2), data = dfFull[analysis_group=="English"]))
summary(lm(log_height ~ predLH + I(predLH^2) + I(predLH^3), data = dfFull[analysis_group=="English"]))





summary(lm(Height ~ Sex + predF3, data = dfFull))
summary(lm(Height ~ Sex + predF3 + I(predF3^2), data = dfFull))
summary(lm(Height ~ Sex + predF3 + I(predF3^2) + I(predF3^3), data = dfFull))

summary(lm(log_height ~ Sex + predF3_log, data = dfFull))
summary(lm(log_height ~ Sex + predF3_log + I(predF3_log^2), data = dfFull))
summary(lm(log_height ~ Sex + predF3_log + I(predF3_log^2) + I(predF3_log^3), data = dfFull))


sink(file = NULL)

jpeg('height.jpg')
plot(dfFull$log_height,dfFull$Height)
abline(a=xxx$coef[1],b=xxx$coef[2], col='green')
abline(lm(dfFull$Height ~ dfFull$log_height), col='green')
dev.off()

