###
# Height normalization
###


library(data.table)
library(lme4)
library(RNOmni)

dfFull = fread(".tsv") # genotypes and phenotypes table
local_path = ""
dim(dfFull)


dfFull[, Centre:=as.factor(Centre)]
dfFull[, Gbatch:=as.factor(Gbatch)]
dfFull$year_visit=relevel(as.factor(dfFull$year_visit), ref = "2008")
dfFull$income=relevel(as.factor(dfFull$Income), ref = "3")

dfFull[, within_sex_norm_score := (score - mean(score)) / sd(score), by = c("analysis_group", "Sex")]
dfFull[, within_sex_std_height := (Height - mean(Height)) / sd(Height), by = c("analysis_group", "Sex")]
dfFull[, within_sex_rank_norm_height := rankNorm(Height), by = c("analysis_group", "Sex")]


# which model to use for RP calculation??

fit_within_sex_std_height = lmer(within_sex_std_height ~ within_sex_norm_score + income + Year_of_birth + year_visit + as.factor(analysis_group) + 
                                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + (1 | Centre) +  (1 | Gbatch), data = dfFull)
dfFull$RP_within_sex_std = scale(
  predict(fit_within_sex_std_height) -
    summary(fit_within_sex_std_height)$coef["within_sex_norm_score", "Estimate"] *
    dfFull$within_sex_norm_score
)



fit_within_sex_rank_norm_height = lmer(within_sex_rank_norm_height ~ within_sex_norm_score + income + Year_of_birth + year_visit + as.factor(analysis_group) + 
                                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + (1 | Centre) +  (1 | Gbatch), data = dfFull)
dfFull$RP_within_sex_std_inv = scale(
  predict(fit_within_sex_rank_norm_height) -
    summary(fit_within_sex_rank_norm_height)$coef["within_sex_norm_score", "Estimate"] *
    dfFull$within_sex_norm_score
)



summary(lm(formula = "within_sex_std_height ~ within_sex_norm_score * RP_within_sex_std", data = dfFull))


summary(lm(formula = "within_sex_rank_norm_height ~ within_sex_norm_score * RP_within_sex_std_inv", data = dfFull[analysis_group=="English"]))




