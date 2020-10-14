###
# analysis of PGHS2 (Yengo et al.)
###

library(data.table)
local_path = ""
gen_and_phen = readRDS( paste(local_path, '.rds', sep = "")) # table with phenotypes and PGHS2

summary(lm(data = gen_and_phen[Sex==1], formula = Height ~ score_meta)) 
summary(lm(data = gen_and_phen[Sex==1], formula = log_height ~ score_meta)) 

set.seed(1)
men_sample = gen_and_phen[analysis_group == "English"][Sex == 1][sample(.N, 10**4)]
set.seed(1)
women_sample = gen_and_phen[analysis_group == "English"][Sex == 0][sample(.N, 10**4)]

mean(men_sample$score_meta)
sd(men_sample$score_meta)

mean(women_sample$score_meta)
sd(women_sample$score_meta)

library(car)
dt = rbind(men_sample, women_sample)

leveneTest(dt$score_meta, group = as.factor(dt$Sex))
t.test(men_sample$score_meta, women_sample$score_meta)

library(ggplot2)
ggplot(data = dt, aes(x = score_meta, fill = as.factor(Sex))) + geom_histogram(position="identity", alpha=0.5, bins = 100) + theme_bw()


#----
FiellerRatioCI_uncorrlated <- function(a, b, se_a, se_b, alpha=0.05){
  theta <- a/b
  v11 <- se_a**2
  v12 <- 0
  v22 <- se_b**2
  
  z <- qnorm(1-alpha/2)
  g <- z*v22/b^2
  #print(g)
  C <- sqrt(v11/a^2 - 2*v12/(a*b) + v22/b^2 - g/a^2*(v11-v12^2/v22))
  #print(C)
  minS <- (1/(1-g))*(theta- g*v12/v22 - z * theta * C)
  maxS <- (1/(1-g))*(theta- g*v12/v22 + z * theta * C)
  #print(z * theta * C)
  return(c(ratio=theta,min=minS,max=maxS))
}

get_beta_ratio = function(model1, model2){
  beta1 = summary(model1)$coefficients['score_meta','Estimate']
  se1 = summary(model1)$coefficients['score_meta','Std. Error']
  beta2 = summary(model2)$coefficients['score_meta','Estimate']
  se2 = summary(model2)$coefficients['score_meta','Std. Error']
  return(FiellerRatioCI_uncorrlated(beta1, beta2, se1, se2))
}
std <- function(x) sd(x)/sqrt(length(x))

m_h = lm(data = gen_and_phen[analysis_group == "English"][Sex == 1], formula = Height ~ score_meta)
w_h = lm(data = gen_and_phen[analysis_group == "English"][Sex == 0], formula = Height ~ score_meta)

get_beta_ratio(m_h, w_h)
FiellerRatioCI_uncorrlated(mean(gen_and_phen[analysis_group == "English"][Sex == 1]$Height), 
                           mean(gen_and_phen[analysis_group == "English"][Sex == 0]$Height), 
                           std(gen_and_phen[analysis_group == "English"][Sex == 1]$Height), 
                           std(gen_and_phen[analysis_group == "English"][Sex == 0]$Height))

m_lh = lm(data = gen_and_phen[analysis_group == "English"][Sex == 1], formula = log_height ~ score_meta)
w_lh = lm(data = gen_and_phen[analysis_group == "English"][Sex == 0], formula = log_height ~ score_meta)

get_beta_ratio(m_lh, w_lh)


# variance explained

summary(m_h)$r.squared
summary(w_h)$r.squared

get_rsq_table = function(basic_formula, dt){
  trained_model = lm(as.formula(paste("Height", basic_formula , sep = " ~ ")), dt)
  trained_model_log = lm(as.formula(paste("log_height", basic_formula, sep = " ~ ")), dt)
  
  dt$pred_h = predict(trained_model, dt)
  dt$pred_logh = predict(trained_model_log, dt)
  rsq_table = dt[, list(
    r_sq_original = 1 - (var(Height - pred_h)) / var(Height),
    r_sq_log_exp = 1 - (var(Height - 10 ** pred_logh)) / var(Height),
    r_sq_log = 1 - (var(log_height - pred_logh)) / var(log_height),
    r_sq_original_log = 1 - (var(log_height - log10(pred_h))) / var(log_height),
    N = .N),
    by = analysis_group]
  rsq_table[, diff_orig_log_exp := r_sq_log_exp - r_sq_original]
  rsq_table[, diff_log_original_log := r_sq_log - r_sq_original_log]
  rsq_table$model = basic_formula
  
  return(rsq_table[order(N, decreasing = T)])
}

get_rsq_table("score_meta",gen_and_phen[analysis_group == "English"][Sex == 1])
get_rsq_table("score_meta",gen_and_phen[analysis_group == "English"][Sex == 0])


#-----
percentiles = c(0.1, 1, 5, 95, 99, 99.9)

get_prediction = function(model, dt){
  model$coefficients['(Intercept)'] + model$coefficients['score_meta'] * quantile(dt$score_meta, percentiles/100)
}

get_prediction(m_h, gen_and_phen[analysis_group == "English"][Sex == 1])
get_prediction(w_h, gen_and_phen[analysis_group == "English"][Sex == 0])

10**get_prediction(m_lh, gen_and_phen[analysis_group == "English"][Sex == 1])
10**get_prediction(w_lh, gen_and_phen[analysis_group == "English"][Sex == 0])

