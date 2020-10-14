"Residual variance as a function of predicted height in the test data of XXX individuals (cm.) 
(A) Residual variance. 
Red: prediction is made using an assumption that height is distributed normally. 
Green: prediction is made using an assumption that log-height is distributed normally. 
(B) Ratio between residual variances of normal and log-normal models."

library(data.table)
library(stringr)
library('R.filesets')
current.user = Sys.info()['user']

local_path = ""
gen_and_phen_path = ".tsv" # Table with phenotypes and genotypes 

gen_and_phen = loadRDS('.rds') # table with phenotypes and PGHS
gen_and_phen[, norm_score_meta := (score_meta - mean(score_meta)) /
       sd(score_meta), by = analysis_group]

RP1 <- fread(".tsv") # table with phentypes and RP
gen_and_phen[, predF3.1 := RP1$predF3]


get_rsq_table = function(basic_formula){
  covariates = basic_formula
  trained_model = lm(as.formula(paste("Height", covariates , sep = " ~ ")), gen_and_phen[analysis_group == "English"])
  trained_model_log = lm(as.formula(paste("log_height", covariates, sep = " ~ ")), gen_and_phen[analysis_group == "English"])
  
  gen_and_phen$pred_h = predict(trained_model, gen_and_phen)
  gen_and_phen$pred_logh = predict(trained_model_log, gen_and_phen)
  rsq_table = gen_and_phen[, list(
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

system.time(t1 <- get_rsq_table("Sex"))
system.time(t2 <- get_rsq_table("Sex + predF3.1"))
system.time(t3 <- get_rsq_table("Sex + predF3.1 + norm_score"))
system.time(t4 <- get_rsq_table("Sex + predF3"))
system.time(t5 <- get_rsq_table("Sex + predF3 + norm_score_meta"))

rsq_table = rbind(t1, t2, t3, t4, t5)

write.table(rsq_table, file = "rsq_table.tsv", quote = F, sep = "\t", row.names = F)


#----
gen_and_phen$pred_h = predict(trained_model, gen_and_phen)
gen_and_phen$pred_log_exp_h = 10**(predict(trained_model_log, gen_and_phen))

set.seed(666)
with(gen_and_phen[sample(.N, 10**4)], plot(x = pred_h, y= (Height - pred_h)))
set.seed(666)
with(gen_and_phen[sample(.N, 10**4)], plot(x = pred_log_exp_h, y= (Height - pred_log_exp_h)))

n_breaks = 50

orig = gen_and_phen[, sd(Height - pred_h), 
             by = cut(pred_h, breaks = quantile(pred_h, probs = seq(0, 1, by = 1/n_breaks)), labels = 1:n_breaks, right = FALSE)][order(cut)]
setnames(orig,"V1", "orig")
log_orig = gen_and_phen[,sd(Height - pred_log_exp_h), 
             by = cut(pred_log_exp_h, breaks = quantile(pred_log_exp_h, probs = seq(0, 1, by = 1/n_breaks)), labels = 1:n_breaks, right = FALSE)][order(cut)]
setnames(log_orig,"V1", "log_orig")

mdt  = merge(orig, log_orig, by = "cut")
plot(x = orig$orig, y = log_orig$log_orig)
summary(lm(log_orig~orig, mdt))
