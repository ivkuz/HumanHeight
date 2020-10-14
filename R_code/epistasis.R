###
# Epistasis estimation
###


library(data.table)
library(stringr)

local_path = ""
gen_and_phen_path = ".tsv" # table with genotypes and phenotypes


gen_and_phen = fread(gen_and_phen_path)


rs_variables = colnames(gen_and_phen)[str_detect(colnames(gen_and_phen), "rs")]
rs_ids = substr(x = rs_variables, start = 1, nchar(rs_variables)-4)
colnames(gen_and_phen)[str_detect(colnames(gen_and_phen), "rs")] = rs_ids


mfi_table = fread(paste(local_path, '.tsv', sep = "")) # summary statistics from GWAS of height
mfi_table = mfi_table[rs_id %in% rs_ids]



get_score = function(rs_ids){
  return( as.matrix(gen_and_phen[, ..rs_ids]) %*% mfi_table[match(rs_ids, rs_id)]$Beta )
}

set.seed(1)


snp_subset1 = sample(rs_ids, 152) 
snp_subset2 = rs_ids[! rs_ids %in% snp_subset1]

gen_and_phen$score = get_score(rs_ids)
gen_and_phen$score1 = get_score(snp_subset1)
gen_and_phen$score2 = get_score(snp_subset2)
gen_and_phen[, score := (score - mean(score)) / sd(score), by = analysis_group]
gen_and_phen[, score1 := (score1 - mean(score1)) / sd(score1), by = analysis_group]
gen_and_phen[, score2 := (score2 - mean(score2)) / sd(score2), by = analysis_group]


summary(lm(Height ~ Sex + score + score*Sex, data = gen_and_phen))$coeff
summary(lm(Height ~ Sex + score1*Sex + score2*Sex + score1*score2, data = gen_and_phen))$coeff#["score1:score2", "Pr(>|t|)"]
summary(lm(log_height ~ Sex + score1*score2, data = gen_and_phen))$coeff


get_p_vals = function(){
  snp_subset1 = sample(rs_ids, 152) 
  snp_subset2 = rs_ids[! rs_ids %in% snp_subset1]
  
  gen_and_phen$score1 = get_score(snp_subset1)
  gen_and_phen$score2 = get_score(snp_subset2)
  gen_and_phen[, score1 := (score1 - mean(score1)) / sd(score1), by = analysis_group]
  gen_and_phen[, score2 := (score2 - mean(score2)) / sd(score2), by = analysis_group]
  p_h = summary(lm(Height ~ Sex + score1*score2, data = gen_and_phen))$coeff["score1:score2", "Pr(>|t|)"]
  p_log = summary(lm(log_height ~ Sex + score1*score2, data = gen_and_phen))$coeff["score1:score2", "Pr(>|t|)"]
  return(c("p_h" = p_h, "p_log" = p_log))
}

p_vals = replicate(n = 100, get_p_vals())

p_dt = as.data.table(t(p_vals))
p_dt = melt(p_dt)
library(ggplot2)


ggplot(p_dt, aes(x = variable, y = log10(value))) + geom_boxplot() + 
  geom_hline(yintercept = log10(0.05)) + theme_bw()

pdf("../Figures/epistasis.R")
ggplot(p_dt, aes(x = variable, y = log10(value))) + geom_boxplot() + 
  geom_hline(yintercept = log10(0.05)) + theme_bw()
dev.off()



