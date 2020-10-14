###
# Basic analysis of height factors
###

#simple analysis 
library(data.table)
gen_and_phen= fread(".tsv") # able with phenotypes, PGHS, RP

s1 = summary(lm(formula = "Height ~ Sex + norm_score + predF3", data = dt))
s1$r.squared

s2 = summary(lm(formula = "Height ~ Sex + norm_score", data = dt))
s2$r.squared

s3 = summary(lm(formula = "Height ~ Sex", data = dt))
s3$r.squared


s = summary(lm(formula = "Height ~ Sex + norm_score + predF3 + 
                 norm_score*Sex + predF3*Sex + norm_score*predF3", data = dt))


modelAdditive <- lm(formula = "Height ~ Sex + norm_score + predF3", data=gen_and_phen)
modelMultiplicative <- lm(formula = "log_height ~ Sex + norm_score + predF3", data=gen_and_phen)

AIC(modelAdditive)
AIC(modelMultiplicative) + 2*sum(gen_and_phen$log_height)


groups = gen_and_phen[,.(Sex =factor(Sex == 1,
                                     level=c(TRUE, FALSE), label=c('Male','Female')),
                         PGHS = factor(norm_score >= median(norm_score), 
                                       level=c(TRUE, FALSE), label=c('high', 'low')),
                         RP = factor(predF3 >= median(predF3),
                                     level=c(TRUE, FALSE), label=c('high', 'low')),
                         ID=ID, Height, log_height),
                      by='analysis_group'] 

gen_and_phen[, PGHS_binary := as.numeric(norm_score >= median(norm_score)), by='analysis_group'] 
gen_and_phen[, RP_binary := as.numeric(predF3 >= median(predF3)), by='analysis_group'] 

gen_and_phen[, Sex_cent := Sex - mean(Sex)]
gen_and_phen[, PGHS_bin := PGHS_binary - mean(PGHS_binary)]
gen_and_phen[, RP_bin := RP_binary - mean(RP_binary)]
gen_and_phen[, SexPGHS_bin := Sex*PGHS_binary - mean(Sex*PGHS_binary)]
gen_and_phen[, SexRP_bin := Sex*RP_binary - mean(Sex*RP_binary)]
gen_and_phen[, PGHS_binRP_bin := PGHS_binary*RP_binary - mean(PGHS_binary*RP_binary)]


get_simpple_table = function(simple_model){
  
  get_dt_with_res_of_model = function(ethn, height_column){
    if (ethn=="ALL"){
      d.t <- gen_and_phen
    }
    else {d.t = gen_and_phen[analysis_group==ethn]}
    
    my_model =  paste(height_column, simple_model, sep="~")
    model_summary = summary(lm(formula = as.formula(my_model), 
                               data = d.t))
    
    res = data.table(model_summary$coefficients, keep.rownames = T)
    res$model_rsq = c(rep(NA, (nrow(res)-1)), model_summary$r.squared)
    res$analysis_group = ethn
    res$height_column = eval(height_column)
    res$mean = mean(d.t[,get(height_column)])
    
    res[rn=="norm_score"]$rn = "PGHS"
    res[rn=="predF3"]$rn = "RP"
    res[rn=="Sex:norm_score"]$rn = "Sex:PGHS"
    res[rn=="(Intercept)"]$rn = "Intercept"
    res[rn=="norm_score:predF3"]$rn = "PGHS:RP"
    res[rn=="Sex:predF3"]$rn = "Sex:RP"
    return(res)
    
  }
  
  simple_table_a = rbindlist(lapply(c(names(sort(table(gen_and_phen$analysis_group), decreasing = T)), "ALL"), 
                                     get_dt_with_res_of_model, height_column = "Height"))
  simple_table_b = rbindlist(lapply(c(names(sort(table(gen_and_phen$analysis_group), decreasing = T)), "ALL"), 
                                    get_dt_with_res_of_model, height_column = "log_height"))
  simple_table = rbind(simple_table_a, simple_table_b)
  return(simple_table)
}


simple_table_1 = get_simpple_table("Sex + norm_score + predF3") 

write.table(simple_table_1, file = "../Tables/simple_table_1.tsv", 
            quote = F, row.names = F, sep = "\t", dec = ",", na = "")



simple_table_3 = get_simpple_table("Sex + norm_score + predF3 + norm_score*Sex")
write.table(simple_table_3, file = "../Tables/simple_table_3.tsv", 
            quote = F, row.names = F, sep = "\t", dec = ",", na = "")


simple_table_4 = get_simpple_table("Sex + norm_score + predF3 + 
                                   norm_score*Sex + predF3*Sex + norm_score*predF3")
write.table(simple_table_4, file = "../Tables/simple_table_4.tsv", 
            quote = F, row.names = F, sep = "\t" , dec = ",", na = "")


simple_table_6= get_simpple_table("Sex + PGHS_binary + RP_binary + 
                                   PGHS_binary*Sex + RP_binary*Sex + PGHS_binary*RP_binary")
write.table(simple_table_6, file = "../Tables/simple_table_6.tsv", 
            quote = F, row.names = F, sep = "\t" , dec = ",", na = "")

simple_table_7 = get_simpple_table("Sex_cent + PGHS_bin + RP_bin + 
                                   SexPGHS_bin + SexRP_bin + PGHS_binRP_bin")
write.table(simple_table_7, file = "../Tables/simple_table_7.tsv", 
            quote = F, row.names = F, sep = "\t" , dec = ",", na = "")


sort(
sapply(unique(gen_and_phen$analysis_group), function(ethn){
  temp = gen_and_phen[analysis_group==ethn][Sex==0]
  nrow(temp[Age < 55])/nrow(temp)
})) 

library(car)

get_stat = function(d.t, height_column="Height"){
  ans = d.t[,list(mean = mean(get(height_column)), sd = sd(get(height_column)), 
                           N = .N, `variation coeff` = sd(get(height_column)) / mean(get(height_column)) * 100 )]
  return(ans)
}

simple_table_2 = rbindlist(lapply(c(names(sort(table(gen_and_phen$analysis_group), decreasing = T)), "ALL_but_Eng", "ALL"),
                                  function(ethn) {
  get_stats = function(d.t1, d.t2, name1, name2, height_column){
    res  = rbind(get_stat(d.t1, height_column), get_stat(d.t2, height_column))
    res$group = c(name1, name2)
    d.t = rbind(d.t1, d.t2)
    res$levine_pval = 
      c(leveneTest(d.t[,get(height_column)], group = as.factor(d.t$ID %in% d.t1$ID))$P[1], NA) 
    res$ethn = ethn
    res$feature = eval(height_column)
    return(res)
  }
  
  if (ethn=="ALL"){
    temp_dt <- gen_and_phen
  } else { 
    if (ethn=="ALL_but_Eng"){
    temp_dt <- gen_and_phen[analysis_group!="English"]
    } 
    else {temp_dt = gen_and_phen[analysis_group==ethn]}
  }
  
  # by sex
  women = temp_dt[Sex==0]
  men  = temp_dt[Sex==1]
  # by predF3 
  low_predF3 = temp_dt[(predF3 < median(predF3))] 
  high_predF3 = temp_dt[(predF3 >= median(predF3))]
  # by score
  low_score = temp_dt[(norm_score < median(norm_score))] 
  high_score = temp_dt[(norm_score >= median(norm_score))] 
  
  res = rbind(
    get_stats(men, women, "men","women", "Height"),
    get_stats(men, women, "men","women", "log_height"),
    get_stats(high_score, low_score, "high_PGHS","low_PGHS", "Height"),
    get_stats(high_score, low_score, "high_PGHS","low_PGHS", "log_height"),
    get_stats(high_predF3, low_predF3, "high_RP","low_RP", "Height"),
    get_stats(high_predF3, low_predF3, "high_RP","low_RP", "log_height")
  )
  
  }))

simple_table_2 = simple_table_2[order(feature)]
write.table(simple_table_2, file = "../Tables/simple_table_2.tsv", 
            quote = F, row.names = F, sep = "\t", dec = ",", na = "")



get_table_zero = function(d.t){
  simple_table_a = rbindlist(lapply(unique(d.t$analysis_group), function(ethn){
    res = get_stat(d.t[analysis_group==ethn], height_column = "Height")
    res$ethn = ethn
    res$feature = "Height"
    return(res)
  }))
  simple_table_b = rbindlist(lapply(unique(d.t$analysis_group), function(ethn){
    res = get_stat(d.t[analysis_group==ethn], height_column = "log_height")
    res$ethn = ethn
    res$feature = "log_height"
    return(res)
  }))
  simple_table = rbind(simple_table_a, simple_table_b)
  return(simple_table)
}

simple_table_0  = get_table_zero(gen_and_phen)
write.table(simple_table_0, file = "../Tables/simple_table_0.tsv", 
            quote = F, row.names = F, sep = "\t", dec = ",")

sapply(unique(simple_table_3$analysis_group), function(x){
  t = simple_table_3[analysis_group==x][height_column=="Height"][rn=="score"]$`t value`
  t_log = simple_table_3[analysis_group==x][height_column=="log_height"][rn=="score"]$`t value`
  ((t**2 - t_log**2)/t**2)*100})


with(simple_table_1[height_column=="Height"][rn=="norm_score"], plot(mean, Estimate))
