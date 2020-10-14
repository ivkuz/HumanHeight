###
# comparing of models in the ability of determining high and short people
###


library('data.table')
library('R.filesets')
library('pROC')


# ===========Funcions=============
get_r2 <- function(real_value, fit_value){
  r_sq_test_set = 1 - var(real_value - fit_value)/var(real_value)
  return(r_sq_test_set)
}

get_model_params <- function(model){
  stat = summary(model)
  return(stat$coefficients[,c('Estimate', 'Std. Error')])
}

get_model_estimation <- function(data_test, data_train,
                                 height_threshold=c(165,185)){
  mod_h = lm("Height ~ score_meta", data=data_train)
  mod_logH = lm("log_height ~ score_meta", data=data_train)
  h_score = predict.lm(mod_h, data_test)
  log_h_score = predict.lm(mod_logH, data_test)
  h_short=as.numeric(data_test$Height < height_threshold[1])
  h_tall=as.numeric(data_test$Height > height_threshold[2])
  params_mod_h = get_model_params(mod_h)
  params_mod_logH = get_model_params(mod_logH)
  delta_tall_h = sqrt(sum((data_test[h_tall == 1]$Height - 
                  h_score[h_tall == 1])^2)/(sum(h_tall) - 1))
  delta_tall_log_h = sqrt(sum((data_test[h_tall == 1]$Height - 
                      10^log_h_score[h_tall == 1])^2)/(sum(h_tall) -1))
  delta_short_h = sqrt(sum((data_test[h_short == 1]$Height - 
                   h_score[h_short == 1])^2)/(sum(h_short) -1))
  delta_short_log_h = sqrt(sum((data_test[h_short == 1]$Height - 
                         10^log_h_score[h_short == 1])^2)/(sum(h_short) -1))
  
  delta_h = sqrt(sum((data_test$Height - 
                             h_score)^2)/(length(h_tall)))
  delta_log_h = sqrt(sum((data_test$Height - 
                                 10^log_h_score)^2)/(length(h_tall)))
  return(data.table("AUC" = c(auc(roc(h_tall, 10**log_h_score)),
                             auc(roc(h_tall, h_score)),
                             auc(roc(h_short, 10**log_h_score)),
                             auc(roc(h_short, h_score))),
                    "Model" = c("Tall Height",
                                "Tall logHeight",
                                "Short Height",
                                "Short logHeight"),
                    "R2" = c(get_r2(data_test$Height, 10 ** log_h_score),
                             get_r2(data_test$Height, h_score),
                             get_r2(data_test$Height, 10 ** log_h_score),
                             get_r2(data_test$Height, h_score)),
                    "R2 tail" = c(get_r2(data_test[h_tall == 1]$Height,
                                         10^log_h_score[h_tall == 1]),
                                  get_r2(data_test[h_tall == 1]$Height,
                                         h_score[h_tall == 1]),
                                  get_r2(data_test[h_short == 1]$Height,
                                         10^log_h_score[h_short == 1]),
                                  get_r2(data_test[h_short == 1]$Height,
                                         h_score[h_short == 1])),
                    "correlation_tail" = c(cor(data_test[h_tall == 1]$log_height,
                                         log_h_score[h_tall == 1]) ^ 2,
                                    cor(data_test[h_tall == 1]$Height,
                                         h_score[h_tall == 1]) ^ 2,
                                    cor(data_test[h_short == 1]$log_height,
                                         log_h_score[h_short == 1])^2,
                                    cor(data_test[h_short == 1]$Height,
                                         h_score[h_short == 1]) ^ 2),
                    "Deviation tail" = c(delta_tall_h, delta_tall_log_h,
                                         delta_short_h ,delta_short_log_h),
                    "Deviation" = c(delta_h, delta_log_h,
                                         delta_h ,delta_log_h),
                    'Model paramters: Intersept' = c(params_mod_h[1,c('Estimate')],
                                          params_mod_logH[1,c('Estimate')],
                                          params_mod_h[1,c('Estimate')],
                                          params_mod_logH[1,c('Estimate')]),
                    'Model paramters: Intersept SD' = c(params_mod_h[1,c('Std. Error')],
                                                     params_mod_logH[1,c('Std. Error')],
                                                     params_mod_h[1,c('Std. Error')],
                                                     params_mod_logH[1,c('Std. Error')]),
                    'Model paramters: Slope' = c(params_mod_h[2,c('Estimate')],
                                                     params_mod_logH[2,c('Estimate')],
                                                     params_mod_h[2,c('Estimate')],
                                                     params_mod_logH[2,c('Estimate')]),
                    'Model paramters: Slope SD' = c(params_mod_h[2,c('Std. Error')],
                                                     params_mod_logH[2,c('Std. Error')],
                                                     params_mod_h[2,c('Std. Error')],
                                                     params_mod_logH[2,c('Std. Error')])
                    ))
}

# ===========Analysis=============
data = loadRDS('.rds') # table with phenotypes and PGHS
data=data.table(data)
data_mean = fread('.csv', # table with mean values of height in assessment centres
                  sep=',')

data_subset = data[Centre %in% c('11005', '11010')]
data_mean_subset = data_mean[Centre %in% c('11005', '11010')]
data_mean_subset$unique_analysis_group =  ifelse(data_mean_subset$n_people > 1000,
                                                 data_mean_subset$analysis_group,
                                                 "Other groups")
data_mean_subset$unique_analysis_group = paste(data_mean_subset$meaning,
                                               data_mean_subset$unique_analysis_group,
                                               data_mean_subset$Sex,
                                               sep="_")

data_mean_subset$Sex = factor(data_mean_subset$Sex,
                       level=c('Male','Female'),
                       label=c(1,0))
results = data.table()
for (sex in c(0,1)){
  for (gr1 in unique(data_mean_subset[Sex == sex]$unique_analysis_group)){
    for (gr2 in unique(data_mean_subset[Sex == sex]$unique_analysis_group)){
      print(gr1)
        train_desc = data_mean_subset[unique_analysis_group == gr1]
        data_train = data[(Centre %in% train_desc$Centre) &
                            (Sex %in% train_desc$Sex) &
                            (analysis_group %in% train_desc$analysis_group),]
        
        test_desc = data_mean_subset[unique_analysis_group == gr2]
        data_test = data[(Centre %in% test_desc$Centre) &
                            (Sex %in% test_desc$Sex) &
                            (analysis_group %in% test_desc$analysis_group),]
        if (sex == 0){
          tmp = get_model_estimation(data_test, data_train,
                               height_threshold=c(152,172))
        }else{
          tmp = get_model_estimation(data_test, data_train,
                               height_threshold=c(165,185))
        }
        tmp$analysis_group_train = gr1
        tmp$analysis_group_test= gr2
        print(tmp)
        tmp$n_people_train = sum(train_desc$n_people)
        tmp$mean_h_train = 
          sum(train_desc$n_people*train_desc$h_mean)/sum(train_desc$n_people)
        tmp$sd_h_train = 
          sum(train_desc$n_people*train_desc$h_sd)/sum(train_desc$n_people)
        tmp$SEM_train = tmp$sd_h_train/sqrt(tmp$n_people_train)
        tmp$n_people_test = sum(test_desc$n_people)
        tmp$mean_h_test = 
          sum(test_desc$n_people*test_desc$h_mean)/sum(test_desc$n_people)
        tmp$sd_h_test = 
          sum(test_desc$n_people*test_desc$h_sd)/sum(test_desc$n_people)
        tmp$SEM_test = tmp$sd_h_test/sqrt(tmp$n_people_test)
        results = rbind(results, tmp)
    }
  }
}

# delta height-predict height
write.table(results[,c("analysis_group_train" ,"analysis_group_test" , "AUC" ,
                    "Model" ,"R2" ,"R2 tail" ,"correlation_tail" ,
                    'Deviation tail',
                    "Model paramters: Intersept" ,
                    "Model paramters: Intersept SD", 
                    "Model paramters: Slope" ,"Model paramters: Slope SD" ,
                    "n_people_train" ,"mean_h_train" ,"sd_h_train" ,
                    "SEM_train" ,"n_people_test" ,"mean_h_test" ,
                    "sd_h_test" ,"SEM_test")],
            'Tables/auc_Leeds_vc_Edinburgh.csv',
            row.names = F, quote = F, sep=',')
tmp = results[analysis_group_train == analysis_group_test,
              c("analysis_group_train", "Model paramters: Intersept" ,
                "Model paramters: Intersept SD", 
                "Model paramters: Slope" ,"Model paramters: Slope SD", "R2",
                'Deviation'
                )]
tmp$model = rep(c('H', 'logH'), 20)
merge(tmp, data_mean_subset[,c('unique_analysis_group', 'Sex', 'meaning')],
      by.x = "analysis_group_train",
      by.y = "unique_analysis_group")
write.table(tmp,
           '.csv',
           row.names = F, quote = F, sep=','
)


