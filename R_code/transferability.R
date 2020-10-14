###
# Predicting height in ethnic groups with English training set
###

library(data.table)
library('R.filesets')

raw_pheno <- loadRDS(".rds") # raw phenotypes
scores <- fread("") # polygenic scores
pheno_scores <- merge(raw_pheno, scores, by.x="ID", by.y="rn")
pheno_scores$norm_score <- (pheno_scores$score-mean(pheno_scores$score))/sd(pheno_scores$score)

nameGroup <- function(eb, pb){
  # eb - Ethnic background
  # pb - Place of birth
  if(is.na(eb)){return("NA")
  } else if(eb == 1001){
    if(is.na(eb)){return("NA")
    } else if(pb == "1"){
      return("English")
    } else if(pb == 2){return("Welsh")
    } else if(pb == 3){return("Scottish")
    } else if(pb %in% 4:6){return("British")
    } 
  } else if(eb == 1002){return("Irish")
  } else if(eb == 1003){return("White")
  } else if(eb %in% 2001:2004){return("Mixed")
  } else if(eb == 3001){return("Indian")
  } else if(eb == 3002){return("Pakistani")
  } else if(eb %in% 3003:3004){return("Asian")
  } else if(eb == 4001){return("Carribean")
  } else if(eb == 4002){return("African")
  } else if(eb == 5){return("Chinese")
  } else if(eb == 6){return("Other")
  } else {return("NA")}
}

pheno_scores[, group := nameGroup(eb=`Ethnic background`, pb=Country_of_birth_UK),
             by=1:nrow(pheno_scores)]
pheno_scores <- pheno_scores[!is.na(Height) & !is.na(Sex) & !is.na(Income) & 
                               !is.na(`Age when attended assessment centre`) & 
                               group != "NA",]


height_lm1 <- lm(Height ~ Sex + norm_score + Income + 
                 `Age when attended assessment centre`,
               data = pheno_scores[group=="English",])
logheight_lm1 <- lm(log10(Height) ~ Sex + norm_score + Income + 
                  `Age when attended assessment centre`,
                data = pheno_scores[group=="English",])
height_lm2 <- lm(Height ~ Sex + norm_score + Income + 
                   `Age when attended assessment centre` +
                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                 data = pheno_scores[group=="English",])
logheight_lm2 <- lm(log10(Height) ~ Sex + norm_score + Income + 
                      `Age when attended assessment centre` +
                      PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                    data = pheno_scores[group=="English",])
pheno_scores[, c("predicted_height1", "predicted_logheight1",
                 "predicted_height2", "predicted_logheight2") := 
               .(predict(height_lm1, pheno_scores), 
                 predict(logheight_lm1, pheno_scores),
                 predict(height_lm2, pheno_scores), 
                 predict(logheight_lm2, pheno_scores))]

summary_table <- pheno_scores[, .(.N, mean = mean(Height), sd = sd(Height), 
                                  cv = sd(Height)/mean(Height), 
                                  r2_additive = cor(Height, predicted_height1)^2,
                                  r2_multip = cor(Height, 10^predicted_logheight1)^2,
                                  r2_additive_log = cor(log10(Height), log10(predicted_height1))^2,
                                  r2_multip_log = cor(log10(Height), predicted_logheight1)^2,
                                  bias_additive = mean(predicted_height1 - Height),
                                  bias_multip = mean(10^predicted_logheight1 - Height),
                                  bias_additive_log = mean(log10(predicted_height1) - log10(Height)),
                                  bias_multip_log = mean(predicted_logheight1 - log10(Height)),
                                  median_bias_additive = median(predicted_height1 - Height),
                                  median_bias_multip = median(10^predicted_logheight1 - Height),
                                  median_bias_additive_log = median(log10(predicted_height1) - log10(Height)),
                                  median_bias_multip_log = median(predicted_logheight1 - log10(Height)),
                                  r2_pc_additive = cor(Height, predicted_height2)^2,
                                  r2_pc_multip = cor(Height, 10^predicted_logheight2)^2,
                                  r2_pc_additive_log = cor(log10(Height), log10(predicted_height2))^2,
                                  r2_pc_multip_log = cor(log10(Height), predicted_logheight2)^2,
                                  bias_pc_additive = mean(predicted_height2 - Height),
                                  bias_pc_multip = mean(10^predicted_logheight2 - Height),
                                  bias_pc_additive_log = mean(log10(predicted_height2) - log10(Height)),
                                  bias_pc_multip_log = mean(predicted_logheight2 - log10(Height)),
                                  median_bias_pc_additive = median(predicted_height2 - Height),
                                  median_bias_pc_multip = median(10^predicted_logheight2 - Height),
                                  median_bias_pc_additive_log = median(log10(predicted_height2) - log10(Height)),
                                  median_bias_pc_multip_log = median(predicted_logheight2 - log10(Height))),
                                  by=group]
write.table(summary_table, "/home/ikuznetsov/height/transferability.tsv",
            row.names = F, quote = F, sep = "\t")
