###
# Quantile-Quantile plot for height in UK Biobank. Supplementary figure 1
###

library(data.table)
library(R.filesets)
library(gridExtra)
library(car)
library(ggplot2)
library(qqplotr)
library(reshape2)
library(plyr)

data = loadRDS('.rds') # table with height, Sex, PGHS, RP

get_dfTarget <- function(tmp){
  model <- lm(value ~ norm_score + predF3, data = tmp)
  tmp$target <- tmp[, value - (norm_score*model$coefficients[2] +
                                 predF3*model$coefficients[3])]
  return(tmp)
}

drawQQPlot <- function(tmp){
  if (nrow(tmp) > 30000){
    tmp <- tmp[sample(30000),]
  }
  model <- lm(value ~ norm_score + predF3, data = tmp)
  tmp$target <- tmp[, value - (norm_score*model$coefficients[2] +
                                      predF3*model$coefficients[3])]
  if ((unique(tmp$Sex)==0) & 
      (unique(tmp$variable) == 'Height')){
    y_lab=unique(tmp$analysis_group)
  }else{y_lab=""}

  gg <- ggplot(data = tmp, mapping = aes(sample = target)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point(colour="royalblue", size=0.2) +
    theme_minimal() +
    labs(x = "", y = y_lab)
  if (unique(tmp$analysis_group) == "English"){
    s <- ifelse(unique(tmp$Sex) == 0, 'Female', 'Male')
    title = unique(paste(s, tmp$variable))
    gg <- gg +  ggtitle(title)
  }
  
  return(gg)
}


labels_order = c("0_English_height" , "0_Scottish_height","0_Welsh_height",
                 "0_English_log(height)", "0_Scottish_log(height)",
                 "0_Welsh_log(height)",
                 "1_English_height" , "1_Scottish_height","1_Welsh_height",
                 "1_English_log(height)", "1_Scottish_log(height)", 
                 "1_Welsh_log(height)")

subset <- data[data$analysis_group %in% c("English", "Scottish", "Welsh")]
t <- melt(subset,
          id.vars=c('Sex', 'analysis_group', 'predF3', 'norm_score'),
          measure.vars = c('Height', 'log_height'))
t_var <- t$variable
t_var <- gsub("Height", "height", t_var)
t_var <- gsub("log_height", "log(height)", t_var)
t$variable <- t_var
t$label = paste(t$Sex, t$analysis_group, t$variable, sep = "_")

list_gg <- lapply(labels_order, function(l) {
  tmp <- t[t$label == l];
  p <- drawQQPlot(tmp);
  return(p)})

pdf("../Figures/multipage1.pdf", width = 10, height = 6)
ml <- marrangeGrob(list_gg, nrow=3, ncol=4,
                   bottom = "Sample Quantiles",
                   left = "Theoretical Quantiles",
                   top='Q-Q plot')
ml
dev.off()
