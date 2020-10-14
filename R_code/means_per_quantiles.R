###
# Plot supplementary figure 2
library('data.table')
library('R.filesets')
library("ggplot2")
library("reshape2")
library("forcats")
library("BSDA")
library(ggpubr)
library(RColorBrewer)

QC <- function(data){
  for(group in unique(data$analysis_group)){
    for(sex in unique(data$Sex)){
      mean_h <- mean(data[Sex == sex & analysis_group == group, Height])
      sd_h <- sd(data[Sex == sex & analysis_group == group, Height])
      mean_lh <- mean(data[Sex == sex & analysis_group == group, log_height])
      sd_lh <- sd(data[Sex == sex & analysis_group == group, log_height])
      data_filt <- data[Sex != sex | analysis_group != group |
                          (abs(Height - mean_h) <= 4.75*sd_h &
                             abs(log_height - mean_lh) <= 4.75*sd_lh),]
      print(paste(group, sex, nrow(data)-nrow(data_filt)))
      data <- data_filt
    }
  }
  return(data)
}

getPredictions <- function(data, RP = T){
  if(RP){
    mod_h <- lm("Height ~ norm_score_meta + predF3", data=data)
    mod_logH <- lm("log_height ~ norm_score_meta + predF3", data=data)
  } else{
    mod_h <- lm("Height ~ norm_score_meta", data=data)
    mod_logH <- lm("log_height ~ norm_score_meta", data=data)
  }
  h_score <- predict.lm(mod_h, data)
  log_h_score <- predict.lm(mod_logH, data)
  data <- data[, c("h_score", "log_h_score") := list(h_score, log_h_score)]
  return(data[order(norm_score_meta),])
}

determGroup <- function(score){
  if(score < qnorm(0.01)){
    return("[0,1]")
  } else if(score < qnorm(0.05)){
    return("[1,5]")
  } else if(score < qnorm(0.20)){
    return("[5,20]")
  } else if(score < qnorm(0.40)){
    return("[20,40]")
  } else if(score < qnorm(0.60)){
    return("[40,60]")
  } else if(score < qnorm(0.80)){
    return("[60,80]")
  } else if(score < qnorm(0.95)){
    return("[80,95]")
  } else if(score < qnorm(0.99)){
    return("[95,99]")
  } else {
    return("[99,100]")
  }
}

groupParams <- function(data){
  tab <- NULL
  for(gr in unique(data$group)){
    for(h in c("Height", "h_score", "log_height", "log_h_score")){
      tmp <- list()
      tmp$Sex <- data[, as.character(Sex[1])]
      tmp$group <- gr
      tmp$variable <- h
      tmp$N <- nrow(data[group == gr,])
      tmp$Mean <- mean(data[group == gr, get(h)])
      tmp$Median <- round(median(data[group == gr, get(h)]), 10)
      tmp$SD <- sd(data[group == gr, get(h)])
      tmp$SEM <- round(tmp$SD/sqrt(tmp$N), 10)
      tmp$Mean <- round(tmp$Mean, 10)
      tmp$SD <- round(tmp$SD, 10)
      if(h == "h_score" | h == "log_h_score"){
        
        if(h == "h_score"){sc_h <- "Height"}else{sc_h <- "log_height"}
        
        tmp$w.test <- formatC(wilcox.test(data[group == gr, get(sc_h)], 
                                          data[group == gr, get(h)], 
                                          paired = TRUE, 
                                          alternative = "two.sided")$p.value, 
                              format = "e", digits = 2)
        tmp$t.test <- formatC(t.test(data[group == gr, get(sc_h)], 
                                     data[group == gr, get(h)], 
                                     paired = TRUE, 
                                     alternative = "two.sided")$p.value, 
                              format = "e", digits = 2)
        tmp$z.test <- formatC(z.test(data[group == gr, get(sc_h)], 
                                     data[group == gr, get(h)], 
                                     alternative = "two.sided",
                                     sigma.x = sd(data[group == gr, get(sc_h)]),
                                     sigma.y = sd(data[group == gr, get(h)])
        )$p.value, 
        format = "e", digits = 2)
        tmp$ks.test <- formatC(ks.test(data[group == gr, get(sc_h)], 
                                       data[group == gr, get(h)], 
                                       alternative = "two.sided")$p.value, 
                               format = "e", digits = 2)
      } else{
        tmp$w.test <- "-"
        tmp$t.test <- "-"
        tmp$z.test <- "-"
        tmp$ks.test <- "-"
        
      }
      tab <- rbind(tab, unlist(tmp))
    }
  }
  tab <- data.table(tab)
  tab$variable <- factor(tab$variable, 
                         levels = c("Height", "h_score", 
                                    "log_height", "log_h_score"))
  return(tab)
}

plotGroups <- function(data, t, analysis_gr, fig_list){
  
  pd <- position_dodge(0.4)
  colors = brewer.pal(3, "Dark2")
  
  t_h <- t[variable == "Height" | variable == "h_score",]
  t_lh <- t[variable == "log_height" | variable == "log_h_score",]
  
  for(t in list(t_h, t_lh)){
    clr <- "red2"
    p <- ggplot(t, aes(x = group, y = as.numeric(Mean),
                       group = variable)) + theme_bw() + 
      scale_x_discrete(limits=c("[0,1]", "[1,5]", "[5,20]", "[20,40]",
                                "[40,60]", "[60,80]", "[80,95]",
                                "[95,99]", "[99,100]")) +
      geom_errorbar(aes(ymin=as.numeric(Mean)-2.58*as.numeric(SEM), 
                        ymax=as.numeric(Mean)+2.58*as.numeric(SEM),
                        color = variable), 
                    width=.1, position=pd) +
      geom_point(position=pd, aes(color = variable), shape = 15)+ 
      scale_color_manual(values=c('royalblue', clr),#c(colors[1], clr),#c("black", clr),
                         labels=c("Observed", "Predicted"),
                         name = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.text = element_text(size = 12))
    if(unique(t$Sex) == "Female"){
      if(identical(t, t_h)){
        p <- p + labs(y = "Mean height (cm)", title = "Female") +
          theme(axis.text.x=element_blank(), axis.title.x=element_blank())
      }else{
        p <- p + labs(y = bquote("Mean " ~ log[10] ~ "(height) (cm)"),
                      x = "Quantile range, %")
      }
    }else{
      p <- p + theme(axis.title.y=element_blank())
      if(identical(t, t_h)){p <- p + labs(title = "Male") +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank())
      } else{
        p <- p + labs(x = "Quantile range, %")
      }
    }
    fig_list <- append(fig_list, list(p))
  }
  
  return(fig_list)
}

PGHS <- 1
analysis_gr <- "English"

data = loadRDS('phen_meta_score.rds') # Table with phenotypes and PGHS2
data <- data[,.(Sex, analysis_group, Height, log_height, predF3,
                norm_score, score_meta)]

if(PGHS == 1){
  RP <- fread('phenotypes_with_pred3.tsv')  # Table with phenotypes and PGHS
  data[, predF3 := RP$predF3]
  data[, norm_score_meta := norm_score]
} else{
  data[, norm_score_meta := (score_meta - mean(score_meta)) /
         sd(score_meta), by = analysis_group]
}
data <- data[analysis_group == analysis_gr,]


data$Sex = factor(data$Sex,
                  level=c(1,0),
                  label=c('Male','Female'))

data_m <- data[Sex == "Male",]
data_w <- data[Sex == "Female",]


data_m <- getPredictions(data_m, RP = T)
data_w <- getPredictions(data_w, RP = T)

data_m$group <- sapply(data_m$norm_score_meta, determGroup) 
data_w$group <- sapply(data_w$norm_score_meta, determGroup)

t_m <- groupParams(data_m)
t_w <- groupParams(data_w)
tab <- rbind(t_m, t_w)
write.table(tab, ".tsv", row.names = F, quote = F, sep = "\t")

fig_list <- list()
fig_list <- plotGroups(data_m, t_m, analysis_gr, fig_list)
fig_list <- plotGroups(data_w, t_w, analysis_gr, fig_list)
fig <- ggarrange(fig_list[[3]], fig_list[[1]], fig_list[[4]], fig_list[[2]], 
                 ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
pdf("groups_boxplots_normgr_newpercent_Eng_PGHS1_3_final.pdf", 
    width = 7, height = 7)
annotate_figure(fig, top = analysis_gr) #, bottom = "Quantile range, %"
dev.off()
