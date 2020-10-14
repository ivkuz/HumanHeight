###
# plot Figure 2
###
library(data.table)
library(RColorBrewer)
phen = fread('.tsv', sep='\t') # table with phenotypes and RP
###split stage
groups = phen[,.(f1_sex=factor(Sex == 1,
                               level=c(TRUE, FALSE), label=c('Male','Female')),
                 f2_prs=factor(norm_score >= median(norm_score), 
                               level=c(TRUE, FALSE), label=c('high', 'low')),
                 f3=factor(predF3 >= median(predF3),
                           level=c(TRUE, FALSE), label=c('high', 'low')),
                 ID=ID),
              by='analysis_group'] 
phen[,.(Med_score = median(norm_score), Med_predF3 = median(predF3)),
     by='analysis_group'] 
groups <- merge(phen[,.(ID, norm_score, predF3, 
                        Sex, Height, log_height),], 
                groups,
                by= c('ID'))
get_model_params<-function(model){
  s_m = summary(model)
  coef = s_m$coefficients
  l <- c(coef[2,'Estimate'], coef[2,'Std. Error'], s_m$adj.r.squared)
  names(l) <- c('beta', 'se', 'r2adj')
  return(l)
}

data_sex = groups[, .(h = mean(Height), 
                      sd_h = sd(Height),
                      log_h = mean(log_height),
                      sd_log_h = sd(log_height),
                      params_h = list(get_model_params(lm(Height ~ Sex))),
                      params_log_h = list(get_model_params(lm(log_height ~ Sex))),
                      N=length(Height)),
                  by=c("analysis_group","f2_prs",
                       "f3")]
data_sex$beta_h = sapply(data_sex$params_h, 
                         function (x) { return(x['beta'])})
data_sex$beta_log_h = sapply(data_sex$params_log_h, 
                             function (x) { return(x['beta'])})

data_prs = groups[, .(h = mean(Height), 
                      sd_h = sd(Height),
                      log_h = mean(log_height),
                      sd_log_h = sd(log_height),
                      params_h = list(get_model_params(lm(Height ~ norm_score))),
                      params_log_h = list(get_model_params(lm(log_height ~ norm_score))),
                      N=length(Height)),
                  by=c("analysis_group","f1_sex",
                       "f3")]
data_prs$beta_h = sapply(data_prs$params_h, 
                         function (x) { return(x['beta'])})
data_prs$beta_log_h = sapply(data_prs$params_log_h, 
                             function (x) { return(x['beta'])})

data_f3 = groups[, .(h = mean(Height), 
                     sd_h = sd(Height),
                     log_h = mean(log_height),
                     sd_log_h = sd(log_height),
                     params_h = list(get_model_params(lm(Height ~ predF3))),
                     params_log_h = list(get_model_params(lm(log_height ~ predF3))),
                     N=length(Height)),
                 by=c("analysis_group","f1_sex",
                      "f2_prs")]
data_f3$beta_h = sapply(data_f3$params_h, 
                        function (x) { return(x['beta'])})



data_f3$beta_log_h = sapply(data_f3$params_log_h, 
                            function (x) { return(x['beta'])})



eth_names = c("English", "Scottish", "Welsh", "Irish",
              "Other_British", "Other_white")

eth_prs_label = paste(eth_names,"high", sep='_')


paste_reg_results_for_h<- function(data, name=""){
  w = sqrt(data$N)
  model = lm(data=data, beta_h ~ h, weights = w^2)
  params = get_model_params(model)
  abline(lm(model),
         col='grey44', lty=2)
  
  if (params['r2adj']>0) {
    mtext(side=3,#x=min(data$log_h), y=max(data$beta_log_h),
          paste0("k=",signif(params['beta'],2),
                 ", SE=",signif(params['se'],2),
                 ", R²=",signif(params['r2adj'],2),
                 " (***)"),
          cex=0.7,adj=0,ps=12,
          at=par('usr')[1])
  } else {
    mtext(side=3,#x=min(data$log_h), y=max(data$beta_log_h),
          paste0("k=",signif(params['beta'],2),
                 ", SE=",signif(params['se'],2),
                 ", R²<0.001"),
          cex=0.7,adj=0,ps=12,
          at=par('usr')[1])
  }
  mtext(side=3,
        text = name,
        cex=1, adj=3, line=1.2,
        at=par('usr')[1])
}

paste_reg_results_for_log_h<- function(data, name=""){
  w = sqrt(data$N)
  model = lm(data=data, beta_log_h ~ log_h, weights = w^2)
  params = get_model_params(model)
  abline(lm(model),
         col='grey44', lty=2)
  if (params['r2adj']>0) {
    mtext(side=3,#x=min(data$log_h), y=max(data$beta_log_h),
          paste0("k=",signif(params['beta'],2),
                 ", SE=",signif(params['se'],2),
                 ", R²=",signif(params['r2adj'],2)),
          cex=0.7,adj=0,
          at=par('usr')[1])
  } else {
    mtext(side=3,#x=min(data$log_h), y=max(data$beta_log_h),
          paste0("k=",signif(params['beta'],2),
                 ", SE=",signif(params['se'],2),
                 ", R²<0.001"),
          cex=0.7,adj=0,
          at=par('usr')[1])
  }
  mtext(side=3,
        text = name,
        cex=1, adj=3, line=1.2,
        at=par('usr')[1])
  
  
}
###PLOT#####

pdf('.pdf', width = 8, height = 6)#, width = 1240, height = 1000, res = 150)
layout(matrix(c(1,2,7,3,4,7,5:7),nrow = 3, byrow = TRUE), 
       widths = c(0.38, 0.38, 0.25))
par(col='black',col.axis='black', bty='n', pch=16,
    cex.lab = 1.2, ps=14)
par(mar = c(5,5,3,2))
size=5
colors = brewer.pal(length(eth_names), "Dark2")
colors = sapply(colors, 
                function(x){paste0(x,'B4')})

#-----Sex-------

colors_names = factor(data_sex$analysis_group,
                      levels= eth_names,
                      labels = colors)

bg_names = factor(paste(data_sex$analysis_group,
                        data_sex$f3, sep="_"),
                  levels= eth_prs_label,
                  labels = colors)

pch_names =  factor(data_sex$f2_prs, 
                    levels= c("low", "high"), 
                    labels = c(21,24))

w = sqrt(data_sex$N)
plot(x=data_sex$h,
     y = data_sex$beta_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     pch=as.numeric(as.vector(pch_names)),
     bg=as.vector(bg_names),
     ylab='Effect size of sex', xlab='',#xlab='Mean Height',
     xlim=c(signif(min(data_f3$h),3), 
            round(max(data_f3$h),3)+ 1)
)

paste_reg_results_for_h(data_sex, "A")

plot(x=data_sex$log_h,
     y = data_sex$beta_log_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     pch=as.numeric(as.vector(pch_names)),
     ylab ='',
     xlab='',#  xlab='Mean log(Height)',
     xlim=c(signif(min(data_sex$log_h),4) - 0.0015, 
            round(max(data_sex$log_h),3) + 0.001),
     bg=as.vector(bg_names))

paste_reg_results_for_log_h(data_sex, "B")

#-----PRS-------
colors_names = factor(data_prs$analysis_group,
                      levels= eth_names,
                      labels = colors)

bg_names = factor(paste(data_prs$analysis_group,
                        data_prs$f3, sep="_"),
                  levels= eth_prs_label,
                  labels = colors)

w = sqrt(data_prs$N)
plot(x=data_prs$h,
     y = data_prs$beta_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     bg=as.vector(bg_names),
     pch=21,
     ylab='Effect size of PGHS', 
     xlab='',# xlab='Mean Height',
     xlim=c(signif(min(data_prs$h),3) - 1, 
            round(max(data_prs$h),3) + 2)
)


paste_reg_results_for_h(data_prs, "C")
text(x=mean(data_prs[data_prs$f1_sex=='Male']$h), y = 2.7,'Male')
text(x=mean(data_prs[data_prs$f1_sex!='Male']$h), y = 2.7,'Female')


plot(x=data_prs$log_h,
     y = data_prs$beta_log_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     bg=as.vector(bg_names),
     pch=21,
     ylab ='',
     xlab='',#xlab='Mean log(Height)',
     xlim=c(signif(min(data_prs$log_h),3) - 0.01 , 
            round(max(data_prs$log_h),3) + 0.003),
     ylim=c(signif(min(data_prs$beta_log_h),2), 
            signif(max(data_prs$beta_log_h),2)))

paste_reg_results_for_log_h(data_prs, "D")
text(x=mean(data_prs[data_prs$f1_sex=='Male']$log_h), y = 0.0067,'Male')
text(x=mean(data_prs[data_prs$f1_sex!='Male']$log_h), y = 0.0067,'Female')

#-----Year-------
colors_names =  factor(data_f3$analysis_group,
                       levels= eth_names,
                       labels = colors)
pch_names =  factor(data_f3$f2_prs, 
                    levels= c("low", "high"), 
                    labels = c(16,17))
w = sqrt(data_f3$N)
plot(x=data_f3$h,
     y = data_f3$beta_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     pch=as.numeric(as.vector(pch_names)),
     ylab='Effect size of RP', xlab='Mean height (cm)',
     xlim=c(signif(min(data_f3$h),3), 
            round(max(data_f3$h),3)+ 1))

paste_reg_results_for_h(data_f3, "E")
text(x=mean(data_f3[data_f3$f1_sex=='Male']$h), y = 2.2,'Male')
text(x=mean(data_f3[data_f3$f1_sex!='Male']$h), y = 2.2,'Female')

plot(x=data_f3$log_h,
     y = data_f3$beta_log_h,
     cex = size*w/max(w),
     col=as.vector(colors_names),
     pch=as.numeric(as.vector(pch_names)),
     ylab ='',
     xlab=expression('Mean log'[10]*'(height) (cm)'),
     ylim=c(signif(min(data_f3$beta_log_h),2)- 0.0003, 
            signif(max(data_f3$beta_log_h),2)),
     xlim=c(signif(min(data_f3$log_h),3) , 
            round(max(data_f3$log_h),3) + 0.003))
paste_reg_results_for_log_h(data_f3, "F")
text(x=mean(data_f3[data_f3$f1_sex=='Male']$log_h), y = 0.0054,'Male')
text(x=mean(data_f3[data_f3$f1_sex!='Male']$log_h), y = 0.0054,'Female')

#----Legend-----
par(mar = c(5,1,2,0))
plot(1, type = "n", axes=F, xlab="", ylab="")
legend(x=0.6, y = 1.2, # x = max(data$mean), y = max(data$sd),
       col = levels(colors_names), 
       legend = eth_names, bty='n',
       pch=16,
       cex = 1,
       # title='Ethnicity',
       title.adj = 1)
legend(x=0.6, y = 1,#"bottomleft", # x = max(data$mean), y = max(data$sd),
       col = 'grey35', 
       legend = c("low PGHS,  low RP",
                  "low PGHS,  high RP",
                  "high PGHS, low RP",
                  "high PGHS, high RP"),
       bty='n',
       pch=c(21,16,24,17),
       bg='grey35',
       title.adj = 0.2,
       cex = 1)

legend(x=0.6, y = 1,#"bottomleft", # x = max(data$mean), y = max(data$sd),
       col = 'grey35', 
       legend = c("",
                  "",
                  "",
                  ""),
       bty='n',
       pch=c(21,21,24,24),
       bg='grey35',
       title.adj = 0.2,
       cex = 1)

dev.off()
