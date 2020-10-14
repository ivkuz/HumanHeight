###
# plot Figure 3
###


library(data.table)
library(RColorBrewer)
phen = fread('.tsv',  sep='\t') # phenotypes with RP

###split stage

groups = phen[,.(f1_sex=factor(Sex == 1,
                               level=c(TRUE, FALSE), label=c('Male','Female')),
                 f2_prs=factor(norm_score >= median(norm_score), 
                               level=c(TRUE, FALSE), label=c('high', 'low')),
                 f3=factor(predF3 >= median(predF3),
                           level=c(TRUE, FALSE), label=c('high', 'low')),
                 ID=ID),
              by='analysis_group']

groups <- merge(phen[,.(ID, norm_score,Year_of_birth, 
                        Sex, Height, log_height),], 
                groups,
                by= c('ID'))

##==============================================================================
height_mean.dt <- groups[, .(h = mean(Height), 
                             sd_h = sd(Height),
                             log_h = mean(log_height),
                             sd_log_h = sd(log_height),
                             N=length(Height)),
                         by=c("analysis_group", "f2_prs",
                              "f3", "f1_sex")]

##Figure 3
eth_names = c("English", "Scottish", "Welsh", "Irish",
              "Other_British", "Other_white")

eth_prs_label = paste(eth_names,"high", sep='_')


colors = brewer.pal(length(eth_names), "Dark2")
colors = sapply(colors, 
                function(x){paste0(x,'B4')})
colors_names = factor(height_mean.dt$analysis_group,
                      levels= eth_names,
                      labels = colors)

bg_names = factor(paste(height_mean.dt$analysis_group,
                        height_mean.dt$f3, sep="_"),
                  levels= eth_prs_label,
                  labels = colors)

pch_names =  factor(height_mean.dt$f2_prs, 
                    levels= c("low", "high"), 
                    labels = c(21,24))
size=5

get_model_params<-function(model){
  s_m = summary(model)
  coef = s_m$coefficients
  l <- c(coef[2,'Estimate'], coef[2,'Std. Error'], s_m$adj.r.squared)
  names(l) <- c('beta', 'se', 'r2adj')
  return(l)
}

paste_reg_results_for_h<- function(data, name=""){
  w = sqrt(data$N)
  model = lm(data=data, sd_h ~ h, weights = w^2)
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
  model = lm(data=data, sd_log_h ~ log_h, weights = w^2)
  params = get_model_params(model)
  abline(lm(model),
         col='grey44', lty=2)
  if (params['r2adj']>0) {
    mtext(side=3,#x=min(data$log_h), y=max(data$beta_log_h),
          paste0("k=",signif(params['beta'],2),
                 ", SE=",signif(params['se'],2),
                 ", R²=",signif(params['r2adj'],2)),
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

pdf('.pdf', width = 8, height = 2.5)#, res = 150)

layout(matrix(c(1:3,4,4,3),nrow = 2, byrow = TRUE),
       widths = c(0.38, 0.38, 0.25),
       heights = c(0.9,0.1))
par(bty='n', ps = 14, cex.lab = 1.2)
par(mar = c(5,5,3,2))
w = sqrt(2*height_mean.dt$N)
plot(x=height_mean.dt$h, y=height_mean.dt$sd_h, 
     cex = size*w/max(w), 
     xlab="Mean height (cm)",
     ylab="SD height (cm)",
     pch= as.numeric(as.vector(pch_names)),
     xlim=c(signif(min(height_mean.dt$h),2) - 1, 
            round(max(height_mean.dt$h),-1)),
     ylim = c(round(min(height_mean.dt$sd_h),1) - 0.05,
              round(max(height_mean.dt$sd_h),1)+ 0.3),
     col = as.vector(colors_names),
     bg = as.vector(bg_names),
     axes = F
)
abline(lm(data=height_mean.dt, sd_h ~ h, weights = w^2), col='grey44', lty=2)

axis(1,
     at = seq(signif(min(height_mean.dt$h),2),
              round(max(height_mean.dt$h),-1), 
              by=5),
     col='gray44')
axis(2,
     at = seq((round(min(height_mean.dt$sd_h),1)),
              round(max(height_mean.dt$sd_h),1)+ 0.3, by=0.5),
     col='gray44')
text(x=mean(height_mean.dt[f1_sex=='Male']$h), y=6.9, 'Male')
text(x=mean(height_mean.dt[f1_sex=='Female']$h), y=6.9, 'Female')
mtext('A', side=3, at=155, line=1)
paste_reg_results_for_h(data=height_mean.dt)


plot(x=height_mean.dt$log_h, y=height_mean.dt$sd_log_h,
     cex = size*w/max(w), 
     bty='n', 
     xlab=expression('Mean log'[10]*'(height) (cm)'),
     ylab=expression('SD log'[10]*'(height) (cm)'),
     xlim=c(signif(min(height_mean.dt$log_h),3), 
            round(max(height_mean.dt$log_h),2)),
     ylim=c(signif(min(height_mean.dt$sd_log_h),2) - 0.0001, 
            signif(max(height_mean.dt$sd_log_h),2) + 0.0005),
     col = as.vector(colors_names),
     pch= as.numeric(as.vector(pch_names)),
     bg = as.vector(bg_names),
     axes=F)

abline(lm(data=height_mean.dt, sd_log_h ~ log_h, 
          weights = w^2), 
       col='grey44', lty=2)
axis(1,
     at = seq(signif(min(height_mean.dt$log_h),3), 
              round(max(height_mean.dt$log_h),2), 
              by=0.01),
     col='gray44')


axis(2,
     at = seq(signif(min(height_mean.dt$sd_log_h) - 0.0002,2), 
              signif(max(height_mean.dt$sd_log_h)+ 0.0002,3), 
              by=0.5*10^-3),
     col='gray44')
text(x=mean(height_mean.dt[f1_sex=='Male']$log_h),
     y=0.0164, 'Male')
text(x=mean(height_mean.dt[f1_sex=='Female']$log_h), 
     y=0.0164, 'Female')
mtext('B', side=3, at=2.19, line=1)
paste_reg_results_for_log_h(data=height_mean.dt)

par(mar = c(5,1,2,0))
plot(1, type = "n", axes=F, xlab="", ylab="")
legend(x=0.6, y = 1.4, # x = max(data$mean), y = max(data$sd),
       col = levels(colors_names), 
       legend = eth_names,
       pch=16,
       bty='n',
       # title='Ethnicity',
       title.adj = 0.1)

legend(x=0.6, y = 0.9,#"bottomleft", # x = max(data$mean), y = max(data$sd),
       col = 'grey35', 
       legend = c("low PGHS,  low RP",
                  "low PGHS,  high RP",
                  "high PGHS, low RP",
                  "high PGHS, high RP"),
       bty='n',
       pch=c(21,16,24,17),
       bg='grey35',
       title.adj = 0.2)

legend(x=0.6, y = 0.9,#"bottomleft", # x = max(data$mean), y = max(data$sd),
       col = 'grey35', 
       legend = c("",
                  "",
                  "",
                  ""),
       bty='n',
       pch=c(21,21,24,24),
       bg='grey35',
       title.adj = 0.2)

dev.off()

### Model describing for Heigh and log(Heigh)
w = 2*height_mean.dt$N
log_h_model = lm(data=height_mean.dt, sd_log_h ~ log_h, weights = w)
summary(log_h_model)
h_model = lm(data=height_mean.dt, sd_h ~ h, weights = w)
summary(h_model)

