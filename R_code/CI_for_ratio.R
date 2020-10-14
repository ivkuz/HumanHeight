###
# Confidence Intervals for mean and SD ratio
###

FiellerRatioCI_basic <- function(a,b,V,alpha=0.05){
  theta <- a/b
  v11 <- V[1,1]
  v12 <- V[1,2]
  v22 <- V[2,2]
  
  z <- qnorm(1-alpha/2)
  g <- z*v22/b^2
  C <- sqrt(v11/a^2 - 2*v12/(a*b) + v22/b^2 - g/a^2*(v11-v12^2/v22))
  minS <- (1/(1-g))*(theta- g*v12/v22 - z * theta * C)
  maxS <- (1/(1-g))*(theta- g*v12/v22 + z * theta * C)
  return(c(ratio=theta,min=minS,max=maxS))
}

delta <-  function(a,b,V,alpha=0.05){
  theta <- a/b
  v11 <- V[1,1]
  v12 <- V[1,2]
  v22 <- V[2,2]
  z <- qnorm(1-alpha/2)
  C = sqrt(v11/a^2 + v22/b^2 - 2*v12/(a*b))
  minS <- theta - z * theta * C
  maxS <- theta + z * theta * C
  print(z * theta * C)
  return(c(ratio=theta,min=minS,max=maxS))
}
  
## from table 3A

model = lm(data=phen, Height ~ Sex + norm_score + predF3 + 
                                   norm_score*Sex + predF3*Sex + norm_score*predF3)
s.model = summary(model)
coef = s.model$coefficients
a = coef['Sex:norm_score','Estimate'] #SEX/PGHS
b = coef['norm_score','Estimate'] #PGHS
N =dim(phen)[1]
SD.m = s.model$cov.unscaled[c('Sex:norm_score','norm_score'),
c('Sex:norm_score','norm_score')] * s.model$sigma^2

FiellerRatioCI_basic(a,b, SD.m)
delta(a,b, SD.m)



get_CI_for_coef_ration = function(x1, x2, s.model) {
  coef = s.model$coefficients
  a = coef[x1,'Estimate'] #SEX/PGHS
  b = coef[x2,'Estimate'] #PGHS
  SD.m = s.model$cov.unscaled[c(x1,x2),
                              c(x1,x2)] * s.model$sigma**2
  FiellerRatioCI_basic(a,b, SD.m)
}

get_CI_for_coef_ration("Sex:norm_score", "norm_score", s.model )
get_CI_for_coef_ration("Sex:predF3", "predF3", s.model )

# ratio of SD
get_CI_for_SD_ratio = function(feature, analysis_gr=unique(gen_and_phen$analysis_group)) {
  d.t = gen_and_phen[analysis_group %in% analysis_gr][,.(sd=sd(Height), N= .N), by=feature]
  a = d.t[get(feature)==1]$sd
  b = d.t[get(feature)==0]$sd
  var_a = a**2 / (2 * (d.t[get(feature)==1]$N - 1))
  var_b = b**2 / (2 * (d.t[get(feature)==0]$N - 1))
  FiellerRatioCI_basic(a,b,  matrix(c(var_a, 0, 0, var_b), nrow = 2))*100-100
}

#ratio of means
get_CI_for_mean_ratio = function(feature, analysis_gr=unique(gen_and_phen$analysis_group)) {
  d.t = gen_and_phen[analysis_group %in% analysis_gr][,.(mean=mean(Height), N = .N, VAR = var(Height)), by=feature]
  a = d.t[get(feature)==1]$mean
  b = d.t[get(feature)==0]$mean
  var_a = d.t[get(feature)==1]$VAR / (d.t[get(feature)==1]$N - 1)
  var_b = d.t[get(feature)==0]$VAR / (d.t[get(feature)==0]$N - 1)
  FiellerRatioCI_basic(a,b,  matrix(c(var_a, 0, 0, var_b), nrow = 2))*100-100
}


gen_and_phen[, PGHS_binary := as.numeric(norm_score >= median(norm_score))] 
gen_and_phen[, RP_binary := as.numeric(predF3 >= median(predF3))] 


get_CI_for_SD_ratio("Sex", analysis_gr = "English")
get_CI_for_SD_ratio("PGHS_binary")
get_CI_for_SD_ratio("RP_binary")

get_CI_for_mean_ratio("Sex", analysis_gr = "English")
get_CI_for_mean_ratio("PGHS_binary")
get_CI_for_mean_ratio("RP_binary")

dt_forest = rbind(get_CI_for_mean_ratio("Sex", analysis_gr = "English"),
                            get_CI_for_mean_ratio("PGHS_binary", analysis_gr = "English"),
                            get_CI_for_mean_ratio("RP_binary", analysis_gr = "English"),
                  get_CI_for_SD_ratio("Sex", analysis_gr = "English"),
                  get_CI_for_SD_ratio("PGHS_binary", analysis_gr = "English"),
                  get_CI_for_SD_ratio("RP_binary", analysis_gr = "English"))
                             
dt_forest = as.data.table(dt_forest)                          
dt_forest$feature = c(rep("mean",3), rep("sd", 3))
dt_forest$pred = rep(c("Sex", "PGHS_binary", "RP_binary"),2)
dt_forest[order(pred)]
dt_forest$color = with(dt_forest, paste(feature, pred, sep = "_"))

library(ggplot2)
p = ggplot(data=dt_forest,
           aes(x = pred,y = ratio, ymin = min, ymax = max ))+
  geom_pointrange(aes(col=feature))+
  xlab('pred')+ ylab("(95% Confidence Interval)")+
  geom_errorbar(aes(ymin=min, ymax=max,col=color),width=0.5,cex=1)+ 
  facet_wrap(~pred,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip() 
p
