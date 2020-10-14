###
# plot figure 1
###

library(data.table)
library(ggplot2)
library(ggpubr)

QC <- function(t, ss = F){
  t1 <- t[abs(mean(`Mean Height`)-`Mean Height`) <= 
            3*sd(`Mean Height`) &
            abs(mean(`SD Height`)-`SD Height`) <= 3*sd(`SD Height`) &
            if(ss){`Sample Size` >= 1000}else{T},]
  return(t1)
}

plotSDCV <- function(t1, weight = F, color = F){
  t1$CV <- t1[, `SD Height`/`Mean Height`]
  if(weight == T){
    sum_sd_w <-summary(t1[,lm(`SD Height` ~ `Mean Height`, weights = `Sample Size`)])
    sum_cv_w <- summary(t1[,lm(CV ~ `Mean Height`, 
                               weight = `Sample Size`)])
    ggSD_w <- ggplot(t1, aes(x = `Mean Height`, y = `SD Height`, 
                             weight = `Sample Size`, size = `Sample Size`)) + 
      theme_bw() + geom_smooth(method='lm',formula=y~x, size = 1,
                               color = "royalblue") + 
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
            title = element_text(size = 11), 
            axis.text.y = element_text(angle = 90, hjust = 0.5)) +
      labs(y = "SD", x = "Mean height (cm)", title = bquote("Female weighted   "~"k="~
                                                              .(formatC(sum_sd_w$coefficients[2,1], 
                                                                        digits = 2))~", "~k[p.value]~"="~
                                                              .(formatC(sum_sd_w$coefficients[2,4], digits = 1, format = "e"))
                                                            ~", "~R^2~"="~.(formatC(sum_sd_w$r.squared, 
                                                                                    digits = 2)))) +
      geom_hline(yintercept=t1[, weighted.mean(`SD Height`, `Sample Size`)], 
                 linetype="dashed")
    ggCV_w <- ggplot(t1, aes(x = `Mean Height`, y = `SD Height`/`Mean Height`, 
                             weight = `Sample Size`, size = `Sample Size`)) + theme_bw()+
      geom_smooth(method='lm',formula=y~x, size = 1, color = "royalblue") + 
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
            title = element_text(size = 11), 
            axis.text.y = element_text(angle = 90, hjust = 0.5)) +
      scale_y_continuous(breaks=c(0.0375, 0.0425)) +
      labs(y = "CV", x = "Mean height (cm)", title = bquote("Female weighted   "~"k="~
                                                              .(formatC(sum_cv_w$coefficients[2,1], 
                                                                        digits = 2))~", "~k[p.value]~"="~
                                                              .(formatC(sum_cv_w$coefficients[2,4], digits = 1, format = "e"))
                                                            ~", "~R^2~"="~.(formatC(sum_cv_w$r.squared, 
                                                                                    digits = 2)))) +
      geom_hline(yintercept=t1[, weighted.mean(`SD Height`/`Mean Height`, 
                                               `Sample Size`)], 
                 linetype="dashed")
    if(color == T){
      ggSD_w <- ggSD_w + geom_point(aes(color = Sample))
      ggCV_w <- ggCV_w + geom_point(aes(color = Sample))
    } else{
      ggSD_w <- ggSD_w + geom_point()
      ggCV_w <- ggCV_w + geom_point()
    }
    
    ggarrange(ggSD_w, ggCV_w, ncol=2, labels = c("A", "B"))
  }
  
  sum_sd <-summary(t1[,lm(`SD Height` ~ `Mean Height`)])
  sum_cv <- summary(t1[,lm(CV ~ `Mean Height`)])
  ggSD <- ggplot(t1, aes(x = `Mean Height`, y = `SD Height`)) + 
    theme_bw() + geom_smooth(method='lm',formula=y~x, size = 1,
                             color = "royalblue") + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
          title = element_text(size = 11), 
          axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    labs(y = "SD", x = "Mean height (cm)", title = bquote("k="~.(formatC(sum_sd$coefficients[2,1], 
                                                                         digits = 2))~","~SE~"="~
                                                            .(formatC(sum_sd$coefficients[2,2], digits = 2))
                                                          ~","~R^2~"="~.(formatC(sum_sd$r.squared, 
                                                                                 digits = 2))~"(***)")) +
    geom_hline(yintercept=t1[, mean(`SD Height`)], linetype="dashed")
  print(sum_sd)
  print(sum_cv)
  ggCV <- ggplot(t1, aes(x = `Mean Height`, y = `SD Height`/`Mean Height`)) + 
    theme_bw() + geom_smooth(method='lm',formula=y~x, size = 1, 
                             color = "royalblue") + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
          title = element_text(size = 11), 
          axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    scale_y_continuous(breaks=c(0.0375, 0.0425)) +
    labs(y = "CV", x = "Mean height (cm)", title = bquote("k="~
                                                            .(formatC(sum_cv$coefficients[2,1], 
                                                                      digits = 2))~","~SE~"="~
                                                            .(formatC(sum_cv$coefficients[2,2], digits = 1))
                                                          ~","~R^2~"="~.(formatC(sum_cv$r.squared, 
                                                                                 digits = 2))~"(ns)")) +
    geom_hline(yintercept=t1[, mean(`SD Height`/`Mean Height`)], 
               linetype="dashed")
  if(color == T){
    ggSD <- ggSD + geom_point(aes(color = Sample))
    ggCV <- ggCV + geom_point(aes(color = Sample))
  } else{
    ggSD <- ggSD + geom_point()
    ggCV <- ggCV + geom_point()
  }
  
  return(ggarrange(ggSD, ggCV, ncol=2, labels = c("A", "B"), 
                   font.label = list(size = 18, face="plain")))
  
}

# t_dev_f <- fread("/home/ikuznetsov/height/Tables/Table1_PlosONE1.txt")
t_dev_f <- fread(".txt") # table with height information for figure 1 A and B
t_dev_f$`Sample Size` <- t_dev_f$`Sample Size`*1000
t_dev_f <- t_dev_f[-1,]
t_dev_f$Sample <- "Developing_female"
t_dev_f_qc <- QC(t_dev_f, ss = T)

ab1 <- plotSDCV(t_dev_f_qc, weight = T)

####################################################
####################################################
####################################################

dw <- fread(".csv") # table with height information for figure 1 C
colnames(dw) <- c('Country', 'Male', 'Female', 'Age')
dw_mf <- dw[Male != "-" & Female != "-" & Age != "-",]
dw_mf[, c("Male", "Female") := list(as.numeric(sapply(strsplit(Male," "), `[`, 1)),
                                    as.numeric(sapply(strsplit(Female," "), `[`, 1)))]
write.table(dw_mf, "m_f_height.tsv", sep = "\t", quote = F, row.names = F)
dw_mf <- dw_mf[abs(Male - mean(Male)) <= 3*sd(Male) & 
                 abs(Female - mean(Female)) <= 3*sd(Female),]
dw_mf <- dw_mf[c(1,2,3,5,8,9,11:13,17,19,22,25,30:33,35:37,39,42,43,
                 48,49,54,55,58,59,61,62,65,67:70,73,74,76:81,
                 83,84,87,88,91,93,95,97,100,102:104,106:114,
                 116,120:125,128,130,132,135:138,140),]

lm_tel <- lm(Male ~ Female, data = dw_mf)

m_pred <- predict(lm_tel)
SE <- sqrt(sum((as.vector(m_pred)-dw_mf$Male)^2)/
             (nrow(dw_mf)-2))/
  sqrt(sum((dw_mf$Female-
              mean(dw_mf$Female))^2))
coef <- summary(lm_tel)$coefficients[2,1]
pvalue <- pt(q = (coef - 1)/SE, df = nrow(dw_mf) - 2, 
             lower.tail = FALSE) * 2
sum_lm <- summary(lm_tel)
print(c(coef, pvalue))

p <- ggplot(dw_mf, aes(x = Female, y = Male)) +
  geom_point() + theme_bw() + theme(axis.text = element_text(size = 12),
                                    axis.title = element_text(size = 14), 
                                    title = element_text(size = 11), 
                                    axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  geom_smooth(method='lm',formula=y~x, show.legend = F, color = "royalblue") +
  labs(x = "Mean female height (cm)", y = "Mean male height (cm)", 
       title = bquote("k="~.(formatC(sum_lm$coefficients[2,1], 
                                     digits = 3))~","~SE~"="~
                        .(formatC(sum_lm$coefficients[2,2], digits = 2))
                      ~","~R^2~"="~.(formatC(sum_lm$r.squared, 
                                             digits = 2))~"(*)"))
c1 <- ggarrange(p, labels = c("C"), font.label = list(size = 18, face="plain"))

pdf("figure1abc.pdf", width = 12, height = 4)
ggarrange(ab1, c1, widths = c(2,1))
dev.off()
