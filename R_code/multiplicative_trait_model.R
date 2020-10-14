###
# multiplicative and additional polygenic height score simulations
###

library(data.table)
library(pROC)

get_r2 <- function(real_value, fit_value){
  r_sq_test_set = 1 - var(real_value - fit_value)/var(real_value)
  return(r_sq_test_set)
}

# Sample generation
generateSample <- function(n_sample, n_snp,  AF, real_coef, 
                           log_mean, log_sd_residual){
  
  # Genotypes generation
  genotypes <- t(matrix(rbinom(n_sample*n_snp, 2, AF), nrow=n_snp, ncol=n_sample))
  colnames(genotypes) <- paste0("SNP", 1:n_snp)
  
  # Genetic component of the trait
  log_trait_gen <- apply(genotypes, 1, function(x) sum(x*real_coef)) + log_mean
  trait_gen <- 10^log_trait_gen

  # Environmental component of the trait addition
  log_trait <- log_trait_gen + rnorm(n_sample, 0, log_sd_residual)
  trait <- 10^log_trait
  print(paste("h2_trait =", cor(trait_gen, trait)^2))
  print(paste("h2_log_trait =", cor(log_trait_gen, log_trait)^2))

  # Mean and SD of the trait
  print(paste("mean =", mean(trait)))
  print(paste("SD =", sd(trait)))
  
  return(list(genotypes = genotypes, trait = trait, log_trait = log_trait,
              trait_gen = trait_gen, log_trait_gen = log_trait_gen))
}


# Model parameters
log_sd <- 0.0164
log_mean = 2.2107
h2 <- 0.4
n_snp <- 300
n_sample <- 25000

# Distribution parameters
log_var <- log_sd^2
log_gen_var <- h2*log_var
log_sd_av_snp <- sqrt(log_gen_var/n_snp)
log_sd_residual <- sqrt((1-h2)*log_var)
AF <- runif(n_snp, min = 0.05, max = 0.95)
real_coef <-  runif(n_snp, -sqrt(3*log_sd_av_snp**2/(2*AF*(1-AF))),
                    sqrt(3*log_sd_av_snp**2/(2*AF*(1-AF))))

sample_train <- generateSample(n_snp = n_snp,
                               n_sample = n_sample, 
                               AF = AF,
                               real_coef = real_coef,
                               log_mean = log_mean,
                               log_sd_residual = log_sd_residual)

sample_test <- generateSample(n_snp = n_snp,
                               n_sample = n_sample, 
                               AF = AF,
                               real_coef = real_coef,
                               log_mean = log_mean+0.01,
                               log_sd_residual = log_sd_residual)
sample_test2 <- generateSample(n_snp = n_snp,
                              n_sample = n_sample, 
                              AF = AF,
                              real_coef = real_coef,
                              log_mean = log_mean,
                              log_sd_residual = log_sd_residual)


# Linear and Log scores calculations
l_lin <- lm(sample_train$trait ~ sample_train$genotypes)
l_log <- lm(sample_train$log_trait ~ sample_train$genotypes)

score_lin <- apply(sample_test$genotypes, 1, function(x) sum(x*l_lin$coefficients[-1]))
score_log <- apply(sample_test$genotypes, 1, function(x) sum(x*l_log$coefficients[-1]))

score_lin2 <- apply(sample_test2$genotypes, 1, function(x) sum(x*l_lin$coefficients[-1]))
score_log2 <- apply(sample_test2$genotypes, 1, function(x) sum(x*l_log$coefficients[-1]))


cal_score_lin_model_lin <- lm(sample_test$trait ~ score_lin)
cal_score_lin_model_log <- lm(sample_test$log_trait ~ score_lin)
cal_score_log_model_lin <- lm(sample_test$trait ~ score_log)
cal_score_log_model_log <- lm(sample_test$log_trait ~ score_log)

cal_score_lin_model_lin2 <- lm(sample_test2$trait ~ score_lin2)
cal_score_lin_model_log2 <- lm(sample_test2$log_trait ~ score_lin2)
cal_score_log_model_lin2 <- lm(sample_test2$trait ~ score_log2)
cal_score_log_model_log2 <- lm(sample_test2$log_trait ~ score_log2)

score_lin_model_lin <- cal_score_lin_model_lin$coefficients[1] +
  score_lin*cal_score_lin_model_lin$coefficients[2]
score_lin_model_log <- cal_score_lin_model_log$coefficients[1] +
  score_lin*cal_score_lin_model_log$coefficients[2]
score_log_model_lin <- cal_score_log_model_lin$coefficients[1] +
  score_log*cal_score_log_model_lin$coefficients[2]
score_log_model_log <- cal_score_log_model_log$coefficients[1] +
  score_log*cal_score_log_model_log$coefficients[2]

score_lin_model_lin2 <- cal_score_lin_model_lin$coefficients[1] +
  score_lin2*cal_score_lin_model_lin$coefficients[2]
score_lin_model_log2 <- cal_score_lin_model_log$coefficients[1] +
  score_lin2*cal_score_lin_model_log$coefficients[2]
score_log_model_lin2 <- cal_score_log_model_lin$coefficients[1] +
  score_log2*cal_score_log_model_lin$coefficients[2]
score_log_model_log2 <- cal_score_log_model_log$coefficients[1] +
  score_log2*cal_score_log_model_log$coefficients[2]

scores <- data.table(sample_test$trait, score_lin_model_lin, 10^score_lin_model_log,
                     score_log_model_lin, 10^score_log_model_log)
colnames(scores)[c(1,3,5)] <- c("trait", "score_lin_model_log", "score_log_model_log")

scores2 <- data.table(sample_test2$trait, score_lin_model_lin2, 10^score_lin_model_log2,
                     score_log_model_lin2, 10^score_log_model_log2)
colnames(scores2) <- c("trait", "score_lin_model_lin", "score_lin_model_log", 
                       "score_log_model_lin", "score_log_model_log")
