###
# Comparing of models in the groups of short and tall individuals
###


library(data.table)
library(pROC)

calibratePGRS2 <- function(genotypes, centre = NA, analysis_group = NA){
  if(is.na(centre)){centre <- unique(genotypes$Centre)}
  if(is.na(analysis_group)){analysis_group <- unique(genotypes$analysis_group)}
  genotypes <- genotypes[Centre %in% centre &
                           analysis_group %in% analysis_group,]
  rs_ids <- colnames(genotypes)[startsWith(colnames(genotypes), "rs")]
  
  genotypes_m <- genotypes[Sex == "1"]
  genotypes_f <- genotypes[Sex == "0"]
  
  fmlaH <- as.formula(paste("Height ~ ", paste(rs_ids, collapse= "+")))
  fmlaLH <- as.formula(paste("log_height ~ ", paste(rs_ids, collapse= "+")))
  
  # female
  lm_height_f <- lm(formula = fmlaH, data = genotypes_f)
  lm_log_height_f <- lm(formula = fmlaLH, data = genotypes_f)
  
  # male
  lm_height_m <- lm(formula = fmlaH, data = genotypes_m)
  lm_log_height_m <- lm(formula = fmlaLH, data = genotypes_m)
  
  lmPGRS2 <- list(h_f = lm_height_f, lh_f = lm_log_height_f,
                  h_m = lm_height_m, lh_m = lm_log_height_m)

  return(lmPGRS2)
  
}


applyPGRS2 <- function(lmPGRS2, genotypes, centre = NA, analysis_group = NA){
  if(is.na(centre)){centre <- unique(genotypes$Centre)}
  if(is.na(analysis_group)){analysis_group <- unique(genotypes$analysis_group)}
  genotypes <- genotypes[Centre %in% centre &
                           analysis_group %in% analysis_group,]
  genotypes_f <- genotypes[Sex == "0"]
  genotypes_m <- genotypes[Sex == "1"]
  
  scores_h_f <- as.matrix(
    data.table(rep(1, nrow(genotypes_f)), 
               genotypes_f[,21:ncol(genotypes_f)])) %*%
    lmPGRS2$h_f$coefficients
  
  scores_lh_f <- as.matrix(
    data.table(rep(1, nrow(genotypes_f)), 
               genotypes_f[,21:ncol(genotypes_f)])) %*%
    lmPGRS2$lh_f$coefficients
  
  scores_h_m <- as.matrix(
    data.table(rep(1, nrow(genotypes_m)), 
               genotypes_m[,21:ncol(genotypes_m)])) %*%
    lmPGRS2$h_m$coefficients
  
  scores_lh_m <- as.matrix(
    data.table(rep(1, nrow(genotypes_m)), 
               genotypes_m[,21:ncol(genotypes_m)])) %*%
    lmPGRS2$lh_m$coefficients
  
  scores <- list(h_f = as.vector(scores_h_f), lh_f = as.vector(scores_lh_f), 
                 h_m = as.vector(scores_h_m), lh_m = as.vector(scores_lh_m))
  
  return(scores)
}


characterizeModel <- function(genotypes, scores, centre = NA, analysis_group = NA){
  if(is.na(centre) | centre == "All"){centre <- unique(genotypes$Centre)}
  if(is.na(analysis_group)){analysis_group <- unique(genotypes$analysis_group)}
  genotypes <- genotypes[Centre %in% centre &
                           analysis_group %in% analysis_group,]
  genotypes_f <- genotypes[Sex == "0"]
  genotypes_m <- genotypes[Sex == "1"]
  
  s <- list(hhf = "h_f", lhf = "h_f",
                     hhm = "h_m", lhm = "h_m",
                     hlhf = "lh_f", llhf = "lh_f",
                     hlhm = "lh_m", llhm = "lh_m")
  genotypes_response <- list(hhf = ifelse(genotypes_f$Height >= 172, 1, 0), 
                             lhf = ifelse(genotypes_f$Height <= 152, 1, 0), 
                             hhm = ifelse(genotypes_m$Height >= 185, 1, 0), 
                             lhm = ifelse(genotypes_m$Height <= 165, 1, 0), 
                             hlhf = ifelse(genotypes_f$log_height >= log10(172), 1, 0), 
                             llhf = ifelse(genotypes_f$log_height <= log10(152), 1, 0), 
                             hlhm = ifelse(genotypes_m$log_height >= log10(185), 1, 0), 
                             llhm = ifelse(genotypes_m$log_height <= log10(165), 1, 0)
                             )
  genotypes_real <- list(hhf = genotypes_f$Height,
                             lhf = genotypes_f$Height,
                             hhm = genotypes_m$Height,
                             lhm = genotypes_m$Height,
                             hlhf = genotypes_f$log_height,
                             llhf = genotypes_f$log_height,
                             hlhm = genotypes_m$log_height,
                             llhm = genotypes_m$log_height
                         )
  genotypes_real_r2 <- list(hhf = genotypes_f$Height[genotypes_f$Height >= 172], 
                         lhf = genotypes_f$Height[genotypes_f$Height <= 152], 
                         hhm = genotypes_m$Height[genotypes_m$Height >= 185], 
                         lhm = genotypes_m$Height[genotypes_m$Height <= 165], 
                         hlhf = genotypes_f$log_height[genotypes_f$Height >= 172], 
                         llhf = genotypes_f$log_height[genotypes_f$Height <= 152], 
                         hlhm = genotypes_m$log_height[genotypes_m$Height >= 185], 
                         llhm = genotypes_m$log_height[genotypes_m$Height <= 165]
  )
  scores_r2 <- list(hhf = scores[[s[["hhf"]]]][genotypes_f$Height >= 172], 
                    lhf = scores[[s[["lhf"]]]][genotypes_f$Height <= 152], 
                    hhm = scores[[s[["hhm"]]]][genotypes_m$Height >= 185], 
                    lhm = scores[[s[["lhm"]]]][genotypes_m$Height <= 165], 
                    hlhf = scores[[s[["hlhf"]]]][genotypes_f$Height >= 172], 
                    llhf = scores[[s[["llhf"]]]][genotypes_f$Height <= 152], 
                    hlhm = scores[[s[["hlhm"]]]][genotypes_m$Height >= 185], 
                    llhm = scores[[s[["llhm"]]]][genotypes_m$Height <= 165]
  )
  
  r_2 <- list()
  rauc <- list()
  
  ###R2 for high or low samples
  for(thr in names(genotypes_response)){
    if(startsWith(s[[thr]], "l")){
      r_2[thr] <- cor(10^(genotypes_real_r2[[thr]]), 10^(scores_r2[[thr]]))^2
    }else{
      r_2[thr] <- cor(genotypes_real_r2[[thr]], scores_r2[[thr]])^2
    }
    rauc[thr] <- auc(roc(genotypes_response[[thr]], scores[[s[[thr]]]], quiet = TRUE))
  }
  
  return(list(r2 = r_2, auc = rauc))
}


genotypes <- fread(".tsv") # table with genotypes


asses_centres <- as.data.table(read.csv("assessment_centre.csv")) # pivot table on assessment sentres data

results <- list()
for(analysis_gr in c(unique(genotypes$analysis_group), "All")){

  print(paste0("analysis_group = ", analysis_gr))
  results[[analysis_gr]] <- list()

  for(centre1 in as.character(unique(genotypes$Centre))){

    print(paste0("centre1 = ", centre1))
    if(asses_centres[Centre == centre1, get(paste0(analysis_gr, "_F"))] < 1000 |
       asses_centres[Centre == centre1, get(paste0(analysis_gr, "_M"))] < 1000){
      next()
    }
    # calibrate PGRS2
    lmpgrs2 <- calibratePGRS2(genotypes = genotypes,
                              centre = centre1, analysis_group = analysis_gr)
    results[[analysis_gr]][[centre1]] <- list()

    for(centre2 in as.character(unique(genotypes$Centre))){

      print(paste0("centre2 = ", centre2))
      if(asses_centres[Centre == centre2, get(paste0(analysis_gr, "_F"))] < 100 |
         asses_centres[Centre == centre2, get(paste0(analysis_gr, "_M"))] < 100){
        next()
      }
      # calculate scores
      scores <- applyPGRS2(lmPGRS2 = lmpgrs2, genotypes = genotypes,
                           centre = centre2, analysis_group = analysis_gr)
      # characterize model
      res <- characterizeModel(genotypes = genotypes, scores = scores,
                               centre = centre2, analysis_group = analysis_gr)
      results[[analysis_gr]][[centre1]][[centre2]] <- res

    }
  }
}

saveRDS(results, "results_r2.2.RDS")

print("Succsess!")
