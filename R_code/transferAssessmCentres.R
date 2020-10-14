###
# Analysis of transferability of models between assessment centres
###

library(data.table)
library(hydroTSM)

getMatrix <- function(results, asc, analysis_gr, type, model){
  nc <- nrow(asc)
  mat <- matrix(ncol = nc, nrow = nc)
  colnames(mat) <- asc[,meaning]
  rownames(mat) <- asc[,meaning]
  mat_nd <- mat
  
  for(i in 1:nc){
    for(j in 1:nc){
      v <- results[[analysis_gr]][[as.character(asc[i, Centre])]][[as.character(asc[j, Centre])]][[type]][[model]]
      if(!is.null(v)){
        mat[i, j] <- v
        if(i != j){
          mat_nd[i, j] <- v
        }
      }
    }
  }
  
  mat <- mat[rowSums(is.na(mat)) != nrow(mat), colSums(is.na(mat)) != nrow(mat)]
  
  return(mat)
}

compareModels <- function(results, asc, analysis_gr, 
                          type = "auc", height = "h", sex = "m", plots = F){
  mat_h <- getMatrix(results, asc, analysis_gr, type = type,
                     model = paste0(height, "h", sex))
  mat_lh <- getMatrix(results, asc, analysis_gr, type = type,
                      model = paste0(height, "lh", sex))
  
  mat <- mat_lh - mat_h
  
  if(plots){
    if(height == "h"){
      h <- "high"
    }else if(height == "l"){
      h <- "low"
    }
    if(sex == "m"){
      s <- "men"
    }else if(sex == "f"){
      s <- "women"
    }
    
    plot(density(mat), main = paste("f(LogH) - f(H) for",
                                    h, analysis_gr, s, type,
                                    sep = " "))
    print(matrixplot(t(mat[nrow(mat):1,]), 
                     main = paste("f(LogH) - f(H) for", 
                                  h, analysis_gr, s, type,
                                  sep = " ")))
    heatmap(mat, main = paste("f(LogH) - f(H) for",
                              h, analysis_gr, s, type,
                              sep = " "))
  }
  
  p <- wilcox.test(mat_h, mat_lh)$p.value
  meanv <- mean(mat)
  medv <- median(mat)
  
  return(formatC(c(p, meanv, medv), format = "e", digits = 1))
}


results <- readRDS(".RDS") # table with R2

asc <- as.data.table(read.csv(".csv")) # table with information on assessment centres
asc <- setorderv(asc, cols = "All_F", order = -1)
asc <- asc[order(match(asc[,country], c("England", "Scotland", "Wales")))]


res_table <- data.frame("P-value" = NA, "Mean" = NA, "Median" = NA)
pdf("compare_models2.pdf")
for(analysis_gr in c("English", "Scottish", "All")){
  for(h in c("h", "l")){
    for(s in c("f", "m")){
      for(type in c("auc", "r2")){
        res <- compareModels(results, asc, analysis_gr, type = type, 
                             height = h, sex = s, plots = T)
        res_table <- rbind(res_table, res)
        rownames(res_table)[nrow(res_table)] <- paste(analysis_gr, 
                                                      h, s, type, 
                                                      sep = "_")
      }
    }
  }
}
dev.off()

res_table <- res_table[-1,]

###
res_table <- res_table[-c(2, 4, 10, 12, 18, 20),]
rownames(res_table)[c(4, 6, 10, 12, 16, 18)] <- c("English_f_r2", "English_m_r2", "Scottish_f_r2", "Scottish_m_r2", "All_f_r2", "All_m_r2")
res_table <- res_table[c(1:3, 5, 4, 6, 7:9, 11, 10, 12, 13:15, 17, 16, 18),]



write.table(res_table, "compare_models2.tsv", quote = F, sep = "\t")
