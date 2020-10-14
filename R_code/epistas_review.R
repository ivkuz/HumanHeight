###
# Epistasis estimation for standard and powerful PGHSes
###

library(data.table)
library(stringr)
library(R.filesets)
library(IRanges)

get_gen_and_phen_data <- function(filename){
  gen_and_phen = loadRDS(filename)
  rs_variables = colnames(gen_and_phen)[str_detect(colnames(gen_and_phen), "rs")]
  rs_ids = gsub("(.+)\\_(.+)\\_(.+)", "\\1",rs_variables)
  colnames(gen_and_phen)[str_detect(colnames(gen_and_phen),
                                    "rs")] = rs_ids
  predF3 = phen[,.(ID, predF3)]
  gen_and_phen <- merge(gen_and_phen, predF3, by = "ID")
  return(gen_and_phen)
}

add_locusID_to_table <- function(chr_subset, window=10^6){
  rs_pos = IRanges(chr_subset$Pos, chr_subset$Pos + window)
  locus_range = reduce(rs_pos)
  rs_overlaps_locus = findOverlaps(rs_pos, locus_range)
  chr=unique(chr_subset$Chr)
  chr_subset$locus_id = paste(chr,rs_overlaps_locus@to, sep = '_') 
  return(chr_subset)
}

format_snp_info_table <- function(snp_info_table){
  snp_format_table <- ldply(unique(snp_info_table$Chr), function(chr, df){
    add_locusID_to_table(df[df$Chr == chr])},
    df = snp_info_table)
  snp_format_table$N_ID <- 1:nrow(snp_info_table)
  return(snp_format_table)
}

get_score <- function(rs_ids, gen_and_phen_matrix, snp_table){
  return(gen_and_phen_matrix[,rs_ids] %*% 
           snp_table[match(rs_ids,  rs_id)]$Beta)
}


get_p_values <- function(gen_and_phen, gen_score_matrix, snp_table){
  l = length(unique(snp_table$locus_id))
  loci_subset = sample(unique(snp_table$locus_id), as.integer(l/2))
  snp_subset = snp_table[locus_id %in% loci_subset]$rs_id
  gen_and_phen$score1 = rowSums(gen_score_matrix[,snp_subset])
  gen_and_phen$score2 = gen_and_phen$score - gen_and_phen$score1
  
  model_data = gen_and_phen[, .(score1 = (score1 - mean(score1)) / sd(score1),
                                score2 = (score2 - mean(score2)) / sd(score2),
                                Sex = Sex,
                                predF3 = predF3, 
                                Height=Height, 
                                log_height=log_height), 
                            by = analysis_group]
  
  p_h = summary(lm(Height ~ Sex + predF3 + score1 + score2 + score1*score2, 
                   data = model_data))$coeff["score1:score2", "Pr(>|t|)"]
  p_log = summary(lm(log_height ~ Sex + predF3 + score1*score2,
                     data = model_data))$coeff["score1:score2", "Pr(>|t|)"]
  return(c("p_H" = p_h, "p_logH" = p_log))
}


runSimulation <- function(gen_and_phen, snp_table, rs_ids, n){
  set.seed(1)
  gen_matrix = as.matrix(gen_and_phen[,..rs_ids])
  gen_score_matrix = gen_matrix * matrix(rep(snp_table[match(rs_ids, 
                                                             rs_id)]$Beta,
                                             each=nrow(gen_matrix)),
                                         nrow = nrow(gen_matrix))
  gen_and_phen$score = rowSums(gen_score_matrix)
  p_vals = replicate(n = n, get_p_values(gen_and_phen, gen_score_matrix,
                                         snp_table)
                     )
  p_dt = as.data.table(t(p_vals))
  return(p_dt)
}


make_test <- function(p_values, p_treschold=0.05){
  n_success <- sum(p_values <  p_treschold)
  bt <- formatC(binom.test(n_success, n=length(p_values), p =  p_treschold,
                           alternative = "two.sided")$p.value,
                format = "e", digits = 2)
  c(p_treschold, n_success, bt)
}


## =========== Start Pipeline ========

#load data
if (PGHS==1){
  gen_and_phen <- get_gen_and_phen_data(filename = "PGRS1.rds")
  snp_info_table <- fread('.tsv') # Table with phenotypes and genotypes for SNPs from Wood et al.
  snp_supp <- fread("../ST1_Wood.tsv",
                  select = c("SNP", "Chr", "Position(bp)"))
  snp_info_table <- merge(snp_info_table, snp_supp, 
                          by.y = 'SNP', by.x = "rs_id")
  snp_info_table <- snp_info_table[,.(rs_id, other_allele, effect_allele,
                                      Beta_GIANT, Chr,
                                      `Position(bp)`)]
if (PGHS==2){
  f_gen_phen = ".rds" # Table with phenotypes and genotypes for SNPs from Yengo et al.
  gen_and_phen <- get_gen_and_phen_data(filename = f_gen_phen)
  snp_info_table <- fread('.txt') # Yengo et al. summary statistics
  snp_info_table <- snp_info_table[,.(SNP, Other_Allele, Tested_Allele,
                                      Beta_СOJO_sign, CHR,
                                      POS)]
}
  

rs_ids = colnames(gen_and_phen)[str_detect(colnames(gen_and_phen), "rs")]
colnames(snp_info_table) <- c("rs_id", "OA", "EA", "Beta", "Chr", "Pos")
snp_info_format_table <- data.table(format_snp_info_table(snp_info_table))
snp_info_format_table <- snp_info_format_table[ snp_info_format_table$rs_id %in% 
                                                  rs_ids]

results <- runSimulation(gen_and_phen[analysis_group == 'English',], 
                         snp_info_format_table, rs_ids,
                         1000)
## output is the resuts of binominal test for height and logheight p-values
apply(results, 2, make_test, p_treschold=0.05)

## сheck significance of square of full score
summary(lm(Height ~ Sex + predF3 + score + I(score^2),
           data = gen_and_phen[analysis_group == 'English',]))$coeff
summary(lm(log_height ~ Sex + predF3 + score + I(score^2),
           data = gen_and_phen[analysis_group == 'English',]))$coeff
