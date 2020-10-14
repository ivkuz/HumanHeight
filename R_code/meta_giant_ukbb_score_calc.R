###
# PGHS2 (Yengo et al.) calculation
library(rbgen)
library(data.table)

get_genotype_matrix <- function(rs.ids, bgen.loaded){
  genotypes = matrix(0, nrow=length(bgen.loaded$samples), ncol=length(rs.ids))
  colnames(genotypes) = rs.ids
  rownames(genotypes) = sample_ids
  for (rs.id in rs.ids){
    genotypes[,rs.id] = bgen.loaded$data[rs.id, ,"g=1"] * 1 + bgen.loaded$data[rs.id, ,"g=2"] * 2
  }
  # add effect allele to col_name
  colnames(genotypes) = sapply(rs.ids, function(rs.id) {
    paste(rs.id, 
          bgen.loaded$variants[rs.id,c("allele0")], 
          bgen.loaded$variants[rs.id,c("allele1")], sep="_") 
  })
  return(genotypes)
}

get_score = function(rs_chunk){
  
  bgen_loaded = bgen.load(filename = paste(local_path, prefix, ".bgen", sep = ''), rsids = rs_chunk)
  print('bgen file has been loaded')
  
  genotype_matrix = get_genotype_matrix(rs_chunk, bgen.loaded = bgen_loaded)
  print('genotype matrix has been calculated')
  
  
  rs_ids_O_E = colnames(genotype_matrix)
  
  eff_all_dt = data.table(SNP = rs_chunk, ef_allele = substring(rs_ids_O_E, nchar(rs_ids_O_E)))
  table_for_beta = merge(eff_all_dt, meta_table, by = "SNP")
  table_for_beta[,BETA_COJO := ifelse(ef_allele == Tested_Allele, BETA_COJO, -BETA_COJO) ]
  
  
  betas = table_for_beta$BETA_COJO
  
  genetic_scores = genotype_matrix %*% betas
  return( genetic_scores[,1])
} 


current.user = Sys.info()['user']
prefix = "meta_3k_snps"

local_path = ""
sample_path = "" # UK Biobank .sample table
sample_ids = fread(sample_path, skip = 2)$V1

meta_table = fread(paste(local_path, "Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt", sep = "")) # Top summary statistics Yengo et al.

rs_ids = meta_table$SNP

rs_ids_list = split(rs_ids, ceiling(seq_along(rs_ids)/500))


test = get_score(rs_ids_list[[15]])


score_table = sapply(rs_ids_list, get_score)

genetic_scores = apply(score_table, 1, sum)
genetic_scores = data.table(ID = sample_ids, score_meta = genetic_scores)



processed_phenotypes = readRDS(paste(local_path, 'phen.rds', sep = ""))
gen_and_phen = merge(processed_phenotypes, genetic_scores, by ="ID")



summary(lm(data = gen_and_phen[Sex==1], formula = Height ~ score_meta)) 
summary(lm(data = gen_and_phen[Sex==1], formula = log_height ~ score_meta)) 
