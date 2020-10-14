###
# get genotype matrix and mfi table
###

library(rbgen)
library(data.table)

local_path = ""
sample_path = "" # UK Biobank sample table

giant_table <- fread(paste(local_path,'.csv', sep = ""), select=c(1:10,17)) # Wood et al. table
giant_table <- giant_table[order(`P-value`, decreasing = F)]
giant_table <- unique(giant_table, by = "Locus ID")

# 423 SNPs top SNPs from 423 loci

# remove 4 SNP which are either absent or duplicated in UKBB
bim_table = fread(paste(local_path, "extracted.bim", sep = ""))
absent_snps = giant_table[!SNP %in%  bim_table$V2]$SNP
duplicated_snps = bim_table[duplicated(bim_table, by = "V2")]$V2
bad_snps = c(absent_snps, duplicated_snps)
giant_table = giant_table[!SNP %in% bad_snps]
# 419 SNPs remain
rm(bim_table, absent_snps, bad_snps, duplicated_snps)

rs_ids = giant_table$SNP
sample_ids = fread(sample_path, skip = 2)$V1

bgen_loaded = bgen.load(filename = paste(local_path, "extracted.bgen", sep = ''), rsids = rs_ids)
print('bgen file has been loaded')

get_genotype_matrix <- function(rs.ids, bgen.loaded=bgen_loaded){
  genotypes = matrix(0, nrow=length(bgen.loaded$samples), ncol=length(rs.ids))
  colnames(genotypes) = rs.ids
  rownames(genotypes) = sample_ids
  for (rs.id in rs.ids){
    genotypes[,rs.id] = bgen.loaded$data[rs.id, ,"g=1"] * 1 + bgen.loaded$data[rs.id, ,"g=2"] * 2
  }
  # add effect allele to col_name
  colnames(genotypes) = sapply(rs.ids, function(rs.id) {
    paste(rs.id, 
          bgen_loaded$variants[rs.id,c("allele0")], 
          bgen_loaded$variants[rs.id,c("allele1")], sep="_") 
    })
  return(genotypes)
}

genotype_matrix = get_genotype_matrix(rs_ids)
print('genotype matrix has been calculated')
saveRDS(genotype_matrix, paste(local_path, "genotype_matrix.rds", sep = ""))

# calculate EAF (based on samples with)
processed_phenotypes = fread(paste(local_path, ".tsv", sep = "")) # table with precessed phenotypes

EAF = apply(genotype_matrix[rownames(genotype_matrix) %in% processed_phenotypes$ID, ], 
            MARGIN = 2, mean)/2

#read table with imputaion quality 
mfi_table = fread(paste(local_path, 'mfi.txt', sep = ""))
mfi_table = mfi_table[V2 %in% rs_ids]
mfi_table[,c("V1","V3","V6", "V7"):=NULL]
setnames(mfi_table,old = c("V2", "V4", "V5", "V8"), 
         new=c("rs_id","other_allele", "effect_allele", "imputation_quality"))


if (sum(mfi_table[match(rownames(bgen_loaded$variants), mfi_table$rs_id)]$effect_allele 
  == bgen_loaded$variants[,"allele1"]) <  length(rs_ids) ){
  stop("Mess withh alleles order")
}
# add EAF from UKB and GIANT
good_part = merge(mfi_table, giant_table[,c("SNP", "Other\nAllele" , "Effect\nAllele", 
                                            "Beta", "Frequency \n(effect allele)")],
                  by.x=c("rs_id", "other_allele", "effect_allele"), 
                  by.y = c("SNP", "Other\nAllele" , "Effect\nAllele"))
reverse_part = merge(mfi_table, giant_table[,c("SNP", "Other\nAllele" , "Effect\nAllele", 
                                               "Beta", "Frequency \n(effect allele)")], 
                     by.x=c("rs_id", "other_allele", "effect_allele"), 
                     by.y = c("SNP", "Effect\nAllele", "Other\nAllele" ))
reverse_part$Beta = reverse_part$Beta*-1
reverse_part$`Frequency \n(effect allele)` = 1 - reverse_part$`Frequency \n(effect allele)`

mfi_table = rbind(good_part, reverse_part)
setnames(mfi_table, old = c("Frequency \n(effect allele)", "Beta"), 
         new=c("EAF_GIANT", "Beta_GIANT"))
mfi_table = mfi_table[match(rs_ids, mfi_table$rs_id)]
mfi_table$EAF_UKBB  = EAF
mfi_table$rs_id_O_E = names(EAF)

write.table(mfi_table, paste(local_path, "imp_and_beta_giant.tsv", sep = ""), 
            quote = F, row.names = F, sep = '\t')

